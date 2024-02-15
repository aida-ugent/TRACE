import openTSNE
import numpy as np
import os
import time
import math
import argparse
import anndata as ad
from utils import normalizeEmbedding
from numpy.random import MT19937
from numpy.random import RandomState, SeedSequence

# call eg.g. with
# python tsne.py ../data/mouse_fibro/mouseCD45neg.h5ad -output ../data/mouse_fibro/embeddings/tsne_skrodzki_0.csv -exaggeration 5 -init "X_totalVI" -use-rep "X_totalVI" -hd-metric cosine
# python tsne.py ../data/mouse_fibro/mouseCD45neg.h5ad -output ../data/mouse_fibro/embeddings/tsne_skrodzki_1.csv -exaggeration 4 -early-exag-iter 0 -init ../data/mouse_fibro/embeddings/tsne_skrodzki_0.csv  -use-rep "X_totalVI" -hd-metric cosine


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("adata", type=str, help="Path to .hd5ad anndata file.")
    parser.add_argument(
        "-output", type=str, help="Path to .csv output filepath.", default=None
    )
    parser.add_argument(
        "-sampling-frac",
        type=float,
        help="Subsampling points for faster embeddings.",
        default=0.1,
    )
    parser.add_argument(
        "-exaggeration",
        type=float,
        help="Exaggeration factor to apply after early exaggeration.",
        default=1,
    )
    parser.add_argument(
        "-early-exag-iter",
        type=int,
        help="Number of iterations with early exaggeration. ",
        default=250,
    )
    parser.add_argument(
        "-init",
        type=str,
        help="Path to .csv file with embedding to initialize the optimization.",
        default=None,
    )
    parser.add_argument(
        "-use-rep",
        type=str,
        help="Name ob adata.obsm key to use for distance calculation",
        default="X_pca",
    )
    parser.add_argument(
        "-hd-metric",
        type=str,
        help="Metric to use for distance calculation",
        default="euclidean",
    )
    args = parser.parse_args()
    return args


def prolongate_embedding(data, coarse_embedding, sampling_indices, coarse_knn_index):
    """Based on an embedding for the points in sampling_indices,
    place the remaining points on top of their nearest neighbor in the dataset.

    Args:
        data (numpy.ndarray): The dataset containing all points.
        coarse_embedding (numpy.ndarray): The embedding for the points in sampling_indices.
        sampling_indices (numpy.ndarray): The indices of the points in the sampling set.
        coarse_knn_index (knn_index): The k-nearest neighbor index for the coarse embedding.

    Returns:
        numpy.ndarray: The new embeddings with the remaining points placed on top of their nearest neighbors.
    """
    neighbors, _ = coarse_knn_index.query(
        data, k=1
    )  # the nearest neighbor is the point itself
    neighbors = neighbors.flatten()
    neighbors[sampling_indices] = np.arange(len(sampling_indices))

    # place the remaining points at the location of their nearest neighbor
    new_embeddings = np.take(coarse_embedding, neighbors, axis=0)
    return new_embeddings


def tsne_skrodzki(
    data,
    init=None,
    sampling_frac=0.1,
    exaggeration=1,
    n_iter=750,
    early_exag_iter=250,
    smoothing_iter=250,
    smoothing_perplexity=30,
    save_fpath=None,
    hd_metric="euclidean",
):
    """
    Compute a single t-SNE embedding using the sampling idea of Skrodzki et al.,
    'Tuning the perplexity for and computing sampling-based t-SNE embeddings' (2023).

    In their approach, a fraction k of the data is embedded using a perplexity
    perp = (N*k)/100. Then, all remaining points are embedded at the location
    of their nearest neighbor from the sample (prolongation). The final smoothing
    step refines the full embedding with a perplexity of 30 for a couple of iterations.

    Args:
        data (array-like): The input data to be embedded.
        init (array-like, optional): The initial embedding. If not provided, PCA initialization will be computed.
        sampling_frac (float, optional): The fraction of data to be embedded using a perplexity perp = (N * sampling_frac) / 100. Default is 0.1.
        exaggeration (float, optional): The exaggeration factor for t-SNE. Default is 1.
        n_iter (int, optional): The total number of iterations for t-SNE. Default is 750.
        early_exag_iter (int, optional): The number of iterations for early exaggeration. Default is 250.
        smoothing_iter (int, optional): The number of iterations for smoothing. Default is 250.
        smoothing_perplexity (float, optional): The perplexity for smoothing. Default is 30.
        save_fpath (str, optional): The file path to save the embedding. Default is None.
        hd_metric (str, optional): The metric to use for computing distances. Default is "euclidean".

    Returns:
        np.ndarray: The computed t-SNE embedding.

    """
    rs = RandomState(MT19937(SeedSequence(123456789)))
    data = np.asarray(data, dtype=np.float32)
    data_size = data.shape[0]
    sampling_size = math.ceil(data_size * sampling_frac)
    sample_ind = np.random.choice(data_size, size=sampling_size, replace=False)
    coarse_perp = math.ceil((data_size * sampling_frac) / 100)
    print(f"Computing affinities with perplexity {coarse_perp}...")
    # computing coarse embedding
    start_aff = time.time()
    aff_coarse = openTSNE.affinity.PerplexityBasedNN(
        data[sample_ind, :],
        perplexity=coarse_perp,
        method="annoy",
        n_jobs=8,
        random_state=rs,
        metric=hd_metric,
        verbose=True,
    )
    print("openTSNE: Coarse NN search", time.time() - start_aff, flush=True)

    # initialization
    if init is None:
        print(f"Computing PCA initialization...")
        init = openTSNE.initialization.pca(data[sample_ind, :])
    else:
        init = openTSNE.initialization.rescale(init[sample_ind, :])

    coarse_embedding = openTSNE.TSNEEmbedding(
        embedding=init,
        affinities=aff_coarse,
        n_jobs=8,
        verbose=True,
        random_state=rs,
        negative_gradient_method="fft",
    )
    
    coarse_embedding.optimize(early_exag_iter, exaggeration=12, inplace=True)
    coarse_embedding.optimize(n_iter=n_iter, exaggeration=exaggeration, inplace=True)
    print("openTSNE: Coarse embedding total", time.time() - start_aff, flush=True)


    # now need affinities for whole dataset
    print(f"Computing affinities for whole dataset with perplexity {smoothing_perplexity}...")
    aff_fine_start = time.time()
    aff_fine = openTSNE.affinity.PerplexityBasedNN(
        data,
        perplexity=smoothing_perplexity,
        n_jobs=8,
        random_state=rs,
        metric=hd_metric,
        method="annoy"
    )
    print("openTSNE: Fine NN search", time.time() - aff_fine_start, flush=True)

    fine_init = prolongate_embedding(
        data, coarse_embedding, sample_ind, aff_coarse.knn_index
    )
    fine_init = openTSNE.initialization.rescale(fine_init)

    smooth_embedding = openTSNE.TSNEEmbedding(
        embedding=fine_init,
        affinities=aff_fine,
        n_jobs=8,
        verbose=True,
        random_state=rs,
        negative_gradient_method="fft"
    )
    smooth_embedding.optimize(smoothing_iter, exaggeration=exaggeration, inplace=True)
    print("openTSNE: Fine embedding total", time.time() - aff_fine_start, flush=True)

    
    smooth_embedding = normalizeEmbedding(smooth_embedding)

    if save_fpath is not None:
        np.savetxt(save_fpath, X=smooth_embedding, delimiter=",")

    return np.asarray(smooth_embedding)


def compute_tsne_series(
    X: np.ndarray,
    max_exaggeration: int,
    fpath_prefix: str,
    hd_metric: str = "euclidean",
    init: np.ndarray = None,
    align_sequence: bool = True,
):
    """
    Compute a series of tSNE embeddings with decreasing exaggeration.

    Args:
        X (np.ndarray): The input data matrix of shape (n_samples, n_features).
        max_exaggeration (int): The maximum exaggeration value for tSNE.
        fpath_prefix (str): The file path prefix for saving the embeddings.
        hd_metric (str, optional): The metric used for high-dimensional space distance calculation (default: "euclidean").
        init (np.ndarray, optional): The initial embedding to use for tSNE (default: None).
        align_sequence (bool, optional): Whether to align (i.e. initialize with previous embedding)
            the tSNE embeddings in the sequence (default: True).

    Returns:
        embeddings (Dict[str, np.ndarray]): A dictionary containing the tSNE embeddings, where the keys contain the exaggeration and the values are the embeddings.
    """
    embeddings = {}
    init = init

    for exag in range(max_exaggeration, 0, -1):
        if exag == max_exaggeration and align_sequence:
            early_exag_iter = 250
        else:
            early_exag_iter = 0
        embedding = tsne_skrodzki(
            data=X,
            init=init,
            early_exag_iter=early_exag_iter,
            exaggeration=exag,
            save_fpath=fpath_prefix + f"_exag_{exag}.csv",
            hd_metric=hd_metric,
        )
        if align_sequence:
            init = embedding
        embeddings[f"tSNE_exag_{exag}"] = embedding
    return embeddings


if __name__ == "__main__":
    args = parse_args()
    adata = ad.read_h5ad(args.adata)

    if args.use_rep in adata.obsm_keys():
        data = adata.obsm[args.use_rep]
        print(f"Using {args.use_rep} as HD data.")
    else:
        data = np.asarray(adata.X)

    if args.init is not None and os.path.isfile(args.init):
        print(f"Using {args.init} as initialization.")
        init = np.loadtxt(args.init, delimiter=",")
    elif args.init is not None and args.init in adata.obsm_keys():
        print(f"Using {args.init} from obsm as initialization.")
        init = adata.obsm[args.init][:, 0:2]
    else:
        init = None

    tsne_skrodzki(
        data=data,
        init=init,
        sampling_frac=args.sampling_frac,
        exaggeration=args.exaggeration,
        early_exag_iter=args.early_exag_iter,
        save_fpath=args.output,
        hd_metric=args.hd_metric,
    )