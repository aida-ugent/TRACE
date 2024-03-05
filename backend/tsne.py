import openTSNE
import numpy as np
import os
import time
import math
import anndata as ad


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


def tsne_exaggeration(
    data: np.ndarray,
    exag_iter=[(10, 200), (5, 200), (3, 200), (1, 200)],
    fpath_prefix: str = None,
    hd_metric: str = "euclidean",
    init: np.ndarray = None,
    perplexity: int = 30,
    random_state: int = 42,
    return_affinities: bool = False,
    **kwargs,
):
    embeddings = {}
    start = time.time()
    knn_index = openTSNE.affinity.get_knn_index(
        data,
        "annoy",
        int(3 * perplexity),
        hd_metric,
        n_jobs=8,
        random_state=None,
        verbose=True,
    )

    print(f"Computing affinities with perplexity {perplexity}...")
    # computing embedding
    affinities = openTSNE.affinity.PerplexityBasedNN(
        perplexity=perplexity,
        method="annoy",
        n_jobs=8,
        random_state=random_state,
        metric=hd_metric,
        verbose=True,
        knn_index=knn_index,
    )

    # initialization
    if init is None:
        print("Computing PCA initialization...")
        init = openTSNE.initialization.pca(data)
    else:
        init = openTSNE.initialization.rescale(init)

    embedding = openTSNE.TSNEEmbedding(
        embedding=init,
        affinities=affinities,
        n_jobs=8,
        verbose=True,
        random_state=random_state,
        **kwargs,
    )

    # optimize embedding
    for exag, n_iter in exag_iter:
        embedding.optimize(n_iter=n_iter, exaggeration=exag, inplace=True)
        embeddings[exag] = np.asarray(embedding).copy()
        if fpath_prefix is not None:
            np.savetxt(
                f"{fpath_prefix}_tsne_exg_{exag}.csv", X=embedding, delimiter=","
            )
    embedding_time = time.time() - start
    print(
        f"Done. ({embedding_time:.2f}s)",
        flush=True,
        end="\n",
    )

    if return_affinities:
        return embeddings, affinities
    else:
        return embeddings


def compute_tsne_series(
    data: np.ndarray,
    coarse_exag_iter=[(12, 200)],
    fine_exag_iter=[(10, 200), (5, 200), (3, 200), (1, 200)],
    fpath_prefix: str = None,
    hd_metric: str = "euclidean",
    init: np.ndarray = None,
    sampling_frac: float = 0.01,
    smoothing_perplexity: int = 30,
    random_state: int = 42,
    **kwargs,
):
    """
    Compute a single t-SNE embedding using the sampling idea of Skrodzki et al.,
    'Tuning the perplexity for and computing sampling-based t-SNE embeddings' (2023).

    In their approach, a fraction k of the data is embedded using a perplexity
    perp = (N*k)/100. Then, all remaining points are embedded at the location
    of their nearest neighbor from the sample (prolongation). The final smoothing
    step refines the full embedding with a low perplexity for a couple of iterations.

    Args:
        X (np.ndarray): The input data matrix of shape (n_samples, n_features).
        coarse_exag_iter (list(Tuple(int, int))): A list of tuples containing the exaggeration and the
            number of iterations for the coarse embedding. Only used if sampling frac is in (0, 1).
        fine_exag_iter (list(Tuple(int, int))): A list of tuples containing the exaggeration and the
            number of iterations for the fine embedding.
        fpath_prefix (str, optional): The file path prefix for saving the embeddings.
        hd_metric (str, optional): The metric used for high-dimensional space distance calculation (default: "euclidean").
        init (np.ndarray, optional): The initial embedding to use for tSNE (default: None).
        smoothing_perplexity (float, optional): The perplexity for smoothing. Default is 30.
        kwargs: Additional keyword arguments for TSNEEmbedding from openTSNE.

    Returns:
        embeddings (Dict[str, np.ndarray]): A dictionary containing the tSNE embeddings,
            where the keys contain the exaggeration (from the fine_exag_iter list) and the values are the embeddings.
    """
    n = data.shape[0]
    fine_init = None

    # only do coarse embedding if sampling frac is in (0, 1)
    if sampling_frac > 0 and sampling_frac < 1:
        sampling_size = math.ceil(n * sampling_frac)
        sample_ind = np.random.choice(n, size=sampling_size, replace=False)
        coarse_perp = math.ceil((n * sampling_frac) / 100)

        coarse_embeddings, coarse_affinities = tsne_exaggeration(
            data=data[sample_ind, :],
            exag_iter=coarse_exag_iter,
            fpath_prefix=None,
            hd_metric=hd_metric,
            init=init[sample_ind, :],
            perplexity=coarse_perp,
            random_state=random_state,
            return_affinities=True,
            **kwargs,
        )

        print(
            "Prolongating embedding by placing each point on their nearest landmark..."
        )
        coarse_embedding = coarse_embeddings[coarse_exag_iter[-1][0]]
        fine_init = prolongate_embedding(
            data, coarse_embedding, sample_ind, coarse_affinities.knn_index
        )
        fine_init = openTSNE.initialization.rescale(fine_init)
    else:
        if init is None:
            print(f"Computing PCA initialization...")
            fine_init = openTSNE.initialization.pca(data)
        else:
            fine_init = openTSNE.initialization.rescale(init)

    fine_embeddings = tsne_exaggeration(
        data=data,
        exag_iter=fine_exag_iter,
        fpath_prefix=fpath_prefix,
        hd_metric=hd_metric,
        init=fine_init,
        perplexity=smoothing_perplexity,
        random_state=random_state,
        **kwargs,
    )
    return fine_embeddings
