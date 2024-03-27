import numpy as np
import scanpy as sc
import sys
import time
import datetime
import pandas as pd
import argparse
from numpy.random import default_rng

sys.path.insert(1, "../../backend/")
from tsne import compute_tsne_series
from dataset import Dataset as TraceData


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "output",
        type=str,
        help="Output file for the benchmarking results.",
    )
    parser.add_argument(
        "--tsne",
        action="store_true",
        default=False,
        help="Compute tSNE embeddings for the benchmarking instead of using random two-dimensional embeddings. Default: False.",
    )
    parser.add_argument(
        "--d",
        type=int,
        nargs="+",
        default=[20],
        help="Dataset dimensionalities to benchmark. Default: 20",
    )
    parser.add_argument(
        "--n",
        type=float,
        nargs="+",
        default=[1e1],
        help="Dataset sizes to benchmark. Default: [1e2, 1e3, 1e4, 1e5, 1e6]",
    )

    return parser.parse_args()


# this function is copied from Böhm et al. (2022)
# https://github.com/berenslab/ne-spectrum/blob/56e7204710258d541fb716033d3542a4fca2705e/jnb_msc/generator/gauss_line.py#L72
def gauss_clusters(
    n_clusters=10,
    dim=20,
    pts_cluster=100,
    random_state=None,
    cov=1,
    stepsize=6,
):
    if random_state is None:
        rng = default_rng()
    else:
        rng = random_state

    s = stepsize / np.sqrt(dim)
    means = np.linspace(np.zeros(dim), n_clusters * s, num=n_clusters, endpoint=False)
    cov = np.eye(dim) * cov

    clusters = np.array(
        [rng.multivariate_normal(m, cov, size=(pts_cluster)) for m in means]
    )

    X = np.reshape(clusters, (-1, dim))
    y = np.repeat(np.arange(n_clusters), pts_cluster)
    return X, y

def random_data(n: int, d: int):
    rng = np.random.default_rng()
    data = rng.standard_normal(size=(n, d))
    return data, np.zeros(n, dtype=int)

def second_to_str(sec):
    str_tm = str(datetime.timedelta(seconds=sec))
    if "day" in str_tm:
        day = str_tm.split(",")[0]
        str_tm = str_tm.split(",")[1]
    else:
        day = "0 day "
    hour, minute, second = str_tm.split(":")
    return f"{day}{hour} hour {minute} min {second} sec"


def embedding_with_quality(size, d, compute_tsne=True):       
    timings = {"n": int(size),
                "d": int(d)}
    size = int(size)
    print(f"Dataset size: n={size}, d={d}")
    if size > 1000 and d > 1000:
        data, labels = random_data(size, d)
    else:
        data, labels = gauss_clusters(pts_cluster=int(size / 10), dim=d)

    start_pca = time.time()
    pca_emb = sc.pp.pca(data, n_comps=2, zero_center=True)
    timings["pca_s"] = time.time() - start_pca

    start_tsne = time.time()
    if compute_tsne:
        print("Computing tSNE embeddings...")
        sampling_frac = 1 if size < 100000 else 1 / (10 * (np.log10(size) - 3))
        tsne_embs = compute_tsne_series(
            data=data,
            coarse_exag_iter=[(12, 200)],
            fine_exag_iter=[(10, 200), (5, 200), (3, 200), (1, 200)],
            hd_metric="euclidean",
            init=pca_emb,
            sampling_frac=sampling_frac,
            smoothing_perplexity=min(30, int(size/3)),
            random_state=42,
            n_jobs=12,
        )
    else:
        tsne_embs = {}
        for t in range(4):
            tsne_embs[t] = default_rng().multivariate_normal(
                size=size, mean=np.zeros(2), cov=np.eye(2)
            )
    timings["tsne_s"] = time.time() - start_tsne

    print("Constructing TraceData...")
    start_trace = time.time()
    trace_data = TraceData(
        hd_data=data,
        name=f"Gaussian Line {size}",
        verbose=True,
        hd_metric="euclidean",
        n_jobs=12,
    )
    trace_data.add_metadata({"labels": labels.astype(int)})
    trace_data.add_embedding(
        name="PCA",
        embedding=pca_emb,
        category="PCA",
    )

    for i, emb in tsne_embs.items():
        if compute_tsne:
            trace_data.add_embedding(
                name=f"tSNE_{i}", embedding=emb, category="tSNE"
            )
        else:
            trace_data.add_embedding(
                name=f"random_{i}", embedding=emb, category="random"
            )
    timings["trace_s"] = time.time() - start_trace
    
    start_quality = time.time()
    trace_data.compute_quality()
    timings["quality_s"] = time.time() - start_quality
    timings["quality_readable"] = second_to_str(
        timings["quality_s"]
    )
    print("")
    return timings
    
    
if __name__ == "__main__":
    print("Running Gaussian Line Benchmarking")

    args = parse_args()
    timings = pd.DataFrame(
        columns=["n", "d", "pca_s", "tsne_s", "trace_s", "quality_s", "quality_readable"]
    )
    
    # numba warmup
    _ = embedding_with_quality(1e1, 20, compute_tsne=True)

    for size in args.n:
        for dim in args.d:
            result_dict = embedding_with_quality(size, d=dim, compute_tsne=args.tsne)
            timings = pd.concat([timings if not timings.empty else None, pd.DataFrame([result_dict])], ignore_index=True)
            
            timings.to_csv(
                args.output,
                float_format="%.2f",
                index=False,
                )
