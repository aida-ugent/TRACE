import numpy as np
import scanpy as sc
import sys
import time
import datetime
import pandas as pd

sys.path.insert(1, "../../backend/")
from tsne import compute_tsne_series
from dataset import Dataset as TraceData


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
        rng = np.random.RandomState()
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


def second_to_str(sec):
    str_tm = str(datetime.timedelta(seconds=sec))
    if "day" in str_tm:
        day = str_tm.split(",")[0]
        str_tm = str_tm.split(",")[1]
    else:
        day = "0 day "
    hour, minute, second = str_tm.split(":")
    return f"{day}{hour} hour {minute} min {second} sec"


if __name__ == "__main__":
    print("Running Gaussian Line Benchmarking")

    data_sizes = [1e3, 1e4]
    timings = pd.DataFrame(
        columns=["pca_s", "tsne_s", "quality_s", "quality_readable"], index=data_sizes
    )

    for size in data_sizes:
        size = int(size)
        print(f"Dataset size: {size}")
        data, labels = gauss_clusters(pts_cluster=int(size / 10))

        start_pca = time.time()
        pca_emb = sc.pp.pca(data, n_comps=2, zero_center=True)
        timings.loc[size]["pca_s"] = time.time() - start_pca

        start_tsne = time.time()
        sampling_frac = 1 if size < 100000 else int(1 / (10 * (np.log10(size) - 3)))
        # tsne_embs = compute_tsne_series(
        #     data=data,
        #     coarse_exag_iter=[(12, 200), (5, 200)],
        #     fine_exag_iter=[(10, 200), (5, 200), (3, 200), (1, 200)],
        #     hd_metric="euclidean",
        #     init=pca_emb,
        #     sampling_frac=sampling_frac,
        #     smoothing_perplexity=30,
        #     random_state=42,
        # )
        tsne_embs = {}
        for t in range(4):
            tsne_embs[t] = np.random.multivariate_normal(
                size=size, mean=np.zeros(2), cov=np.eye(2)
            )
        timings.loc[size]["tsne_s"] = time.time() - start_tsne

        trace_data = TraceData(
            hd_data=data,
            name=f"Gaussian Line {size}",
            verbose=True,
            hd_metric="euclidean",
        )
        trace_data.add_metadata({"labels": labels.astype(int)})
        trace_data.add_embedding(
            name="PCA",
            embedding=pca_emb,
            category="PCA",
        )

        for exag, emb in tsne_embs.items():
            trace_data.add_embedding(
                name=f"tSNE_{exag}", embedding=emb, category="tSNE"
            )
        start_quality = time.time()
        trace_data.compute_quality()
        timings.loc[size]["quality_s"] = time.time() - start_quality
        timings.loc[size]["quality_readable"] = second_to_str(
            timings.loc[size]["quality_s"]
        )
        timings.to_csv(
            "gauss_line_benchmarking.csv",
            float_format="%.2f",
            index_label="dataset_size",
        )
