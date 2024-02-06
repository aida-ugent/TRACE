import numpy as np
import pandas as pd
import sklearn
from scipy import stats
from sklearn.metrics import pairwise_distances
from sklearn.cluster import kmeans_plusplus


def get_cluster_centroids(data: np.ndarray, labels: np.ndarray):
    """
    Calculate the centroids of each cluster.

    Args:
        data (np.ndarray): The data points.
        labels (np.ndarray): The cluster labels for each data point.

    Returns:
        np.ndarray: The centroids of each cluster.
    """
    centroids = []
    for label in labels.unique():
        ind = list(np.nonzero(labels == label)[0])
        centroids.append(np.mean(data[ind, :], axis=0))
    centroids = np.row_stack(centroids)
    return centroids


def compute_centroid_distances(data: np.ndarray, labels: np.ndarray, metric: str):
    """
    Compute the pairwise distances between the cluster centroids.

    Parameters:
        data (np.ndarray): The data space used to compute the centroids.
        labels (np.ndarray): Labels of the data points.
        metric (str, optional): The distance metric to be used.

    Returns:
        ndarray: The pairwise distances between the cluster centroids.
    """
    if metric == "angular":
        metric = "cosine"
    assert metric in sklearn.metrics.pairwise.PAIRWISE_DISTANCE_FUNCTIONS
    centroids = get_cluster_centroids(data=data, labels=labels)
    centroid_distances = pairwise_distances(centroids, metric=metric)
    return centroid_distances


def distance_correlation(hd_centroid_distances, ld_centroid_distances):
    """
    Compute the centroid correlations between high-dimensional and low-dimensional centroid distances.

    Parameters:
    - hd_centroid_distances (np.ndarray): Array of high-dimensional centroid distances.
    - ld_centroid_distances (np.ndarray): Array of low-dimensional centroid distances.

    Returns:
    - corr (np.ndarray): Array of centroid correlations.
    """
    corr = np.empty((hd_centroid_distances.shape[0],))
    for i in range(hd_centroid_distances.shape[0]):
        corr[i] = stats.spearmanr(
            a=hd_centroid_distances[i, :],
            b=ld_centroid_distances[i, :],
            axis=0,
        ).correlation
    return corr


def get_hd_centroid_distance_df(
    labels: pd.Series, label_key: str, hd_centroid_distances: np.ndarray
):
    """
    Store the distances between the centroids of each cluster.

    Args:
        labels (pd.Series): Cluster labels.
        cluster_unique (np.ndarray): Unique cluster labels.
        hd_centroid_distances (np.ndarray): High-dimensional centroid distances.
    """
    labels_unique = labels.unique()
    distance_dict = {label_key: labels_unique}
    distance_dict.update(
        dict(
            (f"{labels_unique[i]}_HD_centroid_dist", hd_centroid_distances[i])
            for i in np.argsort(labels_unique)
        )
    )
    return pd.DataFrame(distance_dict)


def sample_cluster_points(data: np.ndarray, labels: np.ndarray, max_samples: int):
    """Sampling points from the HD data space using the Kmeans++ algorithm.

    Args:
        data (np.ndarray): HD data
        labels (np.ndarray): sampling at least one point from each cluster
        max_samples (int): maximum number of samples.

    Returns:
        np.ndarray: sampled points
    """
    centroids = []
    max_samples = min(max_samples, data.shape[0] * 0.1)

    for label in labels.unique():
        ind = list(np.nonzero(labels == label)[0])
        # how many samples to take from this cluster
        n_samples = max(1, int(max_samples * (len(ind) / data.shape[0])))
        # Calculate seeds from k-means++
        _, indices = kmeans_plusplus(data[ind, :], n_clusters=n_samples, random_state=0)
        orig_indices = np.array(ind)[indices]
        centroids.extend(orig_indices)
    return np.asarray(centroids)


def sample_landmarks(data: np.ndarray, max_samples: int):
    """Sampling points from the HD data space using the Kmeans++ algorithm.

    Args:
        data (np.ndarray): HD data
        max_samples (int): maximum number of samples.

    Returns:
        np.ndarray: sampled points
    """
    max_samples = int(min(max_samples, data.shape[0] * 0.1))
    _, indices = kmeans_plusplus(data, n_clusters=max_samples, random_state=0)
    return np.asarray(indices)


def compute_pairwise_distance(X: np.ndarray, Y: np.ndarray, metric: str):
    if metric == "angular":
        metric = "cosine"
    assert metric in sklearn.metrics.pairwise.PAIRWISE_DISTANCE_FUNCTIONS
    return pairwise_distances(X, Y, metric=metric)
