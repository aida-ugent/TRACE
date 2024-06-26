import numpy as np
from numpy.random import default_rng
import pandas as pd
from scipy import stats
from sklearn.metrics import pairwise_distances
from sklearn.cluster import kmeans_plusplus
from sklearn.neighbors import NearestNeighbors


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


def sample_cluster_points(data: np.ndarray, labels: np.ndarray, num_samples: int):
    """Sampling points from the HD data space using the Kmeans++ algorithm.

    Args:
        data (np.ndarray): HD data
        labels (np.ndarray): sampling at least one point from each cluster
        num_samples (int): maximum number of samples.

    Returns:
        np.ndarray: sampled points
    """
    centroids = []
    num_samples = min(num_samples, data.shape[0] * 0.1)

    for label in labels.unique():
        ind = list(np.nonzero(labels == label)[0])
        # how many samples to take from this cluster
        n_samples = max(1, int(num_samples * (len(ind) / data.shape[0])))
        # Calculate seeds from k-means++
        _, indices = kmeans_plusplus(data[ind, :], n_clusters=n_samples, random_state=0)
        orig_indices = np.array(ind)[indices]
        centroids.extend(orig_indices)
    return np.asarray(centroids)


def sample_landmarks(data: np.ndarray, num_samples: int):
    """Sampling points from the HD data space using the Kmeans++ algorithm.

    Args:
        data (np.ndarray): HD data
        num_samples (int): maximum number of samples.

    Returns:
        np.ndarray: sampled points
    """
    num_samples = min(num_samples, data.shape[0])
    _, indices = kmeans_plusplus(data, n_clusters=num_samples, random_state=0)
    
    indices = np.unique(np.asarray(indices))
    if len(indices) < num_samples:
        print(f"Warning: Not enough unique landmarks from kmeans++. Adding {num_samples - len(indices)} random samples.")
        possible_indices = np.setdiff1d(np.arange(data.shape[0]), indices)
        indices = np.concatenate([indices, default_rng().choice(possible_indices, num_samples - len(indices), replace=False)])
    return np.asarray(indices)


def compute_pairwise_distance(X: np.ndarray, Y: np.ndarray, metric: str):
    if metric == "angular":
        metric = "cosine"
    return pairwise_distances(X, Y, metric=metric)


def compute_landmark_correlation(
    ld_data,
    hd_data,
    landmark_indices,
    hd_landmark_distances,
    LD_landmark_neighbors,
    hd_metric,
    ld_metric,
):
    ld_landmark_distances = compute_pairwise_distance(
        X=ld_data[landmark_indices], Y=None, metric=ld_metric
    )

    corr = distance_correlation(hd_landmark_distances, ld_landmark_distances)

    if LD_landmark_neighbors:
        # find the nearest landmark neighbor in LD space for all points
        # alternative: find the nearest landmark neighbor in HD space for all points.
        landmark_nbr_index = NearestNeighbors(
            n_neighbors=2, algorithm="auto", metric="euclidean"
        ).fit(ld_data[landmark_indices])
        neighbors = landmark_nbr_index.kneighbors(
            ld_data, n_neighbors=1, return_distance=False
        )
    else:
        # nearest neighbor index for HD Data
        landmark_nbr_index = NearestNeighbors(
            n_neighbors=2,
            algorithm="auto",
            metric="cosine" if hd_metric == "angular" else hd_metric,
        ).fit(hd_data[landmark_indices, :])
        neighbors = landmark_nbr_index.kneighbors(
            hd_data, n_neighbors=1, return_distance=False
        )

    neighbors = neighbors.flatten()

    # place the remaining points at the location of their nearest neighbor
    return np.take(corr, neighbors, axis=0)
