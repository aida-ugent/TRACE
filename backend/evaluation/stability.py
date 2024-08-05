import numpy as np
from numba import jit, prange
from numpy.random import default_rng
import warnings

# This stability measure is based on the paper "Dynamic visualization of
# high-dimensional data" from Eric Sun et al. (2023)
# The different is that we are using embeddings of different algorithms or with
# different hyperparameters as basis for the stability calculation while
# they use bootstrapping.
# https://www.nature.com/articles/s43588-022-00380-4


def stability_across_embeddings(embeddings, num_samples=50, alpha=20):
    """
    This follows the definition of point stability from "Dynamic visualization
    of high-dimensional data" from Eric Sun et al. Instead of bootstrapping the data,
    we use different embeddings (from different algorithms or computed with different
    parameters) as basis for the stability calculation.

    Args:
        embeddings (ndarray (e, n, 2)): Two-dimensional embeddings for all data-points as basis for point variance.
        num_samples (int, optional): _description_. Defaults to 50.
        alpha (float, optional): Exponent for the stability calculation.
            Values between 5 and 50 are typically sufficient. Defaults to 20.
    """
    if len(embeddings) < 5:
        # warning instead of value error
        warnings.warn(
            "We recommend at least five embeddings for stability calculation."
        )
    
    if alpha < 0:
        raise ValueError("The alpha parameter must be greater than 0.")

    variances = variance_across_embeddings(np.asarray(embeddings), num_samples)
    stability = stability_from_variance(variances, alpha)
    return stability


def variance_across_embeddings(embeddings, num_samples=50):
    num_embeddings = len(embeddings)
    n = embeddings[0].shape[0]
    embeddings = np.ascontiguousarray(embeddings)
    variance = np.zeros(n)
    # the distances are normalized by the average distances over all points and embeddings
    distances = np.zeros(n)
    rng = default_rng()
    samples = rng.choice(n, (n, num_samples), replace=True)

    compute_variance(embeddings, variance, distances, samples)

    # Normalize the variance scores
    variance /= num_samples
    variance /= np.sum(distances) / (n * num_samples * num_embeddings)
    return variance


@jit(nopython=True, parallel=True, cache=True)
def compute_variance(embeddings, variance, sum_distances, samples):
    n = embeddings[0].shape[0]
    num_embeddings = embeddings.shape[0]

    for i in prange(n):
        # sample = rng.choice(n, num_samples, replace=False)
        distances = np.sqrt(
            np.sum(
                (
                    embeddings[:, samples[i], :]
                    - embeddings[:, np.repeat(i, samples.shape[1]), :]
                )
                ** 2,
                axis=2,
            )
        )

        sum_distances[i] += np.sum(distances)
        distances -= np.sum(distances, axis=0) / samples.shape[1]

        # compute the variance of distances along the samples
        embedding_means = np.sum(distances, axis=0) / num_embeddings
        variance[i] = np.sum(
            np.sum((distances - embedding_means) ** 2, axis=0) / (num_embeddings - 1)
        )


def stability_from_variance(variance, alpha):
    """
    Transforms variance to stability using the alpha parameter.

    Args:
        variance (ndarray): Variance of the distances of points in different embeddings.
        alpha (float, optional): Exponent for the stability calculation.
            Values between 5 and 50 are typically sufficient.
    """
    stability = 1 / ((1 + variance) ** alpha)
    return stability
