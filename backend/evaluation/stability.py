import numpy as np
from tqdm import tqdm
from numba import jit, prange
from numba_progress import ProgressBar

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
        raise ValueError(
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
    variance = np.zeros((n))
    # the distances are normalized by the average distances over all points and embeddings
    distances = np.empty((n,))

    with ProgressBar(total=n, dynamic_ncols=True) as numba_progress:
        variance_progress(embeddings, variance, distances, num_samples, numba_progress) 

    # Normalize the variance scores
    variance /= num_samples
    variance /= np.sum(distances) / (n * num_samples * num_embeddings)
    return variance

@jit(nopython=True, parallel=True, cache=True)
def variance_progress(embeddings, variance, sum_distances, num_samples, progress_hook):
    n = embeddings[0].shape[0]
    
    for i in prange(n):
        sample = np.random.choice(n, num_samples, replace=False)
        distances = np.sqrt(np.sum((embeddings[:, sample, :] - embeddings[:, np.repeat(i, num_samples), :]) ** 2, axis=2))
        
        sum_distances[i] += np.sum(distances)
        distances -= np.sum(distances, axis=0) / num_samples

        # compute the variance of distances along the samples
        embedding_means = np.sum(distances, axis=0) / (len(embeddings) - 1)
        variance[i] = np.sum(np.sum((distances - embedding_means) ** 2, axis=0) / (len(embeddings) - 1))

        if progress_hook is not None:
            progress_hook.update(1)

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
