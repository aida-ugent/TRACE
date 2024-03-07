import numpy as np
from tqdm import tqdm

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
        embeddings (_type_): _description_
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

    variances = variance_across_embeddings(embeddings, num_samples)
    stability = stability_from_variance(variances, alpha)
    return stability


def variance_across_embeddings(embeddings, num_samples=50):
    """
    Computes the variance of points in different embeddings using a random sample of points.

    Args:
        embeddings (list of ndarray): List of embeddings as basis for point variance.
        samples (int): Number of samples to use for the variance calculation. Defaults to 50.
    """

    num_embeddings = len(embeddings)
    n = embeddings[0].shape[0]

    variance = np.zeros((n))
    # the distances are normalized by the average distances over all points and embeddings
    sum_of_distances = 0

    # Iterate through all points
    for i in tqdm(range(n)):
        # sample random points
        sample = np.random.choice(n, num_samples, replace=False)

        # compute distances of i to all points in the sample for all embeddings
        distances = np.zeros((len(embeddings), num_samples))
        for e in range(num_embeddings):
            distances[e] = np.array(
                [np.linalg.norm(embeddings[e][i] - embeddings[e][j]) for j in sample]
            )

        sum_of_distances += np.sum(distances)

        # normalize the distances over embeddings
        distances -= np.mean(distances, axis=0)

        # compute the variance of the distances
        variance[i] = np.sum(np.var(distances, axis=0, ddof=1))

    # Normalize the variance scores
    variance /= num_samples
    variance /= sum_of_distances / (n * num_samples * num_embeddings)
    return variance


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