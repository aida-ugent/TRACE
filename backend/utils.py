import numpy as np


def normalizeEmbedding(arr):
    """
    Normalize the provided embedding to [-1,1].

    Parameters:
    arr (np.ndarray): The embedding to be normalized.

    Returns:
    np.ndarray: The normalized embedding array.
    """
    min = np.min(arr)
    diff = np.max(arr) - min
    arr = (2 * ((arr - min) / diff)) - 1
    return arr
