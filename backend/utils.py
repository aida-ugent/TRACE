import numpy as np
import os
import yaml


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


def get_available_datasets(config_filepath):
    """Removes datasets from the config that do not have a file present.

    Args:
        config (dict): Dictionary containing the dataset configurations.
            The keys are the dataset names and the values are dictionaries

    Returns:
        dict: Cleaned config dictionary
    """
    with open(config_filepath, "r") as file:
        config = yaml.safe_load(file)

    datasets = list(config.keys())

    # check if the file is present
    for name in datasets:
        if not os.path.isfile(config[name]["filepath"]):
            config.pop(name)
    return config
