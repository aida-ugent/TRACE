import numpy as np
from explanation.histogram_intersection import (
    get_histogram_intersection,
    wasserstein_histogram_distance,
)

methods = ["histogram", "wasserstein"]


def get_feature_importance(A: np.ndarray, B: np.ndarray, method: str):
    assert method in methods, f"Method {method} not in {methods}"

    if method == "histogram":
        intersection = get_histogram_intersection(A, B)

        # smaller intersection means more important
        feature_ordering = np.argsort(intersection)
        return (feature_ordering, intersection[feature_ordering])

    elif method == "wasserstein":
        distance = wasserstein_histogram_distance(A, B)

        # larger distance means more important
        feature_ordering = np.argsort(-distance)
        return (feature_ordering, distance[feature_ordering])
