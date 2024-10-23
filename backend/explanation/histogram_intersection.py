import numpy as np
from numba import jit, prange
from numba.typed import List


@jit(nopython=True, parallel=True, cache=True)
def get_histogram_intersection(A: np.ndarray, B: np.ndarray):
    # A and B must have the same dimensionality

    d = A.shape[1]
    hist_intersections = np.empty((d,))
    eps = 1e-15

    for f in prange(d):
        a = A[:, f]
        b = B[:, f]

        min_val = min(np.min(a), np.min(b))
        max_val = max(np.max(a), np.max(b))
        range_val = max(max_val - min_val, eps)
        
        # if there is no variation in the data for feature f, skip
        if (max_val-min_val) == 0:
            hist_intersections[f] = 1
            continue

        tmpA = (a - min_val) / range_val
        tmpB = (b - min_val) / range_val
        tmpAB = np.concatenate((tmpA, tmpB))

        binW = max(3.5 * np.std(tmpAB) / np.cbrt(tmpAB.shape[0]), eps)
        nBins = int((1.0 / binW) + 1)

        counts = np.zeros((2, nBins))

        for val in tmpA:
            binIndex = int(val / binW)
            counts[0, binIndex] += 1

        for val in tmpB:
            binIndex = int(val / binW)
            counts[1, binIndex] += 1
            
        counts[0] = counts[0] / A.shape[0]
        counts[1] = counts[1] / B.shape[0]

        min_counts = np.zeros(nBins)
        for i in range(nBins):
            min_counts[i] = min(counts[0, i], counts[1, i])
            
        hist_intersections[f] = min_counts.sum()
    return hist_intersections


@jit(nopython=True, parallel=True, cache=True)
def wasserstein_histogram_distance(A: np.ndarray, B: np.ndarray):
    d = A.shape[1]
    distances = np.empty((d,))
    eps = 1e-15

    for f in prange(d):
        a = A[:, f]
        b = B[:, f]

        min_val = min(np.min(a), np.min(b))
        max_val = max(np.max(a), np.max(b))
        range_val = max(max_val - min_val, eps)
        
        # if there is no variation in the data for feature f, skip
        if max_val-min_val == 0:
            distances[f] = 0
            continue

        tmpA = (a - min_val) / range_val
        tmpB = (b - min_val) / range_val
        tmpAB = np.concatenate((tmpA, tmpB))

        # if data only contains 0 and 1 values
        if np.unique(tmpAB).shape[0] == 2:
            binW = 0.5+eps
            nBins = 2
        else:
            binW = max(3.5 * np.std(tmpAB) / np.cbrt(tmpAB.shape[0]), eps)
            nBins = int((1.0 / binW) + 1)

        counts = np.zeros((2, nBins))

        for val in tmpA:
            binIndex = int(val / binW)
            counts[0, binIndex] += 1

        for val in tmpB:
            binIndex = int(val / binW)
            counts[1, binIndex] += 1

        counts[0] = counts[0] / A.shape[0]
        counts[1] = counts[1] / B.shape[0]
        deltas = np.full(nBins, 1/nBins)
        A_cdf = np.cumsum(counts[0])
        B_cdf = np.cumsum(counts[1])
        distances[f] = np.sum(np.multiply(np.abs(A_cdf - B_cdf), deltas))
    return distances


# precompute the binwidts for the complete dataset
@jit(nopython=True, parallel=True, cache=True)
def get_binwidths(X: np.ndarray):
    d = X.shape[1]
    eps = 1e-15
    binwidths = np.empty((d,))
    counts = List()

    for f in prange(d):
        a = X[:, f]

        min_val = np.min(a)
        max_val = np.max(a)
        range_val = max(max_val - min_val, eps)
        
        # if there is no variation in the data for feature f, skip
        if (max_val-min_val) == 0:
            binwidths[f] = 1
            counts.append([1.0])
            continue

        # normalize between 0 and 1
        tmpA = (a - min_val) / range_val

        # there is no variance in the data
        if (max_val - min_val) == 0:
            binW = 1
        else:
            binW = max(3.5 * np.std(tmpA) / np.cbrt(tmpA.shape[0]), eps)

        nBins = int((1.0 / binW) + 1)
        f_count = np.zeros(nBins)

        for val in tmpA:
            binIndex = int(val / binW)
            f_count[binIndex] += 1

        # relative frequency
        f_count = f_count / X.shape[0]
        binwidths[f] = binW
        counts.append(f_count)

    return binwidths, counts


@jit(nopython=True, parallel=True, cache=True)
def get_partial_histogram_intersection(
    X: np.ndarray, indices: np.ndarray, binwidths: np.ndarray, counts: list
):
    d = X.shape[1]
    hist_intersections = np.empty((d,))
    eps = 1e-15

    for f in prange(d):
        f_values = X[:, f]
        # normalize data
        min_val = np.min(f_values)
        max_val = np.max(f_values)
        range_val = max(max_val - min_val, eps)
        
        # if there is no variation in the data for feature f, skip
        if (max_val-min_val) == 0:
            hist_intersections[f] = 1
            continue

        f_cluster = (f_values[indices] - min_val) / range_val

        binW = binwidths[f]
        counts_all = counts[f]
        nBins = counts_all.shape[0]
        counts_cluster = np.zeros(nBins)

        for val in f_cluster:
            binIndex = int(val / binW)
            counts_cluster[binIndex] += 1

        # relative frequency
        counts_cluster = counts_cluster / f_cluster.shape[0]

        min_counts = np.zeros(nBins)
        for i in prange(nBins):
            min_counts[i] = min(counts_all[i], counts_cluster[i])

        hist_intersections[f] = 1 + min_counts.sum()

    return hist_intersections