import numpy as np
from numpy import ndarray


# source: https://gist.github.com/nh2/bc4e2981b0e213fefd4aaa33edfb3893#file-rigid-transform-with-scale-py-L20-L39
def umeyama_alignment(P: ndarray, Q: ndarray) -> ndarray:
    """
    Kabsch/Umeyama alignment of P to Q.

    Args:
        P (np.ndarray): (N,D) embedding of points
        Q (np.ndarray): (N,D) embedding of points

    Returns:
        np.ndarray (N, D): Aligned matrix P.
    """
    n, d = P.shape

    # center
    P_centered = P - P.mean(axis=0)
    Q_centered = Q - Q.mean(axis=0)

    # cross-covariance of dimension (D, D)
    C: ndarray = np.dot(np.transpose(P_centered), Q_centered) / n

    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # create Rotation matrix U
    R = np.dot(V, W)

    # scale factor
    varP = np.var(P, axis=0).sum()
    c = 1 / varP * np.sum(S)

    # translation
    t = Q.mean(axis=0) - P.mean(axis=0).dot(c * R)
    aligned_P = P.dot(c * R) + t

    return aligned_P
