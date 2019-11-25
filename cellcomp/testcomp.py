from math import pi
from typing import Tuple

import numpy as np


def abs_cap(val, max_abs_val=1):
    """
    Returns the value with its absolute value capped at max_abs_val.
    Particularly useful in passing values to trignometric functions where
    numerical errors may result in an argument > 1 being passed in.
    """
    return max(min(val, max_abs_val), -max_abs_val)


def lattice_from_parameters(a: float,
                            b: float,
                            c: float,
                            alpha: float,
                            beta: float,
                            gamma: float):
    """
    Create a Lattice using unit cell lengths and angles (in degrees).
    """
    angles_r = np.radians([alpha, beta, gamma])
    cos_alpha, cos_beta, cos_gamma = np.cos(angles_r)
    sin_alpha, sin_beta, sin_gamma = np.sin(angles_r)
    val = (cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta)
    # Sometimes rounding errors result in values slightly > 1.
    val = abs_cap(val)
    gamma_star = np.arccos(val)
    vector_a = [a * sin_beta, 0.0, a * cos_beta]
    vector_b = [
        -b * sin_alpha * np.cos(gamma_star),
        b * sin_alpha * np.sin(gamma_star),
        b * cos_alpha,
    ]
    vector_c = [0.0, 0.0, float(c)]
    return [vector_a, vector_b, vector_c]


def lengths(matrix):
    return np.sqrt(np.sum(matrix ** 2, axis=1)).tolist()


def angles(m) -> Tuple[float]:
    """
    Returns the angles (alpha, beta, gamma) of the lattice.
    """
    lengt = lengths(m)
    angles = np.zeros(3)
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        angles[i] = abs_cap(np.dot(m[j], m[k]) / (lengt[j] * lengt[k]))
    angles = np.arccos(angles) * 180.0 / pi
    return angles.tolist()


if __name__ == '__main__':
    from spglib import spglib
    from cellcomp import ncdist
    c1 = [5.2601, 9.1644, 10.6090, 104.851, 104.324, 100.457]
    c2 = [5.2684, 9.2080, 10.6641, 69.559, 76.132, 79.767]
    c3 = [7.745 ,  11.728  , 21.321  , 78.03  , 90.00 ,  70.72]
    matrix = spglib.niggli_reduce(lattice_from_parameters(*c1))
    matrix2 = spglib.niggli_reduce(lattice_from_parameters(*c2))
    #prim = spglib.standardize_cell((lattice_from_parameters(*c3), [[0, 0, 0]], [0]), to_primitive=True)
    #print('primitive:', lengths(prim[0]), angles(prim[0]))
    #print(matrix)
    print('niggli reduced:', lengths(matrix), angles(matrix))
    ncd = ncdist.ncdist(lengths(matrix)+angles(matrix), lengths(matrix2)+angles(matrix2))
    print(ncd)