from math import pi, cos, sqrt
from typing import Tuple

import numpy as np


def abs_cap(val, max_abs_val=1):
    """
    Returns the value with its absolute value capped at max_abs_val.
    Particularly useful in passing values to trignometric functions where
    numerical errors may result in an argument > 1 being passed in.
    """
    return max(min(val, max_abs_val), -max_abs_val)


def lattice_from_parameters(a: float, b: float, c: float, alpha: float, beta: float, gamma: float):
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


def cell_to_g6(uc):
    """ Take a reduced Niggli Cell, and turn it into the G6 representation """
    a = uc[0] ** 2
    b = uc[1] ** 2
    c = uc[2] ** 2
    d = 2. * uc[1] * uc[2] * cos(uc[3])
    e = 2. * uc[0] * uc[2] * cos(uc[4])
    f = 2. * uc[0] * uc[1] * cos(uc[5])
    return [a, b, c, d, e, f]


if __name__ == '__main__':
    from spglib import spglib
    from cellcomp import ncdist

    c1 = [5.2601, 9.1644, 10.6090, 104.851, 104.324, 100.457]
    c2 = [5.2684, 9.2080, 10.6641, 69.559, 76.132, 79.767]
    c3 = [float(x) for x in "100 100 100 90 90 90".split()]
    c4 = [float(x) for x in "99 99 99 89 89 89".split()]
    
    c5 = [57.98, 57.98, 57.98, 92.02, 92.02, 92.02]
    c6 = [80.36, 80.36, 99.44, 90, 90, 120]
    c7 = [80.949, 80.572, 57.098, 90.0, 90.35, 90.0]
    #
    c8 = [78.961, 82.328, 57.031, 90.00, 93.44, 90.00]
    c9 = [80.36, 80.36, 99.44, 90, 90, 120]
    matrix = spglib.niggli_reduce(lattice_from_parameters(*c3))
    matrix2 = spglib.niggli_reduce(lattice_from_parameters(*c4))
    #prim = spglib.standardize_cell((lattice_from_parameters(*c6), [[0, 0, 0]], [0]), to_primitive=True)
    #print('primitive:', lengths(prim[0]), angles(prim[0]))
    #print(matrix)
    print('niggli reduced:', lengths(matrix), angles(matrix))
    g61 = cell_to_g6(lengths(matrix)+angles(matrix))
    g62 = cell_to_g6(lengths(matrix2)+angles(matrix2))
    print('g6:', g61, g62)
    ncd = ncdist.ncdist(g61, g62)
    print('with g6:', 0.1*sqrt(ncd))
    print('with standard cell:', 0.1*sqrt(ncdist.ncdist(lengths(matrix)+angles(matrix), lengths(matrix2)+angles(matrix2))))
