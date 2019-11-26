from math import pi, cos, sqrt, radians
from typing import Tuple, Union, List

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
    vector_b = [-b * sin_alpha * np.cos(gamma_star),
                b * sin_alpha * np.sin(gamma_star),
                b * cos_alpha, ]
    vector_c = [0.0, 0.0, float(c)]
    return [vector_a, vector_b, vector_c]


def lengths(matrix):
    if not isinstance(matrix, np.ndarray):
        matrix = np.ndarray(matrix, dtype=float)
    return np.sqrt(np.sum(matrix ** 2, axis=1)).tolist()


def angles(matrix) -> Tuple[float]:
    """
    Returns the angles (alpha, beta, gamma) of the lattice.
    """
    lengt = lengths(matrix)
    angles = np.zeros(3)
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        angles[i] = abs_cap(np.dot(matrix[j], matrix[k]) / (lengt[j] * lengt[k]))
    angles = np.arccos(angles) * 180.0 / pi
    return angles.tolist()


def lattice_to_cell(lattice):
    return lengths(lattice) + angles(lattice)


def cell_to_g6(uc):
    """ Take a reduced Niggli Cell, and turn it into the G6 representation """
    a = uc[0] ** 2
    b = uc[1] ** 2
    c = uc[2] ** 2
    d = 2 * uc[1] * uc[2] * cos(radians(uc[3]))
    e = 2 * uc[0] * uc[2] * cos(radians(uc[4]))
    f = 2 * uc[0] * uc[1] * cos(radians(uc[5]))
    return [a, b, c, d, e, f]


def ncdist_fromcell(cell1, cell2):
    """
    Does a niggli reduction, G6 vector and distance calculation from two given unit cells.
    """
    reduced1 = lattice_to_cell(spglib.niggli_reduce(lattice_from_parameters(*cell1)))
    G6a = cell_to_g6(reduced1)
    reduced2 = lattice_to_cell(spglib.niggli_reduce(lattice_from_parameters(*cell2)))
    G6b = cell_to_g6(reduced2)
    return 0.01 * sqrt(ncdist.ncdist(G6a, G6b))


def IsRhombohedralAsHex(v):
    HexPerp = np.identity(6) - np.array(1. / 3.) * np.array([[1, 1, 0, 0, 0, -1],
                                                             [1, 1, 0, 0, 0, -1],
                                                             [0, 0, 3, 0, 0, 0],
                                                             [0, 0, 0, 0, 0, 0],
                                                             [0, 0, 0, 0, 0, 0],
                                                             [-1, -1, 0, 0, 0, 1]])
    RhmPerp = np.identity(6) - np.array(1. / 3.) * np.array([[1, 1, 1, 0, 0, 0],
                                                             [1, 1, 1, 0, 0, 0],
                                                             [1, 1, 1, 0, 0, 0],
                                                             [0, 0, 0, 1, 1, 1],
                                                             [0, 0, 0, 1, 1, 1],
                                                             [0, 0, 0, 1, 1, 1]])
    return np.linalg.norm(HexPerp * v) < np.linalg.norm(RhmPerp * v)


def make_primitive(cell: Union[List, Tuple], latt_type: str):
    """
    Port of C++ code from https://github.com/yayahjb/ncdist
    """
    latt_type = latt_type.upper()
    if latt_type == 'P':
        return np.identity(6)
    elif latt_type == 'I':
        # for monoclinic, assumes b unique:
        return np.array([[1, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0],
                         [0.25, 0.25, 0.25, 0.25, 0.25, 0.25],
                         [0, 1, 0, 0.5, 0, 0.5],
                         [1, 0, 0, 0, .5, .5],
                         [0, 0, 0, 0, 0, 1]])
    elif latt_type == 'A':
        # for monoclinic, assumes b unique:
        return np.array([[1, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0],
                         [0, .25, .25, .25, 0, 0],
                         [0, 1, 0, .5, 0, 0, ],
                         [0, 0, 0, 0, .5, .5],
                         [0, 0, 0, 0, 0, 1]])
    elif latt_type == 'B':
        # for monoclinic, assumes c unique:
        return np.array([[1, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0],
                         [.25, 0, .25, 0, .25, 0],
                         [0, 0, 0, .5, 0, .5],
                         [1, 0, 0, 0, .5, 0],
                         [0, 0, 0, 0, 0, 1]])
    elif latt_type == 'C':
        # for monoclinic, assumes b unique:
        return np.array([[1, 0, 0, 0, 0, 0],
                         [.25, .25, 0, 0, 0, .25],
                         [0, 0, 1, 0, 0, 0],
                         [0, 0, 0, .5, .5, 0],
                         [0, 0, 0, 0, 1, 0],
                         [1, 0, 0, 0, 0, .5]])
    elif latt_type == 'F':
        return np.array([[.25, .25, 0, 0, 0, .25],
                         [.25, 0, .25, 0, .25, 0],
                         [0, .25, .25, .25, 0, 0],
                         [0, 0, .5, .25, .25, .25],
                         [0, .5, 0, .25, .25, .25],
                         [.5, 0, 0, .25, .25, .25]])
    elif latt_type == 'R' or latt_type == 'H' and IsRhombohedralAsHex(cell):
        return (1.0 / 9.0) * np.array([[1, 1, 1, 1, -1, -1],
                                       [4, 1, 1, 1, 2, 2],
                                       [1, 4, 1, -2, -1, 2],
                                       [-4, -4, 2, -1, 1, -5],
                                       [2, -4, 2, -1, -2, 1],
                                       [-4, 2, 2, 2, 1, 1]])


Pc = np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 1.]])
Pr = np.array([[2. / 3., -1./3.0, -1/3.], [1/3, 1/3, -2/3], [1/3, 1/3, 1/3]])


if __name__ == '__main__':
    from spglib import spglib
    from cellcomp import ncdist

c1 = [5.2601, 9.1644, 10.6090, 104.851, 104.324, 100.457]
c2 = [5.2684, 9.2080, 10.6641, 69.559, 76.132, 79.767]
c3 = [float(x) for x in "100 100 100 90 90 90".split()]
c4 = [float(x) for x in "99 99 99 89 89 89".split()]
c3a = [float(x) for x in "7.5675 13.1966 11.3486   90.000  103.608   90.000 ".split()]
c4a = [float(x) for x in "7.6870 13.2020 11.5790   90.000  105.840   90.000 ".split()]

c5 = [57.98, 57.98, 57.98, 92.02, 92.02, 92.02]  # R32  1FE5
c6 = [80.36, 80.36, 99.44, 90, 90, 120]  # R3   1U4J, primitive= 57.02, 57.02, 57.02, 89.605, 89.605, 89.605
c7 = [80.949, 80.572, 57.098, 90.0, 90.35, 90.0]  # C2  1G2X
#
c8 = [78.961, 82.328, 57.031, 90.00, 93.44, 90.00]
c9 = [80.36, 80.36, 99.44, 90, 90, 120]
#
c10 = [3.1457, 3.1457, 3.1451, 60.089, 60.0887, 60.104]
c11 = [3.1456, 3.1458, 3.1451, 90.089, 119.907, 119.898]
#
# primitice cell:
c12 = [10., 10., 10., 112., 112.5, 112.9]  # G6: (100, 100, 100, -74.9, -76.5, -77.8)
# reduced cell:
c13 = [8.41, 10.0, 10.0, 112, 106.3, 106.8]  # reduced G6: (70.7, 100, 100,-74.9,-47-3,-48.5)

reduced = lattice_to_cell(spglib.niggli_reduce(lattice_from_parameters(*c12)))
print('g6 of 12:', cell_to_g6(c12))
print('reduced:', reduced)
G6a = cell_to_g6(reduced)
print('reduced g6 of 12:', G6a)

print(ncdist_fromcell(c1, c2))
print(ncdist_fromcell(c3a, c4a))
print(ncdist_fromcell(c3, c4))
print('c5:')
print(ncdist_fromcell(c6, c5))
print(ncdist_fromcell(c5, c7))
print(lattice_to_cell(spglib.niggli_reduce((lattice_from_parameters(*c7)*Pr))))
