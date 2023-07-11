#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
Utilities for complex beam parameter propagation and string formatting

Note that the ABCD convention used here for beam propagation differs from other references,
e.g. Siegman. Here complex beam parameters q are propagated by

q2 = (A * q1 + B) / (C * q1 + D)

whereas in Siegman

q2/n2 = (AA * (q1/n1) + BB) / (CC * (q1/n1) + DD)

where n1 and n2 are the indices of refraction of the regions where the two beams are propagating.
To convert between the two,
A = AA, B = n1 * BB, C = CC / n2, and D = DD * (n1 / n2)
"""
import numpy as np
from wield.utilities.np import matrix_stack


def np_check_sorted(vals):
    if len(vals.shape) > 1:
        return False
    else:
        return np.all(vals[1:] > vals[:-1])


def matrix_space(L_m):
    return np.matrix(
        [
            [1, L_m],
            [0, 1],
        ]
    )


def matrix_focus(f_m, dL=0):
    if f_m is None:
        mat = np.matrix(
            [
                [1, 0],
                [0, 1],
            ]
        )
    else:
        mat = np.matrix(
            [
                [1, 0],
                [-1 / f_m, 1],
            ]
        )
    if not dL:
        return mat
    return matrix_space(dL) * mat * matrix_space(-dL)


def matrix_telescope(L1, R, dL=0):
    L2 = R * L1
    mat = matrix_focus(-L2) * matrix_space(L1 - L2 - dL) * matrix_focus(L1)
    return mat


def eigen_q(mat):
    pe_A = mat[..., 0, 0]
    pe_B = mat[..., 0, 1]
    pe_C = mat[..., 1, 0]
    pe_D = mat[..., 1, 1]
    q = (pe_A - pe_D) / (2 * pe_C) + 1j * np.sqrt(
        -((pe_D - pe_A) ** 2 + 4 * pe_B * pe_C)
    ) / (2 * abs(pe_C))
    return q


def targets_map_append(targets_map, tname, *targets):
    tlstlst = targets_map.get(tname, None)
    if tlstlst is None:
        tlstlst = []
        targets_map[tname] = tlstlst
    if targets:
        tlstlst.extend(targets)
    return tlstlst


def unit_str(val, unit, d=3, use_c=False, space=True):
    val = float(val)
    v = abs(val)
    suffix = ""
    if v > 1e3:
        prefix = "k"
        div = 1e3
    elif v > 1:
        prefix = ""
        if space:
            suffix = " "
        div = 1
    elif v > 1e-2 and use_c:
        prefix = "c"
        div = 1e2
    elif v > 1e-3:
        prefix = "m"
        div = 1e-3
    elif v > 1e-6:
        # prefix = u'μ'
        prefix = "u"
        div = 1e-6
    elif v > 1e-9:
        prefix = "n"
        div = 1e-9
    elif v > 1e-15:
        prefix = "f"
        div = 1e-15
    elif v > 1e-18:
        prefix = "a"
        div = 1e-18
    elif v > 1e-21:
        prefix = "z"
        div = 1e-21
    elif v == 0:
        div = 1
        prefix = ""
        if space:
            suffix = " "
    nval = val / div
    if abs(nval) >= 100:
        d -= 2
    elif abs(nval) >= 10:
        d -= 1
    if d < 0:
        d = 0
    if space:
        return u"{0: .{1}f}{2}{3}{4}".format(nval, d, prefix, unit, suffix)
    else:
        return u"{0:.{1}f}{2}{3}{4}".format(nval, d, prefix, unit, suffix)


def str_m(val, d=3, use_c=False, space=True):
    return unit_str(val, d=d, unit="m", use_c=use_c, space=space)


def str_D(val, d=3, use_c=False, space=True):
    return unit_str(val, d=d, unit="D", use_c=use_c, space=space)


def interface_ROC(ROC_m, n_from, n_to, neg=False):
    """
    ABCD matrix for transmission through a curved surface at normal incidence

    Parameters
    ----------
    ROC_m : float
        Radius of curvature of the surface [m]
    n_from : float
        Index of refraction for the incident beam
    n_to : float
        Index of refraction for the transmitted beam
    neg : bool; default: False
        If True, reverses the sign of the radius of curvature
    """
    nft = n_from / n_to
    if ROC_m is not None:
        if neg:
            ROC_m = -ROC_m
        return matrix_stack(
            [
                [1, 0],
                [(nft - 1) / ROC_m, nft],
            ]
        )
    else:
        return matrix_stack(
            [
                [1, 0],
                [0, nft],
            ]
        )


def interface_ROC_AOI_Y(ROC_m, n_from, n_to, AOI_rad, neg=False):
    """
    ABCD matrix for transmission through a curved surface at arbitrary angle of incidence in the
    saggital plane

    Parameters
    ----------
    ROC_m : float
        Radius of curvature of the surface [m]
    n_from : float
        Index of refraction for the incident beam
    n_to : float
        Index of refraction for the transmitted beam
    AOI_rad : float
        Angle of incidence [rad]
    neg : bool; default: False
        If True, reverses the sign of the radius of curvature
    """
    nft = n_from / n_to
    if ROC_m is not None:
        adj = (1 - (nft * np.sin(AOI_rad)) ** 2) ** 0.5
        if neg:
            ROC_m = -ROC_m
        return matrix_stack(
            [
                [1, 0],
                [(nft * np.cos(AOI_rad) - adj) / ROC_m, nft],
            ]
        )
    else:
        return matrix_stack(
            [
                [1, 0],
                [0, nft],
            ]
        )


def interface_ROC_AOI_X(ROC_m, n_from, n_to, AOI_rad, neg=False):
    """
    ABCD matrix for transmission through a curved surface at arbitrary angle of incidence in the
    tangential plane

    Parameters
    ----------
    ROC_m : float
        Radius of curvature of the surface [m]
    n_from : float
        Index of refraction for the incident beam
    n_to : float
        Index of refraction for the transmitted beam
    AOI_rad : float
        Angle of incidence [rad]
    neg : bool; default: False
        If True, reverses the sign of the radius of curvature
    """
    nft = n_from / n_to
    if ROC_m is not None:
        adj = (1 - (nft * np.sin(AOI_rad)) ** 2) ** 0.5
        if neg:
            ROC_m = -ROC_m
        return matrix_stack(
            [
                [adj / np.cos(AOI_rad), 0],
                [
                    (nft / adj - 1 / np.cos(AOI_rad)) / ROC_m,
                    nft * np.cos(AOI_rad) / adj,
                ],
            ]
        )
    else:
        return matrix_stack(
            [
                [1, 0],
                [0, nft],
            ]
        )


def REFL_ROC_Y(ROC_m, AOI_rad):
    """
    ABCD matrix for reflection from a curved surface at arbitrary angle of incidence in the
    saggital plane

    Parameters
    ----------
    ROC_m : float
        Radius of curvature of the surface [m]
    AOI_rad : float
        Angle of incidence [rad]
    """
    if ROC_m is not None:
        return matrix_stack(
            [
                [1, 0],
                [2 * np.cos(AOI_rad) / ROC_m, 1],
            ]
        )
    else:
        return matrix_stack(
            [
                [1, 0],
                [0, 1],
            ]
        )


def REFL_ROC_X(ROC_m, AOI_rad):
    """
    ABCD matrix for reflection from a curved surface at arbitrary angle of incidence in the
    tangential plane

    Parameters
    ----------
    ROC_m : float
        Radius of curvature of the surface [m]
    AOI_rad : float
        Angle of incidence [rad]
    """
    if ROC_m is not None:
        return matrix_stack(
            [
                [1, 0],
                [2 / (ROC_m * np.cos(AOI_rad)), 1],
            ]
        )
    else:
        return matrix_stack(
            [
                [1, 0],
                [0, 1],
            ]
        )
