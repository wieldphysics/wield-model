#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import numpy as np
import declarative


class MMTransporter(object):
    def __init__(
        self,
        oLp_path,
        Wk,
        prop,
        inc,
        prop_ol2idx,
        inc_ol2idx,
    ):
        self.oLp_path    = oLp_path
        self.Wk        = Wk
        self.prop      = prop
        self.inc       = inc
        self.prop_ol2idx  = prop_ol2idx
        self.inc_ol2idx   = inc_ol2idx

        mat_inc = np.eye(2)
        len_inc = 0
        build_mat = [mat_inc]
        build_len = [len_inc]
        for (len_m, func, mat) in self.inc:
            mat_inc = mat @ mat_inc
            len_inc += len_m
            build_mat.append(mat_inc)
            build_len.append(len_inc)
        self.inc_build_mat = build_mat
        self.inc_build_len = build_len

    @property
    def full_trip_mat(self):
        if not self.inc_build_mat:
            return np.eye(2)
        return self.inc_build_mat[-1]

    @property
    def full_trip_length(self):
        if not self.inc_build_len:
            return 0
        return self.inc_build_len[-1]

    def ol2mat(self, ol):
        idx = self.inc_ol2idx[ol]
        if idx is None:
            return np.eye(2)
        else:
            return self.inc_build_mat[idx]

    def ol2mat_end(self, ol):
        idx = self.inc_ol2idx[ol]
        if idx is None:
            return np.eye(2)
        mat_inc = np.eye(2)
        for (len_m, func, mat) in self.inc[idx:]:
            mat_inc = mat @ mat_inc
        return mat_inc

    def ol2z(self, ol):
        if isinstance(ol, (tuple, list)):
            d = dict()
            for _ol in ol:
                try:
                    d[_ol] = self.inc_build_len[self.inc_ol2idx[_ol]]
                except KeyError:
                    pass
            return d
        else:
            return self.inc_build_len[self.inc_ol2idx[ol]]


    def z2mat(self, Z):
        Z = np.asarray(Z)
        Zorig = Z
        if Z.shape == ():
            Z = Z.reshape(1)
        mat_by_z  = []
        idx_mat   = 1
        matX_last = np.eye(2)
        z_last    = 0
        z_next    = self.inc_build_len[idx_mat]
        for idx in np.argsort(Z):
            z = Z[idx]
            while z > z_next:
                matX_last = self.inc_build_mat[idx_mat]
                z_last = self.inc_build_len[idx_mat]
                if idx_mat + 1 >= len(self.inc_build_len):
                    z_next = float('infinity')
                else:
                    idx_mat += 1
                    z_next = self.inc_build_len[idx_mat]
            z_diff = z - z_last
            if z_diff < 0:
                NaN = float('NaN')
                mat_by_z.append(np.array([[NaN, NaN], [NaN, NaN]]))
                continue
            if z > self.inc_build_len[-1]:
                NaN = float('NaN')
                mat_by_z.append(np.array([[NaN, NaN], [NaN, NaN]]))
                continue

            len_m, func, mat = self.inc[idx_mat - 1]
            if func is None:
                matX = mat @ matX_last
            else:
                matX = func(z_diff) @ matX_last
            mat_by_z.append(matX)
        if Zorig.shape == ():
            return mat_by_z[0]
        else:
            return np.array(mat_by_z)

