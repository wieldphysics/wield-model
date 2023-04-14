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
from wavestate.bunch import Bunch


class MMTransporter(object):
    def __init__(
        self,
        oLp_path,
        Wk,
        prop,
        inc,
        prop_ol2idx,
        inc_ol2idx,
        shifts,
    ):
        """
        shifts is a list of bunch elements with 

        Bunch(
            idx_inc=len(Xinc),
            idx_prop=len(Xprop),
            name=name,
            obj=obj,
            shift=shift,
        )
        where idxXprop is the index into the prop list that the shift applies.
        name is the shift type, obj is the object it applies to.
        shift is a 2-vector with the shift derivative. The units should be understood
        from the name, such as
        [0; 2] would be the beam angle change in rad from rad of optic/mirror tilt.
        """
        self.oLp_path = oLp_path
        self.Wk = Wk
        self.prop = prop
        self.inc = inc
        self.prop_ol2idx = prop_ol2idx
        self.inc_ol2idx = inc_ol2idx

        mat_inc = np.eye(2)
        len_inc = 0
        build_mat = [mat_inc]
        build_len = [len_inc]
        for (len_m, func, mat) in self.inc:
            mat_inc = mat @ mat_inc
            len_inc += len_m
            build_mat.append(mat_inc)
            build_len.append(len_inc)
        # these are the incrementally built propation matrices
        self.inc_build_mat = build_mat
        self.inc_build_len = build_len
        self.shifts = shifts

    _shifts_in_referred = None
    @property
    def shifts_in_referred(self):
        """
        Collects the list of shifts through the transporter and casts them all to be
        at the input of the transport. If the same shift occurs multiple times,
        it is co-added.
        """
        if self._shifts_in_referred is None:
            _shifts_in_referred = {}
            for sB in self.shifts:
                M = self.inc_build_mat[sB.idx_inc]
                shift = np.linalg.inv(M) @ sB.shift
                shift_prev = _shifts_in_referred.setdefault((sB.obj, sB.name), shift)
                if shift_prev is not shift:
                    _shifts_in_referred[(sB.obj, sB.name)] = shift_prev + shift

            self._shifts_in_referred = _shifts_in_referred
        return self._shifts_in_referred

    _shifts_out_referred = None
    @property
    def shifts_out_referred(self):
        """
        Collects the list of shifts through the transporter and casts them all to be
        at the output of the transport. If the same shift occurs multiple times,
        it is co-added.
        """
        if self._shifts_out_referred is None:
            _shifts_out_referred = {}
            for sB in self.shifts:
                M = self.inc_build_mat[sB.idx_inc]
                shift = self.inc_build_mat[-1] @ np.linalg.inv(M) @ sB.shift
                shift_prev = _shifts_out_referred.setdefault((sB.obj, sB.name), shift)
                if shift_prev is not shift:
                    _shifts_out_referred[(sB.obj, sB.name)] = shift_prev + shift
            self._shifts_out_referred = _shifts_out_referred
        return self._shifts_out_referred

    def shifts_out(self, shifts_in, include_own = True):
        """
        Takes the shift_in dictionary and applies this path transporter, adding its own shifts along the way
        """
        shifts_out = {}
        for shift_key, shift in shifts_in.items():
            shift = self.inc_build_mat[-1] @ shift
            # doesn't need  to check for an alias, since the input can't have any
            shifts_out[shift_key] = shift
        if include_own:
            for sB in self.shifts:
                M = self.inc_build_mat[sB.idx_inc]
                shift = self.inc_build_mat[-1] @ np.linalg.inv(M) @ sB.shift
                # use shift_prev to check for an alias with an existing shift
                shift_prev = shifts_out.setdefault((sB.obj, sB.name), shift)
                if shift_prev is not shift:
                    shifts_out[(sB.obj, sB.name)] = shift_prev + shift
        return shifts_out

    @property
    def full_trip_mat(self):
        """
        Returns the matrix for the entire transport chain
        """
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
        """
        Returns the path length distance z for a given object-link (ol)
        """
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
        """
        Creates the beam propation up to a given Z value
        """
        Z = np.asarray(Z)
        Zorig = Z
        if Z.shape == ():
            Z = Z.reshape(1)
        mat_by_z = []
        idx_mat = 1
        matX_last = np.eye(2)
        z_last = 0
        z_next = self.inc_build_len[idx_mat]
        for idx in np.argsort(Z):
            z = Z[idx]
            while z > z_next:
                matX_last = self.inc_build_mat[idx_mat]
                z_last = self.inc_build_len[idx_mat]
                if idx_mat + 1 >= len(self.inc_build_len):
                    z_next = float("infinity")
                else:
                    idx_mat += 1
                    z_next = self.inc_build_len[idx_mat]
            z_diff = z - z_last
            if z_diff < 0:
                NaN = float("NaN")
                mat_by_z.append(np.array([[NaN, NaN], [NaN, NaN]]))
                continue
            if z > self.inc_build_len[-1]:
                NaN = float("NaN")
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
