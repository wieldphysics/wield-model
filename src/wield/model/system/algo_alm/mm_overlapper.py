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
from .plot_alm import OverlapperPlotter
from wavestate.bunch import Bunch
from . import mm_annotate
from wavestate.utilities.strings import table
from ... import optics


class ModeMatchingOverlapper(object):
    plotter = OverlapperPlotter()

    def __init__(
        self,
        algo_pa,
        algo_mm,
        targetsB_to,
        targetsB_fr,
        oLp_path_center,
        Wk,
        branching,
        pbg=None,
        shifts_use=False,
    ):
        self.pa = algo_pa
        self.mm = algo_mm
        self.targetsB_to = targetsB_to
        self.targetsB_fr = targetsB_fr
        self.oLp_path_center = oLp_path_center
        self.target1 = None
        self.target2 = None
        self.Wk = Wk
        self.branching = branching
        self.shifts_use = shifts_use

        if pbg is None:
            self.pbg = self.mm.pbg
        else:
            self.pbg = pbg

        self.trans_center = self.mm._path_transporters(self.oLp_path_center, Wk=Wk, shifts_use=shifts_use)

        def setup_refer_to_start(tB, direction):
            tB_new = Bunch()
            tB_new.type = tB.tspecB
            if tB.inv_start:
                tB_new.targB = self.mm._target_complete(tB.tspecB, ol=tB.oLp_path[-1], shifts_use=shifts_use)
            else:
                tB_new.targB = self.mm._target_complete(tB.tspecB, ol=tB.oLp_path[0], shifts_use=shifts_use)

            if not tB.inv_start and tB_new.type == 'cavity':
                shiftsX = tB_new.targB.cav_shiftX
                shiftsY = tB_new.targB.cav_shiftY
            else:
                # TODO, allow reverse shifts
                shiftsX = {}
                shiftsY = {}

            tB_new.trans = self.mm._path_transporters(tB.oLp_path, Wk=Wk, shifts_use=shifts_use)
            matXfr = tB_new.trans.X.full_trip_mat
            matYfr = tB_new.trans.Y.full_trip_mat
            if tB.inv_start:
                matXfr = np.linalg.inv(matXfr)
                matYfr = np.linalg.inv(matYfr)
            else:
                shiftsX = tB_new.trans.X.shifts_out(shiftsX)
                shiftsY = tB_new.trans.Y.shifts_out(shiftsY)

            tB_new.inv_start = tB.inv_start
            tB_new.type = direction
            if direction == "from":
                pass
            elif direction == "to":
                matXfr = np.linalg.inv(self.trans_center.X.full_trip_mat) @ matXfr
                matYfr = np.linalg.inv(self.trans_center.Y.full_trip_mat) @ matYfr
            else:
                raise RuntimeError("Bad Direction")

            tB_new.qX = tB_new.targB.qX.propagate_matrix(matXfr)
            tB_new.qY = tB_new.targB.qY.propagate_matrix(matYfr)
            tB_new.qXend = tB_new.qX.propagate_matrix(self.trans_center.X.full_trip_mat)
            tB_new.qYend = tB_new.qY.propagate_matrix(self.trans_center.Y.full_trip_mat)

            if not tB.inv_start:
                shiftsX = self.trans_center.X.shifts_out(shiftsX)
                shiftsY = self.trans_center.Y.shifts_out(shiftsY)

            tB_new.shiftsXend = shiftsX
            tB_new.shiftsYend = shiftsY
            return tB_new

        transB_fr = dict()
        # overlap is referred to the start of the path
        # TODO, collapse the logic of both of these to a single function call
        for t_fr, frB in targetsB_fr.items():
            transB_fr[t_fr] = setup_refer_to_start(frB, direction="from")

        transB_to = dict()
        for t_to, toB in targetsB_to.items():
            transB_to[t_to] = setup_refer_to_start(toB, direction="to")

        self.transB_to = transB_to
        self.transB_fr = transB_fr
        return

    def rebuild(self, pbg=None, overrides=None):
        if pbg is None:
            pbg = self.pbg
        if overrides is not None:
            pbg = pbg.copy()
            for override_param, override_val in overrides.items():
                pbg.override_value(override_param, override_val)
        olap = self.__class__(
            algo_pa=self.pa,
            algo_mm=self.mm,
            targetsB_to=self.targetsB_to,
            targetsB_fr=self.targetsB_fr,
            oLp_path_center=self.oLp_path_center,
            Wk=self.Wk,
            branching=self.branching,
            pbg=pbg,
        )
        olap.target1 = self.target1
        olap.target2 = self.target2
        return olap

    def propagate_reference(self, target, name=None):
        raise NotImplementedError("Not sure this is working as expected")
        if name is None:
            name = target + " reference"
        if target in self.transB_fr:
            self.target1 = target
            self.target2 = name
            frB = self.transB_fr[target]
            toB_new = Bunch()
            # makes a null transfer
            toB_new.trans = self.mm._path_transporters([], Wk=self.Wk, shifts_use=self.shifts_use)
            qX = frB.qX.propagate_matrix(self.trans_center.X.full_trip_mat)
            qY = frB.qY.propagate_matrix(self.trans_center.Y.full_trip_mat)
            toB_new.qX = qX
            toB_new.qY = qY
            toB_new.targB = Bunch()
            toB_new.targB.type = "specified"
            toB_new.targB.qX = qX
            toB_new.targB.qY = qY
            toB_new.type = "to"
            toB_new.inv_start = not frB.inv_start
            self.transB_to[name] = toB_new
        else:
            toB = self.transB_to[target]
            self.target2 = target
            self.target1 = name
            frB_new = Bunch()
            # makes a null transfer
            frB_new.trans = self.mm._path_transporters([], Wk=self.Wk, shifts_use=self.shifts_use)
            qX = toB.qX
            qY = toB.qY
            frB_new.qX = qX
            frB_new.qY = qY
            frB_new.targB = Bunch()
            frB_new.targB.type = "specified"
            frB_new.targB.qX = qX
            frB_new.targB.qY = qY
            frB_new.inv_start = not toB.inv_start
            frB_new.type = "from"
            self.transB_fr[name] = frB_new
        return

    def set_targets(self, target1, target2):
        self.target1 = target1
        self.target2 = target2
        return

    def target_list(self):
        targets_set = set(self.targetsB_fr.keys())
        targets_set.update(self.targetsB_to.keys())

        targets = []
        if self.target1 is not None:
            targets_set.remove(self.target1)
            targets.append(self.target1)
        if self.target2 is not None:
            targets_set.remove(self.target2)
            targets.append(self.target2)
        targets.extend(sorted(targets_set))
        return targets

    def object_z(self, ref, obj=None):
        # TODO
        lset = self.mm.bg.rAp2oLp_set(ref, obj=obj)
        print("lset", lset)
        # TODO, specify X or Y or check that it doesn't matter
        d = self.trans_center.X.ol2z(lset)
        return min(d.values())

    @property
    def length_m(self):
        return self.trans_center.X.full_trip_length

    def __getitem__(self, tname):
        """
        Indexing the overlapper by a target name
        returns the bunch containing the target information
        """
        ret = self.transB_fr.get(tname, None)
        if ret is not None:
            return ret
        return self.transB_to[tname]

    def z2target_qY(self, tname, Z):
        tB = self[tname]
        mat = self.trans_center.Y.z2mat(Z)
        return tB.qY.propagate_matrix(mat)

    def z2target_qX(self, tname, Z):
        tB = self[tname]
        mat = self.trans_center.X.z2mat(Z)
        return tB.qX.propagate_matrix(mat)

    @property
    def overlapX_field(self):
        frB = self.transB_fr[self.target1]
        toB = self.transB_to[self.target2]
        return frB.qX.overlap_HG(toB.qX)

    @property
    def overlapY_field(self):
        frB = self.transB_fr[self.target1]
        toB = self.transB_to[self.target2]
        return frB.qY.overlap_HG(toB.qY)

    @property
    def overlap_field(self):
        try:
            frB = self.transB_fr[self.target1]
        except KeyError:
            frB = self.transB_to[self.target1]
        try:
            toB = self.transB_to[self.target2]
        except KeyError:
            toB = self.transB_fr[self.target2]
        return frB.qX.overlap_HG(toB.qX) * frB.qY.overlap_HG(toB.qY)

    @property
    def overlap(self):
        return abs(self.overlap_field) ** 2

    def plot(self, fname=None, *args, **kwargs):
        return self.plotter.plot(*args, overlapper=self, fname=fname, **kwargs)

    def plot_scan(self, fname=None, *args, **kwargs):
        return self.plotter.plot_scan(*args, overlapper=self, fname=fname, **kwargs)

    def plot_descriptions(
        self,
        axis="Y",
        waists_target=None,
        reverse=False,
        tags=[],
    ):
        """
        tags is a list of dictionaries that look like

        dict(
        obj = '/path/to/obj' or obj,
        to = distance[m], or
        fr = distance[m],
        label = "label W={W}, Gouy={gouy}deg ..."
        )
        """
        descriptions = self.object_descriptions(
            axis=axis,
            waists_target=waists_target,
        )

        descriptions = mm_annotate.annotate_tags(self.pbg, descriptions, tags)
        descriptions = mm_annotate.annotate_waists(descriptions)
        descriptions = mm_annotate.annotate_qspaces(descriptions, reverse=reverse)
        descriptions = mm_annotate.annotate_clearspaces(descriptions)

        self._plot_descriptions = descriptions
        return self._plot_descriptions

    def annotation_target(self, waists_target=None):
        if waists_target is None:
            if self.target1 is not None:
                waists_target = 1
            else:
                waists_target = 2

        if waists_target == 1:
            waists_target = self.target1
        elif waists_target == 2:
            waists_target = self.target2
        return waists_target

    def object_descriptions(
        self,
        axis="Y",
        waists_target=None,
    ):
        waists_target = self.annotation_target(waists_target)
        tB = self[waists_target]

        if axis.lower() == "y":
            descriptions = mm_annotate.transporter2objects(
                self.mm,
                tB.qY,
                self.trans_center.Y,
            )
        elif axis.lower() == "x":
            descriptions = mm_annotate.transporter2objects(
                self.mm,
                tB.qX,
                self.trans_center.X,
            )
        else:
            raise RuntimeError("Unrecognized Axis")

        return descriptions

    def overlap_field_between(self, target1, target2):
        try:
            frB = self.transB_fr[target1]
        except KeyError:
            frB = self.transB_to[target1]
        try:
            toB = self.transB_to[target2]
        except KeyError:
            toB = self.transB_fr[target2]
        return frB.qX.overlap_HG(toB.qX) * frB.qY.overlap_HG(toB.qY)

    def overlap_table(self, first=(), last=()):
        """
        First and last modify the ordering
        """
        targets = list(self.transB_fr.keys()) + list(self.transB_to.keys())
        firsts = []
        mids = []
        lasts = []
        for idx, targ in enumerate(targets):
            for f_ in first:
                if f_ in targ:
                    firsts.append(targ)
                    break
            else:
                for f_ in last:
                    if f_ in targ:
                        lasts.append(targ)
                        break
                else:
                    mids.append(targ)
        targets = firsts[::-1] + mids + lasts

        olaps = []
        for t1 in targets:
            olaps2 = []
            for t2 in targets:
                olaps2.append(abs(self.overlap_field_between(t1, t2)) ** 2)
            olaps.append(olaps2)
        alpha = "ABCDEFGHIJKLMNOPQRSTUV"
        headers = []
        labels = []
        for idx, t1 in enumerate(targets):
            labels.append("({}) {}".format(alpha[idx], t1))
            headers.append("({})".format(alpha[idx]))
        return table(
            olaps,
            headers=headers,
            labels=labels,
            # headers_modify = headers_modify,
            diag=None,
        )

    def gouy_table(self, objects=None, **kwargs):
        descriptions = self.object_descriptions(**kwargs)
        q_start = descriptions[0].q_start

        descB_list = []
        for descB in descriptions:
            if isinstance(descB.object, optics.Space):
                continue
            if objects is not None:
                if descB.object in objects:
                    pass
                elif descB.desc in objects:
                    pass
                else:
                    continue
            descB_list.append(descB)

        obj_list = []
        desc_list = []
        gouy_table = []
        gouy_list = []
        diameters_list = []
        q_list = []
        q2_list = []
        length_table = []
        M_table = []
        A_table = []
        B_table = []
        C_table = []
        D_table = []
        for idx_to, descB_to in enumerate(descB_list):
            obj_list.append(descB_to.object)
            desc_list.append(descB_to.desc)
            q_list.append(descB_to.q_start)
            q2_list.append(descB_to.q_end)
            gouy_diff = []
            length_diff = []
            gouy_list.append(
                np.angle(descB_to.q_start.gouy_phasor / q_start.gouy_phasor, deg=True)
            )
            diameters_list.append(2 * descB_to.q_start.W)
            M_list = []
            A_list = []
            B_list = []
            C_list = []
            D_list = []
            for idx_fr, descB_fr in enumerate(descB_list):
                if idx_fr < idx_to:
                    gouy_diff.append(
                        np.angle(
                            (
                                descB_to.q_start.gouy_phasor
                                / descB_fr.q_start.gouy_phasor
                            ),
                            deg=True,
                        )
                    )
                    length_diff.append(descB_to.z1_m - descB_fr.z2_m)
                    M = descB_to.mat1 @ np.linalg.inv(descB_fr.mat1)
                elif idx_fr == idx_to:
                    gouy_diff.append(
                        np.angle(
                            (descB_to.q_end.gouy_phasor / descB_fr.q_start.gouy_phasor),
                            deg=True,
                        )
                    )
                    length_diff.append(descB_to.z2_m - descB_to.z1_m)
                    M = descB_to.mat2 @ np.linalg.inv(descB_fr.mat1)
                else:
                    gouy_diff.append(
                        np.angle(
                            (descB_to.q_end.gouy_phasor / descB_fr.q_end.gouy_phasor),
                            deg=True,
                        )
                    )
                    length_diff.append(descB_fr.z1_m - descB_to.z2_m)
                    M = descB_to.mat2 @ np.linalg.inv(descB_fr.mat2)
                M_list.append(M)
                A_list.append(M[0, 0])
                B_list.append(M[0, 1])
                C_list.append(M[1, 0])
                D_list.append(M[1, 1])
            gouy_table.append(gouy_diff)
            length_table.append(length_diff)
            M_table.append(M_list)
            A_table.append(A_list)
            B_table.append(B_list)
            C_table.append(C_list)
            D_table.append(D_list)

        gouy_table = np.array(gouy_table)
        gouy_list = np.array(gouy_list)
        length_table = np.array(length_table)
        diameters_list = np.array(diameters_list)
        A_table = np.array(A_table)
        B_table = np.array(B_table)
        C_table = np.array(C_table)
        D_table = np.array(D_table)

        def gouy_table_str(
            headers=desc_list,
            labels=desc_list,
            headers_modify="bind",
            units="deg",
            diag=None,
            **kwargs
        ):
            if units == "deg":
                umult = 1
            elif units == "rad":
                umult = np.pi / 180
            else:
                raise RuntimeError("Unrecognized units")

            if diag is None:
                diag = "Gouy separations [{}]".format(units)

            return table(
                umult * gouy_table,
                headers=headers,
                labels=labels,
                headers_modify=headers_modify,
                diag=diag,
                **kwargs
            )

        def length_table_str(
            headers=desc_list,
            labels=desc_list,
            headers_modify="bind",
            units="in",
            diag=None,
            **kwargs
        ):
            if units == "in":
                umult = 1 / 0.0254
            elif units == "m":
                umult = 1
            else:
                raise RuntimeError("Unrecognized units")

            if diag is None:
                diag = "Length separations [{}]".format(units)

            return table(
                umult * length_table,
                headers=headers,
                labels=labels,
                headers_modify=headers_modify,
                diag=diag,
                **kwargs
            )

        def diameters_str(
            headers=["diameter [um]"], labels=desc_list, units="um", **kwargs
        ):
            if units == "um":
                umult = 1e6
            else:
                raise RuntimeError("Unrecognized units")

            return table(
                umult * diameters_list.reshape(-1, 1),
                headers=headers,
                labels=labels,
                **kwargs
            )

        def Qs_str(headers=None, labels=desc_list, side="input", **kwargs):
            if side == "input":
                if headers is None:
                    headers = ["Beam Q's", "incoming side"]
                qs = np.array([q.string() for q in q_list])
                return table(
                    qs.reshape(-1, 1), headers=headers, labels=list(labels), **kwargs
                )
            elif side == "output":
                if headers is None:
                    headers = ["Beam Q's", "outgoing side"]
                qs = np.array([q.string() for q in q2_list])
                return table(
                    qs.reshape(-1, 1),
                    headers=headers,
                    labels=list(labels) + list(labels),
                    **kwargs
                )

        def A_table_str(
            headers=desc_list,
            labels=desc_list,
            headers_modify="bind",
            diag=None,
            **kwargs
        ):
            if diag is None:
                diag = "Displacement Gain (Abcd) [m/m]"

            return table(
                A_table,
                headers=headers,
                labels=labels,
                headers_modify=headers_modify,
                diag=diag,
                **kwargs
            )

        def D_table_str(
            headers=desc_list,
            labels=desc_list,
            headers_modify="bind",
            diag=None,
            **kwargs
        ):
            if diag is None:
                diag = "Angle Gain (abcD) [rad/rad]"

            return table(
                D_table,
                headers=headers,
                labels=labels,
                headers_modify=headers_modify,
                diag=diag,
                **kwargs
            )

        def B_table_str(
            headers=desc_list,
            labels=desc_list,
            headers_modify="bind",
            diag=None,
            units="mm/mrad",
            **kwargs
        ):
            if units == "in":
                umult = 1 / 0.0254
            elif units == "m/rad":
                umult = 1
            elif units == "mm/mrad":
                umult = 1
            elif units == "mm/rad":
                umult = 1e3
            else:
                raise RuntimeError("Unrecognized units")

            if diag is None:
                diag = "Deflection (aBcd) [{}]".format(units)

            return table(
                umult * B_table,
                headers=headers,
                labels=labels,
                headers_modify=headers_modify,
                diag=diag,
                **kwargs
            )

        def C_table_str(
            headers=desc_list,
            labels=desc_list,
            headers_modify="bind",
            diag=None,
            units="mrad/mm",
            **kwargs
        ):
            if units == "in":
                umult = 1 / 0.0254
            elif units == "rad/m":
                umult = 1
            elif units == "mrad/mm":
                umult = 1
            elif units == "rad/mm":
                umult = 1e-3
            else:
                raise RuntimeError("Unrecognized units")

            if diag is None:
                diag = "Deflection (abCd) [{}]".format(units)

            return table(
                umult * C_table,
                headers=headers,
                labels=labels,
                headers_modify=headers_modify,
                diag=diag,
                **kwargs
            )

        return Bunch(
            A_table_str=A_table_str,
            B_table_str=B_table_str,
            C_table_str=C_table_str,
            D_table_str=D_table_str,
            A_table=A_table,
            B_table=B_table,
            Qs_str=Qs_str,
            q_list=q_list,
            q2_list=q2_list,
            obj_list=obj_list,
            desc_list=desc_list,
            descB_list=descB_list,
            gouy_table=gouy_table,
            gouy_table_str=gouy_table_str,
            gouy_list=gouy_list,
            length_table=length_table,
            length_table_str=length_table_str,
            diameters_list=diameters_list,
            diameters_str=diameters_str,
        )

    def shifts_table(
        self,
        axis="y",
        waists_target=None,
    ):
        waists_target = self.annotation_target(waists_target)
        tB = self[waists_target]

        shifts_pos = {}
        shifts_ang = {}

        if axis.lower() == "y":
            trans_center = self.trans_center.Y
        elif axis.lower() == "x":
            trans_center = self.trans_center.X

        # use the ol2idx mappings to determine which matrices to use based on the object locations within the path
        # to use it, it must first be inverted
        reverse_ol2idx = {}
        for k, v in trans_center.inc_ol2idx.items():
            l = reverse_ol2idx.setdefault(v, [])
            l.append(k)
        idxs = []
        keys = []
        
        for idx, key in sorted(reverse_ol2idx.items()):
            idxs.append(idx)
            keys.append(key)

        if axis.lower() == "y":
            mats = np.array([self.trans_center.Y.inc_build_mat[i] for i in idxs]) @ np.linalg.inv(self.trans_center.Y.inc_build_mat[-1])
            for shift_key, shift in tB.shiftsYend.items():
                shifts = mats @ shift
                shifts_pos[shift_key] = shifts[..., 0, 0]
                shifts_ang[shift_key] = shifts[..., 1, 0]
        elif axis.lower() == "x":
            mats = np.array([self.trans_center.X.inc_build_mat[i] for i in idxs]) @ np.linalg.inv(self.trans_center.X.inc_build_mat[-1])
            for shift_key, shift in tB.shiftsXend.items():
                shifts = mats @ shift
                shifts_pos[shift_key] = shifts[..., 0, 0]
                shifts_ang[shift_key] = shifts[..., 1, 0]
        else:
            raise RuntimeError("Unrecognized Axis")

        return Bunch(
            shifts_keys = keys,
            shifts_pos = shifts_pos,
            shifts_ang = shifts_ang
        )

    def shifts_table_str(
        self,
        axis="y",
        var="pos",
        waists_target=None,
    ):
        """
        Creates a table of the shifts from optic motions to and from cavities.

        Currently it drops spaces in favor of non-space optics when forming the labels
        TODO: make the optic naming better for the labels and headers of the table.
        """
        axis = axis.lower()
        var = var.lower()
        assert(axis in ['x', 'y'])
        assert(var in ['pos', 'ang'])

        sB = self.shifts_table(
            axis = axis,
            waists_target=waists_target,
        )

        if var == 'pos':
            shifts = sB.shifts_pos
        elif var == 'ang':
            shifts = sB.shifts_ang

        keys = sorted(shifts.keys())

        vals = np.array([shifts[k] for k in keys]).T

        self.trans_center.X

        labels = []
        for shift_key_list in sB.shifts_keys:
            for sk in shift_key_list:
                if isinstance(sk[0], optics.Space):
                    continue
                else:
                    break
            labels.append(str(sk))
        return table(
            vals,
            headers=[str(k) for k in keys],
            labels=labels,
            diag=var,
        )



