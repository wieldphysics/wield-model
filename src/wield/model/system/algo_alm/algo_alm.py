#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import numpy as np
import collections
from wield.bunch import Bunch

# from ..matrix import SRE_copy
from ...optics import alm

from . import algo_mm_linkages
from . import mm_optimize
from .. import algo_tups


class ModeMatchingAlgorithm(algo_mm_linkages.ModeMatchingLinkageAlgorithm):
    def __init__(self, pa, previous=None):
        super(ModeMatchingAlgorithm, self).__init__(pa)

        # from target names to ol-loops
        # target names are parameter refs
        self.cavity_targets = dict()
        self.cavity_map = dict()

        self._targets = dict()

        # a dictionary of target names to be constructed upon access
        self._targets_deferred = dict()

        self._populate_targets()

        # This is a bit of a hack to forward the explicit calls into this object
        self.previous = previous
        if self.previous is not None:
            for k, v in self.previous._targets.items():
                self._targets.setdefault(k, v)
            for k, v in self.previous.cavity_targets.items():
                self.cavity_targets.setdefault(k, v)
            for k, v in self.previous.cavity_map.items():
                self.cavity_map.setdefault(k, v)
        return

    def _populate_targets(self):
        """
        Visits objects in the graph looking for ModeMatchingTargetBase objects
        such as Cavity, Target and TargetMeasurement. Those objects then populate
        the target list
        """
        for obj in self.pbg.object_iter():
            try:
                visit_algo = obj.visit_mode_matching_targets
            except AttributeError:
                # don't actually run the algorithm within the try because it can
                # eat exceptions that occur within the visit method
                continue
            else:
                pass

            manip = MMAlgorithmManipulator(
                obj=obj,
                mm_algo=self,
            )

            visit_algo(manip)
        return

    def target_add(
        self,
        target_name,
        waypoints,
        q=None,
        qX=None,
        qY=None,
        z_target=0,
        z_start=0,
        wavelength=None,
        obj=None,
    ):
        if qX is None:
            qX = q
        if qY is None:
            qY = q

        assert qX is not None
        assert qY is not None

        if isinstance(waypoints, str):
            waypoints = [waypoints]

        if wavelength is not None:
            Wk = self.fs.parameter_to_wk(wavelength)
        else:
            Wk = None

        if isinstance(q, alm.ComplexBeamParam):
            if Wk is None:
                Wk = self.fs.parameter_to_wk(q.wavelength_m)
            else:
                assert self.fs.parameter_to_wk(q.wavelength_m) == Wk

        if isinstance(qX, alm.ComplexBeamParam):
            if Wk is None:
                Wk = self.fs.parameter_to_wk(qX.wavelength_m)
            else:
                assert self.fs.parameter_to_wk(qX.wavelength_m) == Wk

        # print(qY, Wk)
        if isinstance(qY, alm.ComplexBeamParam):
            if Wk is None:
                Wk = self.fs.parameter_to_wk(qY.wavelength_m)
            else:
                assert self.fs.parameter_to_wk(qY.wavelength_m) == Wk

        if Wk is None:
            raise RuntimeError("Must specify wavelength for transporting beam targets")
        # use substrates
        oLp_set_seq = []
        for ref in waypoints:
            oLp_set_seq.append(self.bg.rAp2oLp_set(ref, obj=obj))
        if len(oLp_set_seq) == 0:
            raise RuntimeError(
                ("Must specify at least one waypoint port" " to bind the target")
            )
        elif len(oLp_set_seq) == 1:
            if len(oLp_set_seq[0]) > 1:
                with self.pbg.preferred():
                    raise RuntimeError(
                        (
                            "If only one waypoint is specified, it must uniquely define a port. "
                            "waypoint '{}' defines ports {}. Be more specific or add more waypoints"
                        ).format(waypoints[0], oLp_set_seq[0])
                    )
            oLp_path = [next(iter(oLp_set_seq[0]))]
        else:
            oLp_path = self._safe_oLp_path(oLp_set_seq)

        if not isinstance(qX, alm.ComplexBeamParam):
            qX = alm.ComplexBeamParam(qX, wavelength_m=self.fs.wavelength_from_Wk(Wk))
        if not isinstance(qY, alm.ComplexBeamParam):
            qY = alm.ComplexBeamParam(qY, wavelength_m=self.fs.wavelength_from_Wk(Wk))

        trans = self._path_transporters(oLp_path, Wk)
        if z_target != z_start:
            Mx = trans.X.z2mat(z_target - z_start)
            My = trans.Y.z2mat(z_target - z_start)
            qX = qX.propagate_matrix(np.linalg.inv(Mx))
            qY = qY.propagate_matrix(np.linalg.inv(My))
        target_oP = self.pbg.referred_vtup(target_name)
        self._targets[target_oP] = Bunch(
            ol=oLp_path[0],
            oLp_path=oLp_path,
            oltrans=trans,
            qX=qX,
            qY=qY,
            Wk=Wk,
            type="specified",
        )
        return

    def cavity_add(self, target_name, waypoints, obj=None):
        # TODO, what about travelling-wave cavities
        # not sure that I love how this works yet..

        oLp_set_seq = []
        for ref in waypoints:
            oLp_set_seq.append(self.bg.rAp2oLp_set(ref, obj=obj))

        link_seq = self._safe_oLp_path(oLp_set_seq, loop=True)

        target_oP = self.pbg.referred_vtup(target_name)
        # print(oP)
        self.cavity_targets[target_oP] = link_seq
        self.cavity_map[target_name] = target_oP
        # self.cavity_params[target_name] = link_seq

        objpath = self.bg.oLp_path2objpath(link_seq, as_refs=True)
        return objpath

    def _cavity_params(self, target_name, ol, Wk, shifts_use=False):
        target_oP = self.pbg.referred_vtup(target_name)
        # print(oP)
        cav_path = self.cavity_targets[target_oP]
        # print("ol", ol)
        idx = cav_path.index(ol)
        # rotate the path to the given target start, and add the final node back
        # to form the full loop
        path = cav_path[idx:] + cav_path[:idx] + [cav_path[idx]]
        # print("PATH: ", path)

        trans = self._path_transporters(path, Wk, shifts_use=shifts_use)

        matX = trans.X.full_trip_mat
        matY = trans.Y.full_trip_mat

        qX = alm.eigen_q(matX)
        qY = alm.eigen_q(matY)

        eye = np.eye(2)
        cav_shiftX = {}
        for shift_key, shift in trans.X.shifts_out_referred.items():
            shiftX = -np.linalg.inv(matX - eye) @ shift
            # print("SHIFTY: ", shift_key, shift, shiftX)
            cav_shiftX[shift_key] = shiftX
        cav_shiftY = {}
        for shift_key, shift in trans.Y.shifts_out_referred.items():
            shiftY = -np.linalg.inv(matY - eye) @ shift
            cav_shiftY[shift_key] = shiftY

        if not np.isfinite(abs(qX)) or not np.isfinite(abs(qY)):
            # TODO, include matrices in error message?
            mat_full_tot = np.eye(2)
            for ol1, ol2 in algo_mm_linkages.path_pairs(trans.Y.oLp_path):
                # print(ol1, ol2)
                idx1, idx2 = trans.Y.inc_ol2idx.get(ol1, None), trans.Y.inc_ol2idx.get(
                    ol2, None
                )
                if idx1 is not None and idx2 is not None:
                    mat_full = np.eye(2)
                    for l, pfunc, mat in trans.Y.inc[idx1:idx2]:
                        mat_full = mat @ mat_full
                        mat_full_tot = mat @ mat_full_tot
                    if np.any(mat_full != np.eye(2)):
                        print(mat_full)
            # print("Total RT Matrix")
            # print(mat_full_tot)
            raise RuntimeError("Cavity {} is not stable".format(target_name))
        qX = alm.ComplexBeamParam(qX, wavelength_m=self.fs.wavelength_map[Wk])
        qY = alm.ComplexBeamParam(qY, wavelength_m=self.fs.wavelength_map[Wk])

        # TODO, annotate target type to relate it to innate targets
        return Bunch(
            qX=qX,
            qY=qY,
            Wk=Wk,
            ol=ol,
            matX=matX,
            matY=matY,
            type="cavity",
            cavity_path=path,
            cavity_trans=trans,
            cav_shiftX = cav_shiftX,
            cav_shiftY = cav_shiftY,
        )

    def _target_get(self, target, Wk=None):
        """
        Get the basic form of target, but don't complete the target computations if it is a cavity
        """
        target_oP = self.pbg.referred_vtup(target)
        cavB = self.cavity_targets.get(target_oP, None)
        if cavB is not None:
            return Bunch(
                oLp_set=cavB,
                target=target,
                Wk=Wk,
                type="cavity",
            )
        else:
            targB = self._targets[target_oP]
            return Bunch(
                oLp_set=[targB.ol],
                target=target,
                Wk=Wk,
                targB=targB,
                type="specified",
            )

    def _target_complete(
            self,
            targB,
            ol : algo_tups.ObjectLinkageTup,
            shifts_use=False
    ):
        """
        Take the cavity object and complete the remaining computations

        ol is an ObjectLinkage tuple
        """
        if targB.type == "cavity":
            return self._cavity_params(targB.target, ol=ol, Wk=targB.Wk, shifts_use=shifts_use)
        elif targB.type == "specified":
            target_oP = self.pbg.referred_vtup(targB.target)
            return self._targets[target_oP]

    def cavity_parameters(self, cavity_name, waypoint, Wk, obj=None, shifts_use=True):
        Wk = self.fs.parameter_to_wk(Wk)
        if isinstance(waypoint, tuple):
            wp = waypoint
        else:
            wp = self.bg.rAp2oLp_set(waypoint, obj=obj, dir="out").pop()

        params = self._cavity_params(cavity_name, wp, Wk=Wk, shifts_use=shifts_use)

        return Bunch(
            cavB=params,
            qX=params.qX,
            qY=params.qY,
            gouyX_deg=np.angle(
                params.qX.propagate_matrix(params.matX).gouy_phasor
                / params.qX.gouy_phasor,
                deg=True,
            ),
            gouyY_deg=np.angle(
                params.qY.propagate_matrix(params.matY).gouy_phasor
                / params.qY.gouy_phasor,
                deg=True,
            ),
        )

    def cavity_digest(self, Wk, waypoint=None, print=print):
        with self.pbg.preferred():
            for target, target_op in sorted(self.cavity_map.items()):
                cavity_seq = self.cavity_targets[target_op]

                def print_wp(cB, wp):
                    print("-------------")
                    print("-----", target, "  at  ", wp)
                    print(
                        "         Gouy X,Y [deg]  {:.2f}, {:.2f}".format(
                            cB.gouyX_deg, cB.gouyY_deg
                        )
                    )
                    print(
                        "         Gouy Frac {:.4f}, {:.4f}".format(
                            (cB.gouyX_deg / 360 + 0.5) % 1 - 0.5,
                            (cB.gouyY_deg / 360 + 0.5) % 1 - 0.5,
                        )
                    )
                    if waypoint:
                        print("         qX ", cB.qX)
                        print("         qY ", cB.qY)
                        print(
                            "         diameter[m] X, Y ",
                            alm.str_m(cB.qX.W * 2),
                            alm.str_m(2 * cB.qY.W),
                        )

                if waypoint is None:
                    wp = cavity_seq[0]
                    cB = self.cavity_parameters(target, wp, Wk)
                    print_wp(cB, wp)
                elif isinstance(waypoint, list):
                    for wp in waypoint:
                        cB = self.cavity_parameters(target, wp, Wk)
                        print_wp(cB, wp)
                else:
                    wp = waypoint
                    cB = self.cavity_parameters(target, wp, Wk)
                    print_wp(cB, wp)
        return

    def overlap(
        self,
        target_fr,
        target_to,
        targets_fr=None,
        targets_to=None,
        waypoints=None,
        Wk=None,
        obj=None,
        shifts_use=False,
        _just_center=False,
    ):
        """
        Create an overlap objects
        """
        Wk = self.fs.parameter_to_wk(Wk)

        if isinstance(waypoints, str):
            waypoints = [waypoints]

        oLp_set_seq = []
        if waypoints is not None:
            for ref in waypoints:
                ref_oP = self.pbg.referred_vtup(ref)
                if ref_oP in self.cavity_targets:
                    oLp_set_seq.append(self.cavity_targets[ref_oP])
                else:
                    # print("WP: ", self.bg.rAp2oLp_set(ref, obj = obj))
                    oLp_set_seq.append(self.bg.rAp2oLp_set(ref, obj=obj))

            if len(oLp_set_seq) > 1:
                oLp_path_center = self._safe_oLp_path(oLp_set_seq)
                ol_imed = oLp_path_center
            else:
                oLp_path_center = []
                ol_imed = oLp_set_seq[0]
        else:
            oLp_path_center = []
            ol_imed = []

        def targets_normalize(target, targets):
            targets_d = {}
            if target is None:
                if targets is not None:
                    target = next(iter(targets))
            elif isinstance(target, str):
                pass
            else:
                for t in target:
                    targets_d[t] = []
                if not target:
                    target = None
                else:
                    target = next(iter(target))

            if targets is not None:
                if not isinstance(targets, collections.abc.Mapping):
                    for t in targets:
                        targets_d[t] = []
                else:
                    targets_d.update(targets)

            if target not in targets_d:
                if target is not None:
                    targets_d[target] = []
            return target, targets_d

        target_fr, targets_fr = targets_normalize(target_fr, targets_fr)
        target_to, targets_to = targets_normalize(target_to, targets_to)

        targetsB_fr = dict()
        targetsB_to = dict()

        # target only approach
        if not oLp_path_center:
            if _just_center:
                assert(len(ol_imed) == 1)
                oLp_path_center = list(ol_imed)
            else:
                if target_fr is None and target_to is None:
                    raise RuntimeError(
                        "Must specify from and to targets if no waypoint path is provided"
                    )
                oLp_set_seq = []

                # these should be included in the code below, rather than special-cased here
                if target_fr is not None:
                    tspecB_fr = self._target_get(target_fr, Wk=Wk)
                    frB = Bunch()
                    targetsB_fr[target_fr] = frB
                    frB.tspecB = tspecB_fr
                    oLp_set_seq.append(tspecB_fr.oLp_set)
                    for ref in targets_fr[target_fr]:
                        oLp_set_seq.append(self.bg.rAp2oLp_set(ref, obj=obj))
                    targets_fr.pop(target_fr)

                if ol_imed:
                    oLp_set_seq.append(ol_imed)

                # these should be included in the code below, rather than special-cased here
                if target_to is not None:
                    tspecB_to = self._target_get(target_to, Wk=Wk)
                    toB = Bunch()
                    targetsB_to[target_to] = toB
                    toB.tspecB = tspecB_to

                    for ref in targets_to[target_to]:
                        oLp_set_seq.append(self.bg.rAp2oLp_set(ref, obj=obj))
                    targets_to.pop(target_to)
                    oLp_set_seq.append(tspecB_to.oLp_set)

                # print("SEQ", oLp_set_seq)
                oLp_path_center = self._safe_oLp_path(oLp_set_seq)

                # these should be included in the code below, rather than special-cased here
                if target_fr is not None:
                    frB.oLp_path = [oLp_path_center[0]]
                    frB.include_center = False
                    frB.inv_start = False
                if target_to is not None:
                    toB.oLp_path = [oLp_path_center[-1]]
                    toB.include_center = True
                    toB.inv_start = True

        for t_fr, wp_fr in targets_fr.items():
            # TODO, currently from targets can only aim to the start of the path
            # not in the middle of it. ol_imed is supposed to find the nearest
            # point of intersection to the waypoint path. Must add additional handling
            # for the mid-path-intersections.
            tB_fr = self._target_get(t_fr, Wk=Wk)
            oLp_set_seq = [tB_fr.oLp_set]
            for ref in wp_fr:
                oLp_set_seq.append(self.bg.rAp2oLp_set(ref, obj=obj))
            # oLp_set_seq.append(ol_imed)
            # Currently, just aim at the start of the waypoint path
            oLp_set_seq.append(oLp_path_center)
            # TODO, not sure if this loop variable is correct
            oLp_path_fr = self._safe_oLp_path(
                oLp_set_seq, loop=False, allow_non_unique=True
            )

            include_center = False
            if oLp_path_fr[-1] != oLp_path_center[0]:
                if len(oLp_path_fr) == 1:
                    mid = oLp_path_fr[0]
                    idx = oLp_path_center.index(mid)
                    oLp_path_fr = oLp_path_center[: idx + 1]
                    inv_start = True
                else:
                    mid = oLp_path_fr[-1]
                    idx = oLp_path_center.index(mid)

                    oLp_path_fr = oLp_path_fr + oLp_path_center[idx + 1:]
                    #raise NotImplementedError(
                    #    "MM Doesn't support middle-injection beams (target {})".format(
                    #        t_fr
                    #    )
                    #)
                    tspecB_fr = self._target_get(t_fr, Wk=Wk)
                    frB = Bunch()
                    frB.tspecB = tspecB_fr
                    frB.oLp_path = oLp_path_fr
                    frB.include_center = False
                    frB.inv_start = False
                    targetsB_to[t_fr] = frB
                    continue
            else:
                inv_start = False

            tspecB_fr = self._target_get(t_fr, Wk=Wk)
            frB = Bunch()
            targetsB_fr[t_fr] = frB
            frB.tspecB = tspecB_fr
            frB.oLp_path = oLp_path_fr
            frB.include_center = include_center
            frB.inv_start = inv_start

        for t_to, wp_to in targets_to.items():
            # also can only handle targets after the waypoint path
            tB_to = self._target_get(t_to, Wk=Wk)
            # oLp_set_seq = [ol_imed]
            # Currently, just aim at the end of the waypoint path
            oLp_set_seq = [oLp_path_center]
            for ref in wp_to:
                oLp_set_seq.append(self.bg.rAp2oLp_set(ref, obj=obj))
            oLp_set_seq.append(tB_to.oLp_set)

            oLp_path_to = self._safe_oLp_path(
                oLp_set_seq, loop=False, allow_non_unique=True
            )

            if oLp_path_to[0] != oLp_path_center[-1]:
                if len(oLp_path_to) == 1:
                    mid = oLp_path_to[-1]
                    idx = oLp_path_center.index(mid)
                    oLp_path_to = oLp_path_center[idx:]
                    inv_start = False
                else:
                    raise NotImplementedError(
                        "MM Doesn't support middle-injection beams (target {})".format(
                            t_to
                        )
                    )
            else:
                inv_start = True

            tspecB_to = self._target_get(t_to, Wk=Wk)
            toB = Bunch()
            targetsB_to[t_to] = toB
            toB.tspecB = tspecB_to
            toB.oLp_path = oLp_path_to
            toB.inv_start = inv_start
            toB.include_center = True

        # now check if paths are branching
        branching = False

        len_longest = -1
        t_longest_fr = None
        for t_fr, frB in targetsB_fr.items():
            if len(frB.oLp_path) > len_longest:
                t_longest_fr = t_fr
                len_longest = len(frB.oLp_path)
        if t_longest_fr:
            oLp_path_longest_fr = targetsB_fr[t_longest_fr].oLp_path
            # now check that the shorter path must be equal to the longer, or it
            # must be a branching path
            for t_fr, frB in targetsB_fr.items():
                if oLp_path_longest_fr[-len(frB.oLp_path) :] != frB.oLp_path:
                    branching = True

        # now check if paths are branching
        len_longest = -1
        t_longest_to = None
        for t_to, toB in targetsB_to.items():
            if len(toB.oLp_path) > len_longest:
                t_longest_to = t_to
                len_longest = len(toB.oLp_path)
        if t_longest_to:
            oLp_path_longest_to = targetsB_to[t_longest_to].oLp_path
            # now check that the shorter path must be equal to the longer, or it
            # must be a branching path
            for t_to, toB in targetsB_to.items():
                if oLp_path_longest_to[: len(toB.oLp_path)] != toB.oLp_path:
                    branching = True

        overlapper = mm_optimize.ModeMatchingOverlapperOptimizing(
            algo_pa=self.pa,
            algo_mm=self,
            targetsB_to=targetsB_to,
            targetsB_fr=targetsB_fr,
            oLp_path_center=oLp_path_center,
            Wk=Wk,
            branching=branching,
            shifts_use=shifts_use,
        )
        overlapper.set_targets(
            target1=target_fr,
            target2=target_to,
        )
        return overlapper


class MMAlgorithmView(object):
    _mm_algo = None
    _obj = None
    _p = None
    _pbg = None

    def __init__(self, obj, mm_algo, p=None, **kw):
        self._obj = obj
        self._mm_algo = mm_algo
        self._pbg = mm_algo.pbg
        if p is None:
            p = self._pbg.view(obj)
        self._p = p
        self.p = p


class MMAlgorithmManipulator(MMAlgorithmView):
    def path(self):
        return self._pbg.path_str(self._obj)

    def parent(self):
        parent, refP = next(iter(self._pbg.object_paths[self._obj]))
        return parent

    def cavity_add(self, name, waypoints):
        return self._mm_algo.cavity_add(name, waypoints, self.parent())

    def target_add(self, name, waypoints, qX, qY, wavelength):
        return self._mm_algo.target_add(
            name, waypoints, qX=qX, qY=qY, obj=self.parent(), wavelength=wavelength
        )

    def parameter_to_wk(self, wparam):
        return self._pa_algo.fs.parameter_to_wk(wparam)
