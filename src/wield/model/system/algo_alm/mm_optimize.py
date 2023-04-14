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
import collections
import itertools
import warnings
from wield.bunch import Bunch
from . import mm_overlapper
from ...optics import alm


class ModeMatchingOverlapperOptimizing(mm_overlapper.ModeMatchingOverlapper):
    def print_optimizable_parameters(self):
        # TODO, use tabulate and have it align on the last slash
        for ref, deps in sorted(self.optimizable_parameters().usable.items()):
            print(ref)
            if deps:
                for dep in deps:
                    print("\t", dep)

    def optimizable_parameters(self):
        comp = self.compile_overlap_calculation(self.pbg)

        return Bunch(
            usable=comp.usable_refs,
            viewed=comp.viewed_refs,
            locked=comp.locked_refs,
        )

    def optimize(
        self,
        opt_params,
        length_constraints=[],
        nfev=None,
        warn_locked=True,
    ):
        """
        The length constraints is a list of lists with [dist, wp1, wp2,...]
        """
        pbg = self.pbg.copy()

        comp = self.compile_overlap_calculation(pbg)

        ########
        #  Work on optimizable parameters and their bounds constraints
        opt_params = set(opt_params)
        # TODO, better error handling/reporting
        assert opt_params <= set(comp.usable_refs.keys())
        # check that none of the supplied parameters are locked out by the optimizer
        assert not (opt_params & set(comp.locked_refs.keys()))

        # triggers to recompute the q-values
        if opt_params & set(comp.t1_refs.keys()):
            t1_recompute = True
        else:
            t1_recompute = False

        if opt_params & set(comp.t2_refs.keys()):
            t2_recompute = True
        else:
            t2_recompute = False

        req_completion = set()
        opt_params_op = set()
        ref_op_map = dict()

        def op_completion(op):
            req_completion.add(op)
            for op_req in pbg.op_reqs[op]:
                if op_req in req_completion:
                    continue
                req_completion.add(op_req)
                op_completion(op_req)

        for ref in opt_params:
            op = pbg.resolve_reference(ref)
            opt_params_op.add(op)
            ref_op_map[ref] = op
            op_completion(op)

        bounds_op_lbub = dict()
        for op in req_completion:
            obj, vtup = op
            lb = pbg.resolve_parameter_default(obj, vtup + ("lower_bound",), None)
            ub = pbg.resolve_parameter_default(obj, vtup + ("upper_bound",), None)
            if lb is not None or ub is not None:
                bounds_op_lbub[op] = (lb, ub)

        ZpropB = build_length_constraints(
            pbg=pbg,
            bg=self.pa.bg,
            trans_center=self.trans_center,
            length_constraints=length_constraints,
        )
        mats = compile_length_constraints(
            ZpropB=ZpropB,
            opt_params_op=opt_params_op,
            warn_locked=warn_locked,
        )
        ##Zprop_deps.update(ops_)
        ##Perform a completion on the dependencies
        deps_all = set()
        for depset in ZpropB.deps:
            deps_all.update(depset)

        Zprop_deptree = pbg.opset2ref_deptree(deps_all)
        Zprop_deps = set(Zprop_deptree.keys())

        # now work out length constraints matrix
        # TODO, make this a separate function
        dep_params = opt_params_op & Zprop_deps
        if dep_params:
            raise RuntimeError(
                (
                    "Length constraints cannot be computed while"
                    " these parameters are being optimized {}"
                ).format(dep_params)
            )

        # generate the ops vector for fitting, removing the linear constraints
        ops_fit = list()
        x0 = []
        for ref in sorted(opt_params):
            op = ref_op_map[ref]
            if mats is not None and op in mats.ops_constrain:
                continue
            ops_fit.append(op)
            val = pbg._resolve_parameter(*op)
            x0.append(val)

        first_value = 0
        lowest_value = 0

        def pbg_set_x(x):
            for op, val in zip(ops_fit, x):
                pbg.override_value_op(op, val)
                # val2 = pbg._resolve_parameter(*op)

            if mats is not None:
                # include the linear constraints
                indep_vec = list()
                for op in mats.ops_indep:
                    val = pbg._resolve_parameter(*op)
                    indep_vec.append(val)
                indep_vec = np.asarray(indep_vec)
                dep_vec = mats.mat_transfer @ indep_vec
                for len_m, op, val in zip(
                    mats.lengths_inv, mats.ops_constrain, dep_vec
                ):
                    pbg.override_value_op(op, len_m - val)
            return

        def f_eval(x):
            nonlocal lowest_value
            pbg_set_x(x)
            bounds_diff_max = 0
            # now check that bounds have not been violated
            for op, (lb, ub) in bounds_op_lbub.items():
                val = pbg._resolve_parameter(*op)
                if lb is not None and lb > val:
                    bdiff = lb - val
                    if bdiff > bounds_diff_max:
                        bounds_diff_max = bdiff
                if ub is not None and ub < val:
                    bdiff = val - ub
                    if bdiff > bounds_diff_max:
                        bounds_diff_max = bdiff

            # maybe a smarter policy exists for this...
            if bounds_diff_max > 0:
                return lowest_value + 100 * bounds_diff_max

            overlap = comp.calculate_overlap(
                pbg,
                t1_recompute=t1_recompute,
                t2_recompute=t2_recompute,
            )
            if overlap is None:
                return first_value

            if -overlap < lowest_value:
                lowest_value = -overlap

            return -overlap

        first_value = f_eval(x0)
        lowest_value = first_value

        import scipy.optimize

        out = scipy.optimize.minimize(
            fun=f_eval,
            x0=x0,
            method="Powell",
            options=dict(
                maxfev=nfev,
                xtol=1e-14,
                ftol=1e-14,
            ),
        )
        # out = scipy.optimize.minimize(
        #    fun = f_eval,
        #    x0 = x0,
        #    method = 'Nelder-Mead',
        #    options = dict(
        #        maxfev = nfev,
        #        xatol = 1e-14,
        #        fatol = 1e-14,
        #        adaptive = True,
        #    ),
        # )
        # print(out)
        # pbg.override_value_op(op, 2)
        # print("OLAP", lowest_value)
        # print("OLAP", out)

        if len(x0) == 1:
            # needs special case since optimize.minimize seems to unwrap the
            # parameter list when it is length 1
            pbg_set_x([out.x])
        else:
            pbg_set_x(out.x)

        pa = self.pa.regenerate(pbg=pbg)

        overlapper = self.__class__(
            algo_pa=pa,
            algo_mm=pa.mm,
            targetsB_to=self.targetsB_to,
            targetsB_fr=self.targetsB_fr,
            oLp_path_center=self.oLp_path_center,
            Wk=self.Wk,
            branching=self.branching,
        )
        overlapper.set_targets(
            target1=self.target1,
            target2=self.target2,
        )
        return overlapper

    def compile_overlap_calculation(
            self,
            pbg=None,
            loc=None,
    ):
        # TODO, also allow 'from' or 'fr'
        assert(loc in ['fr', 'to', None])
        if pbg is None:
            pbg = self.pbg
        t1B = self[self.target1]
        t2B = self[self.target2]

        comp = Bunch()
        # standardize direction
        if t1B.type == "to" and t2B.type == "from":
            t1B, t2B = t2B, t1B

        # the target q-values are the original ones, not the normalized ones
        # stored directly in the bunches, but inside the inner targB bunches
        comp.t1qX = t1B.targB.qX
        comp.t1qY = t1B.targB.qY
        comp.t2qX = t2B.targB.qX
        comp.t2qY = t2B.targB.qY
        if t1B.type == "from" and t2B.type == "to":
            if t1B.inv_start or not t2B.inv_start:
                raise NotImplementedError("mid-path modes not yet supported")
            comp.type = "ft"
            # this one is always in loc == 'to'
            if loc == 'fr':
                comp.t1pX = t1B.trans.X.prop
                comp.t1pY = t1B.trans.Y.prop
                comp.t2pX = self.trans_center.X.prop + t2B.trans.X.prop
                comp.t2pY = self.trans_center.Y.prop + t2B.trans.Y.prop
            else:
                comp.t1pX = t1B.trans.X.prop + self.trans_center.X.prop
                comp.t1pY = t1B.trans.Y.prop + self.trans_center.Y.prop
                comp.t2pX = t2B.trans.X.prop
                comp.t2pY = t2B.trans.Y.prop
        elif t1B.type == "from" and t2B.type == "from":
            if t1B.inv_start or t2B.inv_start:
                raise NotImplementedError("mid-path modes not yet supported")
            comp.type = "ff"
            if loc == 'to':
                comp.t1pX = t1B.trans.X.prop + self.trans_center.X.prop
                comp.t1pY = t1B.trans.Y.prop + self.trans_center.Y.prop
                comp.t2pX = t2B.trans.X.prop + self.trans_center.X.prop
                comp.t2pY = t2B.trans.Y.prop + self.trans_center.Y.prop
            else:
                comp.t1pX = t1B.trans.X.prop
                comp.t1pY = t1B.trans.Y.prop
                comp.t2pX = t2B.trans.X.prop
                comp.t2pY = t2B.trans.Y.prop
        elif t1B.type == "to" and t2B.type == "to":
            if not t1B.inv_start or not t2B.inv_start:
                raise NotImplementedError("mid-path modes not yet supported")
            comp.type = "tt"
            if loc == 'fr':
                comp.t1pX = self.trans_center.X.prop + t1B.trans.X.prop
                comp.t1pY = self.trans_center.Y.prop + t1B.trans.Y.prop
                comp.t2pX = self.trans_center.X.prop + t2B.trans.X.prop
                comp.t2pY = self.trans_center.Y.prop + t2B.trans.Y.prop
            else:
                comp.t1pX = t1B.trans.X.prop
                comp.t1pY = t1B.trans.Y.prop
                comp.t2pX = t2B.trans.X.prop
                comp.t2pY = t2B.trans.Y.prop
        else:
            raise RuntimeError("Unrecognized directions (bug)")

        # qX.propagate_matrix(np.linalg.inv(matXto))
        # qY.propagate_matrix(np.linalg.inv(matYto))
        path_set = set()
        locked_set = set()
        t1_set = set()
        t2_set = set()

        def gen_mat(prop, vset):
            mat_inc = np.eye(2)
            for (obj, prop) in prop:
                p = pbg.view(obj, vset)
                mat = prop(p)
                mat_inc = mat @ mat_inc
            return mat_inc

        mat1X = gen_mat(comp.t1pX, path_set)
        mat1Y = gen_mat(comp.t1pY, path_set)
        mat2X = gen_mat(comp.t2pX, path_set)
        mat2Y = gen_mat(comp.t2pY, path_set)
        mat1X, mat1Y, mat2X, mat2Y  # refer to them just to get lint to stop complaining

        def lock_paths(tB, viewset):
            if tB.targB.type == "cavity":
                cavity_trans = tB.targB.cavity_trans
                gen_mat(cavity_trans.X.prop, viewset)
                gen_mat(cavity_trans.Y.prop, viewset)
            elif tB.targB.type == "specified":
                pass
            else:
                raise RuntimeError("Unrecognized target type (bug)")

        lock_paths(t1B, t1_set)
        lock_paths(t2B, t2_set)
        # print("L", locked_set)
        comp.used_ops = path_set
        comp.locked_ops = locked_set
        path_refs = pbg.opset2ref_deptree(path_set)
        t1_refs = pbg.opset2ref_deptree(t1_set)
        t2_refs = pbg.opset2ref_deptree(t2_set)
        locked_refs = pbg.opset2ref_deptree(locked_set)
        viewed_refs = dict(path_refs)
        viewed_refs.update(t1_refs)
        viewed_refs.update(t2_refs)

        usable_refs = dict()
        # removed the locked values from the usable values
        for k, v in viewed_refs.items():
            if k not in locked_refs:
                usable_refs[k] = v
        comp.viewed_refs = viewed_refs
        comp.path_refs = path_refs
        comp.t1_refs = t1_refs
        comp.t2_refs = t2_refs
        comp.locked_refs = locked_refs
        comp.usable_refs = usable_refs

        wavelength_m = self.mm.fs.wavelength_map[self.Wk]

        def calculate_Qs(
            pbg=pbg,
            t1_recompute=True,
            t2_recompute=True,
        ):
            def gen_mat(prop):
                mat_inc = np.eye(2)
                for (obj, propagator) in prop:
                    p = pbg.view(obj)
                    mat = propagator(p)
                    mat_inc = mat @ mat_inc
                return mat_inc

            mat1X = gen_mat(comp.t1pX)
            mat1Y = gen_mat(comp.t1pY)
            mat2X = gen_mat(comp.t2pX)
            mat2Y = gen_mat(comp.t2pY)
            bad = False

            if not t1_recompute or (t1B.targB.type != "cavity"):
                pre_t1qX = comp.t1qX
                pre_t1qY = comp.t1qY
            else:
                q1Xmat = gen_mat(t1B.targB.cavity_trans.X.prop)
                q1Ymat = gen_mat(t1B.targB.cavity_trans.Y.prop)
                qX = alm.eigen_q(q1Xmat)
                qY = alm.eigen_q(q1Ymat)
                if np.any(~np.isfinite(abs(qX))) or np.any(~np.isfinite(abs(qY))):
                    bad = True
                pre_t1qX = alm.ComplexBeamParam(qX, wavelength_m=wavelength_m)
                pre_t1qY = alm.ComplexBeamParam(qY, wavelength_m=wavelength_m)

            if not t2_recompute or (t2B.targB.type != "cavity"):
                pre_t2qX = comp.t2qX
                pre_t2qY = comp.t2qY
            else:
                q2Xmat = gen_mat(t2B.targB.cavity_trans.X.prop)
                q2Ymat = gen_mat(t2B.targB.cavity_trans.Y.prop)
                qX = alm.eigen_q(q2Xmat)
                qY = alm.eigen_q(q2Ymat)
                if np.any(~np.isfinite(abs(qX))) or np.any(~np.isfinite(abs(qY))):
                    bad = True
                pre_t2qX = alm.ComplexBeamParam(qX, wavelength_m=wavelength_m)
                pre_t2qY = alm.ComplexBeamParam(qY, wavelength_m=wavelength_m)

            if comp.type == "ft":
                t1qX = pre_t1qX.propagate_matrix(mat1X)
                t1qY = pre_t1qY.propagate_matrix(mat1Y)
                t2qX = pre_t2qX.propagate_matrix(np.linalg.inv(mat2X))
                t2qY = pre_t2qY.propagate_matrix(np.linalg.inv(mat2Y))
            elif comp.type == "ff":
                t1qX = pre_t1qX.propagate_matrix(mat1X)
                t1qY = pre_t1qY.propagate_matrix(mat1Y)
                t2qX = pre_t2qX.propagate_matrix(mat2X)
                t2qY = pre_t2qY.propagate_matrix(mat2Y)
            elif comp.type == "tt":
                t1qX = pre_t1qX.propagate_matrix(np.linalg.inv(mat1X))
                t1qY = pre_t1qY.propagate_matrix(np.linalg.inv(mat1Y))
                t2qX = pre_t2qX.propagate_matrix(np.linalg.inv(mat2X))
                t2qY = pre_t2qY.propagate_matrix(np.linalg.inv(mat2Y))
            return Bunch(
                t1qX=t1qX,
                t2qX=t2qX,
                t1qY=t1qY,
                t2qY=t2qY,
                bad=bad,
            )

        def calculate_overlap(
            pbg=pbg,
            t1_recompute=True,
            t2_recompute=True,
        ):
            qB = calculate_Qs(
                pbg=pbg,
                t1_recompute=t1_recompute,
                t2_recompute=t2_recompute,
            )
            overlap = (
                abs(qB.t1qX.overlap_HG(qB.t2qX) * qB.t1qY.overlap_HG(qB.t2qY)) ** 2
            )
            if qB.bad:
                return None
            return overlap

        comp.calculate_overlap = calculate_overlap
        comp.calculate_Qs = calculate_Qs
        overlapA = calculate_overlap(t1_recompute=False, t2_recompute=False)
        # overlapB = calculate_overlap(t1_recompute = True, t2_recompute = True)
        # assert(abs(overlapA - overlapB) < 1e-8)
        comp.overlap = overlapA
        return comp


def build_length_constraints(
    pbg,
    bg,
    trans_center,
    length_constraints,
):
    ##########
    #  Length constraints
    Zprop_deps = []
    Zprop_lengths = []
    Zprop_constraints = []
    for lc in length_constraints:
        # print('length_constraints', lc)
        length_m = lc[0]
        # from icecream import ic; ic(self.trans_center.Zprop_ol2idx)
        path_idx = []
        for wp in lc[1:]:
            olset = bg.rAp2oLp_set(wp)
            # print(wp, olset)
            idxs = set()
            for ol in olset:
                idx = trans_center.Zprop_ol2idx.get(ol, None)
                if idx is not None:
                    idxs.add(idx)
            if len(idxs) == 0:
                raise RuntimeError(
                    "Length constraint not within matching path for waypoint '{}'".format(
                        wp
                    )
                )
            if len(idxs) > 1:
                raise RuntimeError(
                    "Length constraint not unique for waypoint '{}' (be more specific)".format(
                        wp
                    )
                )
            idx = next(iter(idxs))
            path_idx.append(idx)
        Zprop_scales = dict()
        for Zscale in trans_center.Zprop_scales[path_idx[0] : path_idx[-1]]:
            for op, scale in Zscale.items():
                if op in Zprop_scales:
                    Zprop_scales[op] += scale
                else:
                    Zprop_scales[op] = scale
        deps_collect = set()
        for Zdeps in trans_center.Zprop_deps[path_idx[0] : path_idx[-1]]:
            for op, deps in Zdeps.items():
                deps_collect.update(deps)
        Zprop_deps.append(deps_collect)
        Zprop_constraints.append(Zprop_scales)
        Zprop_lengths.append(length_m)

    def check(pbg):
        for len_m, cmap in zip(Zprop_lengths, Zprop_constraints):
            len_m_sum = 0
            for op, scale in cmap.items():
                len_m_sum += pbg._resolve_parameter(*op) * scale
            # print(len_m, len_m_sum)

    return Bunch(
        lengths=Zprop_lengths,
        deps=Zprop_deps,
        constraints=Zprop_constraints,
        check=check,
    )


def compile_length_constraints(
    ZpropB,
    opt_params_op,
    warn_locked,
):
    N_constraints = len(ZpropB.constraints)
    if N_constraints == 0:
        return None
    op_constr = collections.defaultdict(lambda: 0)
    op_len_constr = []
    ops_full_constr = set()
    for cmap in ZpropB.constraints:
        cmap_k = set(cmap.keys())
        ops_full_constr.update(cmap_k)
        op_len_constr.append(cmap_k & opt_params_op)
        if warn_locked and len(op_len_constr[-1]) == 1:
            warnings.warn(
                (
                    "Parameter {} locked by length constraint,"
                    "will require additional parameters to free the lock"
                ).format(next(iter(op_len_constr[-1])))
            )
        for op in opt_params_op:
            if op in cmap:
                op_constr[op] += 1
    ops_possible_constr = set()
    for op, N in op_constr.items():
        if N > 0:
            ops_possible_constr.add(op)

    # must choose orthogonal parameters within the constraints
    # this is technically this problem
    # https://doi.org/10.1016/j.tcs.2009.06.018
    # and this Q/R decomposition is related to the "greedy" algorithm
    # TODO, actually add in the bits to choose largest remainders
    ops_possible_constr = list(ops_possible_constr)
    mat_possible = {}
    for op in ops_possible_constr:
        vec_possible = []
        for cmap in ZpropB.constraints:
            vec_possible.append(cmap.get(op, 0))
        mat_possible[op] = np.asarray(vec_possible)
    ops_chosen = set()

    try:
        for N in range(N_constraints):
            op, vec = mat_possible.popitem()
            ops_chosen.add(op)
            vec_sqsum = np.dot(vec, vec)
            newmat = dict()
            for op2, vec2 in mat_possible.items():
                vec2 = vec2 - np.dot(vec2, vec) * vec / vec_sqsum
                vec2_sqsum = np.dot(vec2, vec2)
                if vec2_sqsum < 1e-10:
                    # don't add it, since it is too small
                    continue
                newmat[op2] = vec2

            mat_possible = newmat
    except KeyError:
        # this only throws if mat_possible has run out of orthogonal variables
        raise RuntimeError(
            (
                "No collection of parameters can satisfy all length constraints."
                "They are degenerate."
            )
        )

    def try_chosen_constraints(ops_constrain):
        ops_constrain = list(ops_constrain)
        ops_indep = list(ops_full_constr - set(ops_constrain))
        mat_constrain = []
        mat_indep = []
        for cmap in ZpropB.constraints:
            vec_constr = []
            for op in ops_constrain:
                vec_constr.append(cmap.get(op, 0))
            mat_constrain.append(vec_constr)
            vec_constr = []
            for op in ops_indep:
                vec_constr.append(cmap.get(op, 0))
            mat_indep.append(vec_constr)

        mat_constrain = np.asarray(mat_constrain)
        mat_indep = np.asarray(mat_indep)

        try:
            mat_indep_inv = np.linalg.inv(mat_constrain)
            mat_transfer = mat_indep_inv @ mat_indep
            lengths_inv = mat_indep_inv @ np.asarray(ZpropB.lengths)
        except np.linalg.linalg.LinAlgError:
            return None
        return Bunch(
            ops_constrain=ops_constrain,
            ops_indep=ops_indep,
            mat_constrain=mat_constrain,
            mat_indep=mat_indep,
            mat_indep_inv=mat_indep_inv,
            mat_transfer=mat_transfer,
            lengths_inv=lengths_inv,
        )

    mats = try_chosen_constraints(ops_chosen)

    if mats is None:
        # relies on combinations emitting in lexicographic sort order to
        # have highest success likelihood on first iteration.
        # removes the ops with the most dependencies first
        # TODO, this should never run since the algorithm above finds a
        # non-degenerate set
        for ops_constrain in itertools.combinations(ops_possible_constr, N_constraints):
            mats = try_chosen_constraints(ops_constrain)
            if mats is not None:
                break
        else:
            # this segment only called if all matricies have linalg error
            raise RuntimeError(
                (
                    "No collection of parameters can satisfy all length constraints."
                    "They are likely degenerate."
                )
            )

    # from icecream import ic; ic(mat_transfer)
    # ic(ops_indep)
    # ic(ops_constrain)
    # ic(opt_params_op)
    # ic(opt_params)
    # ic(bounds_op_lbub)
    return mats
