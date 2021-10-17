#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""
import collections
import numpy as np
from .. import base
from . import algo_phys


class PhysicsDCAlgorithm(object):
    N_tol_rel = 1e-10
    N_tol_abs_thresh = 1e-12

    def __init__(self, pa):
        self.pa = pa
        self.log = pa.log
        self.pbg = pa.pbg
        #pg.print_parameters_eval()
        self.bg = pa.bg
        self.fs = pa.fs
        self.check_build = pa.check_build

        #(obj, lport) -> kvector
        self._sources = dict()
        #obj-> (lportRow, lportCol) -> kmatrix
        self._object_edges = collections.defaultdict(dict)

        #indicates which outputs to monitor for sensitivity
        #views = {(self._obj, lport)}
        self._views = set()

        #mapping from (obj,link) -> (link_fr, link_to)
        #where the link_fr/link_to have the obj assumed
        self._nl_deps = collections.defaultdict(set)

        #indicates which outputs to monitor for sensitivity
        #views = {(self._obj, lport) : kmatrix}
        self._noise = dict()

        self._solutions_DC_prev = dict()
        self._solutions_DC      = dict()

        solve_again = True
        while solve_again:
            self.log.info("solving DC")
            self._build_matrix_DC()
            self._solve_matrix_DC()
            solve_again = self._check_matrix_DC_NL()
            self._solutions_DC_prev = self._solutions_DC

        return

    def _optical_frequency_allowed(self, fk):
        return fk in self.fs.freq_set_optical

    def __call__(
        self,
        port_to,
        obj_to = None,
        dir_to = 'out',
        demod = {},
    ):
        """
        Gets the power or current at a photodiode.

        units, defaults to 'W' for total power, but may be
        'A', or 'Q' for quanta-based counts. For multiple-wavelengths these
        can look inconsistent.

        demod takes a frequency dictionary.
        """
        oLp_to = self.bg.rBp2oLp(port_to, dir = dir_to)
        demod = self.fs.parameter_to_fk(demod)

        solvect = self._solutions_DC[oLp_to]

        datavec = solvect.kmatrix[(demod, '+')][()][..., 0, 0]
        #from icecream import ic; ic(datavec)
        #TODO, check realness on DC-demod
        if demod.DC_is():
            return datavec.real
        else:
            return datavec

    def PDsum(
        self,
        ref,
        obj                 = None,
        demod               = {},
        units               = 'W',
        raw_vector          = False,
        itemize_frequencies = False,
        itemize_wavelengths = False,
        itemize_HOMs        = False,
    ):
        """
        Gets the power or current at a photodiode.

        units, defaults to 'W' for total power, but may be
        'A', or 'Q' for quanta-based counts. For multiple-wavelengths these
        can look inconsistent.

        demod takes a frequency dictionary.
        """
        if obj is None:
            obj = self.pbg.root
        obj = self.pbg.referred_path(obj, ref)

        demod = self.fs.parameter_to_fk(demod)

        try:
            visit_method = obj.visit_photodiode_linkage
        except AttributeError:
            raise RuntimeError((
                "Object or path not a photodiode type, it is instead a {}"
            ).format(repr(obj)))

        view = PhysicsAlgorithmDCView(obj = obj, dc_algo = self)
        #simple method just reports the port
        lport = visit_method(view)

        solvect = self._solutions_DC[(obj, lport)]

        if raw_vector:
            return solvect

        if demod.DC_is():
            inner_prod = base.KeyMatrixSame(
                st = (base.kg_twophoton_pm, ),
                kmatrix = {('-',) : {('+', ) : [[1]]}},
                build = True,
            )
            datavec = (solvect.T @ inner_prod @ solvect).kmatrix[()][()][..., 0, 0]
            #TODO, check realness
            return datavec.real
        else:
            #vdict carries the indexs for the columns (from)
            km = dict()
            for fd in self.pa.optical_basis_DC['frequency'].enumerated:
                vdict = km.setdefault((fd, '-'), dict())
                #the lower sideband moving up
                fd_from = fd + demod
                if self._optical_frequency_allowed(fd_from):
                    vdict[(fd_from, '+')] = [[1]]

                inner_prod = matrix.KeyMatrixSame(
                    st      = (
                        self.pa.optical_basis_DC['frequency'],
                        base.kg_twophoton_pm,
                    ),
                    kmatrix = km,
                    build   = True,
                    check   = self.check_build,
                )

            datavec = (solvect.T @ inner_prod @ solvect).kmatrix[()][()][..., 0, 0]
            return datavec

    def _check_matrix_DC_NL(self):
        must_iterate = False
        for sK, sV in self._solutions_DC.items():
            sVp = self._solutions_DC_prev.get(sK, 0)
            if sVp == 0:
                if sV != 0:
                    must_iterate = True
                    continue
                else:
                    continue
            sVdiff = sVp - sV
            diff_pow = (sVdiff.A @ sVdiff).kmatrix[()][()][..., 0, 0]
            sVpow = (sV.A @ sV).kmatrix[()][()][..., 0, 0]
            sVPpow = (sVp.A @ sVp).kmatrix[()][()][..., 0, 0]

            relerr = diff_pow * (sVpow**-1 + sVPpow**-1)
            #from icecream import ic
            #if sV != 0:
            #    ic(sK, relerr)
            if np.any(relerr > self.N_tol_rel):
                must_iterate = True
        return must_iterate

    def _solve_matrix_DC(self):
        SREIO = self.SREIO_DC(
            map_nodes = True,
            subtract_1 = True,
        )
        (seq, req, edges, inputs, outputs) = SREIO
        del SREIO
        with self.log.heading('DC_inversion'):
            seq, req, edges = matrix.SREkmatrix_inverse(
                seq, req, edges,
                outputs_set = outputs,
                inputs_set  = inputs,
                verbose     = False,
                log         = self.log,
            )
        solutions_DC = dict()
        for output in outputs:
            val = None
            op_to = output
            for n_fr in req[output]:
                op_fr = n_fr
                edge = edges[n_fr, output]
                source = self._sources[op_fr]
                if val is None:
                    val = edge @ source
                else:
                    val = val + (edge @ source)
            if val is None:
                solutions_DC[op_to] = 0
            else:
                solutions_DC[op_to] = val

        self._solutions_DC = solutions_DC

    def _build_matrix_DC(self):
        for obj in self.pbg.object_iter():
            try:
                visit_algo = obj.visit_matrix_algorithm_DCAC
            except AttributeError:
                continue
            else:
                #TODO verbose option for found objects?
                #print(obj)
                pass

            manip = PhysicsAlgorithmDCManipulator(
                obj     = obj,
                dc_algo = self,
            )

            visit_algo(manip)
        return

    def SREIO_DC(
            self,
            map_nodes = None,
            subtract_1 = False
    ):
        seq = collections.defaultdict(set)
        req = collections.defaultdict(set)
        edges = dict()

        for oLp_fr, eset in self.bg.link_seq.items():
            m_fr = (oLp_fr)
            #print("FT: ", m_fr, [map_nodes(oLp_to) for oLp_to in eset])
            for oLp_to in eset:
                m_to = (oLp_to)
                edges[m_fr, m_to] = 1
                seq[m_fr].add(m_to)
                req[m_to].add(m_fr)

        for obj, edict in self._object_edges.items():
            for (l_fr, l_to), e_kmat in edict.items():
                m_to = map_nodes((obj, l_to))
                m_fr = map_nodes((obj, l_fr))
                edges[m_fr, m_to] = e_kmat
                seq[m_fr].add(m_to)
                req[m_to].add(m_fr)

        if subtract_1:
            for e_key, e_val in edges.items():
                edges[e_key] = -e_val

            nodes = set(seq.keys())
            nodes.update(req.keys())
            for node in nodes:
                seq[node].add(node)
                req[node].add(node)
                kdiag = edges.get((node, node), 0)
                kdiag = 1 + kdiag
                edges[node, node] = kdiag

        inputs = set()
        for k_in in self._sources.keys():
            inputs.add(map_nodes(k_in))

        outputs = set()
        for k_in in self._views:
            outputs.add(map_nodes(k_in))
        for k_in in self._nl_deps:
            outputs.add(map_nodes(k_in))
        return (seq, req, edges, inputs, outputs)


class PhysicsAlgorithmDCView(algo_phys.PhysicsAlgorithmView):
    _dc_algo = None
    def __init__(self, dc_algo, **kw):
        super(PhysicsAlgorithmDCView, self).__init__(
            bg_algo = dc_algo.bg,
            pbg     = dc_algo.pbg,
            pa_algo = dc_algo.pa,
            **kw
        )
        self._dc_algo = dc_algo

    def configure_optical_wavenumber(self, wdict):
        #TODO, check that the wavenumber exists
        #assert(len(self._pa_algo.fs.freq_set_wavelengths) == 1)
        return base.FrequencyKey(wdict)

    def parameter_to_fk(self, fparam):
        return self._dc_algo.fs.parameter_to_fk(fparam)

    def optical_frequency_allowed(self, fk):
        return self._dc_algo._optical_frequency_allowed(fk)

    def basis_frequencies(self, with_keys = False):
        #print("FVs: ", self._pa_algo._basis_frequency_values_DC_optical)
        if with_keys:
            return zip(
                self._pa_algo.optical_basis_DC['frequency'].enumerated,
                self._pa_algo._basis_frequency_values_DC_optical
            )
        else:
            return self._pa_algo._basis_frequency_values_DC_optical

    def basis_WFQ_optical_pm(self):
        """
        This method iterates over all wavenumber, frequencies and quantum,
        with conjugation or negative-frequency shifting for negative quantum

        usage:
        for (Wk, Fk, Qk), (wnval, fval, conj) in manip.basis_WFQ_pm():

        while this method doesn't do a ton for DC analysis. It is overloaded
        for AC analysis to make things simpler
        """
        for Wk, wnval in self.basis_wavenumbers(with_keys = True):
            for Fk, fval in self.basis_frequencies(with_keys = True):
                for Qk, conj in (('+', False), ('-', True)):
                    yield (Wk, Fk, Qk), (wnval, fval, conj)
        return

    #is a dictionary of known basis elements as of this linkage
    def link_basis(self, lport):
        op = (self._obj, lport)
        btype = self._bg_algo.link_basis_types[op]
        return {
            'optical' :    self._pa_algo.optical_basis_DC,
            'mechanical' : self._pa_algo.mechanical_basis_DC,
            'signal' :     self._pa_algo.signal_basis_DC,
            'electrical' : self._pa_algo.electrical_basis_DC,
        }[btype]


class PhysicsAlgorithmDCManipulator(PhysicsAlgorithmDCView):
    is_DC = True
    is_AC = False

    def add_source(self, lport, svect):
        basis = set(self.link_basis(lport).values())
        inj = set(svect.stR + svect.dtR)
        assert(basis == inj)
        #TODO, check vector consistency
        self._dc_algo._sources[self._obj, lport] = svect

    def add_view(self, lport):
        self._dc_algo._views.add((self._obj, lport))

    def add_view_DC(self, lport):
        self._dc_algo._views.add((self._obj, lport))

    def add_view_AC(self, lport):
        """
        Adds views only necessary for AC analysis
        """
        self._dc_algo._views.add((self._obj, lport))

    def get_field_DC(self, lport):
        return self._dc_algo._solutions_DC.get((self._obj, lport), 0)

    def add_nl_link(self, lport_fr, lport_to, lport_deps, update_func):
        self._dc_algo._nl_links[self._obj][lport_fr, lport_to] = update_func
        if not isinstance(lport_deps, (list, tuple)):
            lport_deps = [lport_deps]
        for lpdep in lport_deps:
            self._dc_algo._nl_deps[(self._obj, lpdep)].add((lport_fr, lport_to))

    def add_link(self, lport_fr, lport_to, kmatrix, lowering_only = False):
        """
        Adds a link to the system matrix

        The lowering_only indicates that only the lowering (non conjugated)
        operator is supplied, and so needs some additional completion
        """
        if lowering_only:
            km_new = dict()
            for krow, kvect in kmatrix.kmatrix.items():
                kvect_new_p = km_new.setdefault(krow + ('+',), dict())
                kvect_new_n = km_new.setdefault(krow + ('-',), dict())
                for kcol, kdm in kvect.items():
                    kvect_new_p[kcol + ('+', )] = kdm
                    kvect_new_n[kcol + ('-', )] = kdm.conjugate()
            #TODO, not always using the optical basis for this
            kmatrix = kmatrix.__class__(
                stR = kmatrix.stR + (self._pa_algo.optical_basis_DC['quantum'],),
                stC = kmatrix.stC + (self._pa_algo.optical_basis_DC['quantum'],),
                dtR = kmatrix.dtR,
                dtC = kmatrix.dtC,
                kmatrix = km_new,
                build = False,
                check = self.check_build,
            )

        if isinstance(kmatrix, base.KeyMatrixBase):
            basisR = set(self.link_basis(lport_to).values())
            assert(set(kmatrix.stR + kmatrix.dtR) <= basisR)
            basisC = set(self.link_basis(lport_fr).values())
            assert(set(kmatrix.stC + kmatrix.dtC) <= basisC)
            self._dc_algo._object_edges[self._obj][lport_fr, lport_to] = kmatrix
        else:
            #is probably just a number or array
            self._dc_algo._object_edges[self._obj][lport_fr, lport_to] = kmatrix
