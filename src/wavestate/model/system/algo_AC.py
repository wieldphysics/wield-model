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
from .. import base
from . import algo_phys


class PhysicsACAlgorithm(object):
    def __init__(self, pa, dc):
        self.pa = pa
        self.log = pa.log
        self.pbg = pa.pbg
        self.dc = dc
        # pg.print_parameters_eval()
        self.bg = pa.bg
        self.fs = pa.fs
        self.check_build = pa.check_build

        # {(obj, lport)}
        self._drives = set()
        # obj-> (lportRow, lportCol) -> kmatrix
        self._object_edges = collections.defaultdict(dict)

        # indicates which outputs to monitor for sensitivity
        # views = {(self._obj, lport)}
        self._views = set()

        # indicates which outputs to monitor for sensitivity
        # views = {(self._obj, lport) : kmatrix}
        self._noise = dict()

        self._build_matrix_AC()
        self._solve_matrix_AC()
        return

    def _parameter_to_fdict(self, fparam):
        if isinstance(fparam, str):
            # TODO, check that the key exists
            return base.FrequencyKey({fparam: 1})
        elif isinstance(fparam, collections.Mapping):
            return base.FrequencyKey(fparam)

    def _optical_frequency_allowed(self, fk):
        return fk in self.fs.freq_set_optical

    def __call__(
        self,
        port_fr,
        port_to,
        obj_fr=None,
        obj_to=None,
        dir_fr="in",
        dir_to="out",
        demod={},
        quadrature=1,
    ):
        """
        Gets the power or current at a photodiode.

        units, defaults to 'W' for total power, but may be
        'A', or 'Q' for quanta-based counts. For multiple-wavelengths these
        can look inconsistent.

        demod takes a frequency dictionary.
        """
        oLp_fr = self.bg.rBp2oLp(port_fr, dir=dir_fr)
        oLp_to = self.bg.rBp2oLp(port_to, dir=dir_to)

        demod = self.fs.parameter_to_fk(demod)

        seq, req, edges = self._solutions_AC_SRE
        solvect = edges[oLp_fr, oLp_to]

        if quadrature == "I":
            quadrature = 1
        elif quadrature == "Q":
            quadrature = 1j
        else:
            quadrature = quadrature / abs(quadrature)

        # TODO, make a upperphoton_pmAC quantum keyset
        # from icecream import ic
        # ic(solvect.kmatrix)
        datavec_p = solvect.kmatrix[(demod, "+AC")][(base.FrequencyKey({}), "+AC")][
            ..., 0, 0
        ]
        datavec_m = solvect.kmatrix[(demod, "-AC")][(base.FrequencyKey({}), "+AC")][
            ..., 0, 0
        ]
        # TODO, check realness

        return quadrature * datavec_p + quadrature.conjugate() * datavec_m

    def _solve_matrix_AC(self):
        SREIO = self.SREIO_AC(
            subtract_1=True,
        )
        (seq, req, edges, inputs, outputs) = SREIO
        del SREIO
        with self.log.heading("AC_inversion"):
            seq, req, edges = matrix.SREkmatrix_inverse(
                seq,
                req,
                edges,
                outputs_set=outputs,
                inputs_set=inputs,
                verbose=False,
                log=self.log,
            )
            edges2 = dict()
            for (n_fr, n_to), edge in edges.items():
                oLp_fr = n_fr
                oLp_to = n_to
                edges2[oLp_fr, oLp_to] = edge

        self._solutions_AC_SRE = seq, req, edges2

    def _build_matrix_AC(self):
        for obj in self.pbg.object_iter():
            try:
                visit_algo = obj.visit_matrix_algorithm_DCAC
            except AttributeError:
                continue
            else:
                # TODO verbose option for found objects?
                # print(obj)
                pass

            manip = PhysicsAlgorithmACManipulator(
                obj=obj,
                ac_algo=self,
            )

            visit_algo(manip)
        return

    def SREIO_AC(self, subtract_1=False):
        seq = collections.defaultdict(set)
        req = collections.defaultdict(set)
        edges = dict()

        for oLp_fr, eset in self.bg.link_seq.items():
            m_fr = oLp_fr
            # print("FT: ", m_fr, [map_nodes(oLp_to) for oLp_to in eset])
            for oLp_to in eset:
                m_to = oLp_to
                edges[m_fr, m_to] = 1
                seq[m_fr].add(m_to)
                req[m_to].add(m_fr)

        for obj, edict in self._object_edges.items():
            for (l_fr, l_to), e_kmat in edict.items():
                m_to = (obj, l_to)
                m_fr = (obj, l_fr)
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
        for k_in in self._drives:
            inputs.add(k_in)

        outputs = set()
        for k_in in self._views:
            outputs.add(k_in)
        return (seq, req, edges, inputs, outputs)


class PhysicsAlgorithmACView(algo_phys.PhysicsAlgorithmView):
    _ac_algo = None

    def __init__(self, ac_algo, **kw):
        super(PhysicsAlgorithmACView, self).__init__(
            bg_algo=ac_algo.bg, pbg=ac_algo.pbg, pa_algo=ac_algo.pa, **kw
        )
        self._ac_algo = ac_algo
        self._dc_algo = ac_algo.dc

    # is a dictionary of known basis elements as of this linkage
    def link_basis(self, lport):
        op = (self._obj, lport)
        btype = self._bg_algo.link_basis_types[op]
        return {
            "optical": self._pa_algo.optical_basis_AC,
            "mechanical": self._pa_algo.mechanical_basis_AC,
            "signal": self._pa_algo.signal_basis_AC,
            "electrical": self._pa_algo.electrical_basis_AC,
        }[btype]

    def configure_optical_wavelength(self, wdict):
        # TODO, check that the wavenumber exists
        # assert(len(self._pa_algo.fs.freq_set_wavelengths) == 1)
        return base.FrequencyKey(wdict)

    def parameter_to_fdict(self, fparam):
        return self._ac_algo._parameter_to_fdict(fparam)

    def optical_frequency_allowed(self, fk):
        return self._ac_algo._optical_frequency_allowed(fk)

    def basis_frequencies_DC(self, with_keys=False):
        # print("FVs: ", self._pa_algo._basis_frequency_values_DC_optical)
        if with_keys:
            return zip(
                self._pa_algo.optical_basis_DC["frequency"].enumerated,
                self._pa_algo._basis_frequency_values_DC_optical,
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
        for Wk, wnval in self.basis_wavenumbers(with_keys=True):
            for Fk, fval in self.basis_frequencies_DC(with_keys=True):
                for Qk, conj in (("+AC", False), ("-AC", True)):
                    if not conj:
                        yield (Wk, Fk, Qk), (
                            wnval,
                            fval + self._ac_algo.fs.AC_fspan,
                            conj,
                        )
                    else:
                        yield (Wk, Fk, Qk), (
                            wnval,
                            fval - self._ac_algo.fs.AC_fspan,
                            conj,
                        )
        return


class PhysicsAlgorithmACManipulator(PhysicsAlgorithmACView):
    is_DC = False
    is_AC = True

    def get_field_DC(self, lport):
        return self._dc_algo._solutions_DC[(self._obj, lport)]

    def add_drive(self, lport):
        self._ac_algo._drives.add((self._obj, lport))

    def add_view(self, lport):
        self._ac_algo._views.add((self._obj, lport))

    def add_noise(self, lport, nmatrix):
        self._ac_algo._noise[(self._obj, lport)] = nmatrix

    def add_link(self, lport_fr, lport_to, kmatrix, lowering_only=False):
        """
        Adds a link to the system matrix

        The lowering_only indicates that only the lowering (non conjugated)
        operator is supplied, and so needs some additional completion
        """
        if lowering_only:
            km_new = dict()
            for krow, kvect in kmatrix.kmatrix.items():
                kvect_new_p = km_new.setdefault(krow + ("+AC",), dict())
                kvect_new_n = km_new.setdefault(krow + ("-AC",), dict())
                for kcol, kdm in kvect.items():
                    kvect_new_p[kcol + ("+AC",)] = kdm
                    kvect_new_n[kcol + ("-AC",)] = kdm.conjugate()
            # TODO, not always using the optical basis for this
            kmatrix = kmatrix.__class__(
                stR=kmatrix.stR + (self._pa_algo.optical_basis_AC["quantum"],),
                stC=kmatrix.stC + (self._pa_algo.optical_basis_AC["quantum"],),
                dtR=kmatrix.dtR,
                dtC=kmatrix.dtC,
                kmatrix=km_new,
                build=False,
                check=self.check_build,
            )

        if isinstance(kmatrix, base.KeyMatrixBase):
            basisR = set(self.link_basis(lport_to).values())
            assert set(kmatrix.stR + kmatrix.dtR) <= basisR
            basisC = set(self.link_basis(lport_fr).values())
            assert set(kmatrix.stC + kmatrix.dtC) <= basisC
            self._ac_algo._object_edges[self._obj][lport_fr, lport_to] = kmatrix
        else:
            # is probably just a number or array
            self._ac_algo._object_edges[self._obj][lport_fr, lport_to] = kmatrix
