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

from wavestate.utilities.np import matrix_stack

from .. import base


class Space(base.OpticalObject):
    annotate_as_space = True

    def __init__(self):
        super(Space, self).__init__()
        with self._internal():
            self['length[m]'] = None
            self['length[m].lower_bound'] = 0
            self['length_scan[m]'] = 0
            self['gouy_phase[rad]'] = None

        self.values_settable([
            'length[m]',
            'length[m].lower_bound',
            'length_scan[m]',
            'gouy_phase[rad]',
        ])

    def port_chain(self, p, pname):
        bmap = {
            '+A-t' :  (None, '+B'),
            '+B-t' :  (None, '+A'),
        }.get(pname, None)

        if bmap is not None:
            return bmap

        return super(Space, self).port_chain(p, pname)

    @classmethod
    def visit_port_information(cls, manip):
        manip.gen_optical_port('+A', 'A')
        manip.gen_optical_port('+B', 'B')
        return

    def visit_matrix_algorithm_DCAC(self, manip):
        settings = manip.optical_settings()

        #expects basis to be the same in A!i, A!o, B!i, B!o,
        basis = manip.link_basis('A!i')

        if settings.transverse is None:
            kmatrix_TM = 1
        elif settings.transverse == 'HGXY':
            gpX = manip.get_gouy_X()
            gpY = manip.get_gouy_Y()
            GX = np.diag(gpX * np.array(basis['HGX'].enumerated))
            GY = np.diag(gpY * np.array(basis['HGY'].enumerated))

            #TODO, get the algorithm to cache these
            kmatrix_TMX = matrix.KeyMatrixSame(
                dt      = (basis['HGX'],),
                kmatrix = GX,
                build   = True,
                check   = manip.check_build,
            )

            #TODO, get the algorithm to cache these
            kmatrix_TMY = matrix.KeyMatrixSame(
                dt      = (basis['HGY'],),
                kmatrix = GY,
                build   = True,
                check   = manip.check_build,
            )

            kmatrix_TM = kmatrix_TMX @ kmatrix_TMY

        elif settings.transverse == 'HG':
            gp = manip.get_gouy()
            kdm_G = np.diag(gp * np.array(basis['HG'].enumerated))

            kmatrix_TM = matrix.KeyMatrixSame(
                dt      = (basis['HG'],),
                kmatrix = {() : {(): kdm_G}},
                build   = False,
                check   = manip.check_build,
            )
        else:
            raise RuntimeError("Unsupported Transverse mode setting")

        #print("BFREQ", freq)

        #this length is to act on the wavelengths only
        length_scan_m = manip.p['length_scan[m]']
        #this length is for acting on the frequencies
        length_m = manip.p['length[m]'] + length_scan_m
        if length_m is None:
            raise RuntimeError("Must specify length of space")

        km = collections.defaultdict(dict)
        for (Wk, Fk, Qk), (wnval, fval, conj) in manip.basis_WFQ_optical_pm():
            if not conj:
                km[(Wk, Fk, Qk)][(Wk, Fk, Qk)] = [[
                    np.exp(
                        -2j * np.pi * fval * length_m / manip.constants.c_m_s
                        - 1j*wnval*length_scan_m
                    )]]
            else:
                km[(Wk, Fk, Qk)][(Wk, Fk, Qk)] = [[
                    np.exp(
                        +2j * np.pi * fval * length_m / manip.constants.c_m_s
                        + 1j*wnval*length_scan_m
                    )]]
        kmatrix_F = matrix.KeyMatrixSame(
            st      = (
                basis['wavenumber'],
                basis['frequency'],
                basis['quantum'],
            ),
            kmatrix = km,
            build   = True,
            check   = manip.check_build,
        )

        kmatrix = kmatrix_TM @ kmatrix_F

        manip.add_link('B!i', 'A!o', kmatrix)
        manip.add_link('A!i', 'B!o', kmatrix)
        return

    def visit_mode_matching_linkage(self, manip):
        manip.add_link('B!i', 'A!o', manip.p['length[m]'])
        manip.add_link('A!i', 'B!o', manip.p['length[m]'])

    def visit_mode_matching_transport(self, manip):
        length_m = manip.p['length[m]']
        #no need to check these since the space can only be called on proper
        #links and the two directions are identical
        #manip.lport_fr
        #manip.lport_to

        #the P-builders are for fast optimization solving
        def p_builderXY(p):
            length_m = p['length[m]']
            return matrix_stack([[1, length_m], [0, 1]])

        manip.set_XYpropagator(p_builderXY)

        manip.set_Zpropagator({'length[m]' : 1})

        matrix = p_builderXY(manip.p)

        def inc_builder(z):
            return matrix_stack([[1, z], [0, 1]])

        manip.set_XYincremental([
            (length_m, inc_builder, matrix)
        ])
        return


#def _build_dense():
#    km = collections.defaultdict(dict)
#    freq = manip.basis_frequencies()
#    for Wk, wnval in manip.basis_wavenumbers(with_keys = True):
#        kdm_freq = []
#        for idx, fval in enumerate(freq):
#            frow = [0] * len(freq)
#            frow[idx] = np.exp(
#                2j * np.pi * fval * length_m / manip.constants.c_m_s
#                + wnval*length_scan_m
#            )
#            kdm_freq.append(frow)
#        km[(Wk,)][(Wk,)] = kdm_freq
#
#    kmatrix_F = matrix.KeyMatrixSame(
#        st      = (basis['wavenumber'],),
#        dt      = (basis['frequency'],),
#        kmatrix = km,
#        build   = True,
#        check   = manip.check_build,
#    )
