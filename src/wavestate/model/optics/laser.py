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
from . import vacuum


class Laser(base.OpticalObject):
    def __init__(self):
        super(Laser, self).__init__()
        with self._internal():
            self['power[W]'] = None
            self['polarization'] = 'P'
            self['wavelength'] = None
            self['frequency'] = None

    @classmethod
    def visit_port_information(cls, manip):
        manip.gen_optical_port('+A', 'A')
        return

    @classmethod
    def visit_matrix_algorithm_DCAC(cls, manip):
        if manip.is_DC:
            settings = manip.optical_settings()
            basis = manip.link_basis('A!o')

            assert(settings.quantum == 'pm')

            dk = {}

            if settings.polarization:
                dk['polarization'] = manip.p['polarization']

            if settings.transverse is None:
                pass
            elif settings.transverse == 'HGXY':
                dk['HGX'] = 0
                dk['HGY'] = 0
            elif settings.transverse == 'HG':
                dk['HG'] = 0
            else:
                raise RuntimeError("Unsupported Transverse mode setting")

            param_F = manip.p['frequency']
            if param_F is None:
                Fkey = manip.default_optical_freqkey()
            elif isinstance(param_F, collections.Mapping):
                Fkey = base.FrequencyKey(param_F)
            elif isinstance(param_F, str):
                Fkey = base.FrequencyKey({param_F : 1})
            else:
                raise RuntimeError("Unrecognized type")
            dk['frequency'] = Fkey

            param_W = manip.p['wavelength']
            if param_W is None:
                Wkey = manip.default_optical_wavekey()
            elif isinstance(param_W, collections.Mapping):
                Wkey = manip.configure_optical_wavenumber(param_W)
            elif isinstance(param_W, str):
                #TODO, check existence
                Wkey = manip.configure_optical_wavenumber({param_W : 1})
            else:
                raise RuntimeError("Unrecognized type")
            dk['wavenumber'] = Wkey

            #TODO, simplify notation
            kmatrix = {
                base.DictKey(dk, quantum = '+') : {() : [[1]]},
                base.DictKey(dk, quantum = '-') : {() : [[1]]},
            }

            svect = matrix.KeyMatrixDiff(
                stR     = basis,
                kmatrix = kmatrix,
                build   = True,
                check   = manip.check_build,
            )

            manip.add_source('A!o', svect)

        if manip.is_AC:
            kmatrix = vacuum.optical_quantum_noise_matrix(manip)
            manip.add_noise('A!o', kmatrix)

