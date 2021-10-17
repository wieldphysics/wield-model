#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


from wavestate.model import optics
from wavestate.model import base

def system1064(SHG = True):
    obj_sys = base.SimulationObject()
    freqs = base.SimulationObject()

    obj_sys['frequencies/'] = freqs

    freqs['Nd1064/'] = optics.OpticalFrequency()
    freqs['Nd1064/wavelength[m]'] = 1064e-9
    if SHG:
        freqs['Nd1064/order'] = 2
    else:
        freqs['Nd1064/order'] = 1

    aliases1064 = freqs['aliases_1064/'] = optics.OpticalFrequencyAliases()
    aliases1064['to'] = {'Nd1064' : 1}
    aliases1064['names'] = ['1064', 1064, '1064nm', 1064e-9]

    aliases532 = freqs['aliases_532/'] = optics.OpticalFrequencyAliases()
    aliases532['to'] = {'Nd1064' : 2}
    aliases532['names'] = ['532', 532, '532nm', 532e-9]

    return obj_sys


def system1550(SHG = True):
    obj_sys = base.SimulationObject()
    freqs = base.SimulationObject()

    obj_sys['frequencies/'] = freqs

    freqs['1550/'] = optics.OpticalFrequency()
    freqs['1550/wavelength[m]'] = 1550e-9
    if SHG:
        freqs['1550/order'] = 2
    else:
        freqs['1550/order'] = 1

    aliases1550 = freqs['aliases_1550/'] = optics.OpticalFrequencyAliases()
    aliases1550['to'] = {'1550' : 1}
    aliases1550['names'] = ['1550', 1550, '1550nm', 1550e-9]

    aliases775 = freqs['aliases_775/'] = optics.OpticalFrequencyAliases()
    aliases775['to'] = {'1550' : 2}
    aliases775['names'] = ['775', 775, '775nm', 775e-9]

    return obj_sys
