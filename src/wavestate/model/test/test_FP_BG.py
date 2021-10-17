#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

from os import path
import numpy as np
import pytest

from wavestate.model import pgraph
from wavestate.model import optics
from wavestate.model import base
from wavestate import model
from wavestate.model.system import algo_bg
from wavestate.model.system import algo_freq
from wavestate.model.system import algo_phys
from wavestate.model.system import algo_graphs

from wavestate.model.system import algo_log
from wavestate.pytest import Timer

pytestmark = pytest.mark.xfail(reason = "all tests still WIP")

def system_build():
    obj_sys = model.system1064()
    freqs = obj_sys['frequencies/']
    freqs['Nd1064/order'] = 2

    freqs['RF9/'] = base.Frequency()
    freqs['RF9/frequency[Hz]'] = 9e6
    freqs['RF9/order_optical'] = 4
    freqs['RF6/'] = base.Frequency()
    freqs['RF6/frequency[Hz]'] = 6e6

    freqs['span96/'] = base.FrequencySpan()
    freqs['span96/frequencies'] = ['RF9', 'RF6']
    freqs['span96/order'] = 3

    freqs['suppress96_1/'] = base.FrequencySuppress()
    freqs['suppress96_1/suppress'] = {'RF9' : 1, 'RF6' : 1}

    #freqs['suppress96_2/'] = base.FrequencySuppress()
    #freqs['suppress96_2/suppress'] = {'RF9' : [1, 2], 'RF6' : [2, 1]}

    obj_sys['L1/']    = optics.Laser()
    obj_sys['L1/wavelength'] = 'Nd1064'

    obj_sys['EOM1/']  = optics.SimplePhaseModulator()
    obj_sys['EOM1/frequency'] = 'RF9'
    obj_sys['EOM1/index']     = .01

    obj_sys['REFL/']  = optics.PhotodiodeUnphysical()
    obj_sys['M1/']    = optics.Mirror()
    obj_sys['M2/']    = optics.Mirror()
    obj_sys['M1_M2/'] = optics.Space()
    obj_sys['M1_M2/length[m]'] = 1

    obj_sys['L1/power[W]'] = 1
    obj_sys['M1/T']   = .001
    obj_sys['M2/T']   = .001

    obj_sys.bond_add('L1+A | EOM1+A-t | REFL+B-t | M1+A-t | M1_M2+A-t | M2+A')
    return obj_sys

def test_fabry_perot_BG(tpath):
    print()
    log = algo_log.LoggingAlgorithm(
        log_level = 9,
        log_folder = tpath,
        filters = {
            r'digraph' : dict(investigate = True),
        }
    )

    obj_sys = system_build()

    pg = pgraph.ParameterGraph(obj_sys)
    #pg.print_parameters_eval()
    bg = algo_bg.BondGraphAlgorithm(pg)
    print(bg.link_seq)

    ap = algo_freq.FrequenciesAlgorithm(pg)
    print("OPTICAL: ", ap.freq_set_optical)
    print("SIGNAL: ", ap.freq_set_signal)
    print("MECH: ", ap.freq_set_mechanical)
    print("WLEN: ", ap.freq_set_wavenumbers)
    return


#@pytest.mark.skip()
def test_fabry_perot_DC(tpath):
    print()
    log = algo_log.LoggingAlgorithm(
        log_level = 9,
        log_folder = tpath,
        filters = {
            r'digraph' : dict(investigate = True),
        }
    )

    obj_sys = system_build()
    pa = algo_phys.PhysicsAlgorithm(obj_sys, log = log)
    ga = algo_graphs.GraphPlottingAlgorithm(pa)

    ga.graph_bond_graph(path.join(tpath, 'test_bg.pdf'))
    ga.graph_DC_link_matrix(path.join(tpath, 'test_DC_LM.pdf'))
    #from icecream import ic
    #ic(edges)
    nmap, SREIO = pa.dc.SREIO_DC(
        map_nodes = True,
        subtract_1 = True,
    )
    (seq, req, edges, inputs, outputs) = SREIO

    from wavestate.utilities.np.SRE.semidense import SREkmatrix_inverse

    SREkmatrix_inverse(
        seq, req, edges,
        outputs_set = outputs,
        inputs_set  = inputs,
        verbose     = False,
        log         = log,
    )

    #assert(False)

def test_fabry_perot_DC_matrix(tpath):
    print()
    log = algo_log.LoggingAlgorithm(
        log_level = 9,
        log_folder = tpath,
        filters = {
            r'digraph' : dict(investigate = True),
        }
    )

    obj_sys = system_build()
    freqs = obj_sys['frequencies/']
    freqs['RF6/frequency_span[Hz]'] = np.linspace(5.5e6, 6.5e6, 3)
    pa = algo_phys.PhysicsAlgorithm(obj_sys)
    #ga = algo_graphs.GraphPlottingAlgorithm(pa)

    #from icecream import ic
    #ic(edges)

    from wavestate.utilities.np.SRE.semidense import SREkmatrix_inverse

    nmap, SREIO = pa.dc.SREIO_DC(
        map_nodes = True,
        subtract_1 = True,
    )
    (seq, req, edges, inputs, outputs) = SREIO
    with Timer() as t:
        SREkmatrix_inverse(
            seq, req, edges,
            outputs_set = outputs,
            inputs_set  = inputs,
            verbose     = False,
        )
    print("TIME: ", t.interval)

    #assert(False)

@pytest.mark.skip
def test_FP_DC_noise(ic, tpath):
    print()
    log = algo_log.LoggingAlgorithm(
        log_level = 9,
        log_folder = tpath,
        filters = {
            r'digraph' : dict(investigate = True),
        }
    )

    obj_sys = system_build()
    freqs = obj_sys['frequencies/']
    freqs['RF6/frequency_span[Hz]'] = np.linspace(5.5e6, 6.5e6, 3)
    pa = algo_phys.PhysicsAlgorithm(obj_sys)
    #ga = algo_graphs.GraphPlottingAlgorithm(pa)

    ic(pa._noise)
