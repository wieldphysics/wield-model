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
import pytest
from icecream import ic


from wavestate.model import pgraph
from wavestate.model import optics
from wavestate.model import base
from wavestate.model import model
from wavestate.model.system import algo_bg
from wavestate.model.system import algo_freq
from wavestate.model.system import algo_phys
from wavestate.model.system import algo_graphs

from wavestate.model.system import algo_log
from wavestate.pytest import relfile_test

pytestmark = pytest.mark.xfail(reason="all tests still WIP")


def system_build():
    obj_sys = model.system1064()
    freqs = obj_sys["frequencies/"]
    freqs["Nd1064/order"] = 1

    freqs["RF9/"] = base.Frequency()
    freqs["RF9/frequency[Hz]"] = 9e6
    freqs["RF9/order_optical"] = 2

    obj_sys["L1/"] = optics.Laser()
    obj_sys["L1/wavelength"] = "Nd1064"

    obj_sys["EOM1/"] = optics.SimplePhaseModulator()
    obj_sys["EOM1/frequency"] = "RF9"
    obj_sys["EOM1/index"] = 0.01

    obj_sys["REFL/"] = optics.PhotodiodeUnphysical()
    obj_sys["M1/"] = optics.Mirror()
    obj_sys["M2/"] = optics.Mirror()
    obj_sys["M1_M2/"] = optics.Space()
    obj_sys["M1_M2/length[m]"] = 1

    obj_sys["L1/power[W]"] = 1
    obj_sys["M1/T"] = 0.001
    obj_sys["M2/T"] = 0.001

    obj_sys.bond_add("L1+A | EOM1+A-t | REFL+B-t | M1+B-t | M1_M2+A-t | M2+A")
    return obj_sys


# @pytest.mark.skip()
def test_FP_DC(algo_log):
    obj_sys = system_build()
    pa = algo_phys.PhysicsAlgorithm(obj_sys, log=algo_log)

    from icecream import ic

    ic(pa.dc._solutions_DC)
    ic(pa.dc("REFL/+Wpd"))
    ic(pa.dc("REFL/+Wpd", demod="RF9"))


# @pytest.mark.skip()
def test_FP_DC_scan(plot):
    fpath = relfile_test(__file__, "plots")
    log = algo_log.LoggingAlgorithm(
        log_level=9,
        log_folder=fpath,
        filters={
            # r'digraph' : dict(investigate = True),
        },
    )
    obj_sys = system_build()
    obj_sys["M1_M2/length_scan[m]"] = np.linspace(-1e-9, 1e-9, 300)
    obj_sys["EOM1/"] = optics.SimpleSSBUpperModulator()
    obj_sys["M2/T"] = 0.001
    pa = algo_phys.PhysicsAlgorithm(obj_sys, log=log)

    # from icecream import ic
    # ic(pa.dc._solutions_DC)
    # ic(pa.dc('REFL/'))
    # ic(pa.dc('REFL/', demod = 'RF9'))

    if plot:
        from wavestate.utilities.mpl import mplfigB

        axB = mplfigB(Nrows=2)
        axB.ax0.plot(
            obj_sys["M1_M2/length_scan[m]"],
            pa.dc("REFL/+Wpd"),
        )
        axB.ax1.plot(
            obj_sys["M1_M2/length_scan[m]"],
            pa.dc("REFL/+Wpd", demod="RF9").real,
        )
        axB.ax1.plot(
            obj_sys["M1_M2/length_scan[m]"],
            pa.dc("REFL/+Wpd", demod="RF9").imag,
        )
        # axB.ax3.matshow(ikm_dense != 0)
        from os import path

        axB.save(path.join(fpath, "testFP_PDH.pdf"))
    else:
        print("use --plot to plot")


def test_FP_AC(plot):
    fpath = relfile_test(__file__, "plots")
    log = algo_log.LoggingAlgorithm(
        log_level=9,
        log_folder=fpath,
        filters={
            r"digraph": dict(investigate=True),
        },
    )
    obj_sys = system_build()
    F_AC_Hz = obj_sys["AC/frequency_span[Hz]"] = np.logspace(1, 7, 1000)
    obj_sys["M1_M2/length_scan[m]"] = 0.50e-9
    # obj_sys['EOM1/']  = optics.SimpleSSBUpperModulator()
    pa = algo_phys.PhysicsAlgorithm(obj_sys, log=log)

    # from icecream import ic
    # ic(pa.dc._solutions_DC)
    # ic(pa.dc('REFL/'))
    # ic(pa.dc('REFL/', demod = 'RF9'))

    if plot:
        from wavestate.utilities.mpl import mplfigB

        axB = mplfigB(Nrows=3)
        axB.ax0.loglog(
            F_AC_Hz,
            abs(pa.ac("M1+Dl", "REFL+Wpd")),
        )
        axB.ax1.loglog(
            F_AC_Hz,
            abs(pa.ac("M1+Dl", "REFL+Wpd", demod="RF9", quadrature="I")),
        )
        axB.ax1.loglog(
            F_AC_Hz,
            abs(pa.ac("M1+Dl", "REFL+Wpd", demod="RF9", quadrature="Q")),
        )
        axB.ax2.semilogx(
            F_AC_Hz,
            np.angle(pa.ac("M1+Dl", "REFL+Wpd", demod="RF9", quadrature="I"), deg=True),
        )
        axB.ax2.semilogx(
            F_AC_Hz,
            np.angle(pa.ac("M1+Dl", "REFL+Wpd", demod="RF9", quadrature="Q"), deg=True),
        )
        # axB.ax3.matshow(ikm_dense != 0)
        from os import path

        axB.save(path.join(fpath, "testFP_AC.pdf"))
    else:
        print("use --plot to plot")


# def test_FP_DC_matrix():
#    obj_sys = system_build()
#    freqs = obj_sys['frequencies/']
#    pa = algo_phys.PhysicsAlgorithm(obj_sys)
#    #ga = algo_graphs.GraphPlottingAlgorithm(pa)
#
#    #from icecream import ic
#    #ic(edges)
#
#    from wavestate.utilities.np.SRE.semidense import SREkmatrix_inverse
#
#    nmap, SREIO = pa.dc.SREIO_DC(
#        map_nodes = True,
#        subtract_1 = True,
#    )
#    (seq, req, edges, inputs, outputs) = SREIO
#    with Timer() as t:
#        SREkmatrix_inverse(
#            seq, req, edges,
#            outputs_set = outputs,
#            inputs_set  = inputs,
#            verbose     = False,
#        )
#    print("TIME: ", t.interval)
#
#    #assert(False)
#
# def test_FP_DC_noise():
#    obj_sys = system_build()
#    freqs = obj_sys['frequencies/']
#    pa = algo_phys.PhysicsAlgorithm(obj_sys)
#    #ga = algo_graphs.GraphPlottingAlgorithm(pa)
#
#    from icecream import ic
#    ic(pa.dc._noise)
#
# import time
# class Timer(object):
#    def __enter__(self):
#        self.start = time.clock()
#        return self
#
#    def __exit__(self, *args):
#        self.end = time.clock()
#        self.interval = self.end - self.start
