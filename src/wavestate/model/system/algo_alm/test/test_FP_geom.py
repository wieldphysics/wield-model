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
from os import path

from wavestate.model import optics
from wavestate import model
from wavestate.model.system import algo_phys


from wavestate.pytest import tpath_join, dprint, plot  # noqa: F401

def linear_cavity(sys=None):
    if sys is None:
        sys = model.base.SimulationObject()

    sys["M0_M1/"] = optics.Space()
    sys["M1_M2/"] = optics.Space()

    sys["M0/"] = optics.Mirror()
    sys["M1/"] = optics.Mirror()
    sys["M2/"] = optics.Mirror()

    sys["M0_M1/length[m]"] = 1
    sys["M1_M2/length[m]"] = 1

    sys["M1/AOI[deg]"] = 0
    sys["M2/AOI[deg]"] = 0

    sys["M1/ROC[m]"] = None
    sys["M2/ROC[m]"] = -2

    sys["M1/depth[m]"] = 0.00
    sys["M2/depth[m]"] = 0.00

    sys["M0/T"] = 1
    sys["M1/T"] = 7600e-6
    sys["M2/T"] = 7540e-6

    sys["cav/"] = optics.Cavity()
    sys["cav/waypoints"] = ["M1/", "M2/"]
    # add the cavity loop
    sys.bond_add(
        "M0/+A1 | M0_M1/+A-t | M1/+B1-t | M1_M2/+A-t | M2/+A1"
    )

    sys.port_forward_add(pfrom="+IC", pto="M1/+B1")

    return sys

def T_linear_cavity_view(tpath_join, algo_log, plot, dprint):
    sys = model.system1064()
    sys["cavity/"] = linear_cavity()
    pa = algo_phys.PhysicsAlgorithm(sys, log=algo_log)

    def print_things(cavity, params):
        dprint("--------------")
        dprint(
            "Gouy: X={:.2f}deg, {:.2f}pi, Y={:.2f}deg, {:.2f}pi".format(
                params.gouyX_deg,
                params.gouyX_deg / 180,
                params.gouyY_deg,
                params.gouyY_deg / 180,
            )
        )


    with pa.pbg.preferred():
        print('------- M1+A in')
        params = pa.mm.cavity_parameters("cavity/cav", "cavity/M1+A1!i", Wk=1064, shifts_use=True)
        print(params.cavB.cavity_trans.X.shifts)
        print(params.cavB.cavity_trans.X.shifts_out_referred)

        print('------- M1+A out')
        params = pa.mm.cavity_parameters("cavity/cav", "cavity/M1+A1!o", Wk=1064, shifts_use=True)
        print(params.cavB.cavity_trans.X.shifts)
        print(params.cavB.cavity_trans.X.shifts_out_referred)

        pa.mm.cavity_digest(Wk=1064)

        overlapper = pa.mm.overlap(
            target_fr=[
                "cavity/cav",
            ],
            target_to=None,
            waypoints=[
                "cavity/M1+A!i",
                "cavity/M0/",
            ],
            Wk=1064,
            shifts_use=True,
        )

        cavB = overlapper['cavity/cav']
        #print(cavB.shiftsXend)
        plB = overlapper.plot(
            tpath_join("linearM1.pdf"),
            # reverse = True,
            # self_overlap = True,
        )
        print(overlapper.shifts_table_str())

        overlapper = pa.mm.overlap(
            target_fr=[
                "cavity/cav",
            ],
            target_to=None,
            waypoints=[
                "cavity/M2/",
                "cavity/M0/",
            ],
            Wk=1064,
            shifts_use=True,
        )

        cavB = overlapper['cavity/cav']
        #print(cavB.shiftsXend)
        plB = overlapper.plot(
            tpath_join("linearM2.pdf"),
            # reverse = True,
            # self_overlap = True,
        )
        #print(overlapper.shifts_table())
        print(overlapper.shifts_table_str(var='pos'))
        print(overlapper.shifts_table_str(var='ang'))
    return

