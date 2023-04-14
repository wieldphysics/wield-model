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
from wield.model import optics
from wield.model import base


def linear_cavity(sys=None):
    if sys is None:
        sys = base.SimulationObject()

    sys["M1_M2/"] = optics.Space()

    sys["M1/"] = optics.Mirror()
    sys["M2/"] = optics.Mirror()

    sys["M1_M2/length[m]"] = 1

    sys["M1/AOI[deg]"] = 0
    sys["M2/AOI[deg]"] = 0

    sys["M1/ROC[m]"] = 2
    sys["M2/ROC[m]"] = 2

    sys["M1/depth[m]"] = 0.01
    sys["M2/depth[m]"] = 0.01

    sys["M1/T"] = 7600e-6
    sys["M2/T"] = 7540e-6

    sys["cav/"] = optics.Cavity()
    sys["cav/waypoints"] = ["M1/", "M2/"]
    # add the cavity loop
    sys.bond_add(
        "M1/+A1 | M1_M2/+A-t | M2/+A1"
    )

    sys.port_forward_add(pfrom="+IC", pto="M1/+B1")

    return sys
