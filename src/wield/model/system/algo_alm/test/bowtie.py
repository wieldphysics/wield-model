#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import numpy as np
from wield.model import optics
from wield.model import base


def dist_m_from_XY(pm_A, pm_B):
    V = pm_A - pm_B
    return np.dot(V, V) ** 0.5


def bowtie(
    sys=None,
    M1loc=np.array([0,       0]),
    M2loc=np.array([.2815,   0]),
    M3loc=np.array([0,       -0.0396023]),
    M4loc=np.array([0.2815,  -0.0396023]),
):
    if sys is None:
        sys = base.SimulationObject()

    sys["1_2/"] = optics.Space()
    sys["2_3/"] = optics.Space()
    sys["3_4/"] = optics.Space()
    sys["4_1/"] = optics.Space()

    sys["M1/"] = optics.Mirror()
    sys["M2/"] = optics.Mirror()
    sys["M3/"] = optics.Mirror()
    sys["M4/"] = optics.Mirror()

    sys["1_2/length[m]"] = dist_m_from_XY(M1loc, M2loc)
    sys["2_3/length[m]"] = dist_m_from_XY(M2loc, M3loc)
    sys["3_4/length[m]"] = dist_m_from_XY(M3loc, M4loc)
    sys["4_1/length[m]"] = dist_m_from_XY(M4loc, M1loc)

    sys["M1/AOI[deg]"] = 4.004
    sys["M2/AOI[deg]"] = 4.004
    sys["M3/AOI[deg]"] = 4.004
    sys["M4/AOI[deg]"] = 4.004

    sys["M3/ROC[m]"] = -2.57321
    sys["M4/ROC[m]"] = -2.57369

    sys["M1/depth[m]"] = 0.01
    sys["M2/depth[m]"] = 0.01
    sys["M3/depth[m]"] = 0.00635
    sys["M4/depth[m]"] = 0.00635

    sys["M1/T"] = 7600e-6
    sys["M2/T"] = 7540e-6
    sys["M3/T"] = 35.9e-6
    sys["M4/T"] = 36.0e-6

    sys["M1M2/"] = optics.Cavity()
    sys["M1M2/waypoints"] = ["M1/", "M2/", "M4/"]
    sys["M1M4/"] = optics.Cavity()
    sys["M1M4/waypoints"] = ["M1/", "M4/", "M2/"]
    # add the cavity loop
    sys.bond_add(
        "M1/+A1 | 1_2/+A-t | M2/+A1-r | 2_3/+A-t | M3/+A1-r | 3_4/+A-t | M4/+A1-r | 4_1/+A-t | M1/+A2"
    )

    sys.port_forward_add(pfrom="+IC", pto="M1/+B1")
    sys.port_forward_add(pfrom="+OC", pto="M2/+B1")

    return sys
