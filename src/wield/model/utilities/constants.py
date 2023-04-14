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
from wield.bunch import Bunch


constants_floats = Bunch(
    c_m_s=299792458,
    kB_J_K=1.380658e-23,
    h_Js=6.6260700408e-34,
    hbar_Js=1.0545718001e-34,
    pi=np.pi,
    i=1j,
    i2pi=np.pi * 2j,
    temp_K=299,
)
