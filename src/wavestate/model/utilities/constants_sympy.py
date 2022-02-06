#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""
from wavestate.bunch import Bunch
import sympy


constants_sympy = Bunch(
    c_m_s=299792458,
    kB_J_K=1.380658e-23,
    h_Js=6.6260700408e-34,
    hbar_Js=1.0545718001e-34,
    pi=sympy.pi,
    i=sympy.I,
    i2pi=2 * sympy.pi * sympy.I,
    temp_K=299,
)
