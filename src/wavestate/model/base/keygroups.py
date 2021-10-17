#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

from wavestate.utilities.np.semidense import KeyGroup


kg_polarization = KeyGroup(
    'polarization',
    ('S', 'P'),
)

kg_twophoton_qp = KeyGroup(
    'quantum',
    ('q', 'p'),
    identifier = 'qp',
)

kg_twophoton_pm = KeyGroup(
    'quantum',
    ('+', '-'),
    identifier = 'pm',
)

kg_twophoton_qpAC = KeyGroup(
    'quantum',
    ('qAC', 'pAC'),
    identifier = 'qpAC',
)

kg_twophoton_pmAC = KeyGroup(
    'quantum',
    ('+AC', '-AC'),
    identifier = 'pmAC',
)

