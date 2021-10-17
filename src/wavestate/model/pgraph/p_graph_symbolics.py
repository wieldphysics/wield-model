#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

from wavestate.utilities import ref_p_split, ref_port_split
import sympy
import collections

from .p_graph import ParameterGraph


class ParameterGraphSymbolics(ParameterGraph):
    """
    For running some computations sybolically
    """
    def __init__(self, root):
        super(ParameterGraphSymbolics, self).__init__(root)
        if not isinstance(root, ParameterGraphSymbolics):
            self.op2symbols = collections.defaultdict(dict)
            self.symbols2op = dict()
        else:
            raise NotImplementedError()

