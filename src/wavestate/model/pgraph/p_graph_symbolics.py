# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals, absolute_import
from transient.utilities import ref_p_split, ref_port_split
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

