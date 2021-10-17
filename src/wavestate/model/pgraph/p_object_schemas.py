#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import collections

class ParameterObjectBaseSchema(object):
    """
    This object isn't for use or importing, it just has code to help document
    the dictionary schema used for the actual object
    """
    def __init__(self):
        self._value_dict = collections.defaultdict(dict)
        #schema
        {
            ("rtup", "vtup") : {
                "default_value" : 123,
                #if this is a string, then it is a default assignment
                #otherwise it is found in _default_functions
                #the function must get all arguments and make all assignments
                #through the passed single argument p
                "default_func" : lambda p : None,
                #if units are known, then they are checked and converted
                "units" :          None,
                "allow_sympy" :    False,
                "allow_numpy" :    False,
                "allow_numrep" :   False,
                "allow_casadiMX" : False,
                "allow_casadiSX" : False,
            },
        }

        self._reference_dict = collections.defaultdict(dict)
        #schema
        {
            "rtup" : {
                "default_value" : 123,
                "default_assign" : "rtup",
                #if this is a string, then it is a default assignment
                #otherwise it is found in _default_functions
                "default_func" : lambda x : x,
                "default_func_deps" : (('rtup', 'vtup'), ...)
            },
        }

        #also used for reference defaults
        self._default_functions = collections.defaultdict(dict)
        {
            'func' : {  # func is the function object (not a string)
                "deps" : (),  # list of parameter deps. None if unknown
                "assigns" : (),  # list of parameter assignments. None if unknown
            }
        }

        self._bond_generators = set()
        return
