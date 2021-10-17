#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import re

re_split_refs = re.compile(r'\.|/')

def ref_2_rtup(key):
    #TODO, assert on bond symbols
    idx = key.rfind('/')
    is_param = False
    if idx == -1:
        is_param = True
        params = key
        refs = ()
    else:
        refs = key[:idx]
        params = key[idx+1:]
        params = params.strip()
        refs = re_split_refs.split(refs)
        if params:
            is_param = True
        if refs[0] == '':
            refs = refs[1:]
    params = params.split('.')
    if is_param:
        raise RuntimeError("was only expecting a reference, not a reference and parameter")
    else:
        return tuple(refs)

def ref_value_split(key):
    #TODO, assert on bond symbols
    idx = key.rfind('/')
    is_param = False
    if idx == -1:
        is_param = True
        params = key
        refs = ()
    else:
        refs = key[:idx]
        params = key[idx+1:]
        params = params.strip()
        refs = re_split_refs.split(refs)
        if params:
            is_param = True
        if refs[0] == '':
            refs = refs[1:]
    params = params.split('.')
    if is_param:
        return tuple(refs), tuple(params)
    else:
        return tuple(refs), None


def ref_port_split(key):
    idx = key.rfind('+')
    if idx == -1:
        raise RuntimeError("Expecting a bond specifier (requires + in name)")
    else:
        refs = key[:idx]
        bond = key[idx:]
        bond = bond.strip()
        refs = re_split_refs.split(refs)
        #clear the last one if it is specified as A/B/C/+bond,
        #rather than A/B/C+bond
        if refs[-1] == '':
            refs = refs[:-1]
    return tuple(refs), bond
