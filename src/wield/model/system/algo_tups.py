#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""
from collections import namedtuple

ObjectLinkageTup = namedtuple('ObjectLinkageTup', ['obj', 'linkage'])

ObjectPortTup = namedtuple('ObjectLinkageTup', ['obj', 'port'])

ObjectParamTup = namedtuple('ObjectLinkageTup', ['obj', 'parameter'])


