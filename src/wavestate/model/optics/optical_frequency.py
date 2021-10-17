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
from .. import base

class OpticalFrequency(base.SimulationObject):
    def __init__(self):
        super(OpticalFrequency, self).__init__()
        with self._internal():
            self['wavelength[m]'] = None
            self['order'] = 1

            self.set_assign(
                pfunc = lambda _lambda : 2 * np.pi/_lambda,
                kto   = 'wavenumber[1_m]',
                kfrom = 'wavelength[m]',
            )


