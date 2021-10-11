
# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import numpy as np
from .. import base

class OpticalFrequency(base.ParameterObject):
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


class OpticalFrequencyAliases(base.ParameterObject):
    pass
