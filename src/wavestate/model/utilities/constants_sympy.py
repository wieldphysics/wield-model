# -*- coding: utf-8 -*-
"""
"""
from wavestate.bunch import Bunch
import sympy


constants_sympy = Bunch(
    c_m_s                 = 299792458,
    kB_J_K                = 1.380658e-23,
    h_Js                  = 6.6260700408e-34,
    hbar_Js               = 1.0545718001e-34,
    pi                    = sympy.pi,
    i                     = sympy.I,
    i2pi                  = 2 * sympy.pi * sympy.I,
    temp_K                = 299,
)

