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
import scipy.optimize
import declarative

from .beam_param import (
    ComplexBeamParam
)

from wavestate.utilities.mpl import (
    mplfigB,
)

from . import utils

class QFit(declarative.OverridableObject):
    @declarative.dproperty
    def wavelength_m(self, val):
        return val

    m2 = 1.00

    @declarative.mproperty
    def D_um(self, arg):
        arg = np.array(arg)
        return arg

    @declarative.mproperty
    def R_m(self, arg = declarative.NOARG):
        if arg is declarative.NOARG:
            arg = self.D_um * 1e-6 / 2
        else:
            arg = np.array(arg)
        return arg

    @declarative.mproperty
    def Z_in(self, arg):
        arg = np.array(arg)
        return arg

    @declarative.mproperty
    def Z_m(self, arg = declarative.NOARG):
        if arg is declarative.NOARG:
            arg = self.Z_in * .0254
        else:
            arg = np.array(arg)
        return arg

    @declarative.mproperty
    def Z0_ZR_init(self, arg = declarative.NOARG):
        if arg is declarative.NOARG:
            idx_W0 = np.argsort(self.R_m)
            W0 = self.R_m[idx_W0[0]]  * 1

            ZR = np.pi*W0**2/(self.wavelength_m) / (self.m2)
            Z0 = -np.mean(self.Z_m[idx_W0[:4]])
            arg = (Z0, ZR)
        return arg

    def waist_func(self, z, z_0, z_R):
        return (self.m2 * self.wavelength_m / (np.pi * z_R) * ((z + z_0)**2 + z_R**2))**.5

    no_prefit = False
    @declarative.mproperty
    def Z0_ZR_fit(self):
        idx_W0 = np.argmin(self.R_m)
        init = self.Z0_ZR_init
        #do a prefit to try and find tiny waists using a subset of the data
        if idx_W0 > 1 and idx_W0 < len(self.R_m) - 1 and not self.no_prefit and len(self.R_m) > 3:
            #don't include the point
            if idx_W0 < len(self.R_m) / 2:
                idx_W0 += 1
                #ignore the actual point itself as it may be across a gap
                init, hess = scipy.optimize.curve_fit(
                    self.waist_func,
                    self.Z_m[idx_W0:],
                    self.R_m[idx_W0:],
                    p0 = self.Z0_ZR_init
                )
            else:
                init, hess = scipy.optimize.curve_fit(
                    self.waist_func,
                    self.Z_m[:idx_W0],
                    self.R_m[:idx_W0],
                    p0 = self.Z0_ZR_init
                )
        (z0, zR), hess = scipy.optimize.curve_fit(self.waist_func, self.Z_m, self.R_m, p0 = init)
        return (z0, zR)

    def waist_func_m2(self, z, z_0, z_R, m2):
        return (abs(m2) * self.wavelength_m / (np.pi * z_R) * ((z + z_0)**2 + z_R**2))**.5

    @declarative.mproperty
    def Z0_ZR_m2_fit(self):
        idx_W0 = np.argmin(self.R_m)
        init = self.Z0_ZR_init + (self.m2,)
        #do a prefit to try and find tiny waists using a subset of the data
        if idx_W0 > 1 and idx_W0 < len(self.R_m) - 1 and not self.no_prefit and len(self.R_m) > 3:
            #don't include the point
            if idx_W0 < len(self.R_m) / 2:
                idx_W0 += 1
                #ignore the actual point itself as it may be across a gap
                init, hess = scipy.optimize.curve_fit(
                    self.waist_func_m2,
                    self.Z_m[idx_W0:],
                    self.R_m[idx_W0:],
                    p0 = self.Z0_ZR_fit + (self.m2,)
                )
            else:
                init, hess = scipy.optimize.curve_fit(
                    self.waist_func_m2,
                    self.Z_m[:idx_W0],
                    self.R_m[:idx_W0],
                    p0 = self.Z0_ZR_fit + (self.m2,)
                )
        (z0, zR, m2), hess = scipy.optimize.curve_fit(self.waist_func_m2, self.Z_m, self.R_m, p0 = init)
        return (z0, zR, m2)

    @declarative.mproperty
    def q_fit(self):
        return ComplexBeamParam.from_Z_ZR(
            self.Z0_ZR_fit[0],
            self.Z0_ZR_fit[1],
            wavelength_m = self.wavelength_m,
        )

    @declarative.mproperty
    def q_fit_m2(self):
        return ComplexBeamParam.from_Z_ZR(
            self.Z0_ZR_m2_fit[0],
            self.Z0_ZR_m2_fit[1],
            wavelength_m = self.wavelength_m,
        )

    @declarative.mproperty
    def q_init(self, initval = None):
        if initval is None:
            return ComplexBeamParam.from_Z_ZR(
                self.Z0_ZR_init[0],
                self.Z0_ZR_init[1],
                wavelength_m = self.wavelength_m,
            )
        else:
            return initval

    def rep(self, place_in = 0):
        print(self.Z0_ZR_fit)
        try:
            print("ComplexBeamParam.from_Z_ZR({0}, {1}, wavelen = {2})".format(
                self.Z0_ZR_fit[0] + place_in * .0254,
                self.Z0_ZR_fit[1],
                self.wavelength_m * 1e9,
            ))
        except Exception as e:
            print(e)

    def plot(
        self,
        with_init = False,
        with_m2fit = True,
        with_fit = True,
    ):
        F = mplfigB()
        diff = max(self.Z_m) - min(self.Z_m)
        Z_pts = np.linspace(min(self.Z_m) - diff/8, max(self.Z_m) + diff/8, 100)

        if int(self.wavelength_m) == 1064:
            color_pts = 'red'
            color_fit = 'orange'
            color_init = 'purple'
        elif int(self.wavelength_m) == 532:
            color_pts = 'blue'
            color_fit = 'green'
            color_init = 'purple'
        else:
            color_pts = 'blue'
            color_fit = 'black'
            color_init = 'purple'

        F.ax0.scatter(
            self.Z_in,
            self.D_um,
            color = color_pts,
            label = 'data',
        )

        fit_label = (u"Fit: $Z_0$ = {Zm} = {Zin:.1f}in\nW0={W0} D0={D0}\nZR={ZR}\n $M^2$={m2:.2f} (L<{MML:.1f}%)".format(
            Zm      = utils.str_m(-self.q_fit.Z, d = 3),
            Zin     = -self.q_fit.Z / .0254,
            ZR      = utils.str_m(self.q_fit.ZR, d = 4),
            W0      = utils.str_m(self.q_fit.W0, d = 4),
            D0      = utils.str_m(2 * self.q_fit.W0, d = 4),
            m2      = self.m2,
            MML = (self.m2 - 1)/7 * 100,
        ))

        if with_fit:
            F.ax0.plot(
                Z_pts / .0254,
                self.m2**0.5 * 2 * 1e6 * self.q_fit.propagate_distance(Z_pts).W,
                color = color_fit,
                label = fit_label,
            )

        if with_m2fit:
            fit_label_m2 = (u"Fit: $Z_0$ = {Zm} = {Zin:.1f}in\nW0={W0} D0={D0}\nZR={ZR}\n $M^2$={m2:.2f} (L<{MML:.1f}%, $L^*$={XL:.1f}%)".format(
                Zm      = utils.str_m(-self.q_fit_m2.Z, d = 3),
                Zin     = -self.q_fit_m2.Z / .0254,
                ZR      = utils.str_m(self.q_fit_m2.ZR, d = 4),
                W0      = utils.str_m(self.q_fit_m2.W0, d = 4),
                D0      = utils.str_m(2 * self.q_fit.W0, d = 4),
                m2      = self.Z0_ZR_m2_fit[2],
                MML = (1 - (1 - (self.Z0_ZR_m2_fit[2] - 1)/4)**2) * 100,
                XL  = (1 - abs(self.q_fit_m2.overlap_LG(self.q_fit))**2) * 100,
            ))
            #https://www.dataray.com/blogs/dataray-blog/m-sup2-and-high-order-modes
            F.ax0.plot(
                Z_pts / .0254,
                self.Z0_ZR_m2_fit[2]**0.5 * 2 * 1e6 * self.q_fit_m2.propagate_distance(Z_pts).W,
                color = color_fit,
                label = fit_label_m2,
                ls = '--'
            )
            F.ax0.plot(
                Z_pts / .0254,
                2 * 1e6 * self.q_fit_m2.propagate_distance(Z_pts).W,
                color = color_fit,
                label = 'embedded HG00 in m2 fit',
                ls = ':'
            )

        if with_init:
            F.ax0.plot(
                Z_pts / .0254,
                2 * 1e6 * self.q_init.propagate_distance(Z_pts).W,
                color = color_init,
                label = 'Initial',
            )
        F.ax0.set_xlabel('Inches from reference')
        F.ax0.set_ylabel('2σ intensity\ndiameter[μm]')
        F.ax0.set_title('Beam Parameter Fit (at {0:.0f}nm)'.format(self.wavelength_m * 1e9))
        F.ax0.legend(loc = 'best')
        return F

