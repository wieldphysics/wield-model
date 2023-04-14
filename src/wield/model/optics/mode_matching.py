#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import numpy as np
from .. import base
from .alm.beam_fit import QFit


class ModeMatchingTargetBase(base.OpticalObject):
    # disallow reference aliasing so that the naming scheme for waypoints
    # is unique
    _allow_multiple_parents = False
    pass


class Cavity(ModeMatchingTargetBase):
    def __init__(self):
        super(Cavity, self).__init__()
        with self._internal():
            # if name is None, then use the object path name as the name
            self["name"] = None

            self["waypoints"] = None

    def visit_mode_matching_targets(self, manip):
        name = manip.p["name"]
        if name is None:
            name = manip.path()
            # remove trailing path symbol
            if name[-1] == "/":
                name = name[:-1]
        waypoints = manip.p["waypoints"]
        manip.cavity_add(name, waypoints)
        return


class Target(ModeMatchingTargetBase):
    def __init__(self):
        super(Target, self).__init__()
        with self._internal():
            # if name is None, then use the object path name as the name
            self["name"] = None
            self["waypoints"] = None

            self["wavelength"] = None

            self["distance[m]"] = None
            self["distance[in]"] = None

            self["q"] = None

            @self.deco_one_one("q")
            def qX(q):
                return q

            @self.deco_one_one("q")
            def qY(q):
                return q

        return

    def visit_mode_matching_targets(self, manip):
        name = manip.p["name"]
        if name is None:
            name = manip.path()
            # remove trailing path symbol
            if name[-1] == "/":
                name = name[:-1]
        waypoints = manip.p["waypoints"]
        wavelength = manip.p["wavelength"]
        qX = manip.p["qX"]
        qY = manip.p["qY"]
        # print('qX', qX)
        manip.target_add(
            name,
            waypoints,
            qX=qX,
            qY=qY,
            wavelength=wavelength,
        )
        return


class TargetMeasurement(ModeMatchingTargetBase):
    def __init__(self):
        super(TargetMeasurement, self).__init__()
        with self._internal():
            # if name is None, then use the object path name as the name
            self["name"] = None
            self["waypoints"] = None

            self["distance[in]"] = None
            self["distance[m]"] = None

            self["diameter[in]"] = None
            self["diameterX[in]"] = None
            self["diameterY[in]"] = None

            self["reversed"] = False

            @self.deco_many_one(
                dict(
                    wavelength_m="wavelength",
                    Z_m="distance[m]",
                    D_m="diameterX[m]",
                )
            )
            def qfitX(wavelength_m, Z_m, D_m):
                return QFit(
                    wavelength_m=wavelength_m,
                    Z_in=Z_m / 0.0254,
                    D_um=1e6 * D_m,
                )

            @self.deco_many_one(
                dict(
                    wavelength_m="wavelength",
                    Z_m="distance[m]",
                    D_m="diameterY[m]",
                )
            )
            def qfitY(wavelength_m, Z_m, D_m):
                return QFit(
                    wavelength_m=wavelength_m,
                    Z_in=Z_m / 0.0254,
                    D_um=1e6 * D_m,
                )

    def visit_mode_matching_targets(self, manip):
        # TODO, need to do some additional lookups to allow abstract wavelength typing
        name = manip.p["name"]
        if name is None:
            name = manip.path()
            # remove trailing path symbol
            if name[-1] == "/":
                name = name[:-1]
        waypoints = manip.p["waypoints"]
        wavelength = manip.p["wavelength"]

        if not manip.p["reversed"]:
            manip.target_add(
                name,
                waypoints,
                qX=manip.p["qfitX"].q_fit,
                qY=manip.p["qfitY"].q_fit,
                wavelength=wavelength,
            )
        else:
            manip.target_add(
                name,
                waypoints,
                qX=manip.p["qfitX"].q_fit.reversed(),
                qY=manip.p["qfitY"].q_fit.reversed(),
                wavelength=wavelength,
            )
        return
