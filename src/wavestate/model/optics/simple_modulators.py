#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""
from .. import base


class SimpleModulator(base.OpticalObject):
    mod_type = None

    def __init__(self):
        super(SimpleModulator, self).__init__()
        with self._internal():
            # this should either be a name or a F-dictionary
            self["frequency"] = None
            self["index"] = None

    def port_chain(self, p, pname):
        bmap = {
            "+A-t": (None, "+B"),
            "+B-t": (None, "+A"),
        }.get(pname, None)

        if bmap is not None:
            return bmap

        return super(SimplePhaseModulator, self).port_chain(p, pname)

    @classmethod
    def visit_port_information(cls, manip):
        manip.gen_optical_port("+A", "A")
        manip.gen_optical_port("+B", "B")
        return

    def visit_matrix_algorithm_DCAC(self, manip):
        settings = manip.optical_settings()

        # expects basis to be the same in A!i, A!o, B!i, B!o,
        basis = manip.link_basis("A!i")

        index = manip.p["index"]
        if index is None:
            raise RuntimeError("Must specify modulation index")
        frequency = manip.p["frequency"]

        if frequency is None:
            raise RuntimeError("Must specify modulation frequency")

        fK_mod = manip.parameter_to_fk(frequency)

        if self.mod_type == "phase":
            mod_factor = 1j
            upper = True
            lower = True
        elif self.mod_type == "amplitude":
            mod_factor = 1
            upper = True
            lower = True
        elif self.mod_type == "SSB upper":
            mod_factor = 2
            upper = True
            lower = False
        elif self.mod_type == "SSB lower":
            mod_factor = 2
            upper = False
            lower = True
        else:
            raise RuntimeError("Unrecognized modulation type")

        # TODO, use bessel functions or otherwise conserve energy
        # ideally numeric integration..
        km = dict()
        for fK in basis["frequency"].enumerated:
            # vdict carries the indexs for the columns (from)
            vdict = km.setdefault((fK,), dict())

            vdict[(fK,)] = [[1 - abs(index) ** 2 / 4]]
            # the lower sideband moving up
            if upper:
                fd_from = fK - fK_mod
                if manip.optical_frequency_allowed(fd_from):
                    vdict[(fd_from,)] = [[mod_factor * index / 2]]

            # the upper sideband moving down
            if lower:
                fd_from = fK + fK_mod
                if manip.optical_frequency_allowed(fd_from):
                    vdict[(fd_from,)] = [[mod_factor * index.conjugate() / 2]]

        kmatrix = matrix.KeyMatrixSame(
            st=(basis["frequency"],),
            kmatrix=km,
            build=True,
            check=manip.check_build,
        )

        manip.add_link("B!i", "A!o", kmatrix, lowering_only=True)
        manip.add_link("A!i", "B!o", kmatrix, lowering_only=True)
        return

    def visit_mode_matching_linkage(self, manip):
        manip.add_link("B!i", "A!o", None)
        manip.add_link("A!i", "B!o", None)


class SimplePhaseModulator(SimpleModulator):
    mod_type = "phase"


class SimpleAmplitudeModulator(SimpleModulator):
    mod_type = "amplitude"


class SimpleSSBUpperModulator(SimpleModulator):
    mod_type = "SSB upper"


class SimpleSSBLowerModulator(SimpleModulator):
    mod_type = "SSB lower"
