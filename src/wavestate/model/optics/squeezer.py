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

from .. import base


class SimpleSqueezer(base.OpticalObject):
    mod_type = None

    def __init__(self):
        super(SimpleSqueezer, self).__init__()
        with self._internal():
            # this should either be a name or a F-dictionary
            self["frequency"] = 0
            self["gain"] = 1j
            self["wavelength"] = None
            self["index"] = None

    def port_chain(self, p, pname):
        bmap = {
            "+A-t": (None, "+B"),
            "+B-t": (None, "+A"),
        }.get(pname, None)

        if bmap is not None:
            return bmap

        return super(SimpleSqueezer, self).port_chain(p, pname)

    @classmethod
    def visit_port_information(cls, manip):
        manip.gen_optical_port("+A", "A")
        manip.gen_optical_port("+B", "B")
        return

    def visit_matrix_algorithm_DCAC(self, manip):
        settings = manip.optical_settings()

        # expects basis to be the same in A!i, A!o, B!i, B!o,
        basis = manip.link_basis("A!i")

        gain = manip.p["gain"]
        gain_c = gain.conjugate()
        gain_d = (1 + gain * gain_c) ** 0.5

        fK_center = manip.parameter_to_fk(manip.p["frequency"])
        wK_center = manip.parameter_to_fk(manip.p["wavelength"], allow_default=True)

        if manip.is_DC:
            q_p = "+"
            q_n = "-"
        elif manip.is_AC:
            q_p = "+AC"
            q_n = "-AC"

        km = collections.defaultdict(dict)
        for wK in basis["wavelength"].enumerated:
            for fK in basis["frequency"].enumerated:
                # build the diagonal

                if wK != wK_center:
                    km[(wK, fK, q_p)][(wK, fK, q_p)] = [[1]]
                    km[(wK, fK, q_n)][(wK, fK, q_n)] = [[1]]
                    continue

                fK_mirror = fK - fK_center
                if not manip.optical_frequency_allowed(fK_mirror):
                    km[(wK, fK, q_p)][(wK, fK, q_p)] = [[1]]
                    km[(wK, fK, q_n)][(wK, fK, q_n)] = [[1]]
                    print("SQUEEZER WARNING")
                    continue

                km[(wK, fK, q_p)][(wK, fK, q_p)] = [[gain_d]]
                km[(wK, fK, q_n)][(wK, fK, q_n)] = [[gain_d]]

                km[(wK, fK, q_p)][(wK, fK_mirror, q_n)] = [[gain]]
                km[(wK, fK_mirror, q_n)][(wK, fK, q_p)] = [[gain_c]]

        kmatrix = matrix.KeyMatrixSame(
            st=(basis["wavelength"], basis["frequency"], basis["quantum"]),
            kmatrix=km,
            build=True,
            check=manip.check_build,
        )

        manip.add_link("B!i", "A!o", kmatrix)
        manip.add_link("A!i", "B!o", kmatrix)
        return

    def visit_mode_matching_linkage(self, manip):
        manip.add_link("B!i", "A!o", None)
        manip.add_link("A!i", "B!o", None)
