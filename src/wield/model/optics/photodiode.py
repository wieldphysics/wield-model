#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""
import collections
from .. import base

from . import vacuum


class Photodiode(base.OpticalObject):
    def __init__(self):
        super(Photodiode, self).__init__()

    @classmethod
    def visit_port_information(cls, manip):
        manip.gen_optical_port("+A", "A")
        return

    def visit_matrix_algorithm_DCAC(self, manip):
        manip.add_view("A!i")

        if manip.is_AC:
            kmatrix = vacuum.optical_quantum_noise_matrix(manip)
            manip.add_noise("A!o", kmatrix)
        return

    def visit_photodiode_linkage(self, view):
        return "A!i"


class PhotodiodeUnphysical(base.OpticalObject):
    def __init__(self):
        super(PhotodiodeUnphysical, self).__init__()

    def port_chain(self, p, pname):
        bmap = {
            "+A-t": (None, "+B"),
            "+B-t": (None, "+A"),
        }.get(pname, None)

        if bmap is not None:
            return bmap

        return super(PhotodiodeUnphysical, self).port_chain(p, pname)

    @classmethod
    def visit_port_information(cls, manip):
        manip.gen_optical_port("+A", "A")
        manip.gen_optical_port("+B", "B")
        manip.gen_signal_out_port("+Wpd", "Wpd")
        # manip.gen_signal_out_port('+WpdX', 'WpdX')
        # manip.gen_signal_out_port('+WpdY', 'WpdY')
        # manip.gen_signal_out_port('+Apd',  'Apd')
        # manip.gen_signal_out_port('+ApdX', 'ApdX')
        # manip.gen_signal_out_port('+ApdY', 'ApdY')
        # manip.gen_signal_out_port('+Qpd',  'Qpd')
        # manip.gen_signal_out_port('+QpdX', 'QpdX')
        # manip.gen_signal_out_port('+QpdY', 'QpdY')
        return

    def visit_photodiode_linkage(self, view):
        return "A!i"

    def visit_matrix_algorithm_DCAC(self, manip):
        manip.add_link("B!i", "A!o", 1)
        manip.add_link("A!i", "B!o", 1)
        manip.add_view("A!i")
        manip.add_view("Wpd!o")
        if manip.is_DC:
            manip.add_view_DC("A!i")
            manip.add_view_AC("A!i")
        self.build_lport_DCAC(manip, "Wpd!o")
        return

    def build_lport_DCAC(self, manip, lport_to):
        if lport_to == "Wpd!o":
            basis_fr = manip.link_basis("A!i")
            field_vec = manip.get_field_DC("A!i")
            if field_vec != 0:
                basis_to = manip.link_basis("Wpd!o")
                if manip.is_DC:
                    projection_full = 0
                    for fK in basis_to["frequency"].enumerated:
                        field_to_p = base.KeyMatrixDiff(
                            stR=(
                                basis_to["frequency"],
                                base.kg_twophoton_pm,
                            ),
                            kmatrix={
                                (
                                    fK,
                                    "+",
                                ): {(): [[1]]}
                            },
                            build=True,
                        )
                        field_to_m = base.KeyMatrixDiff(
                            stR=(
                                basis_to["frequency"],
                                base.kg_twophoton_pm,
                            ),
                            kmatrix={
                                (
                                    fK,
                                    "-",
                                ): {(): [[1]]}
                            },
                            build=True,
                        )
                        km_p = collections.defaultdict(dict)
                        km_m = collections.defaultdict(dict)
                        for fK_fr in basis_fr["frequency"].enumerated:
                            fK_diff = fK_fr - fK
                            if not manip.optical_frequency_allowed(fK_diff):
                                continue
                            # TODO, demod phasings
                            # demod complex output
                            km_p[(fK_fr, "+")][(fK_diff, "-")] = [[1]]
                            km_m[(fK_fr, "-")][(fK_diff, "+")] = [[1]]
                            ##demod I-phase
                            # km[(fK_fr, '+')][(fK_diff, '-')] = [[1/2]]
                            # km[(fK_diff, '+')][(fK_fr, '-')] = [[1/2]]
                            ##demod Q-phase
                            # km[(fK_fr, '+')][(fK_diff, '-')] = [[1j/2]]
                            # km[(fK_diff, '+')][(fK_fr, '-')] = [[-1j/2]]

                        inner_prod_p = base.KeyMatrixSame(
                            stR=(
                                basis_fr["frequency"],
                                base.kg_twophoton_pm,
                            ),
                            stC=(
                                basis_fr["frequency"],
                                base.kg_twophoton_pm,
                            ),
                            kmatrix=km_p,
                            build=True,
                        )
                        inner_prod_m = base.KeyMatrixSame(
                            stR=(
                                basis_fr["frequency"],
                                base.kg_twophoton_pm,
                            ),
                            stC=(
                                basis_fr["frequency"],
                                base.kg_twophoton_pm,
                            ),
                            kmatrix=km_m,
                            build=True,
                        )
                        projection_p = field_to_p @ field_vec.T @ inner_prod_p
                        projection_m = field_to_m @ field_vec.T @ inner_prod_m
                        # TODO, this incremental sum is N**2 efficiency
                        if projection_full == 0:
                            projection_full = projection_p + projection_m
                        else:
                            projection_full = (
                                projection_p + projection_m
                            ) + projection_full
                    projection = projection_full
                elif manip.is_AC:
                    projection_full = 0
                    for fK in basis_to["frequency"].enumerated:
                        field_to_p = base.KeyMatrixDiff(
                            stR=(
                                basis_to["frequency"],
                                base.kg_twophoton_pmAC,
                            ),
                            kmatrix={
                                (
                                    fK,
                                    "+AC",
                                ): {(): [[1]]}
                            },
                            build=True,
                        )
                        field_to_m = base.KeyMatrixDiff(
                            stR=(
                                basis_to["frequency"],
                                base.kg_twophoton_pmAC,
                            ),
                            kmatrix={
                                (
                                    fK,
                                    "-AC",
                                ): {(): [[1]]}
                            },
                            build=True,
                        )
                        km_p = collections.defaultdict(dict)
                        km_m = collections.defaultdict(dict)
                        for fK_fr in basis_fr["frequency"].enumerated:
                            fK_diff = fK_fr - fK
                            if not manip.optical_frequency_allowed(fK_diff):
                                continue
                            km_p[(fK_diff, "+")][(fK_fr, "-AC")] = [[1 / 2]]
                            # km_p[(fK_fr, '+')][(fK_diff, '-AC')] = [[1/2]]
                            km_m[(fK_diff, "-")][(fK_fr, "+AC")] = [[1 / 2]]
                            # km_m[(fK_fr, '-')][(fK_diff, '+AC')] = [[1/2]]

                        inner_prod_m = base.KeyMatrixSame(
                            stR=(
                                basis_fr["frequency"],
                                base.kg_twophoton_pm,
                            ),
                            stC=(
                                basis_fr["frequency"],
                                base.kg_twophoton_pmAC,
                            ),
                            kmatrix=km_m,
                            build=True,
                        )
                        inner_prod_p = base.KeyMatrixSame(
                            stR=(
                                basis_fr["frequency"],
                                base.kg_twophoton_pm,
                            ),
                            stC=(
                                basis_fr["frequency"],
                                base.kg_twophoton_pmAC,
                            ),
                            kmatrix=km_p,
                            build=True,
                        )
                        projection_p = field_to_p @ field_vec.T @ inner_prod_p
                        projection_m = field_to_m @ field_vec.T @ inner_prod_m
                        # TODO, this incremental sum is N**2 efficiency
                        if projection_full == 0:
                            projection_full = projection_p + projection_m
                        else:
                            projection_full = (
                                projection_p + projection_m
                            ) + projection_full
                    projection = projection_full
                manip.add_link("A!i", "Wpd!o", projection)
        else:
            raise RuntimeError("Boo")

    def visit_mode_matching_linkage(self, manip):
        manip.add_link("B!i", "A!o", None)
        manip.add_link("A!i", "B!o", None)
