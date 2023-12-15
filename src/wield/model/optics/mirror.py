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
import collections

from wield.utilities.np import matrix_stack

from .. import base
from . import alm


class MirrorBase(base.OpticalObject):
    def __init__(self):
        super(MirrorBase, self).__init__()
        with self._internal():
            self["AC_active"] = True
            self["beam_diameter[m]"] = None
            self["ROC[m]"] = None
            self["ROC_B[m]"] = None
            self["AOI[deg]"] = 0
            self["defocus[D]"] = 0
            self["defocusX[D]"] = 0
            self["defocusY[D]"] = 0

            self["annotate"] = "Mirror"

            self["substrate"] = "fused_silica"
            self["depth[m]"] = 0
            self["depth[m].lower_bound"] = 0

            @self.deco_one_one("AOI[deg]")
            def beamsplitter_ports(AOI):
                if AOI != 0:
                    return True
                else:
                    return False

        self.values_settable(
            [
                "AC_active",
                "substrate",
                "beam_diameter[m]",
                "ROC[m]",
                "ROC_B[m]",
                "depth[m]",
                "AOI[deg]",
                "defocus[D]",
                "defocusX[D]",
                "defocusY[D]",
            ]
        )

    def port_chain(self, p, pname):
        # being a beamsplitter doesn't change the chains
        bmap = {
            "+A-t": (None, "+B"),
            "+B-t": (None, "+A"),
            "+A-r": (None, "+A"),
            "+B-r": (None, "+B"),
            "+A1-t": (None, "+B1"),
            "+B1-t": (None, "+A1"),
            "+A2-t": (None, "+B2"),
            "+B2-t": (None, "+A2"),
            "+A1-r": (None, "+A2"),
            "+B1-r": (None, "+B2"),
            "+A2-r": (None, "+A1"),
            "+B2-r": (None, "+B1"),
        }.get(pname, None)

        if bmap is not None:
            return bmap

        return super().port_chain(p, pname)

    def port_forward(self, p, pname):
        if not p["beamsplitter_ports"]:
            bmap = {
                "+A": None,
                "+B": None,
                "+A1": "+A",
                "+B1": "+B",
                "+A2": "+A",
                "+B2": "+B",
            }.get(pname, None)
        else:
            bmap = {
                "+A": False,
                "+B": False,
                "+A1": None,
                "+B1": None,
                "+A2": None,
                "+B2": None,
            }.get(pname, None)
        if bmap is not None:
            return bmap
        return super().port_forward(p, pname)

    @classmethod
    def visit_port_information(cls, manip):
        # TODO, only generate the ports if the mechanical connection
        # is not generated?
        # (linkages may always be generated or be implicit)

        if manip.p["AC_active"]:
            # the AC optomechanical ports
            manip.gen_signal_out_port("+Fl", "Fl")
            manip.gen_signal_out_port("+Fp", "Fp")
            manip.gen_signal_out_port("+Fy", "Fy")
            manip.gen_signal_io_port("+Dl", "Dl")
            manip.gen_signal_io_port("+Dp", "Dp")
            manip.gen_signal_io_port("+Dy", "Dy")

        if not manip.p["beamsplitter_ports"]:
            manip.gen_optical_port("+A", "A")
            manip.gen_optical_port("+B", "B")
            if manip.p["has_loss"]:
                manip.gen_optical_port("+A_loss", "A_loss")
                manip.gen_optical_port("+B_loss", "B_loss")
        else:
            manip.gen_optical_port("+A1", "A1")
            manip.gen_optical_port("+B1", "B1")
            manip.gen_optical_port("+A2", "A2")
            manip.gen_optical_port("+B2", "B2")
            if manip.p["has_loss"]:
                manip.gen_optical_port("+A1_loss", "A1_loss")
                manip.gen_optical_port("+B1_loss", "B1_loss")
                manip.gen_optical_port("+A2_loss", "A2_loss")
                manip.gen_optical_port("+B2_loss", "B2_loss")
        return

    def matrix_RTL(self, manip):
        raise NotImplementedError()

    def visit_matrix_algorithm_DCAC(self, manip):
        kmatrix_R, kmatrix_T, kmatrix_L = self.matrix_RTL(manip)

        if not manip.p["beamsplitter_ports"]:
            manip.add_link("B!i", "A!o", kmatrix_T, lowering_only=True)
            manip.add_link("A!i", "A!o", kmatrix_R, lowering_only=True)
            manip.add_link("A!i", "B!o", kmatrix_T, lowering_only=True)
            manip.add_link("B!i", "B!o", -kmatrix_R, lowering_only=True)
            if kmatrix_L is not None:
                manip.add_link("A_loss!i", "A!o", kmatrix_L, lowering_only=True)
                manip.add_link("B_loss!i", "B!o", kmatrix_L, lowering_only=True)
                # TODO? this requires t and r from the inputs to be correct
                # manip.add_link('A!i', 'A_loss!o', -kmatrix_L, lowering_only = True)
                # manip.add_link('B!i', 'B_loss!o', -kmatrix_L, lowering_only = True)

            if manip.p["AC_active"]:
                # the AC optomechanical ports
                if manip.is_DC:
                    manip.add_view_AC("A!i")
                    manip.add_view_AC("A!o")
                    manip.add_view_AC("B!i")
                    manip.add_view_AC("B!o")
                if manip.is_AC:
                    manip.add_link("Dl!i", "Dl!o", 1)
                    manip.add_link("Dp!i", "Dp!o", 1)
                    manip.add_link("Dy!i", "Dy!o", 1)
                    manip.add_drive("Dl!i")
                    manip.add_drive("Dp!i")
                    manip.add_drive("Dy!i")

                    self.build_lport_DCAC(manip, "Dl!i")
                    self.build_lport_DCAC(manip, "Fl!o")
        else:
            manip.add_link("B1!i", "A1!o", kmatrix_T, lowering_only=True)
            manip.add_link("A2!i", "A1!o", -kmatrix_R, lowering_only=True)
            manip.add_link("A1!i", "B1!o", kmatrix_T, lowering_only=True)
            manip.add_link("B2!i", "B1!o", kmatrix_R, lowering_only=True)
            manip.add_link("B2!i", "A2!o", kmatrix_T, lowering_only=True)
            manip.add_link("A1!i", "A2!o", -kmatrix_R, lowering_only=True)
            manip.add_link("A2!i", "B2!o", kmatrix_T, lowering_only=True)
            manip.add_link("B1!i", "B2!o", kmatrix_R, lowering_only=True)
            if kmatrix_L is not None:
                manip.add_link("A1_loss!i", "A1!o", kmatrix_L, lowering_only=True)
                manip.add_link("B1_loss!i", "B1!o", kmatrix_L, lowering_only=True)
                manip.add_link("A2_loss!i", "A2!o", kmatrix_L, lowering_only=True)
                manip.add_link("B2_loss!i", "B2!o", kmatrix_L, lowering_only=True)
                # TODO? this requires t and r from the inputs to be correct
                # manip.add_link('A1!i', 'A1_loss!o', -kmatrix_L, lowering_only = True)
                # manip.add_link('B1!i', 'B1_loss!o', -kmatrix_L, lowering_only = True)
                # manip.add_link('A2!i', 'A2_loss!o', -kmatrix_L, lowering_only = True)
                # manip.add_link('B2!i', 'B2_loss!o', -kmatrix_L, lowering_only = True)
            if manip.p["AC_active"]:
                # the AC optomechanical ports
                if manip.is_DC:
                    manip.add_view_AC("A1!i")
                    manip.add_view_AC("A1!o")
                    manip.add_view_AC("B1!i")
                    manip.add_view_AC("B1!o")
                    manip.add_view_AC("A2!i")
                    manip.add_view_AC("A2!o")
                    manip.add_view_AC("B2!i")
                    manip.add_view_AC("B2!o")
        return

    def build_lport_DCAC(self, manip, lport_build):
        if manip.is_DC:
            # TODO
            # DC forces not yet supported
            return
        # expects basis to be the same in all links
        basis = manip.link_basis("A!i")
        sgn_map = {
            "A!i": 1,
            "B!i": -1,
            "A!o": 1,
            "B!o": -1,
        }
        for lport_fr in ["A!i", "B!i", "A!o", "B!o"]:
            sgn = sgn_map[lport_fr]
            field_vec = manip.get_field_DC(lport_fr)
            if field_vec == 0:
                continue
            if lport_build == "Dl!i":
                basis_fr = manip.link_basis("Dl!i")
                for Wk, wnval in manip.basis_wavenumbers(with_keys=True):
                    # displacement drive
                    Dl_prod_left = base.KeyMatrixSame(
                        stR=(
                            basis["wavenumber"],
                            base.kg_twophoton_pmAC,
                        ),
                        stC=(
                            basis["wavenumber"],
                            base.kg_twophoton_pm,
                        ),
                        kmatrix={
                            (
                                Wk,
                                "+AC",
                            ): {(Wk, "+"): [[+1j * wnval / 2]]},
                            (
                                Wk,
                                "-AC",
                            ): {(Wk, "-"): [[-1j * wnval / 2]]},
                        },
                        build=True,
                    )
                    # TODO, use DictKey construction to ensure no missing basis
                    Dl_prod_right = base.KeyMatrixDiff(
                        stR=(),
                        stC=(
                            basis_fr["frequency"],
                            base.kg_twophoton_pmAC,
                        ),
                        kmatrix={
                            (): {
                                (
                                    base.FrequencyKey({}),
                                    "+AC",
                                ): [[1]]
                            },
                        },
                        build=True,
                    )
                    outer_gen = Dl_prod_left @ field_vec @ Dl_prod_right
                    from icecream import ic

                    ic(outer_gen.kmatrix)
                    manip.add_link("Dl!i", lport_fr, sgn * outer_gen)

            elif lport_build == "Fl!o":
                # force signal
                Fl_prod_right = base.KeyMatrixSame(
                    stR=(base.kg_twophoton_pm,),
                    stC=(base.kg_twophoton_pmAC,),
                    kmatrix={
                        ("-",): {("+AC",): [[1 / manip.constants.c_m_s]]},
                        ("+",): {("-AC",): [[1 / manip.constants.c_m_s]]},
                    },
                    build=True,
                )

                # TODO, use DictKey construction to ensure no missing basis
                basis_to = manip.link_basis("Fl!o")
                Fl_prod_left = base.KeyMatrixDiff(
                    stC=(),
                    stR=(
                        basis_to["frequency"],
                        base.kg_twophoton_pmAC,
                    ),
                    kmatrix={
                        (
                            base.FrequencyKey({}),
                            "+AC",
                        ): {(): [[1]]},
                    },
                    build=True,
                )
                outer_gen = Fl_prod_left @ field_vec.T @ Fl_prod_right
                manip.add_link(lport_fr, "Fl!o", sgn * outer_gen)
            else:
                assert False
        return

    def visit_mode_matching_linkage(self, manip):
        if not manip.p["beamsplitter_ports"]:
            manip.add_link("A!i", "A!o", 0)
            manip.add_link("A!i", "B!o", manip.p["depth[m]"])
            manip.add_link("B!i", "A!o", manip.p["depth[m]"])
            manip.add_link("B!i", "B!o", 2 * manip.p["depth[m]"])
        else:
            manip.add_link("A1!i", "A2!o", 0)
            manip.add_link("A1!i", "B1!o", manip.p["depth[m]"])
            manip.add_link("B1!i", "A1!o", manip.p["depth[m]"])
            manip.add_link("B1!i", "B2!o", 2 * manip.p["depth[m]"])
            manip.add_link("A2!i", "A1!o", 0)
            manip.add_link("A2!i", "B2!o", manip.p["depth[m]"])
            manip.add_link("B2!i", "A2!o", manip.p["depth[m]"])
            manip.add_link("B2!i", "B1!o", 2 * manip.p["depth[m]"])

    def visit_mode_matching_transport(self, manip):
        n_A = 1
        n_B = 1
        n_I = manip.IOR_n(manip.p["substrate"])

        ROC_A_m = manip.p["ROC[m]"]
        depth_m = manip.p["depth[m]"]
        ROC_B_m = manip.p["ROC_B[m]"]
        AOI_rad = manip.p["AOI[deg]"] * np.pi / 180

        # TODO, linear approximation
        depth_m = depth_m / np.cos(AOI_rad * n_A / n_I)

        lport_fr_to = manip.lport_fr, manip.lport_to

        if lport_fr_to in [("A!i", "B!o"), ("A1!i", "B1!o"), ("A2!i", "B2!o")]:
            # AB transmission
            matrixAI_X = alm.interface_ROC_AOI_X(
                ROC_A_m, n_from=n_A, n_to=n_I, AOI_rad=AOI_rad
            )
            matrixAI_Y = alm.interface_ROC_AOI_Y(
                ROC_A_m, n_from=n_A, n_to=n_I, AOI_rad=AOI_rad
            )
            matrixII = matrix_stack([[1, depth_m], [0, 1]])
            matrixIB_X = alm.interface_ROC_AOI_X(
                ROC_B_m, n_from=n_I, n_to=n_B, AOI_rad=AOI_rad, neg=True
            )
            matrixIB_Y = alm.interface_ROC_AOI_Y(
                ROC_B_m, n_from=n_I, n_to=n_B, AOI_rad=AOI_rad, neg=True
            )

            def inc_builder(z):
                return matrix_stack([[1, z], [0, 1]])

            manip.set_Xincremental(
                [
                    (0, None, matrixAI_X),
                    (depth_m, inc_builder, matrixII),
                    (0, None, matrixIB_X),
                ]
            )

            manip.set_Yincremental(
                [
                    (0, None, matrixAI_Y),
                    (depth_m, inc_builder, matrixII),
                    (0, None, matrixIB_Y),
                ]
            )
            # the P-builders are for fast optimization solving
            def p_builder_X(p):
                ROC_A_m = p["ROC[m]"]
                depth_m = p["depth[m]"]
                ROC_B_m = p["ROC_B[m]"]
                AOI_rad = p["AOI[deg]"] * np.pi / 180
                depth_m = depth_m / np.cos(AOI_rad * n_A / n_I)

                matrixAI = alm.interface_ROC_AOI_X(
                    ROC_A_m, n_from=n_A, n_to=n_I, AOI_rad=AOI_rad
                )
                matrixII = matrix_stack([[1, depth_m], [0, 1]])
                matrixIB = alm.interface_ROC_AOI_X(
                    ROC_B_m, n_from=n_I, n_to=n_B, AOI_rad=AOI_rad, neg=True
                )
                return matrixIB @ matrixII @ matrixAI

            manip.set_Xpropagator(p_builder_X)

            def p_builder_Y(p):
                ROC_A_m = p["ROC[m]"]
                depth_m = p["depth[m]"]
                ROC_B_m = p["ROC_B[m]"]
                AOI_rad = p["AOI[deg]"] * np.pi / 180
                depth_m = depth_m / np.cos(AOI_rad * n_A / n_I)

                matrixAI = alm.interface_ROC_AOI_Y(
                    ROC_A_m, n_from=n_A, n_to=n_I, AOI_rad=AOI_rad
                )
                matrixII = matrix_stack([[1, depth_m], [0, 1]])
                matrixIB = alm.interface_ROC_AOI_Y(
                    ROC_B_m, n_from=n_I, n_to=n_B, AOI_rad=AOI_rad, neg=True
                )
                return matrixIB @ matrixII @ matrixAI

            manip.set_Ypropagator(p_builder_Y)
            manip.set_Zpropagator(
                {"depth[m]": (1 / np.cos(AOI_rad * n_A / n_I), "AOI[deg]")}
            )
        elif lport_fr_to in [("B!i", "A!o"), ("B1!i", "A1!o"), ("B2!i", "A2!o")]:
            # BA transmission
            matrixBI_X = alm.interface_ROC_AOI_X(
                ROC_B_m, n_from=n_B, n_to=n_I, AOI_rad=AOI_rad
            )
            matrixBI_Y = alm.interface_ROC_AOI_Y(
                ROC_B_m, n_from=n_B, n_to=n_I, AOI_rad=AOI_rad
            )
            matrixII = matrix_stack([[1, depth_m], [0, 1]])
            matrixIA_X = alm.interface_ROC_AOI_X(
                ROC_A_m, n_from=n_I, n_to=n_A, AOI_rad=AOI_rad, neg=True
            )
            matrixIA_Y = alm.interface_ROC_AOI_Y(
                ROC_A_m, n_from=n_I, n_to=n_A, AOI_rad=AOI_rad, neg=True
            )

            def inc_builder(z):
                return matrix_stack([[1, z], [0, 1]])

            manip.set_Xincremental(
                [
                    (0, None, matrixBI_X),
                    (depth_m, inc_builder, matrixII),
                    (0, None, matrixIA_X),
                ]
            )

            manip.set_Yincremental(
                [
                    (0, None, matrixBI_Y),
                    (depth_m, inc_builder, matrixII),
                    (0, None, matrixIA_Y),
                ]
            )

            # the P-builders are for fast optimization solving
            def p_builder_X(p):
                ROC_A_m = p["ROC[m]"]
                depth_m = p["depth[m]"]
                ROC_B_m = p["ROC_B[m]"]
                AOI_rad = p["AOI[deg]"] * np.pi / 180
                depth_m = depth_m / np.cos(AOI_rad * n_A / n_I)

                matrixBI = alm.interface_ROC_AOI_X(
                    ROC_B_m, n_from=n_B, n_to=n_I, AOI_rad=AOI_rad
                )
                matrixII = matrix_stack([[1, depth_m], [0, 1]])
                matrixIA = alm.interface_ROC_AOI_X(
                    ROC_A_m, n_from=n_I, n_to=n_A, AOI_rad=AOI_rad, neg=True
                )
                return matrixIA @ matrixII @ matrixBI

            def p_builder_Y(p):
                ROC_A_m = p["ROC[m]"]
                depth_m = p["depth[m]"]
                ROC_B_m = p["ROC_B[m]"]
                AOI_rad = p["AOI[deg]"] * np.pi / 180
                depth_m = depth_m / np.cos(AOI_rad * n_A / n_I)

                matrixBI = alm.interface_ROC_AOI_Y(
                    ROC_B_m, n_from=n_B, n_to=n_I, AOI_rad=AOI_rad
                )
                matrixII = matrix_stack([[1, depth_m], [0, 1]])
                matrixIA = alm.interface_ROC_AOI_Y(
                    ROC_A_m, n_from=n_I, n_to=n_A, AOI_rad=AOI_rad, neg=True
                )
                return matrixIA @ matrixII @ matrixBI

            manip.set_Xpropagator(p_builder_X)
            manip.set_Ypropagator(p_builder_Y)

            def p_builderZ(p):
                depth_m = p["depth[m]"]
                AOI_rad = p["AOI[deg]"] * np.pi / 180
                depth_m = depth_m / np.cos(AOI_rad * n_A / n_I)
                return depth_m

            manip.set_Zpropagator(
                {"depth[m]": (1 / np.cos(AOI_rad * n_A / n_I), "AOI[deg]")}
            )
        elif lport_fr_to in [("A!i", "A!o"), ("A1!i", "A2!o"), ("A2!i", "A1!o")]:
            # AA reflection
            def p_builder_X(p):
                roc_m = p["ROC[m]"]
                aoi_rad = p["AOI[deg]"] / 180 * np.pi
                defocus_D = p["defocus[D]"] + p["defocusX[D]"]
                M_D = matrix_stack([[1, 0], [defocus_D, 1]])
                return M_D @ alm.REFL_ROC_X(roc_m, aoi_rad)

            def p_builder_Y(p):
                roc_m = p["ROC[m]"]
                aoi_rad = p["AOI[deg]"] / 180 * np.pi
                defocus_D = p["defocus[D]"] + p["defocusY[D]"]
                M_D = matrix_stack([[1, 0], [defocus_D, 1]])
                return M_D @ alm.REFL_ROC_Y(roc_m, aoi_rad)

            manip.set_Xpropagator(p_builder_X)
            manip.set_Ypropagator(p_builder_Y)
            manip.set_Zpropagator()

            manip.set_Xincremental([(0, None, p_builder_X(manip.p))])
            manip.set_Yincremental([(0, None, p_builder_Y(manip.p))])

            # TODO, determine signs
            manip.set_Xshifts('yaw[rad]', matrix_stack([[0], [2]]))
            manip.set_Yshifts('pitch[rad]', matrix_stack([[0], [2]]))
        elif lport_fr_to in [("B!i", "B!o"), ("B1!i", "B2!o"), ("B2!i", "B1!o")]:
            # BB reflection
            # TODO, include astigmatism and defocus for X and Y
            defocus_D = manip.p["defocus[D]"]
            M_D = matrix_stack([[1, 0], [defocus_D, 1]])
            matrixBI = alm.interface_ROC(ROC_B_m, n_from=n_B, n_to=n_I)
            matrixII = matrix_stack([[1, depth_m], [0, 1]])
            if ROC_A_m is not None:
                # TODO add AOI effects
                matrixR = matrix_stack([[1, 0], [-2 / ROC_A_m, 1]])
            else:
                matrixR = matrix_stack([[1, 0], [0, 1]])
            matrixIB = alm.interface_ROC(ROC_B_m, n_from=n_I, n_to=n_B, neg=True)

            def inc_builder(z):
                return matrix_stack([[1, z], [0, 1]])

            manip.set_XYincremental(
                [
                    (0, None, matrixBI),
                    (depth_m, inc_builder, matrixII),
                    (0, None, matrixR),
                    (depth_m, inc_builder, matrixII),
                    (0, None, M_D @ matrixIB),
                ]
            )

            # the P-builders are for fast optimization solving
            def p_builder(p):
                ROC_A_m = manip.p["ROC[m]"]
                depth_m = manip.p["depth[m]"]
                ROC_B_m = manip.p["ROC_B[m]"]
                defocus_D = p["defocus[D]"]
                depth_m = depth_m / np.cos(AOI_rad * n_A / n_I)

                M_D = matrix_stack([[1, 0], [defocus_D, 1]])
                matrixBI = alm.interface_ROC(ROC_B_m, n_from=n_B, n_to=n_I)
                matrixII = matrix_stack([[1, depth_m], [0, 1]])
                if ROC_A_m is not None:
                    # TODO add AOI effects
                    matrixR = matrix_stack([[1, 0], [-2 / ROC_A_m, 1]])
                else:
                    matrixR = matrix_stack([[1, 0], [0, 1]])
                matrixIB = alm.interface_ROC(ROC_B_m, n_from=n_I, n_to=n_B, neg=True)
                return M_D @ matrixIB @ matrixII @ matrixR @ matrixII @ matrixBI

            manip.set_XYpropagator(p_builder)

            def p_builderZ(p):
                depth_m = p["depth[m]"]
                AOI_rad = p["AOI[deg]"] * np.pi / 180
                depth_m = depth_m / np.cos(AOI_rad * n_A / n_I)
                return 2 * depth_m

            manip.set_Zpropagator(
                {"depth[m]": (2 / np.cos(AOI_rad * n_A / n_I), "AOI[deg]")}
            )

            # TODO, determine signs
            manip.set_Xshifts('yaw[rad]', matrix_stack([[0], [-2]]))
            manip.set_Yshifts('pitch[rad]', matrix_stack([[0], [-2]]))

        return

    def visit_mm_anno_description(self, pbg, view, descB):
        desc = []

        def view_add(name, default, name2=None, transform=lambda x: "{:.3f}".format(x)):
            val = view[name]
            if val == default:
                return
            if name2 is None:
                name2 = name
            desc.append("{}={}".format(name2, transform(val)))

        view_add("ROC[m]", None, "ROC", lambda x: alm.str_m(x, space=False))
        view_add("ROC_B[m]", None, "ROC_AR", lambda x: alm.str_m(x, space=False))
        view_add("defocus[D]", 0, "defocus", lambda x: alm.str_D(x, space=False))
        if desc:
            # only add if there is some focus
            view_add(
                "AOI[deg]",
                0,
                "AOI",
                lambda x: alm.unit_str(x, d=2, unit="deg", space=False),
            )
        anno = view["annotate"]
        if anno is None:
            return ", ".join(desc)
        else:
            return anno + " " + ", ".join(desc)


class Mirror(MirrorBase):
    def __init__(self):
        super().__init__()
        with self._internal():
            self["R"] = None
            self["T"] = None
            self["L"] = None

            @self.deco_many_many(
                assignments=["R_eff", "T_eff", "L_eff"],
                dependencies=["R", "T", "L"],
            )
            def param_many(R, T, L):
                if R is None:
                    if T is None:
                        return None, None, None
                        raise RuntimeError("Must specify at least R or T for {}".format(self))
                    if L is None:
                        L = 0
                    R = 1 - T - L
                else:
                    if T is None:
                        if L is None:
                            L = 0
                        T = 1 - R - L
                    else:
                        if L is None:
                            L = 1 - T - L
                        else:
                            # check consistency
                            # TODO, make this handle floats, symbolics
                            assert T == 1 - R - T
                return R, T, L

            @self.deco_one_one("L_eff")
            def has_loss(L):
                if np.any(L != 0):
                    return True
                else:
                    return False

        self.values_settable(
            [
                "R",
                "T",
                "L",
            ]
        )

    def matrix_RTL(self, manip):
        # expects basis to be the same in all links
        basis = manip.link_basis("A!i")

        # mirrors don't care about transverse modes!
        # (any system algorithm will apply I/O basis changes)
        # (They do care about them for torque reasons, once mechanics is added)
        has_loss = manip.p["has_loss"]

        kmatrix_R = matrix.KeyMatrixSame(
            kmatrix={(): {(): [[manip.p["R_eff"] ** 0.5]]}},
            build=True,
            check=manip.check_build,
        )
        kmatrix_T = matrix.KeyMatrixSame(
            kmatrix={(): {(): [[manip.p["T_eff"] ** 0.5]]}},
            build=True,
            check=manip.check_build,
        )
        if has_loss:
            kmatrix_L = matrix.KeyMatrixSame(
                kmatrix={(): {(): [[manip.p["L_eff"] ** 0.5]]}},
                build=True,
                check=manip.check_build,
            )
            return kmatrix_R, kmatrix_T, kmatrix_L
        else:
            return kmatrix_R, kmatrix_T, None


class DichroicMirror(MirrorBase):
    def __init__(self):
        super().__init__()
        with self._internal():
            self["R_default"] = None
            self["T_default"] = None
            self["L_default"] = None

            self["wavelength"] = None

            self["R_dichroic"] = None
            self["T_dichroic"] = None
            self["L_dichroic"] = None

            def param_many(R, T, L):
                if R is None:
                    if T is None:
                        raise RuntimeError("Must specify at least R or T")
                    if L is None:
                        L = 0
                    R = 1 - T - L
                else:
                    if T is None:
                        if L is None:
                            L = 0
                        T = 1 - R - L
                    else:
                        if L is None:
                            L = 1 - T - L
                        else:
                            # check consistency
                            # TODO, make this handle floats, symbolics
                            assert T == 1 - R - T
                return R, T, L

            self.deco_many_many(
                assignments=["R_default_eff", "T_default_eff", "L_default_eff"],
                dependencies=["R_default", "T_default", "L_default"],
                func=param_many,
            )

            self.deco_many_many(
                assignments=["R_dichroic_eff", "T_dichroic_eff", "L_dichroic_eff"],
                dependencies=["R_dichroic", "T_dichroic", "L_dichroic"],
                func=param_many,
            )

            @self.deco_one_one("L_eff")
            def has_loss(L):
                if np.any(L != 0):
                    return True
                else:
                    return False

        self.values_settable(
            [
                "R_default",
                "R_default",
                "R_default",
                "R_dichroic",
                "R_dichroic",
                "R_dichroic",
            ]
        )

    def matrix_RTL(self, manip):
        # expects basis to be the same in all links
        basis = manip.link_basis("A!i")
        has_loss = manip.p["has_loss"]

        Wk_dichroic = manip.parameter_to_wk(manip.p["wavelength"])
        km_R = (collections.defaultdict(dict),)
        km_T = (collections.defaultdict(dict),)
        if has_loss:
            km_L = (collections.defaultdict(dict),)
        else:
            km_L = None

        for Wk in basis["wavenumber"].enumerated:
            if Wk == Wk_dichroic:
                km_R[(Wk,)][(Wk,)] = [[manip.p["R_dichroic_eff"] ** 0.5]]
                km_T[(Wk,)][(Wk,)] = [[manip.p["T_dichroic_eff"] ** 0.5]]
                km_L[(Wk,)][(Wk,)] = [[manip.p["L_dichroic_eff"] ** 0.5]]
            else:
                km_R[(Wk,)][(Wk,)] = [[manip.p["R_default_eff"] ** 0.5]]
                km_T[(Wk,)][(Wk,)] = [[manip.p["T_default_eff"] ** 0.5]]
                km_L[(Wk,)][(Wk,)] = [[manip.p["L_default_eff"] ** 0.5]]

        # mirrors don't care about transverse modes!
        # (any system algorithm will apply I/O basis changes)
        # (They do care about them for torque reasons, once mechanics is added)

        kmatrix_R = matrix.KeyMatrixSame(
            st=(basis["wavenumber"],),
            kmatrix=km_R,
            build=True,
            check=manip.check_build,
        )
        kmatrix_T = matrix.KeyMatrixSame(
            st=(basis["wavenumber"],),
            kmatrix=km_T,
            build=True,
            check=manip.check_build,
        )

        if km_L is not None:
            kmatrix_L = matrix.KeyMatrixSame(
                st=(basis["wavenumber"],),
                kmatrix=km_L,
                build=True,
                check=manip.check_build,
            )
        else:
            kmatrix_L = None

        return kmatrix_R, kmatrix_T, kmatrix_L
