# -*- coding: utf-8 -*-
"""
"""

import numpy as np

from wavestate.utilities.np import matrix_stack
from .. import base
from . import alm


class SubstrateLens(base.OpticalObject):
    def __init__(self):
        super().__init__()
        with self._internal():
            self['ROC[m]']               = None
            self['ROC_B[m]']             = None
            self['substrate']            = 'fused_silica'
            self['depth[m]']             = 0
            self['depth[m].lower_bound'] = 0
            self['AOI[deg]']             = 0
            self['translation[m]']       = 0
            self['defocus[D]']           = 0
            self['annotate']             = 'Lens'
        #TODO, disallow depth=None

    def port_chain(self, p, pname):
        bmap = {
            '+A-t' :  (None, '+B'),
            '+B-t' :  (None, '+A'),
        }.get(pname, None)

        if bmap is not None:
            return bmap

        return super().port_chain(p, pname)

    def port_forward(self, p, pname):
        bmap = {
            '+A'  : None,
            '+B'  : None,
        }.get(pname, None)
        if bmap is not None:
            return bmap
        return super().port_forward(p, pname)

    @classmethod
    def visit_port_information(cls, manip):
        manip.gen_optical_port('+A', 'A')
        manip.gen_optical_port('+B', 'B')
        return

    def visit_matrix_algorithm_ACDC(self, manip):
        manip.add_link('B!i', 'A!o', 1)
        manip.add_link('A!i', 'B!o', 1)

    def visit_mode_matching_linkage(self, manip):
        manip.add_link('B!i', 'A!o', manip.p['depth[m]'])
        manip.add_link('A!i', 'B!o', manip.p['depth[m]'])

    def visit_mode_matching_transport(self, manip):
        n_A = 1
        n_B = 1
        n_I = manip.IOR_n(manip.p['substrate'])

        ROC_A_m = manip.p['ROC[m]']
        depth_m = manip.p['depth[m]']
        ROC_B_m = manip.p['ROC_B[m]']
        AOI_rad = manip.p['AOI[deg]'] * np.pi / 180

        def M_side(p):
            defocus_D = p['defocus[D]']
            M_D = matrix_stack([[1, 0], [-defocus_D, 1]])

            translation_m = p['translation[m]']
            M_t = matrix_stack([[1, translation_m], [0, 1]])
            M_t_i = matrix_stack([[1, -translation_m], [0, 1]])
            return M_D, M_t, M_t_i
        M_D, M_t, M_t_i = M_side(manip.p)

        if manip.lport_to == 'B!o':
            matrixAI_X = alm.interface_ROC_AOI_X(ROC_A_m, n_from = n_A, n_to = n_I, AOI_rad = AOI_rad) @ M_t
            matrixAI_Y = alm.interface_ROC_AOI_Y(ROC_A_m, n_from = n_A, n_to = n_I, AOI_rad = AOI_rad) @ M_t
            matrixII = matrix_stack([[1, depth_m], [0, 1]])
            matrixIB_X = M_t_i @ M_D @ alm.interface_ROC_AOI_X(ROC_B_m, n_from = n_I, n_to = n_B, AOI_rad = AOI_rad, neg = True)
            matrixIB_Y = M_t_i @ M_D @ alm.interface_ROC_AOI_Y(ROC_B_m, n_from = n_I, n_to = n_B, AOI_rad = AOI_rad, neg = True)

            def inc_builder(z):
                return matrix_stack([[1, z], [0, 1]])

            manip.set_Xincremental([
                (0,       None,        matrixAI_X),
                (depth_m, inc_builder, matrixII),
                (0,       None,        matrixIB_X),
            ])

            manip.set_Yincremental([
                (0,       None,        matrixAI_Y),
                (depth_m, inc_builder, matrixII),
                (0,       None,        matrixIB_Y),
            ])

            #the P-builders are for fast optimization solving
            def p_builder_X(p):
                ROC_A_m = p['ROC[m]']
                depth_m = p['depth[m]']
                ROC_B_m = p['ROC_B[m]']
                AOI_rad = p['AOI[deg]'] * np.pi / 180
                M_D, M_t, M_t_i = M_side(p)
                depth_m = depth_m / np.cos(AOI_rad * n_A / n_I)

                matrixAI = alm.interface_ROC_AOI_X(ROC_A_m, n_from = n_A, n_to = n_I, AOI_rad = AOI_rad)
                matrixII = matrix_stack([[1, depth_m], [0, 1]])
                matrixIB = alm.interface_ROC_AOI_X(ROC_B_m, n_from = n_I, n_to = n_B, AOI_rad = AOI_rad, neg = True)
                return M_t_i @ M_D @ matrixIB @ matrixII @ matrixAI @ M_t

            manip.set_Xpropagator(p_builder_X)

            def p_builder_Y(p):
                ROC_A_m = p['ROC[m]']
                depth_m = p['depth[m]']
                ROC_B_m = p['ROC_B[m]']
                AOI_rad = p['AOI[deg]'] * np.pi / 180
                M_D, M_t, M_t_i = M_side(p)
                depth_m = depth_m / np.cos(AOI_rad * n_A / n_I)

                matrixAI = alm.interface_ROC_AOI_Y(ROC_A_m, n_from = n_A, n_to = n_I, AOI_rad = AOI_rad)
                matrixII = matrix_stack([[1, depth_m], [0, 1]])
                matrixIB = alm.interface_ROC_AOI_Y(ROC_B_m, n_from = n_I, n_to = n_B, AOI_rad = AOI_rad, neg = True)
                return M_t_i @ M_D @ matrixIB @ matrixII @ matrixAI @ M_t

            manip.set_Ypropagator(p_builder_Y)
            manip.set_Zpropagator({'depth[m]' : (1 / np.cos(AOI_rad * n_A / n_I), 'AOI[deg]')})
        else:
            matrixBI_X = alm.interface_ROC_AOI_X(ROC_B_m, n_from = n_B, n_to = n_I, AOI_rad = AOI_rad) @ M_t
            matrixBI_Y = alm.interface_ROC_AOI_Y(ROC_B_m, n_from = n_B, n_to = n_I, AOI_rad = AOI_rad) @ M_t
            matrixII = matrix_stack([[1, depth_m], [0, 1]])
            matrixIA_X = M_t_i @ M_D @ alm.interface_ROC_AOI_X(ROC_A_m, n_from = n_I, n_to = n_A, AOI_rad = AOI_rad, neg = True)
            matrixIA_Y = M_t_i @ M_D @ alm.interface_ROC_AOI_Y(ROC_A_m, n_from = n_I, n_to = n_A, AOI_rad = AOI_rad, neg = True)

            def inc_builder(z):
                return matrix_stack([[1, z], [0, 1]])

            manip.set_Xincremental([
                (0,       None,        matrixBI_X),
                (depth_m, inc_builder, matrixII),
                (0,       None,        matrixIA_X),
            ])

            manip.set_Yincremental([
                (0,       None,        matrixBI_Y),
                (depth_m, inc_builder, matrixII),
                (0,       None,        matrixIA_Y),
            ])

            #the P-builders are for fast optimization solving
            def p_builder_X(p):
                ROC_A_m = p['ROC[m]']
                depth_m = p['depth[m]']
                ROC_B_m = p['ROC_B[m]']
                AOI_rad = p['AOI[deg]'] * np.pi / 180
                M_D, M_t, M_t_i = M_side(p)
                depth_m = depth_m / np.cos(AOI_rad * n_A / n_I)

                matrixBI = alm.interface_ROC_AOI_X(ROC_B_m, n_from = n_B, n_to = n_I, AOI_rad = AOI_rad)
                matrixII = matrix_stack([[1, depth_m], [0, 1]])
                matrixIA = alm.interface_ROC_AOI_X(ROC_A_m, n_from = n_I, n_to = n_A, AOI_rad = AOI_rad, neg = True)
                return M_t_i @ M_D @ matrixIA @ matrixII @ matrixBI @ M_t

            def p_builder_Y(p):
                ROC_A_m = p['ROC[m]']
                depth_m = p['depth[m]']
                ROC_B_m = p['ROC_B[m]']
                AOI_rad = p['AOI[deg]'] * np.pi / 180
                depth_m = depth_m / np.cos(AOI_rad * n_A / n_I)

                matrixBI = alm.interface_ROC_AOI_Y(ROC_B_m, n_from = n_B, n_to = n_I, AOI_rad = AOI_rad)
                matrixII = matrix_stack([[1, depth_m], [0, 1]])
                matrixIA = alm.interface_ROC_AOI_Y(ROC_A_m, n_from = n_I, n_to = n_A, AOI_rad = AOI_rad, neg = True)
                return M_t_i @ M_D @ matrixIA @ matrixII @ matrixBI @ M_t

            manip.set_Xpropagator(p_builder_X)
            manip.set_Ypropagator(p_builder_Y)

            manip.set_Zpropagator({'depth[m]' : (1 / np.cos(AOI_rad * n_A / n_I), 'AOI[deg]')})

        manip.set_annotations([
            dict(
                z = 0,
                z_end = depth_m,
                type = 'lens',
            )
        ])
        return

    def visit_mm_anno_description(self, pbg, view, descB):
        desc = []

        def view_add(name, default, name2 = None, transform = lambda x: '{:.3f}'.format(x)):
            val = view[name]
            if val == default:
                return
            if name2 is None:
                name2 = name
            desc.append(
                "{}={}".format(name2, transform(val))
            )
        view_add('ROC[m]', None, 'ROC', lambda x : alm.str_m(x, space=False))
        view_add('ROC_B[m]', None, 'ROCB', lambda x : alm.str_m(x, space=False))
        view_add('defocus[D]', 0, 'defocus', lambda x : alm.str_D(x, space=False))
        #TODO: include substrate
        if desc:
            #only add if there is some focus
            view_add('AOI[deg]', 0, 'AOI', lambda x: alm.unit_str(x, d=2, unit = 'deg', space = False))
        view_add('translation[m]', 0, 'shift', lambda x : alm.str_m(x, space=False))
        anno = view['annotate']
        if anno is None:
            return ', '.join(desc)
        else:
            return anno + ' ' + ', '.join(desc)

