# -*- coding: utf-8 -*-
"""
"""

import numpy as np

from wavestate.utilities.np import matrix_stack
from .. import base
from . import alm


class ThinLens(base.OpticalObject):
    def __init__(self):
        super(ThinLens, self).__init__()
        with self._internal():
            self['focal_length[m]'] = None
            self['defocus[D]']      = 0
            self['annotate']        = 'Lens'

    def port_chain(self, p, pname):
        bmap = {
            '+A-t' :  (None, '+B'),
            '+B-t' :  (None, '+A'),
        }.get(pname, None)

        if bmap is not None:
            return bmap

        return super(ThinLens, self).port_chain(p, pname)

    @classmethod
    def visit_port_information(cls, manip):
        manip.gen_optical_port('+A', 'A')
        manip.gen_optical_port('+B', 'B')
        return

    def visit_matrix_algorithm_ACDC(self, manip):
        manip.add_link('B!i', 'A!o', 1)
        manip.add_link('A!i', 'B!o', 1)
        return

    def visit_mode_matching_linkage(self, manip):
        manip.add_link('B!i', 'A!o', 0)
        manip.add_link('A!i', 'B!o', 0)
        return

    def visit_mode_matching_transport(self, manip):
        #the P-builders are for fast optimization solving
        def p_builderXY(p):
            fl_m = p['focal_length[m]']
            D = p['defocus[D]']
            if fl_m is None:
                M = matrix_stack([[1, 0], [-D, 1]])
            else:
                M = matrix_stack([[1, 0], [-1/fl_m - D, 1]])
            return M

        manip.set_XYpropagator(p_builderXY)
        manip.set_Zpropagator()
        #manip.link_Xpropagator(p_builder)
        #manip.link_Ypropagator(p_builder)
        matrix = p_builderXY(manip.p)

        manip.set_XYincremental([
            (0, None, matrix)
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
        view_add('focal_length[m]', None, 'f', lambda x : alm.str_m(x, space=False))
        view_add('defocus[D]', 0, 'defocus', lambda x : alm.str_D(x, space=False))
        anno = view['annotate']
        if anno is None:
            return ', '.join(desc)
        else:
            return anno + ' ' + ', '.join(desc)


class ThinMirror(base.OpticalObject):
    def __init__(self):
        super(ThinMirror, self).__init__()
        with self._internal():
            self['ROC[m]'] = None
            self['AOI[deg]'] = 0
            self['defocus[D]'] = 0
            self['defocusX[D]'] = 0
            self['defocusY[D]'] = 0
            self['annotate'] = 'HRMirror'

    def port_chain(self, p, pname):
        bmap = {
            '+A1-r' :  (None, '+A2'),
            '+A2-r' :  (None, '+A1'),
        }.get(pname, None)

        if bmap is not None:
            return bmap

        return super(ThinMirror, self).port_chain(p, pname)

    @classmethod
    def visit_port_information(cls, manip):
        manip.gen_optical_port('+A1', 'A1')
        manip.gen_optical_port('+A2', 'A2')
        return

    def visit_matrix_algorithm_ACDC(self, manip):
        manip.add_link('A1!i', 'A2!o', 1)
        manip.add_link('A2!i', 'A1!o', 1)
        return

    def visit_mode_matching_linkage(self, manip):
        manip.add_link('A1!i', 'A2!o', 0)
        manip.add_link('A2!i', 'A1!o', 0)
        return

    def visit_mode_matching_transport(self, manip):
        #the P-builders are for fast optimization solving
        def p_builder_X(p):
            roc_m = p['ROC[m]']
            AOI_rad = p['AOI[deg]'] / 180 * np.pi
            M = alm.REFL_ROC_X(roc_m, AOI_rad)
            D = p['defocus[D]'] + p['defocusX[D]']
            MD = matrix_stack([[1, 0], [D, 1]])
            return MD @ M

        def p_builder_Y(p):
            roc_m = p['ROC[m]']
            AOI_rad = p['AOI[deg]'] / 180 * np.pi
            M = alm.REFL_ROC_Y(roc_m, AOI_rad)
            D = p['defocus[D]'] + p['defocusY[D]']
            MD = matrix_stack([[1, 0], [D, 1]])
            return MD @ M

        manip.set_Xpropagator(p_builder_X)
        manip.set_Ypropagator(p_builder_Y)
        manip.set_Zpropagator()

        manip.set_Xincremental([(0, None, p_builder_X(manip.p))])
        manip.set_Yincremental([(0, None, p_builder_Y(manip.p))])
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
        view_add('defocus[D]', 0, 'defocus', lambda x : alm.str_D(x, space=False))
        if desc:
            #only add if there is some focus
            view_add('AOI[deg]', 0, 'AOI', lambda x: alm.unit_str(x, d=2, unit = 'deg', space = False))
        anno = view['annotate']
        if anno is None:
            return ', '.join(desc)
        else:
            return anno + ' ' + ', '.join(desc)


class ThinLensTranslation(base.OpticalObject):
    """
    This Element represents a translation correction to a thin lens element.
    Can be used to show translation adjustments to existing lenses or collimators
    """
    def __init__(self):
        super(ThinLensTranslation, self).__init__()
        with self._internal():
            self['focal_length[m]'] = None
            self['shift[m]']        = None
            self['annotate']        = 'ShiftLens'

    def port_chain(self, p, pname):
        bmap = {
            '+A-t' :  (None, '+B'),
            '+B-t' :  (None, '+A'),
        }.get(pname, None)

        if bmap is not None:
            return bmap

        return super(ThinLensTranslation, self).port_chain(p, pname)

    @classmethod
    def visit_port_information(cls, manip):
        manip.gen_optical_port('+A', 'A')
        manip.gen_optical_port('+B', 'B')
        return

    def visit_matrix_algorithm_ACDC(self, manip):
        manip.add_link('B!i', 'A!o', 1)
        manip.add_link('A!i', 'B!o', 1)
        return

    def visit_mode_matching_linkage(self, manip):
        manip.add_link('B!i', 'A!o', 0)
        manip.add_link('A!i', 'B!o', 0)
        return

    def visit_mode_matching_transport(self, manip):
        #the P-builders are for fast optimization solving
        def p_builderXY(p):
            fl_m = p['focal_length[m]']
            shift_m = p['shift[m]']
            if fl_m is None:
                return matrix_stack([[1, 0], [0, 1]])
            else:
                Mf  = matrix_stack([[1, 0], [-1/fl_m, 1]])
                Mfi = matrix_stack([[1, 0], [1/fl_m, 1]])
                Ms  = matrix_stack([[1, shift_m], [0, 1]])
                Msi = matrix_stack([[1, -shift_m], [0, 1]])

                return Msi @ Mf @ Ms @ Mfi

        manip.set_XYpropagator(p_builderXY)
        manip.set_Zpropagator()
        #manip.link_Xpropagator(p_builder)
        #manip.link_Ypropagator(p_builder)
        matrix = p_builderXY(manip.p)

        manip.set_XYincremental([
            (0, None, matrix)
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
        view_add('focal_length[m]', None, 'f', lambda x : alm.str_m(x, space=False))
        view_add('shift[m]', 0, 'shift', lambda x : alm.str_m(x, space=False))
        anno = view['annotate']
        if anno is None:
            return ', '.join(desc)
        else:
            return anno + ' ' + ', '.join(desc)


