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
from . import space


class Marker(base.OpticalObject):
    annotate_as_space = False

    def port_chain(self, p, pname):
        bmap = {
            '+A-t' :  (None, '+B'),
            '+B-t' :  (None, '+A'),
        }.get(pname, None)

        if bmap is not None:
            return bmap

        return super().port_chain(p, pname)

    @classmethod
    def visit_port_information(cls, manip):
        manip.gen_optical_port('+A', 'A')
        manip.gen_optical_port('+B', 'B')
        return

    def visit_mode_matching_linkage(self, manip):
        manip.add_link('B!i', 'A!o', None)
        manip.add_link('A!i', 'B!o', None)

    #def visit_mode_matching_transport(self, manip):
    #    length_m = manip.p['length[m]']
    #    #no need to check these since the space can only be called on proper
    #    #links and the two directions are identical
    #    #manip.lport_fr
    #    #manip.lport_to

    #    #the P-builders are for fast optimization solving
    #    def p_builderXY(p):
    #        length_m = p['length[m]']
    #        return matrix_stack([[1, length_m], [0, 1]])

    #    manip.set_XYpropagator(p_builderXY)

    #    manip.set_Zpropagator({'length[m]' : 1})

    #    matrix = p_builderXY(manip.p)

    #    def inc_builder(z):
    #        return matrix_stack([[1, z], [0, 1]])

    #    manip.set_XYincremental([
    #        (length_m, inc_builder, matrix)
    #    ])
    #    return


class MaterialMarker(space.Space):
    annotate_as_space = False
    #TODO, need to apply substrate to propagation
