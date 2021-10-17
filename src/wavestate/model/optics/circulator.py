# -*- coding: utf-8 -*-
"""
"""

from .. import base


class Circulator4(base.OpticalObject):
    """
    A perfect 4-port optical circulator. Cycles the optical ports A->B->C->D->A
    """
    def port_chain(self, p, pname):
        bmap = {
            '+A-t' :  (None, '+B'),
            '+B-t' :  (None, '+C'),
            '+C-t' :  (None, '+D'),
            '+D-t' :  (None, '+A'),
        }.get(pname, None)

        if bmap is not None:
            return bmap

        return super(Circulator4, self).port_chain(p, pname)

    @classmethod
    def visit_port_information(cls, manip):
        manip.gen_optical_port('+A', 'A')
        manip.gen_optical_port('+B', 'B')
        manip.gen_optical_port('+C', 'C')
        manip.gen_optical_port('+D', 'D')
        return

    def visit_matrix_algorithm_ACDC(self, manip):
        manip.add_link('A!i', 'B!o', 1)
        manip.add_link('B!i', 'C!o', 1)
        manip.add_link('C!i', 'D!o', 1)
        manip.add_link('D!i', 'A!o', 1)
        return

    def visit_mode_matching_linkage(self, manip):
        manip.add_link('A!i', 'B!o', None)
        manip.add_link('B!i', 'C!o', None)
        manip.add_link('C!i', 'D!o', None)
        manip.add_link('D!i', 'A!o', None)
        return
