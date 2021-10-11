# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals, absolute_import
import collections

from transient import matrix
from .. import base


class Vacuum(base.OpticalObject):
    def __init__(self):
        super(Vacuum, self).__init__()

    @classmethod
    def visit_port_information(cls, manip):
        manip.gen_optical_port('+A', 'A')
        return

    @classmethod
    def visit_matrix_algorithm_DCAC(cls, manip):
        if manip.is_AC:
            #TODO
            kmatrix = optical_quantum_noise_matrix(manip)
            manip.add_noise('A!o', kmatrix)


def optical_quantum_noise_matrix(manip):
    settings = manip.optical_settings()
    basis = manip.link_basis('A!o')
    assert(settings.quantum == 'pm')

    #this is an array of the frequencies
    km = collections.defaultdict(dict)
    for Fk in basis['frequency'].enumerated:
        for Wk, iwval in manip.basis_wavenumbers(with_keys = True):
            km[(Fk, Wk, '+AC')][(Fk, Wk, '-AC')] = [[manip.constants.c_m_s * manip.constants.h_Js * iwval / 2]]
    kmatrix = matrix.KeyMatrixSame(
        st = (
            basis['frequency'],
            basis['wavenumber'],
            basis['quantum'],
        ),
        kmatrix = km,
        build   = True,
        check   = manip.check_build,
    )
    return kmatrix


