#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

#import numpy as np
#import pytest
from os import path


from wavestate.model import optics
from wavestate import model
from wavestate.model.system import algo_phys
from wavestate.model.system import algo_graphs
from wavestate.model.system import algo_log

from wavestate.pytest import (
    tpath, ic, pprint
)


def system_build():
    obj_sys = model.system1064()
    freqs = obj_sys['frequencies/']
    freqs['Nd1064/order'] = 2

    obj_sys['SRC1/']    = optics.Laser()
    obj_sys['SRC1/wavelength'] = 'Nd1064'

    obj_sys['M1/']    = optics.Mirror()
    obj_sys['M2/']    = optics.Mirror()
    obj_sys['M1_M2/'] = optics.Space()
    obj_sys['M1_M2/length[m]'] = 16
    obj_sys['M1/ROC[m]'] = -18
    obj_sys['M2/ROC[m]'] = -18

    obj_sys['SRC1/power[W]'] = 1
    obj_sys['M1/T']   = .001
    obj_sys['M2/T']   = .001

    obj_sys['SRC_L1/'] = optics.Space()
    obj_sys['L1_L2/']  = optics.Space()
    obj_sys['L2_M1/']  = optics.Space()
    obj_sys['L1/'] = optics.ThinLens()
    obj_sys['L2/'] = optics.ThinLens()

    obj_sys['M1/depth[m]'] = .1
    obj_sys['SRC_L1/length[m]'] = 1
    obj_sys['L1_L2/length[m]'] = 1
    #obj_sys['L2_M1/length[m]'] = 1

    obj_sys['L1_M1_length[m]'] = 2.5
    #add a more complex parameter mapping
    @obj_sys.deco_many_many(
        assignments = (
            'L2_M1/length[m]',
        ),
        dependencies = (
            'L1_M1_length[m]',
            'L1_L2/length[m]',
        ))
    def length_conserve(total, L1L2):
        return (total - L1L2),

    obj_sys['L1/focal_length[m]'] = -1/0.74
    obj_sys['L2/focal_length[m]'] = 2

    obj_sys.bond_add([
        'SRC1+A | SRC_L1+A-t | L1+A-t | L1_L2+A-t',
        'L2+A-t | L2_M1+A-t | M1+B-t | M1_M2+A-t | M2+A'
    ])
    return obj_sys


def test_FP_MM(plot, tpath):
    print()
    log = algo_log.LoggingAlgorithm(
        log_level = 9,
        log_folder = tpath,
        filters = {
            r'digraph' : dict(investigate = True),
        }
    )

    obj_sys = system_build()

    pa = algo_phys.PhysicsAlgorithm(obj_sys)

    #d = pa.mm._dijkstra2out('SRC1/+A', 'M1/+A')
    #print(d.path_shortest)
    #print(d.weight_shortest)
    cavity_loop = pa.mm.cavity_add('cavity', ['M1/', 'M2/'])
    print(cavity_loop)
    cavity_loop = pa.mm.target_add(
        'src',
        waypoints = ['SRC1/+A!o'],
        q = optics.alm.ComplexBeamParam.from_W_R(1.500e-3, 100, wavelength_m = 1064e-9),
    )
    pa.ga.graph_MM_link_matrix(
        path.join(tpath,  'MM_links.pdf'),
        color_cavities = {'cavity' : 'red'},
    )

    overlapper = pa.mm.overlap(
        target_fr = 'src',
        target_to = 'cavity',
        Wk = 1064
    )
    print(overlapper.overlap)
    if plot:
        overlapper.plot(path.join(tpath,  'MM_FP_plot.pdf'))
    return


def test_FP_MM_obj(plot, tpath):
    print()
    obj_sys = system_build()
    obj_sys['SRC1_beam/'] = optics.Target()
    obj_sys['SRC1_beam/waypoints'] = ['SRC1/+A!o']
    obj_sys['SRC1_beam/q'] = optics.ComplexBeamParam.from_W_R(1.500e-3, 100, wavelength_m = 1064e-9)

    obj_sys['cavity/'] = optics.Cavity()
    obj_sys['cavity/waypoints'] = ['M1/', 'M2/']

    pa = algo_phys.PhysicsAlgorithm(obj_sys)

    pa.ga.graph_MM_link_matrix(
        path.join(tpath,  'MM_links.pdf'),
        color_cavities = {'cavity' : 'red'},
    )

    overlapper = pa.mm.overlap(
        target_fr = 'SRC1_beam',
        target_to = 'cavity',
        Wk = 1064
    )
    print(overlapper.overlap)
    if plot:
        overlapper.plot(path.join(tpath,  'MM_FP_plot.pdf'))
    return


def test_FP_MM_optimize(tpath, ic):
    print()
    log = algo_log.LoggingAlgorithm(
        log_level = 9,
        log_folder = tpath,
        filters = {
            r'digraph' : dict(investigate = True),
        }
    )

    obj_sys = system_build()
    obj_sys['L1_L2/length[m]'] = .3
    #obj_sys['SRC_L1/length[m]'] = 1

    pa = algo_phys.PhysicsAlgorithm(obj_sys)

    pa.ga.graph_MM_link_matrix
    #d = pa.mm._dijkstra2out('SRC1/+A', 'M1/+A')
    #print(d.path_shortest)
    #print(d.weight_shortest)
    cavity_loop = pa.mm.cavity_add('cavity', ['M1/', 'M2/'])
    print(cavity_loop)
    cavity_loop = pa.mm.target_add(
        'src',
        waypoints = ['SRC1/+A!o'],
        q = optics.alm.ComplexBeamParam.from_W_R(1.500e-3, 100, wavelength_m = 1064e-9),
    )
    ic(pa.mm._targets)
    pa.ga.graph_MM_link_matrix(
        path.join(tpath,  'MM_links.pdf'),
        color_cavities = {'cavity' : 'red'}
    )

    overlapper = pa.mm.overlap(
        target_fr = 'src',
        target_to = 'cavity',
        Wk = 1064
    )
    print("OVERLAP: ", overlapper.overlap)
    comp = overlapper.compile_overlap_calculation()
    print("OVERLAP2: ", comp.calculate_overlap(t1_recompute = False, t2_recompute = False))

    #for ref in sorted(overlapper.optimizable_parameters().optimizable):
    #    print(ref)

    overlapper.print_optimizable_parameters()
    overlapper.optimize(['L1_L2/length[m]', 'SRC_L1/length[m]'])
    overlapper.plot(path.join(tpath,  'MM_FP_plot.pdf'))
    return
