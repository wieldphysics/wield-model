#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

# from wield import declarative
# import collections
# import numpy as np

try:
    import pygraphviz as pgv
except ImportError:
    pass


class GraphPlottingAlgorithm(object):
    def __init__(self, pa=None, pbg=None, bg=None):
        self.pa = pa
        if bg is None:
            self.bg = self.pa.bg
        if pbg is None:
            self.pbg = self.bg.pbg
        return

    def graph_bond_graph(self, fname):
        G = pgv.AGraph()
        # G.graph_attr['newrank'] = "True"
        edge_str = []
        # TODO, sort to make the plotting stable
        for m_fr, to_set in self.bg.bond_seq.items():
            m_fr = self.pbg.path_str(m_fr[0]) + m_fr[1]
            for m_to in to_set:
                m_to = self.pbg.path_str(m_to[0]) + m_to[1]
                G.add_edge(m_fr, m_to)
                edge_str.append((m_fr, m_to))
        for obj, pset in self.bg.object_ports.items():
            m_nodes = []
            oname = self.pbg.path_str(obj)
            for pname in pset:
                m_nodes.append(oname + pname)
            G.add_subgraph(m_nodes, "cluster_" + oname)
        edge_str.sort()
        # ic(edge_str)
        G.draw(fname, prog="dot")

    def graph_DC_link_matrix(self, fname, cluster_links=False):
        nmap, SREIO = self.pa.dc.SREIO_DC(
            map_nodes=True,
            subtract_1=True,
        )
        (seq, req, edges, inputs, outputs) = SREIO

        # TODO, sort to make the plotting stable
        import pygraphviz as pgv

        edge_str = []
        G = pgv.AGraph(directed=True)
        G.graph_attr["newrank"] = "True"
        for (m_fr, m_to), kmat in edges.items():
            if m_fr == m_to:
                # suppress self-edges
                continue
            G.add_edge(m_fr, m_to)
            edge_str.append((m_fr, m_to, kmat))
        # ic(pa.bg.object_linkages)
        if cluster_links:
            for obj, pset in self.bg.object_linkages.items():
                m_nodes = []
                oname = self.pbg.path_str(obj)
                for pname in pset:
                    m_nodes.append(oname + "+" + pname)
                G.add_subgraph(m_nodes, "cluster_" + oname)
        edge_str.sort()
        # ic(edge_str)
        G.draw(fname, prog="dot")

    def graph_MM_link_matrix(
        self,
        fname,
        color_chains={},
        color_cavities={},
    ):
        SRE = self.pa.mm.SRE_linkage()
        seq, req, edges = SRE

        # TODO, sort to make the plotting stable
        import pygraphviz as pgv

        G = pgv.AGraph(directed=True)
        G.graph_attr["newrank"] = "True"
        for (m_fr, m_to), weight in edges.items():
            if m_fr == m_to:
                # suppress self-edges
                continue
            if weight is not None:
                G.add_edge(m_fr, m_to, label=str(weight))
            else:
                G.add_edge(m_fr, m_to)

        for color, chain in color_chains.items():
            n_fr = chain[0]
            for n_to in chain[1:]:
                G.add_edge(n_fr, n_to, color=color)
                n_fr = n_to

        # TODO, don't use internal API for the cavity targets?
        for target, color in color_cavities.items():
            target_oP = self.pbg.referred_vtup(target)
            chain = list(self.pa.mm.cavity_targets[target_oP])
            # form a loop
            chain = chain + [chain[0]]

            n_fr = chain[0]
            for n_to in chain[1:]:
                G.add_edge(n_fr, n_to, color=color)
                n_fr = n_to

        # ic(pa.bg.object_linkages)
        G.draw(fname, prog="dot")
