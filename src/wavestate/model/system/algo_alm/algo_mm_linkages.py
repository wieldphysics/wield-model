# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals, absolute_import
#import numpy as np
import collections
from collections import abc
import declarative
from transient.utilities.priority_queue import HeapPriorityQueue
from ...optics import alm

from .. import algo_phys

from . import mm_transporter


class ModeMatchingLinkageAlgorithm(object):
    """
    """

    def __init__(self, pa):
        self.log = pa.log
        self.pbg = pa.pbg
        #pg.print_parameters_eval()
        self.fs = pa.fs
        self.bg = pa.bg
        self.pa = pa

        self._object_edges = collections.defaultdict(dict)
        self._build_MM_link_graph()

        self.SRE = self.SRE_linkage()
        return

    def _build_MM_link_graph(self):
        for obj in self.pbg.object_iter():
            try:
                visit_algo = obj.visit_mode_matching_linkage
            except AttributeError:
                continue
            else:
                #TODO verbose option for found objects?
                #print(obj)
                pass

            manip = MMAlgorithmLinkManipulator(
                obj     = obj,
                mm_algo = self,
            )

            visit_algo(manip)
        return

    def SRE_linkage(self):
        seq = collections.defaultdict(set)
        req = collections.defaultdict(set)
        edges = dict()

        for oLp_fr, eset in self.bg.link_seq.items():
            m_fr = oLp_fr
            #print("FT: ", m_fr, eset)
            for oLp_to in eset:
                m_to = oLp_to
                edges[m_fr, m_to] = None
                seq[m_fr].add(m_to)
                req[m_to].add(m_fr)

        for obj, edict in self._object_edges.items():
            for (l_fr, l_to), e_kmat in edict.items():
                m_to = (obj, l_to)
                m_fr = (obj, l_fr)
                edges[m_fr, m_to] = e_kmat
                seq[m_fr].add(m_to)
                req[m_to].add(m_fr)
        return dict(seq), dict(req), edges

    def _dijkstra2(
        self,
        oLp1_group,
        oLp2_group,
    ):
        """
        """
        seq, req, edges = self.SRE

        best_weights = dict()
        unchecked = HeapPriorityQueue()
        parents = dict()
        visited = set()

        if isinstance(oLp1_group, abc.Mapping):
            for oLp1, weight in oLp1_group.items():
                name1 = oLp1
                best_weights[name1] = weight
                parents[name1] = None
                unchecked.push((weight, name1))
            oLp1_group = set(oLp1_group.keys())
        else:
            for oLp1 in oLp1_group:
                name1 = oLp1
                best_weights[name1] = (0, 0)
                parents[name1] = None
                unchecked.push(((0, 0), name1))

        oLp2s_remaining = set(oLp2_group)

        while unchecked:
            weight, node = unchecked.pop()
            weight_len, weight_edgecnt = weight
            if node in visited:
                continue

            for n_to in seq.get(node, ()):
                edge_len = edges[node, n_to]
                if edge_len is None:
                    edge_len = 0
                new_weight = (weight_len + edge_len, weight_edgecnt + 1)
                best_weight = best_weights.setdefault(n_to, new_weight)
                if new_weight < best_weight:
                    best_weights[n_to] = new_weight
                    parents[n_to] = node
                    unchecked.push((new_weight, n_to))
                if new_weight == best_weight:
                    parents[n_to] = node
                    unchecked.push((new_weight, n_to))
            visited.add(node)
            if node in oLp2s_remaining:
                oLp2s_remaining.remove(node)
                if not oLp2s_remaining:
                    break
        if len(oLp2s_remaining) >= len(oLp2_group):
            raise RuntimeError(
                "No path between groups {} and {}".format(
                    [self.pbg.op_str(oLp) for oLp in oLp1_group],
                    [self.pbg.op_str(oLp) for oLp in oLp2_group],
                )
            )

        nonunique_shortest = False
        oLp2_best = None
        weight_best = (float('inf'), 0)
        for oLp2 in oLp2_group:
            try:
                weight = best_weights[oLp2]
            except KeyError:
                continue
            if weight < weight_best:
                oLp2_best = oLp2
                weight_best = weight
            elif weight == weight_best:
                nonunique_shortest = True

        #TODO, check for path-overlap uniqueness
        #this should determine if cavities are linear or travelling-wave

        shortest_path = []
        node = oLp2_best
        while node is not None:
            shortest_path.append(node)
            node = parents.get(node)

        #gives shortest path between the groups
        #remap through the name-mapping to use oLp
        return wavestate.bunch.Bunch(
            path_shortest      = shortest_path[::-1],
            weight_shortest    = weight_best,
            parents            = {k: v if v is not None else None for k, v in parents.items()},
            best_weights       = {k: v for k, v in best_weights.items()},
            oLp1_group         = oLp1_group,
            oLp2_group         = oLp2_group,
            nonunique_shortest = nonunique_shortest,
        )

    def _path_transporters(self, oLp_path, Wk):
        assert(Wk is not None)
        Wk = self.fs.parameter_to_wk(Wk)

        if not oLp_path:
            Xtransporter = mm_transporter.MMTransporter(
                oLp_path  = oLp_path,
                Wk       = Wk,
                prop     = [],
                inc      = [],
                prop_ol2idx = dict(),
                inc_ol2idx  = dict(),
            )
            return wavestate.bunch.Bunch(
                oLp_path     = oLp_path,
                Wk           = Wk,
                X            = Xtransporter,
                Y            = Xtransporter,
            )

        #these will hold (obj, param-function) pairs
        #the object storage allows param views to be made for the function
        Xprop = []
        Yprop = []
        Zprop_scales = []
        Zprop_deps = []

        #these will hold (length, transport-function, matrix) tuples
        #the lengths are how far the transport-function acts and the matrix
        #is the final value of the transport-function
        Xinc   = []
        Yinc   = []

        Xprop_ol2idx = dict()
        Yprop_ol2idx = dict()
        Zprop_ol2idx = dict()
        Xinc_ol2idx  = dict()
        Yinc_ol2idx  = dict()

        #holds the lengths used for path formation, or None for trivial linkage
        edges = self.SRE[2]

        Xprop_ol2idx[oLp_path[0]] = 0
        Yprop_ol2idx[oLp_path[0]] = 0
        Zprop_ol2idx[oLp_path[0]] = 0
        Xinc_ol2idx[oLp_path[0]]  = 0
        Yinc_ol2idx[oLp_path[0]]  = 0

        for oLp_fr, oLp_to in path_pairs(oLp_path):

            edge = edges[oLp_fr, oLp_to]
            if edge is None:
                continue
            obj_fr = oLp_fr[0]
            obj_to = oLp_to[0]
            #only internal linkages can be nontrivial
            assert(obj_fr is obj_to)
            manip = MMAlgorithmTransportManipulator(
                obj      = obj_fr,
                mm_algo  = self,
                lport_fr = oLp_fr[1],
                lport_to = oLp_to[1],
                Wk       = Wk,
            )
            obj_fr.visit_mode_matching_transport(manip)

            Xprop_ol2idx[oLp_fr] = len(Xprop)
            Yprop_ol2idx[oLp_fr] = len(Yprop)
            Zprop_ol2idx[oLp_fr] = len(Zprop_scales)
            Xinc_ol2idx[oLp_fr]  = len(Xinc)
            Yinc_ol2idx[oLp_fr]  = len(Yinc)

            #now commit the manipulator settings into the lists
            Xprop.extend((obj_fr, p_trans) for p_trans in manip._Xprop)
            Yprop.extend((obj_fr, p_trans) for p_trans in manip._Yprop)
            Zp_scale = dict()
            Zp_deps = dict()
            for p_trans, clist in manip._Zprop.items():
                #TODO need a pbg.pref2vtup method rather than this split nonsense
                vtup = tuple(p_trans.split('.'))
                if isinstance(clist, (list, tuple)):
                    Zp_scale[(obj_fr, vtup)] = clist[0]
                    Zp_deps[(obj_fr, vtup)] = [(obj_fr, tuple(c.split('.'))) for c in clist[1:]]
                else:
                    Zp_scale[(obj_fr, vtup)] = clist
            Zprop_scales.append(Zp_scale)
            Zprop_deps.append(Zp_deps)
            Xinc.extend(manip._Xinc)
            Yinc.extend(manip._Yinc)

            Xprop_ol2idx[oLp_to] = len(Xprop)
            Yprop_ol2idx[oLp_to] = len(Yprop)
            Zprop_ol2idx[oLp_to] = len(Zprop_scales)
            Xinc_ol2idx[oLp_to] = len(Xinc)
            Yinc_ol2idx[oLp_to] = len(Yinc)
            #TODO, check that the total inc length equals the edge length

        Xtransporter = mm_transporter.MMTransporter(
            oLp_path  = oLp_path,
            Wk       = Wk,
            prop     = Xprop,
            inc      = Xinc,
            prop_ol2idx = Xprop_ol2idx,
            inc_ol2idx  = Xinc_ol2idx,
        )

        if Xinc is Yinc:
            Ytransporter = Xtransporter
        else:
            Ytransporter = mm_transporter.MMTransporter(
                oLp_path  = oLp_path,
                Wk       = Wk,
                prop     = Yprop,
                inc      = Yinc,
                prop_ol2idx = Yprop_ol2idx,
                inc_ol2idx  = Yinc_ol2idx,
            )

        return wavestate.bunch.Bunch(
            oLp_path     = oLp_path,
            Wk           = Wk,
            X            = Xtransporter,
            Y            = Ytransporter,
            Zprop_scales = Zprop_scales,
            Zprop_deps   = Zprop_deps,
            Zprop_ol2idx = Zprop_ol2idx,
        )

    def _safe_oLp_path(self, oLp_set_seq, loop = False, allow_non_unique = False):
        dijkstra_seq = [self._dijkstra2(oLp_set_seq[0], oLp_set_seq[1])]
        last_oLp_set = oLp_set_seq[1]

        for oLp_set in oLp_set_seq[2:]:
            last_oLp_set_weighted = dict()
            for oLp in last_oLp_set:
                weight = dijkstra_seq[-1].best_weights.get(oLp, None)
                if weight is not None:
                    last_oLp_set_weighted[oLp] = weight
            dijkstra_seq.append(self._dijkstra2(last_oLp_set_weighted, oLp_set))
            last_oLp_set = oLp_set

        if loop:
            starters = {oLp : oLp for oLp in dijkstra_seq[-1].oLp2_group}
            for dijkstra in reversed(dijkstra_seq):
                starters_next = dict()
                for oLp_fr, oLp_to in starters.items():
                    node = oLp_to
                    while node is not None:
                        last_node = node
                        node = dijkstra.parents.get(node)
                    if last_node is not oLp_to:
                        starters_next[oLp_fr] = last_node
                starters = starters_next

            starters_inv = collections.defaultdict(set)
            for oLp_fr, oLp_to in starters.items():
                starters_inv[oLp_to].add(oLp_fr)

            ambiguous_loop = False
            weight_shortest = (float('inf'), 0)
            for oLp_to, oLp_fr_set in starters_inv.items():
                last_oLp_set_weighted = dict()
                for oLp in oLp_fr_set:
                    weight = dijkstra_seq[-1].best_weights.get(oLp, None)
                    if weight is not None:
                        last_oLp_set_weighted[oLp] = weight
                #use the weighted set to dijkstra back the loop
                dijkstra = self._dijkstra2(last_oLp_set_weighted, [oLp_to])
                if dijkstra.nonunique_shortest:
                    ambiguous_loop = True
                if dijkstra.weight_shortest < weight_shortest:
                    link_seq = list(dijkstra.path_shortest[::-1])
                    weight_shortest = dijkstra.weight_shortest
                elif dijkstra.weight_shortest == weight_shortest:
                    ambiguous_loop = True
            if not allow_non_unique and ambiguous_loop:
                raise RuntimeError(
                    "Formulating the unique shortest loop is ambiguous."
                    " Consider adding one more waypoint."
                )

            for dijkstra in reversed(dijkstra_seq):
                node = link_seq[-1]
                while True:
                    node = dijkstra.parents.get(node)
                    if node is not None:
                        link_seq.append(node)
                    else:
                        break
            #reverse and cut off new last element, so that the path is a loop
            link_seq = link_seq[:0:-1]
        else:
            last_oLp_set_weighted = dict()
            for oLp in dijkstra_seq[-1].oLp2_group:
                weight = dijkstra_seq[-1].best_weights.get(oLp, None)
                if weight is not None:
                    last_oLp_set_weighted[oLp] = weight
            sorted_links = sorted(last_oLp_set_weighted, key = lambda k : last_oLp_set_weighted[k])
            best = sorted_links[0]
            if not allow_non_unique and len(sorted_links) > 1 and last_oLp_set_weighted[sorted_links[1]] == last_oLp_set_weighted[best]:
                #then the best path is ambiguous
                raise RuntimeError(
                    "Formulating the unique shortest path is ambiguous."
                    " Consider adding more waypoints."
                )

            link_seq = [best]
            for dijkstra in reversed(dijkstra_seq):
                node = link_seq[-1]
                while True:
                    node = dijkstra.parents.get(node)
                    if node is not None:
                        link_seq.append(node)
                    else:
                        break
            #reverse sequence
            link_seq = link_seq[::-1]
        if not link_seq:
            print("LSEQ ", link_seq, oLp_set_seq)
            assert(False)
        return link_seq


class MMAlgorithmView(algo_phys.PhysicsAlgorithmView):
    _mm_algo = None
    def __init__(self, mm_algo, **kw):
        super(MMAlgorithmView, self).__init__(
            bg_algo = mm_algo.bg,
            pbg     = mm_algo.pbg,
            pa_algo = mm_algo.pa,
            **kw
        )
        self._mm_algo = mm_algo


class MMAlgorithmLinkManipulator(MMAlgorithmView):
    def add_link(self, lport_fr, lport_to, length_m):
        self._mm_algo._object_edges[self._obj][lport_fr, lport_to] = length_m


class MMAlgorithmTransportManipulator(MMAlgorithmView):
    _Zprop       = None
    _Xprop       = None
    _Yprop       = None
    _Xinc        = None
    _Yinc        = None
    _annotations = None

    def __init__(self, mm_algo, lport_fr, lport_to, Wk, **kw):
        self.lport_fr = lport_fr
        self.lport_to = lport_to
        self.Wk = Wk
        super(MMAlgorithmTransportManipulator, self).__init__(
            mm_algo,
            **kw
        )

    def add_link(self, lport_fr, lport_to, length_m):
        self._mm_algo._object_edges[self._obj][lport_fr, lport_to] = length_m

    def IOR_n(self, substrate):
        return alm.substrates[substrate][self.Wk]

    def set_Zpropagator(self, prop = None):
        if prop is None:
            self._Zprop = {}
        else:
            self._Zprop = prop

    def set_XYpropagator(self, prop):
        self.set_Xpropagator(prop)
        self.set_Ypropagator(prop)

    def set_Xpropagator(self, prop):
        self._Xprop = [prop]
        return

    def set_Ypropagator(self, prop):
        self._Yprop = [prop]
        return

    def set_XYincremental(self, inc):
        self.set_Xincremental(inc)
        self.set_Yincremental(inc)
        return

    def set_Xincremental(self, inc):
        #print('inc', self._obj, self.lport_fr, self.lport_to)
        #for (l_m, func, mat) in inc:
        #    print('inc2', mat)
        self._Xinc = inc
        return

    def set_Yincremental(self, inc):
        self._Yinc = inc
        return

    def set_annotations(self, anno):
        self._annotations = anno
        return


def path_pairs(link_path):
    last_link = link_path[0]
    for next_link in link_path[1:]:
        yield last_link, next_link
        last_link = next_link
    return


def path_pairs_loop(link_path):
    last_link = link_path[0]
    for next_link in link_path[1:]:
        yield last_link, next_link
        last_link = next_link
    yield last_link, link_path[0]
    return
