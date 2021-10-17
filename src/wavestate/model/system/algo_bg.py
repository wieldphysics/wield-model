#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


import collections

from ..pgraph.utilities import ref_port_split
from ..base import port_types

from .. import optics
from .. import _base
from .. import base

PORTS_AUTOTERMINATE = {
    port_types.BGT_OPTICS_PHYSICAL : optics.Vacuum,
}

class BondGraphAlgorithm(object):
    def __init__(self, pbg, log = _base.log_default):
        self.pbg = pbg
        self.log = log

        with self.log.heading("Bond Graph"):
            #initially this will hold all BG edges, without yet reducing them to
            #fully port-forwarded form

            #they are in (obj1, pstr1), (obj2, pstr2) form in the edge_list
            #the edge list is "flat" in the sense that all port chains are
            #already resolved
            edge_list = []
            anti_edge_list = []
            for obj in pbg.object_iter():
                #print(pbg.path_str(obj))
                #print(obj.bond_lists(pbg.view(obj)))
                for bl in obj.bond_lists(pbg.view(obj)):
                    is_antilist, edges = self._parse_bond_list(obj, bl)
                    if not is_antilist:
                        edge_list.extend(edges)
                    else:
                        anti_edge_list.extend(edges)

            self.port_forward_mapping = dict()
            bond_edges_map, obj_ports = self._fully_forward_edges(
                edge_list,
                anti_edge_list = anti_edge_list
            )
            self.bond_edges_map = bond_edges_map
            self.object2ports_map = obj_ports

            #info contains info for all ports of all objects, including ports
            #which do not have bonds. It is used to generate the autoterminations
            self.info, self.object_ports_full = self._get_port_information()

            #TODO
            with self.log.heading("Simplify Bond Globs"):
                self._simplify_globs()
            with self.log.heading("Check Bond Pairings"):
                self._check_pairings()
            with self.log.heading("Check Physical Connections"):
                self._check_physical()
            with self.log.heading("Auto Terminate"):
                self._auto_terminate()

            #must be done after autotermination
            self.obj_lport2bport, self.link_basis_types = self._link_map()

            with self.log.heading("Generate Bond Linkages"):
                link_seq, object_linkages = self._gen_linkage_graph()
            self.link_seq = link_seq
            self.object_linkages = object_linkages

        #now with info, some graph simplifications and audits may be performed
        #electronics globbing
        #physical check of optics (only single physical linkages)
        #autotermination of graphs

    def oLp2oBp(self, ol):
        """
        TODO, rename to oLp2oBp
        """
        obj, lport = ol
        bport = self.obj_lport2bport[obj][lport]
        return obj, bport

    def rBp2oLp(self, bpref, obj = None, dir = 'out'):
        """
        rename to rBp2oLp.
        """
        op1 = self.pbg.port_ref(bpref, obj = obj)
        info = self.info[op1]
        lin, lout = info['inout']
        if dir.lower() == 'in':
            ol1 = op1[0], lin
        elif dir.lower() == 'out':
            ol1 = op1[0], lout
        return ol1

    def ref2oLp_set(self, refstr, dir = None, obj = None):
        """
        Returns all links to an object.
        """
        obj = self.pbg.resolve_object(refstr, obj = obj)
        pset = self.object_ports_full[obj]
        ol_set = set()
        for pstr in pset:
            op = obj, pstr
            info = self.info[op]
            lin, lout = info['inout']
            if dir == 'in':
                if lin is not None:
                    ol_set.add((obj, lin))
            elif dir == 'out':
                if lout is not None:
                    ol_set.add((obj, lout))
            elif dir is None:
                if lin is not None:
                    ol_set.add((obj, lin))
                if lout is not None:
                    ol_set.add((obj, lout))
        return ol_set

    def oLp_path2objpath(self, oLp_path, as_refs = False):
        objlist = []
        objlist_reduced = []
        obj_prev = None
        for ol in oLp_path:
            obj, bport = self.oLp2oBp(ol)
            objlist.append(obj)
            if obj_prev is not obj:
                objlist_reduced.append(obj)
            obj_prev = obj
        if not as_refs:
            return objlist_reduced
        else:
            objrefs_reduced = []
            for obj in objlist_reduced:
                objrefs_reduced.append(self.pbg.path_str(obj))
            return objrefs_reduced

    def rAp2oLp_set(self, refstr, dir = None, obj = None):
        """
        This is an "any-port" reference to oLp_set converter. Any port includes
        object-references,
        bond-port references,
        link-port references

        In the case of object and bond-port references, the returned set of link
        ports has more than 1 element.
        """
        if '!' in refstr:
            #link port
            psplit, lsplit = refstr.split('!')
            #TODO, clean up these naming schemes
            #this split allows port fowarding to work
            if lsplit == 'o':
                return self.rAp2oLp_set(psplit, dir = 'out', obj = obj)
            elif lsplit == 'i':
                return self.rAp2oLp_set(psplit, dir = 'in', obj = obj)
            else:
                raise RuntimeError("Unrecognized linkage port direction")
        elif '+' in refstr:
            op1 = self.pbg.port_ref(refstr, obj = obj)
            try:
                info = self.info[op1]
            except KeyError:
                self.port_forward_mapping = self._fully_forward_ports([op1], pmap = self.port_forward_mapping)
                #must be a forwarded port, so then cast through the forwarding
                op1 = self.port_forward_mapping[op1]
                info = self.info[op1]
            lin, lout = info['inout']
            if dir is None:
                ol1 = op1[0], lin
                ol2 = op1[0], lout
                return [ol1, ol2]
            if dir.lower() == 'in':
                ol1 = op1[0], lin
                return [ol1]
            elif dir.lower() == 'out':
                ol1 = op1[0], lout
                return [ol1]
            else:
                raise RuntimeError("Unrecognized direction")
        elif refstr[-1] == '/':
            #object reference, must
            return self.ref2oLp_set(refstr, dir = dir, obj = obj)
        else:
            raise RuntimeError((
                "Unrecognized port '{}', did you forget a trailing '/'?"
            ).format(refstr))

    def _link_map(self):
        """
        Generate the mappings for link-ports back to bond-ports as well as
        generate the type-map for bond type-sets.

        Must be called after autotermination
        """
        obj_lport2bport = collections.defaultdict(dict)
        link_basis_types = dict()
        for op1, info in self.info.items():
            obj, bport = op1
            btype = port_types.BASIS_TYPES[info['type']]
            in1, out1 = info['inout']
            #TODO, testing bug with objects that are unlinked
            #print('linkage: ', obj, in1, out1)
            if in1 is not None:
                obj_lport2bport[obj][in1] = bport
                link_basis_types[obj, in1] = btype
            if out1 is not None:
                obj_lport2bport[obj][out1] = bport
                link_basis_types[obj, out1] = btype
        return obj_lport2bport, link_basis_types

    def _simplify_globs(self):
        return

    def _check_pairings(self):
        return

    def _check_physical(self):
        return

    def _auto_terminate(self):
        terminate_map = dict()
        for op1, info in self.info.items():
            #print("OP", self.pbg.path_str(op1[0]), op1[1])
            edge_set = self.bond_edges_map.get(op1, ())
            if edge_set:
                continue

            obj1, pstr1 = op1
            at = PORTS_AUTOTERMINATE.get(info['type'], None)
            #if at is not None:
            #    print("AUTOTERMINATE: ", self.pbg.path_str(op1[0]), op1[1], at)
            if at is not None:
                terminate_map[op1] = at
                if self.log.hint('lists', default = False):
                    self.log.debug(8, "autoterm:", self.pbg.op_str(op1), at)

        for op1, terminator in terminate_map.items():
            term_name = '<{}{}>'.format(self.pbg.path_str(op1[0]), op1[1])
            term_obj = terminator()
            self.pbg.object_insert((term_name,), term_obj)
            op2 = (term_obj, '+A')
            self.bond_edges_map.setdefault(op1, set()).add(op2)
            self.bond_edges_map.setdefault(op2, set()).add(op1)
            #there are and cannot be any port forwarding in autotermination objects
            self._get_port_information_obj(term_obj, self.info)
        return

    def _gen_linkage_graph(self):
        #for obj, ports in self.object2ports_map.items():

        #using this more like a set, except it can hold some edge information
        link_seq = collections.defaultdict(set)
        object_linkages = collections.defaultdict(set)

        #this will visit all edges twice, once in each direction, this
        #is relied up for symmetry of incoming/outgoing edges
        for op1, edge_set in self.bond_edges_map.items():
            idict1 = self.info[op1]
            obj1, pstr1 = op1

            if idict1['type'] in port_types.INOUT_TYPES:
                #then it is assumed that all op2 are also in INOUT_TYPES since
                #checking was performed
                for op2 in edge_set:
                    obj2, pstr2 = op2
                    idict2 = self.info[op2]
                    in1, out1 = idict1['inout']
                    in2, out2 = idict2['inout']
                    if out1 is not None:
                        link_seq[(obj1, out1)].add((obj2, in2))
                    object_linkages[obj1].add(out1)
                    object_linkages[obj2].add(in2)
        return dict(link_seq), dict(object_linkages)

    def _get_port_information(self):
        full_info = dict()
        object_ports_full = dict()
        for obj in self.object2ports_map:
            pstrs = self._get_port_information_obj(obj, d_populate = full_info)
            object_ports_full[obj] = pstrs
        return full_info, object_ports_full

    def _get_port_information_obj(self, obj, d_populate):
        imanip = PortInformationAlgorithmManipulator(
            obj = obj,
            pbg = self.pbg,
            bg_algo = self,
        )
        obj.visit_port_information(imanip)
        port_info = imanip.info
        pstrs = set()
        for pstr, info in port_info.items():
            d_populate[obj, pstr] = info
            pstrs.add(pstr)
        return pstrs

    ## unused?
    #def flag_map(self, flag_visit_func):
    #    """
    #    generates a dictionary of (obj, port) -> flags using the
    #    flag_visit_func, which must have the signature func(obj, p)->[dict(), ...]

    #    This should be relatively simple, but forwarded port names complicate it
    #    """
    #    flag_mapping = dict()
    #    for obj in self.pbg.object_iter():
    #        #TODO, make port flags be selective about delivered necessary flags
    #        for pstr, pflags in flag_visit_func(obj, self.pbg.view(obj)).items():
    #            rtup, pstr = ref_port_split(pstr)
    #            subobj = self.pbg.referred_tup(obj, rtup)
    #            flag_mapping[(subobj, pstr)] = pflags
    #    pmap = self._fully_forward_ports(flag_mapping.keys(), pmap = self.port_forward_mapping)

    #    fgroup = collections.defaultdict(list)
    #    for op, flags in flag_mapping.items():
    #        op = pmap[op]
    #        fgroup[op].extend(flags)

    #    return dict(fgroup)

    def _fully_forward_ports(self, to_check, pmap = None):
        """
        pmap is an already established port mapping. It is modified by this call
        """
        to_check = set(to_check)
        if pmap is None:
            pmap = dict()
        #this generates the entire mapping of all ports to their forwards.
        #It eagerly consumes edges
        while to_check:
            op1 = to_check.pop()
            many_map = []
            while True:
                oplink = pmap.get(op1, None)
                if oplink is not None:
                    break

                many_map.append(op1)
                obj, pstr = op1
                pret = obj.port_forward(self.pbg.view(obj), pstr)
                if pret is None:
                    oplink = op1
                    break
                elif pret is False:
                    raise RuntimeError("Port name {} of {} not acceptable".format(pstr, self.pbg.path_str(obj)))

                rtup, pstr = ref_port_split(pret)
                obj = self.pbg.referred_tup(obj, rtup)
                op1 = obj, pstr

            for op1 in many_map:
                pmap[op1] = oplink
        return pmap

    def _fully_forward_edges(self, edge_list, anti_edge_list = []):
        """
        This takes a list of edges and generates a dictionary of the edges.
        The dictionary is symmetric, as the bond edges form an undirected graph.

        the first return value bond_edges_map is the edge dictionary

        TODO, should the bond_edges_map include multiplicity? some error checking
        is possible if so.

        the second return value is a mapping of all objects to their used ports
        """
        to_check = set()
        for op1, op2 in edge_list:
            to_check.add(op1)
            to_check.add(op2)
        for op1, op2 in anti_edge_list:
            to_check.add(op1)
            to_check.add(op2)

        pmap = self._fully_forward_ports(to_check, pmap = self.port_forward_mapping)

        #TODO, at this section, certain checks and manipulations could be performed
        #replacement rules and the like

        object2ports_map = collections.defaultdict(set)
        bond_edges_map = collections.defaultdict(set)
        anti_bond_edges_map = collections.defaultdict(lambda : 0)

        for op1, op2 in anti_edge_list:
            op1 = pmap[op1]
            op2 = pmap[op2]
            if (op2, op1) in anti_bond_edges_map:
                anti_bond_edges_map[op2, op1] += 1
            else:
                anti_bond_edges_map[op1, op2] += 1

        for op1, op2 in edge_list:
            op1 = pmap[op1]
            op2 = pmap[op2]
            antiVal = anti_bond_edges_map.get((op1, op2), None)
            if antiVal is None:
                antiVal = anti_bond_edges_map.get((op2, op1), None)
                #swap the ports in the bond since the order is symmetrized
                #below, but isn't symmetric in anti_bond_edges_map
                if antiVal is not None:
                    op2, op1 = op1, op2
            do_add = True
            if antiVal is not None:
                if antiVal > 0:
                    antiVal -= 1
                    anti_bond_edges_map[op1, op2] = antiVal
                    do_add = False
            if do_add:
                bond_edges_map[op1].add(op2)
                bond_edges_map[op2].add(op1)
                object2ports_map[op1[0]].add(op1[1])
                object2ports_map[op2[0]].add(op2[1])

        for (op1, op2), antiVal in anti_bond_edges_map.items():
            if antiVal > 0:
                raise RuntimeError("More anti-bonds than bonds for edge ({}, {}).format(op1, op2)")
        return dict(bond_edges_map), dict(object2ports_map)

    def _parse_bond_list(self, obj, bl):
        """
        This takes a bond list bl centered on object obj and resolves all of the
        sub-references, ports and chains. A list of 2 elements will generate
        a single bond pair [((first_obj, port1), (last_obj, port2))].

        Any middle items in the list must be chained ports that contain a '-'
        and are resolved using subobj.port_chain into two separate ports. Those
        will then be resolved and will expand the list.
        """
        true_lst = []

        def process_str(s):
            for port_chain_str in s.split('|'):
                port_chain_str = port_chain_str.strip()
                yield port_chain_str

        if isinstance(bl, str):
            true_lst.extend(process_str(bl))
        else:  # isinstance(bl, (list, tuple)):
            for b_item in bl:
                if isinstance(b_item, str):
                    true_lst.extend(process_str(b_item))
                else:
                    for b_item2 in b_item:
                        true_lst.extend(process_str(b_item2))

        assert(len(true_lst) > 1)
        idx_firstmid = 1
        first = true_lst[0]
        is_antilist = False

        #check if this is a list of anti-edges which are denoted
        #with a bond list where the first element is "remove"
        if first == 'remove':
            idx_firstmid = 2
            first = true_lst[1]
            is_antilist = True

        last  = true_lst[-1]

        #first and last in sequence are not allowed to be chained ports
        #TODO, needs more error reporting
        if '-' in first:
            raise RuntimeError((
                "First port in a sequence cannot be a chained port {}"
            ).format(first))
        if '-' in last:
            raise RuntimeError((
                "Last port in a sequence cannot be a chained port {}"
            ).format(last))

        #this first object will be applied as a ling
        rtup, pstr = ref_port_split(first)
        subobj = self.pbg.referred_tup(obj, rtup)
        next_obj_port = subobj, pstr

        mid_pairs = []
        for mid_chain in true_lst[idx_firstmid:-1]:
            if '-' not in mid_chain:
                raise RuntimeError((
                    "port in the middle of sequence must be a chained port {}"
                ).format(mid_chain))

            rtup, pstr = ref_port_split(mid_chain)
            subobj = self.pbg.referred_tup(obj, rtup)

            #print(mid_chain, self.pbg.path_str(subobj), pstr)
            ppair = subobj.port_chain(self.pbg.view(subobj), pstr)
            if ppair is None:
                #TODO, more error reporting
                raise RuntimeError("Port Chain not recognized: {}".format(mid_chain))
            port1, port2 = ppair

            #apply convenience rules to the return values
            if port1 is None:
                port1 = pstr.split('-')[0]
            if port2 is None:
                port2 = pstr.split('-')[0]

            mid_pairs.append((next_obj_port, (subobj, port1)))
            next_obj_port = (subobj, port2)

        rtup, pstr = ref_port_split(last)
        subobj = self.pbg.referred_tup(obj, rtup)
        mid_pairs.append((next_obj_port, (subobj, pstr)))

        #TODO, make this a generator?
        return is_antilist, mid_pairs

    def pgv_draw_bg(self):
        """
        """
        import pygraphviz as pgv
        G = pgv.AGraph()
        edge_str = []
        for m_fr, to_set in self.bond_edges_map.items():
            m_fr = self.pbg.path_str(m_fr[0]) + m_fr[1]
            for m_to in to_set:
                m_to = self.pbg.path_str(m_to[0]) + m_to[1]
                G.add_edge(m_fr, m_to)
                edge_str.append((m_fr, m_to))
        for obj, pset in self.object2ports_map.items():
            m_nodes = []
            oname = self.pbg.path_str(obj)
            for pname in pset:
                m_nodes.append(oname + pname)
            G.add_subgraph(m_nodes, 'cluster_' + oname)

        edge_str.sort()
        G.draw('test2.pdf', prog = 'dot')


class BondGraphAlgorithmView(object):
    _bg_algo = None
    _obj     = None
    _p       = None
    _pbg     = None

    def __init__(self, obj, pbg, bg_algo, p = None, **kw):
        self._obj = obj
        self._bg_algo = bg_algo
        self._pbg = pbg
        if p is None:
            p = pbg.view(obj)
        self._p = p
        self.p = p

    #returns the port which were supplied to generate linkages
    def ports_used(self):
        return self._bg_algo.object2ports_map[self._obj]


class LinkageAlgorithmView(BondGraphAlgorithmView):
    def __init__(self, **kw):
        super(LinkageAlgorithmView, self).__init__(**kw)

    #returns the port which were supplied to generate linkages
    def links_used(self):
        return self._bg_algo.object_linkages[self._obj]

    def link_port(self, link):
        raise NotImplementedError()
        return self._bg_algo.object_linkages[self._obj]

    def link_port_flags(self, link):
        """
        Returns the bond-port and associated flags relevant for the linkage
        supplied
        """
        raise NotImplementedError()
        return self._bg_algo.object_linkages[self._obj]


class PortInformationAlgorithmManipulator(BondGraphAlgorithmView):
    def __init__(self, **kw):
        super(PortInformationAlgorithmManipulator, self).__init__(**kw)
        self.info = dict()

    def gen_optical_port(self, pname, lprefix):
        self.port_inform(
            pname, base.port_types.BGT_OPTICS_PHYSICAL,
            inout = (lprefix + '!i', lprefix + '!o'),
        )

    def gen_signal_in_port(self, pname, lprefix):
        self.port_inform(
            pname, base.port_types.BGT_SIGNAL_IN,
            inout = (lprefix + '!i', None),
        )

    def gen_signal_out_port(self, pname, lprefix):
        self.port_inform(
            pname, base.port_types.BGT_SIGNAL_OUT,
            inout = (None, lprefix + '!o'),
        )

    def gen_signal_io_port(self, pname, lprefix):
        self.port_inform(
            pname, base.port_types.BGT_SIGNAL_OUT,
            inout = (lprefix + '!i', lprefix + '!o'),
        )

    def port_inform(self, pname, ptype, **kw):
        idict = dict()
        idict['type'] = ptype
        idict.update(kw)

        idict_test = self.info.setdefault(pname, idict)
        if idict_test is not idict:
            raise RuntimeError("Can't inform a port multiple times")
        return

