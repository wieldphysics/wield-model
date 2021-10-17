#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


#from .utilities import ref_value_split, ref_port_split, ref_2_rtup
from . import utilities
from .p_graph_impl import ParameterGraphImpl


class ParameterGraph(ParameterGraphImpl):

    def override_value(self, p_ref, val, obj = None):
        if obj is None:
            obj = self.root
        rtup, vtup = utilities.ref_value_split(p_ref)
        ref_obj = self.referred_tup(obj, rtup)
        return self.override_value_op((ref_obj, vtup), val)

    def override_value_op(self, op, val):
        self.op_reqs[op].clear()
        obj, vtup = op

        self.object_parameter_values[obj][vtup] = val

        for op_dep in self.op_deps[op]:
            self.clear_value_op(op_dep)
        return

    def resolve_object(self, key, obj = None):
        """
        Takes a reference key and resolves the specific object
        """
        if obj is None:
            obj = self.root
        rtup = utilities.ref_2_rtup(key)
        for ref in rtup:
            obj = self.object_references[obj][ref]
        #TODO, make better error handling
        return obj

    def resolve_reference(self, key, obj = None):
        if obj is None:
            obj = self.root
        rtup, vtup = utilities.ref_value_split(key)
        for ref in rtup:
            obj = self.object_references[obj][ref]
        return obj, vtup

    def get_parameter(self, key, obj = None):
        if obj is None:
            obj = self.root
        return self.resolve_parameter(obj, key)

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.get_parameter(key)
        elif isinstance(key, tuple):
            obj, key = key
            return self.get_parameter(key, obj = obj)
        elif isinstance(key, slice):
            obj, key = key
            return self.get_parameter(key, obj = obj)
        else:
            raise RuntimeError('Unrecognized parameter index format')

    def resolve_parameter(self, obj, key):
        """
        This resolves a local parameter where obj is known
        """
        rtup, vtup = utilities.ref_value_split(key)
        for ref in rtup:
            obj = self.object_references[obj][ref]
        val = self._resolve_parameter(obj, vtup)
        return val

    def clear_value_op(self, op):
        self.op_reqs[op].clear()
        obj, vtup = op

        #safe delete
        self.object_parameter_values[obj].pop(vtup)

        for op_dep in self.op_deps[op]:
            self.clear_value_op(op_dep)
        return

    def port_ref(self, bport_ref, obj = None):
        if obj is None:
            obj = self.root
        rtup, bpname = utilities.ref_port_split(bport_ref)
        obj = self.referred_tup(obj, rtup)
        return obj, bpname

    def object_iter(self):
        return self.object_path_built.keys()

    def path_str(self, obj, vtup = None):
        opath = self.object_path_built[obj]
        if vtup is not None:
            return ('/'.join(opath) + '/' + '.'.join(vtup))
        else:
            return ('/'.join(opath) + '/')

    _path_suffix_trie = None
    def path_str_short(self, obj, pre = 1, vtup = None):
        if self._path_suffix_trie is None:
            self._path_suffix_trie = self._build_suffix_trie()
        opath = self.object_path_built[obj]
        trie = self._path_suffix_trie
        ptup = []
        for p in opath[::-1]:
            ptup.append(p)
            trie = trie[p]
            if len(trie) == 1:
                if pre <= 0:
                    break
                else:
                    pre -= 1

        if len(ptup) < len(opath):
            opath = tuple(['..'] + ptup[::-1])
        else:
            opath = tuple(ptup[::-1])
        if vtup is not None:
            return ('/'.join(opath) + '/' + '.'.join(vtup))
        else:
            return ('/'.join(opath) + '/')

    def _build_suffix_trie(self):
        trie_orig = dict()
        for obj, ptup in self.object_path_built.items():
            trie = trie_orig
            for p in ptup[::-1]:
                trie = trie.setdefault(p, dict())
        return trie_orig

    def _build_prefix_trie(self):
        trie_orig = dict()
        for obj, ptup in self.object_path_built.items():
            trie = trie_orig
            for p in ptup:
                trie = trie.setdefault(p, dict())
        return trie_orig

    def op_str(self, op):
        #TODO, split this between the parameter and port-tup version
        obj, path_or_port = op
        if isinstance(path_or_port, (list, tuple)):
            path_or_port = '.'.join(path_or_port)
        ostr = self.path_str(obj)
        return ostr + path_or_port

    def opset2ref_deptree(self, opset):
        """
        assumes obj, vtup input as a set
        """
        deptree = dict()
        ops_remaining = set(opset)
        ops_done = set()
        while ops_remaining:
            op = ops_remaining.pop()
            if op in ops_done:
                continue
            opstr = self.op_str(op)
            depstrs = set()
            for opdep in self.op_deps[op]:
                depstr = self.op_str(opdep)
                depstrs.add(depstr)
                if opdep not in ops_done:
                    ops_remaining.add(opdep)
            deptree[opstr] = depstrs
        return deptree

    def referred_path(self, obj, key):
        rtup, vtup = utilities.ref_value_split(key)
        assert(vtup is None)
        for ref in rtup:
            obj = self.object_references[obj][ref]
        return obj

    def referred_vtup(self, pkey, obj = None):
        rtup, vtup = utilities.ref_value_split(pkey)
        if obj is None:
            obj = self.root
        for ref in rtup:
            obj = self.object_references[obj][ref]
        return obj, vtup

    def referred_tup(self, obj, rtup):
        for ref in rtup:
            obj = self.object_references[obj][ref]
        return obj

    def resolve_parameter_default(self, obj, vtup, default):
        try:
            return self._resolve_parameter(obj, vtup)
        except KeyError:
            return default

    def print_values(self):
        for obj, vtup_dict in self.object_values.items():
            opath = self.object_path_built[obj]
            for vtup, cdict in vtup_dict.items():
                print('/'.join(opath) + '/' + '.'.join(vtup))

    def print_values_eval(self, obj = None):
        def oprint(obj):
            opath = self.object_path_built[obj]
            for vtup, cdict in vtup_dict.items():
                print('/'.join(opath) + '/' + '.'.join(vtup), self._resolve_parameter(obj, vtup))
        if obj is None:
            for obj, vtup_dict in self.object_values.items():
                oprint(obj)
        else:
            oprint(obj)

    def print_tree(self, print=print):
        ptree = self._build_prefix_trie()
        import yaml
        s = yaml.dump(ptree)
        print(s)
        return s

    def dict_values_eval(self, obj = None):
        pdict = dict()

        def odict(obj):
            opath = self.object_path_built[obj]
            for vtup, cdict in vtup_dict.items():
                key = '/'.join(opath) + '/' + '.'.join(vtup)
                pdict[key] = self._resolve_parameter(obj, vtup)

        if obj is None:
            for obj, vtup_dict in self.object_values.items():
                odict(obj)
        else:
            odict(obj)
        return pdict

