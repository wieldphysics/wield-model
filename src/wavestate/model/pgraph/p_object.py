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
import contextlib
import warnings
import functools

try:
    from collections import abc
except ImportError:
    import collections as abc

from . import utilities
from . import p_graph_impl


@functools.total_ordering
class ParameterObjectBase(object):
    """
    base object for the parameter graph system.
    Does not include bond graph methods
    """
    _allow_multiple_parents = True
    def __init__(self):
        self._value_dict = collections.defaultdict(dict)
        self._reference_dict = collections.defaultdict(dict)
        #also used for reference defaults
        self._default_functions = collections.defaultdict(dict)
        #if None, then the settable values is unrestricted
        self._values_settable = None

        #TODO, what does the internal flag do again?
        self.internal = False
        return

    @contextlib.contextmanager
    def _internal(self):
        self.internal = True
        yield
        self.internal = False

    def values_settable(self, params):
        if self._values_settable is None:
            self._values_settable = set()
            self.internal = False
        vtups = []
        for pkey in params:
            rtup, vtup = utilities.ref_value_split(pkey)
            assert(vtup is not None)
            assert(rtup == ())
            vtups.append(vtup)
        self._values_settable.update(vtups)
        return

    def __getitem__(self, key):
        rtup, vtup = utilities.ref_value_split(key)
        if vtup is not None:
            pdict = self._value_dict[rtup, vtup]
            #TODO, error message
            return pdict['default_value']
        else:
            rdict = self._reference_dict[rtup]
            #TODO, error message
            return rdict['default_value']
        return

    def __setitem__(self, key, val):
        rtup, vtup = utilities.ref_value_split(key)
        if vtup is not None:
            if isinstance(val, ParameterObjectBase):
                warnings.warn((
                    "Parameter assignment {} to object {} is a ParameterObject"
                ).format(key, self))
            pdict = self._value_dict[rtup, vtup]
            pdict['default_value'] = val
            if self.internal:
                pdict['internal'] = True
        else:
            if isinstance(val, ParameterObjectBase):
                rdict = self._reference_dict[rtup]
                rdict['default_value'] = val
                if self.internal:
                    pdict['internal'] = True
            elif isinstance(val, str):
                if val[-1] != '/':
                    raise RuntimeError(
                        "reference assignment must be to another reference"
                    )
                self.set_assign(kto = key, kfrom = val)
            else:
                warnings.warn((
                    "Reference assignment {} to object {} is not a ParameterObject"
                ).format(key, self))
        return

    def set_value(self, key, val):
        """
        Sets a parameter or reference value directly. This is the same as the
        item assignment, except allows additional annotations.
        """
        rtup, vtup = utilities.ref_value_split(key)
        if vtup is not None:
            pdict = self._value_dict[rtup, vtup]
            pdict['default_value'] = val
            if self.internal:
                pdict['internal'] = True
        else:
            rdict = self._reference_dict[rtup]
            rdict['default_value'] = val
            if self.internal:
                pdict['internal'] = True

    def set_function(self, key, func, deps = None, dependencies = None):
        """
        Sets a parameter or reference value directly. This is the same as the
        item assignment, except allows additional annotations.
        """
        rtup, vtup = utilities.ref_value_split(key)

        if dependencies is None:
            dependencies = deps

        if vtup is not None:
            def fwrap(p):
                p[key] = func(p)

            pdict = self._value_dict[rtup, vtup]
            pdict['default_func'] = fwrap
            if self.internal:
                pdict['internal'] = True
            fdict = self._default_functions[fwrap]
            if dependencies is not None:
                fdict['deps'] = tuple(utilities.ref_value_split(key) for key in dependencies)
            fdict['assigns'] = ((rtup, vtup,),)
        else:
            rdict = self._reference_dict[rtup]
            rdict['default_func'] = func
            if dependencies is not None:
                rdict['default_func_deps'] = tuple(utilities.ref_value_split(key) for key in dependencies)
        return

    def set_assign(self, kto, kfrom, pfunc = None):
        rto, pto = utilities.ref_value_split(kto)
        rfrom, pfrom = utilities.ref_value_split(kfrom)
        if (pto is not None):
            #it is a parameter path
            assert(pfrom is not None)

            pdict_to = self._value_dict[rto, pto]
            if self.internal:
                pdict_to['internal'] = True

            if pfunc is None:
                def pwrap(p):
                    p[kto] = p[kfrom]
            else:
                #TODO, use faster lookup/assign methods that use the fact
                #that keys have already been parsed
                def pwrap(p):
                    p[kto] = pfunc(p[kfrom])

            pdict_to['default_func'] = pwrap

            #remove a value if it was assigned previously
            pdict_to.pop('default_value', None)

            fdict = self._default_functions[pwrap]
            fdict['deps'] = ((rfrom, pfrom),)
            fdict['assigns'] = ((rto, pto,),)
        else:
            #it is a reference path
            assert(pfrom is None)
            assert(pfunc is None)

            rdict_to = self._reference_dict[rto]
            rdict_to['default_assign'] = rfrom
        return

    def set_multi_function(self, assignments, func, deps = None, dependencies = None):
        """
        decorator for a parameter generating function. The name of the function
        determines where the value is assigned
        """
        if dependencies is None:
            dependencies = deps

        fdict = self._default_functions[func]
        if dependencies is not None:
            fdict['deps'] = tuple(utilities.ref_value_split(key) for key in dependencies)
        assigns = []
        for key in assignments:
            rtup, vtup = utilities.ref_value_split(key)
            if vtup is not None:
                pdict = self._value_dict[rtup, vtup]
                #remove a value if it was assigned previously
                pdict.pop('default_value', None)

                if self.internal:
                    pdict['internal'] = True
                pdict['default_func'] = func
                assigns.append((rtup, vtup))
            else:
                raise RuntimeError("Cannot multi-assign to references")
        fdict['assigns'] = ((rtup, vtup,),)
        return

    def deco_parameter(self, name = None, func = None, dependencies = None):
        """
        decorator for a parameter generating function. The name of the function
        determines where the value is assigned
        """
        if func is None:
            if name is not None:
                if not isinstance(name, str):
                    func = name
                    name = func.__name__
                return self.set_function(name, func, dependencies = dependencies)

        def deco_parameter_func(func):
            if name is None:
                name_int = func.__name__
            else:
                name_int = name

            self.set_function(
                name = name_int,
                func = func,
                dependencies = dependencies
            )
            return func
        return deco_parameter_func

    def deco_reference(self, name = None, func = None, dependencies = None):
        """
        decorator for a parameter generating function. The name of the function
        determines where the value is assigned
        """
        if func is None:
            if name is not None:
                if not isinstance(name, str):
                    func = name
                    name = func.__name__
                return self.set_function(name + ':', func, dependencies = dependencies)

        def deco_parameter_func(func):
            if name is None:
                name_int = func.__name__ + ':'
            else:
                name_int = name + ':'
            self.set_function(
                name = name_int,
                func = func,
                dependencies = dependencies
            )
            return func
        return deco_parameter_func

    def deco_multi_parameter(self, assignments, func = None, dependencies = None):
        """
        decorator for a multi-parameter generating function
        """
        if func is not None:
            return self.set_multi_function(
                assignments = assignments,
                func = func,
                dependencies = dependencies
            )
        else:
            def deco_multi_parameter_func(func):
                self.set_multi_function(
                    assignments = assignments,
                    func = func,
                    dependencies = dependencies
                )
                return func
            return deco_multi_parameter_func

    def deco_many_many(self, assignments, dependencies, func = None):
        """
        decorator for a multi-parameter generating function
        """
        #TODO, pre-build the ktup/vtups and use faster p-accessors
        def setup_function(func):
            if isinstance(dependencies, abc.Mapping):
                if isinstance(assignments, abc.Mapping):
                    def fwrap(p):
                        kw = dict()
                        for k, pkey in dependencies.items():
                            kw[k] = p[pkey]
                        ret = func(**kw)
                        for k, pkey in assignments.items():
                            p[pkey] = ret[k]
                    self.set_multi_function(
                        func = fwrap,
                        assignments = assignments.value(),
                        dependencies = dependencies.values(),
                    )
                else:
                    def fwrap(p):
                        kw = dict()
                        for k, pkey in dependencies.items():
                            kw[k] = p[pkey]
                        ret = func(**kw)
                        for i, pkey in enumerate(assignments):
                            p[pkey] = ret[i]
                    self.set_multi_function(
                        assignments = assignments,
                        dependencies = dependencies.values(),
                        func = fwrap,
                    )
            else:
                if isinstance(assignments, abc.Mapping):
                    def fwrap(p):
                        args = tuple(p[pkey] for pkey in dependencies)
                        ret = func(*args)
                        for k, pkey in assignments.items():
                            p[pkey] = ret[k]
                    self.set_multi_function(
                        func = fwrap,
                        assignments = assignments.value(),
                        dependencies = dependencies,
                    )
                else:
                    def fwrap(p):
                        args = tuple(p[pkey] for pkey in dependencies)
                        ret = func(*args)
                        for i, pkey in enumerate(assignments):
                            p[pkey] = ret[i]
                    self.set_multi_function(
                        assignments = assignments,
                        dependencies = dependencies,
                        func = fwrap,
                    )
            return func

        if func is not None:
            return setup_function(func)
        else:
            def deco_many_many_func(func):
                return setup_function(func)
            return deco_many_many_func

    def deco_many_one(self, dependencies, assignment = None, func = None):
        """
        decorator for a multi-parameter generating function
        """
        #TODO, pre-build the ktup/vtups and use faster p-accessors
        def setup_function(func):
            if assignment is None:
                name = func.__name__
            else:
                name = assignment

            if isinstance(dependencies, abc.Mapping):
                def fwrap(p):
                    kw = dict()
                    for k, pkey in dependencies.items():
                        kw[k] = p[pkey]
                    p[name] = func(**kw)
                self.set_multi_function(
                    func = fwrap,
                    assignments = [name],
                    dependencies = dependencies.values(),
                )
            else:
                def fwrap(p):
                    args = tuple(p[pkey] for pkey in dependencies)
                    p[name] = func(*args)
                self.set_multi_function(
                    func = fwrap,
                    assignments = [name],
                    dependencies = dependencies.values(),
                )

        if func is not None:
            return setup_function(func)
        else:
            def deco_many_one_func(func):
                return setup_function(func)
            return deco_many_one_func

    def deco_one_one(self, kfrom, func = None):
        """
        decorator for a multi-parameter generating function
        """
        #TODO, pre-build the ktup/vtups and use faster p-accessors
        def setup_function(func):
            kto = func.__name__

            self.set_assign(
                pfunc = func,
                kto = kto,
                kfrom = kfrom,
            )

        if func is not None:
            return setup_function(func)
        else:
            def deco_many_one_func(func):
                return setup_function(func)
            return deco_many_one_func

    def __hash__(self):
        return id(self)

    def __lt__(self, other):
        return id(self) < id(other)

    def __eq__(self, other):
        return id(self) == id(other)


class ParameterObject(ParameterObjectBase):
    def __init__(self):
        super(ParameterObject, self).__init__()
        self._bond_generators = set()
        self._bond_lists      = []
        self._port_forwards   = dict()
        self._port_chains     = dict()

    def port_chain(self, p, pname):
        return self._port_chains.get(pname, None)

    def port_forward_add(self, pfrom, pto):
        self._port_forwards[pfrom] = pto

    def port_chain_add(self, pfrom, pto, pchain):
        self._port_chains[pfrom] = (pto, pchain)

    def port_forward(self, p, pname):
        return self._port_forwards.get(pname, None)

    def port_flags(self, p):
        """
        Return flags that get forwarded into subobjects. The flags are
        to recognize objects within the bond-graph (they annotate the graph)
        """
        return dict()

    @classmethod
    def visit_port_information(cls, algo):
        """
        Subobjects must implement this
        """
        if algo.ports_used():
            raise RuntimeError(
                ("Unrecognized ports, likely a naming error in port forwarding {}"
                 ).format(algo.ports_used())
            )

    def bond_add(self, *blists):
        """
        Bond lists are...
        """
        for blist in blists:
            self._bond_lists.append(blist)

    def bond_lists(self, p):
        #requires the parameter graph to evaluate values
        lst = list(self._bond_lists)
        for gen in self._bond_generators:
            lst.extend(gen(p))
        return lst

    def add_bond_generator(self, func = None):
        """
        decorator for a bond generating function
        """
        self._bond_generators.add(func)

    def deco_bond_generator(self, func = None):
        """
        decorator for a bond generating function
        """
        if func is not None:
            self.add_bond_generator(func)
        else:
            def deco_bond_generator_func(func):
                self.add_bond_generator(func)
                return func
            return deco_bond_generator_func

    def __repr__(self):
        pbg = p_graph_impl._principle_graph.pbg
        if pbg is None:
            return super(ParameterObjectBase, self).__repr__()

        try:
            refstr = pbg.path_str(self)
            return ":{}:".format(refstr)
        except Exception:
            return super(ParameterObjectBase, self).__repr__()


