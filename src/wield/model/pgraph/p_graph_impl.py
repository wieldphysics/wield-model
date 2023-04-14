#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


import numpy as np
import collections
import threading
import contextlib

from wavestate.model.utilities import constants

from . import utilities


# this is used for specifying a preferred graph to inspect during debugging
# knowing the graph allows objects to be named
_principle_graph = threading.local()
_principle_graph.pbg = None


class ParameterGraphImpl(object):
    def __init__(self, root):
        # if root is a ParameterGraphImpl, then do a copy constructor
        self.root = root
        self.object_values = dict()
        self.object_references = dict()
        self.object_parameter_values = collections.defaultdict(dict)
        self.op_deps = collections.defaultdict(set)
        self.op_reqs = collections.defaultdict(set)
        self.object_functions = dict()
        self._path_suffix_trie = None

        # mapping of each object to its multiple paths within root
        self.object_paths = dict()

        # mapping of object to the unique path it was built on (not guaranteed stable)
        self.object_path_built = dict()

        # full of (obj, rtup) tuples
        objects = self._resolve_reference_tree((), root)
        self._resolve_parameter_tree(objects)
        return

    def copy(self):
        # Invokes a copy constructor by calling __new__ on its own
        newself = self.__class__.__new__(self.__class__)
        other = self
        # copy constructor
        newself.root = other.root
        newself.object_values = dict(other.object_values)
        newself.object_references = dict(other.object_references)
        newself.object_parameter_values = collections.defaultdict(dict)
        newself.op_deps = collections.defaultdict(set)
        newself.op_reqs = collections.defaultdict(set)
        newself.object_functions = dict(other.object_functions)

        for k, v in other.object_parameter_values.items():
            newself.object_parameter_values[k].update(v)
        for k, v in other.op_deps.items():
            newself.op_deps[k].update(v)
        for k, v in other.op_reqs.items():
            newself.op_reqs[k].update(v)

        newself.object_paths = dict(other.object_paths)
        newself.object_path_built = dict(other.object_path_built)
        return newself

    def _invalidate(self):
        self._path_suffix_trie = None
        return

    def object_insert(self, rtup, obj):
        """
        Inserts the object at rtup and obj, then resolves it
        """
        objects = self._resolve_reference_tree(rtup, obj)
        self._resolve_parameter_tree(objects)
        return

    def _resolve_reference_tree(self, rtup, root):
        """ """
        objects_inserted = []
        obj_ref_stack = collections.deque()

        # an RI dep is a (object, (rtup)) tuple that can efficiently index
        # object_references and refererences_intermediate
        references_intermediate = collections.defaultdict(
            lambda: collections.defaultdict(lambda: dict())
        )

        Pid_to_obj = dict()
        Pid_to_RI_path = dict()
        Pid_to_cdict = dict()
        Pid_to_Pid_deps = dict()
        Pid_to_Pid_pdeps = dict()
        Pid_to_Pid_req = dict()
        # schema
        # references_intermediate = {'obj': {'ref' : {('r1', 'r2', ...) : dict or from ref_intermed}}}

        def insert(obj):
            # print("INSERT: ", self.object_path_built[obj])
            rdict = {}
            rdict_ret = self.object_references.setdefault(obj, rdict)
            # could test if object already inserted
            if rdict_ret is rdict:
                objects_inserted.append(obj)
                # newly inserted, so generate mappings
                for rtup, cdict in obj._reference_dict.items():
                    dassign = cdict.get("default_assign", None)
                    obj, rtup, Pid = resolve_Pid(obj, rtup)
                    Pid_to_RI_path[Pid] = (obj, rtup)
                    cdict_prev = Pid_to_cdict.setdefault(Pid, cdict)
                    if cdict_prev is not cdict:
                        # this means that it had already been assigned a
                        # construction or alias from a higher system
                        print(cdict_prev, cdict)
                        print("overridden: ", self.object_path_built[obj], rtup)
                    else:
                        # print("New: ", self.object_path_built[obj], rtup)

                        # must be a new reference
                        # check if it is an assignment type
                        if dassign is not None:
                            rtup_ref = dassign
                            obj_ref, rtup_ref, Pid_ref = resolve_Pid(obj, rtup_ref)
                            Pid_to_Pid_req[Pid] = Pid_ref

                            # It's possible that deps were assigned before a
                            # construction plan, so get them and remove the set
                            # but we should add to the pool, removing our own
                            local_deps = Pid_to_Pid_deps.pop(Pid, ())
                            req_deps = Pid_to_Pid_deps.setdefault(Pid_ref, set())
                            req_deps.add(Pid)
                            req_deps.update(local_deps)

                            # to the same for the parent deps (removing them in the process)
                            local_pdeps = Pid_to_Pid_pdeps.pop(Pid, ())
                            req_pdeps = Pid_to_Pid_pdeps.setdefault(Pid_ref, set())
                            req_pdeps.update(local_pdeps)

                            # TODO, there is some referencing logic that isn't quite right here
                            if len(rtup) == 0:
                                pass
                            elif len(rtup) == 1:
                                # TODO, may not need to insert this
                                # req_pdeps.add(Pid)
                                pass
                            else:
                                _, _, Pid_p = resolve_Pid(obj, rtup[:-1])
                                req_pdeps.add(Pid_p)

                            # all of the local deps must be update to point
                            # to the root of the assignments
                            for lPid in local_deps:
                                Pid_to_Pid_req[lPid] = Pid_ref
                        else:
                            # if it was not linked, then if it is short, add it
                            # to the resolution stack for evaluation
                            if len(rtup) == 1:
                                obj_ref_stack.append((obj, rtup[0]))

                pass
            else:
                # TODO, ensure that objects cannot be multiply inserted unless
                # through a multiple assignment
                raise RuntimeError(
                    "Cannot have the same object inserted during construction"
                    " must be bound through linkages"
                )
            return

        def resolve_Pid(obj, rtup):
            while True:
                obj2 = self.object_references[obj].get(rtup[0], None)
                if obj2 is None:
                    break
                rtup = rtup[1:]
                obj = obj2
                if not rtup:
                    break
            Pid = references_intermediate[obj][rtup[0]].get(rtup[1:], None)
            if Pid is None:
                Pid = object()
                references_intermediate[obj][rtup[0]][rtup[1:]] = Pid
            return obj, rtup, Pid

        self.object_paths[root] = set()
        self.object_path_built[root] = rtup
        insert(root)
        while obj_ref_stack:
            obj, ref = obj_ref_stack.pop()
            ridict = references_intermediate[obj][ref]
            Pid = ridict.get((), None)
            if Pid is None:
                # print("HMMMM (no Pid?): ", self.object_path_built[obj], ref)
                # this should be OK, and simply means that the object ref
                # was created already
                continue
            cdict = Pid_to_cdict[Pid]
            if cdict is None:
                print("(no cdict?): ", self.object_path_built[obj], ref)

            # print("PREX_INSERT: ", self.object_path_built[obj], ref)

            req = Pid_to_Pid_req.get(Pid, None)
            if req is not None:
                # then this is a dependent path waiting for req to build
                # it will be resolved on the req and should not be in the
                # obj_ref_stack (so don't return it there)
                # this is the most efficient place to clear these
                continue

            pdeps = Pid_to_Pid_pdeps.get(Pid, None)
            if pdeps is not None:
                # check if the pdep has been built
                for pdep_Pid in list(pdeps):
                    if Pid_to_obj.get(pdep_Pid, None) is not None:
                        pdeps.remove(pdep_Pid)
            if pdeps:
                # need to get more parents before this one is ready
                # should go back in the queue though
                obj_ref_stack.appendleft((obj, ref))
                continue

                # then this is path waiting for all req parents to build

            # print("PRE_INSERT: ", self.object_path_built[obj], ref)

            # now assumes that the build is a dictionary
            dfunc = cdict.get("default_func", None)
            if dfunc is not None:
                # TODO, refererences during build
                raise NotImplementedError()
            else:
                robject = cdict["default_value"]

            self.object_references[obj][ref] = robject
            self.object_path_built[robject] = self.object_path_built[obj] + (ref,)
            self.object_paths[robject] = set([(obj, ref)])

            # resolve all dependencies at once!
            Pid_to_obj[Pid] = robject

            # include current ridict
            deps_RIdicts = [ridict]
            # delete for the current object the way that "pop" does for the deps
            del references_intermediate[obj][ref]
            for dep_Pid in Pid_to_Pid_deps.get(Pid, ()):
                dep_obj, dep_rtup = Pid_to_RI_path[dep_Pid]
                dep_obj2, dep_rtup2, dep_Pid2 = resolve_Pid(dep_obj, dep_rtup)
                # HMMM, maybe can't be confirmed in-general
                # assert(dep_Pid2 is dep_Pid)
                # must be since all parents are resolved
                # print("DEPRTUP", dep_rtup, dep_rtup2)
                assert len(dep_rtup2) == 1
                # assign the reference for the dependency
                if not robject._allow_multiple_parents:
                    # TODO, improve error message
                    raise RuntimeError(
                        "Object {} does not allow multiple"
                        " parents (reference aliasing)".format(
                            self.object_path_built[robject]
                        )
                    )
                self.object_references[dep_obj2][dep_rtup2[0]] = robject
                self.object_paths[robject].add((dep_obj2, dep_rtup2[0]))
                # should remove the old reference with a pop
                deps_RIdicts.append(references_intermediate[dep_obj2].pop(dep_rtup2[0]))
                Pid_to_obj[dep_Pid] = robject

            # now have to synthesize a new one from the others
            synth_ridict = dict()
            for ridict_src in deps_RIdicts:
                for rtup_new, Pid_new in ridict_src.items():
                    if rtup_new == ():
                        # this is the currently, assigned object, so don't analyze it
                        continue
                    # the logic in here must replicate much of what insert does
                    Pid_test = synth_ridict.setdefault(rtup_new, Pid_new)
                    if Pid_test is not Pid_new:
                        # ok, must alias over then
                        cdict_new = Pid_to_cdict.get(Pid_new, None)
                        cdict_test = Pid_to_cdict.get(Pid_test, None)
                        if cdict_new is not None:
                            if cdict_test is None:
                                # so forward the aliased sub-item
                                Pid_to_cdict[Pid_test] = cdict_test
                                Pid_to_Pid_req[Pid_test] = Pid_to_Pid_req[Pid_new]
                            else:
                                raise RuntimeError("Multi-aliasing")
                        # no need to check Pid_req since it only exists if cdict
                        # does and that was just tested
                        Pid_req_test = Pid_to_Pid_req.get(Pid_test, None)
                        if Pid_req_test is not None:
                            Pid_deps_new = Pid_to_Pid_deps.pop(Pid_new, None)
                            Pid_deps_test = Pid_to_Pid_deps.get(Pid_test, None)
                            assert Pid_deps_test is None
                            if Pid_deps_new is not None:
                                Pid_to_Pid_deps.setdefault(Pid_req_test, set()).update(
                                    Pid_deps_new
                                )
                        else:
                            Pid_deps_new = Pid_to_Pid_deps.get(Pid_new, None)
                            Pid_deps_test = Pid_to_Pid_deps.get(Pid_test, None)
                            if Pid_deps_new is not None:
                                if Pid_deps_test is None:
                                    Pid_to_Pid_deps[Pid_test] = Pid_deps_test
                                    Pid_to_Pid_pdeps[Pid_test] = Pid_to_Pid_pdeps[
                                        Pid_new
                                    ]
                                else:
                                    Pid_deps_test.update(Pid_deps_new)
                                    Pid_to_Pid_pdeps[Pid_test].update(
                                        Pid_to_Pid_pdeps[Pid_new]
                                    )

            # null tuple already removed so all rtups have elements
            for rtup, Pid in synth_ridict.items():
                # print(rtup)
                references_intermediate[robject][rtup[0]][rtup[1:]] = Pid
                if len(rtup) == 1:
                    obj_ref_stack.append((robject, rtup[0]))

            insert(robject)

        # TODO, do tests that the references have been depleted
        do_assert = False
        for obj, refdict in references_intermediate.items():
            for ref, rtupdict in refdict.items():
                for rtup, Pid in rtupdict.items():
                    print("LEFT BEHIND:", obj, ref, rtup, Pid)
                    do_assert = True
        if do_assert:
            raise RuntimeError(
                "Remaining Items left in the references dictionary (check for cycles/bugs/missing constructions?)"
            )
        return objects_inserted

    def _resolve_parameter_tree(self, objects):
        """ """
        object_values_sources = collections.defaultdict(
            lambda: collections.defaultdict(list)
        )

        def resolve_obj_rseq(obj, rtup):
            rseq = [obj]
            for ref in rtup:
                obj = self.object_references[obj][ref]
                rseq.append(obj)
            return rseq

        def resolve_obj(obj, rtup):
            for ref in rtup:
                obj = self.object_references[obj][ref]
            return obj

        object_parameter_values = collections.defaultdict(dict)
        object_values = collections.defaultdict(dict)
        object_functions = dict()

        def transfer_cdict(
            active_obj,
            context_obj_rtup_vtup,
            context_rseq,
            vtup,
            original_cdict,
        ):
            cdict = dict(original_cdict)
            context_obj, context_rtup, context_vtup = context_obj_rtup_vtup

            psettable = active_obj._values_settable
            if not cdict.get("internal", False) and psettable is not None:
                if vtup not in psettable:
                    raise RuntimeError(
                        (
                            "Parameter '{}' of object '{}' is not in the assignable set."
                            " most likely it has bad casing or spelling."
                            " Assignable values are: {}"
                        ).format(
                            ".".join(vtup),
                            self.path_str(active_obj),
                            [".".join(_pt) for _pt in active_obj._values_settable],
                        )
                    )
            cdict["context_obj"] = context_obj
            object_values[obj][vtup] = cdict

            ofunc = original_cdict.get("default_func", None)
            if ofunc is not None:
                ofunc_d_test = dict()
                ofunc_d = object_functions.setdefault(
                    (context_obj, ofunc), ofunc_d_test
                )
                if ofunc_d_test is ofunc_d:
                    # must set up the dictionary then
                    ofunc_d.update(context_obj._default_functions[ofunc])
                    ofunc_d["active_assignments"] = dict()
                ofunc_d["active_assignments"][
                    (context_rtup, context_vtup)
                ] = context_rseq
            UNIQUE = original_cdict
            oval = original_cdict.get("default_value", UNIQUE)
            if oval is not UNIQUE:
                object_parameter_values[active_obj][vtup] = oval

        t_seq_cdict = collections.namedtuple(
            "seq_cdict", ["rseq", "cdict", "origin_obj_rtup_vtup"]
        )
        for obj in objects:
            # print("B: ", obj)
            for (rtup, vtup), cdict in obj._value_dict.items():
                # print('A', obj, rtup, vtup)
                obj_rseq = tuple(resolve_obj_rseq(obj, rtup))
                # could/should also store the source object and key command for debugging..
                object_values_sources[obj_rseq[-1]][vtup].append(
                    t_seq_cdict(obj_rseq, cdict, (obj, rtup, vtup))
                )

        for obj, pdict in object_values_sources.items():
            for vtup, seq_cdict_lst in pdict.items():
                # now, must determine the principle cdict
                # cdict_obj carrying None are alias cancelations
                seq_cdict_lst.sort(key=lambda seq_cdict_tup: -len(seq_cdict_tup.rseq))
                # shortest last, for popping off
                did_cover = set()
                not_covered = list()
                not_covered_none = list()
                while seq_cdict_lst:
                    scd = seq_cdict_lst.pop()
                    # print("rseq", [self.object_path_built[_o] for _o in scd.rseq])
                    # needs to be cancelled or overridden
                    for scd_test in seq_cdict_lst:
                        # print("rseq_test", [self.object_path_built[_o] for _o in scd_test.rseq])
                        if scd_test.rseq[-len(scd.rseq) :] == scd.rseq:
                            did_cover.add(scd_test.rseq)
                            # was masked by something
                            break
                    else:
                        if cdict is None:
                            if scd.rseq not in did_cover:
                                raise RuntimeError("Parameter Cover was not utilized")
                            not_covered_none.append(scd)
                        else:
                            not_covered.append(scd)

                if len(not_covered) == 0:
                    raise RuntimeError("Parameter missing a cdict!")
                elif len(not_covered) > 1:
                    raise RuntimeError("Parameter needs cancellation")
                scd = not_covered[0]
                # should maybe test these..
                # for scd_test in not_covered_none:
                #    if cdict is None:
                #        #TODO, improve error handling
                #        raise RuntimeError("Bad Alias Cancellation")
                transfer_cdict(
                    active_obj=obj,
                    context_obj_rtup_vtup=scd.origin_obj_rtup_vtup,
                    context_rseq=scd.rseq,
                    vtup=vtup,
                    original_cdict=scd.cdict,
                )

        self.object_parameter_values.update(object_parameter_values)
        self.object_values.update(object_values)
        self.object_functions.update(object_functions)
        return

    def _resolve_parameter(self, obj, vtup):
        pdict = self.object_parameter_values[obj]
        try:
            return pdict[vtup]
        except KeyError:
            # ok, must build it
            pass

        cdict = self.object_values[obj][vtup]
        # print(self.path_str(obj, vtup = vtup), cdict)
        ofunc = cdict["default_func"]
        ofunc_obj = cdict["context_obj"]
        fdict = self.object_functions[ofunc_obj, ofunc]
        fview = ParameterGraphFunctionView(self, ofunc_obj, fdict["active_assignments"])
        ofunc(fview)

        ret_val = fview.assign_dict[obj, vtup]
        # transfer the values over
        for (aobj, vtup), val in fview.assign_dict.items():
            op_param = aobj, vtup
            self.object_parameter_values[aobj][vtup] = val
            for op_acc in fview.access_dict:
                self.op_deps[op_param].add(op_acc)
                self.op_reqs[op_acc].add(op_param)
        return ret_val

    def view(self, obj, access_opset=None):
        """
        Returns a view of the parameter graph suitable for function evaluations
        """
        return ParameterGraphView(self, obj, access_opset=access_opset)

    @contextlib.contextmanager
    def preferred(self):
        prev = _principle_graph.pbg
        _principle_graph.pbg = self
        yield
        _principle_graph.pbg = prev
        return


class ParameterGraphView(object):
    """
    View object for ParameterGraph traversals
    """

    def __init__(self, pg, context_obj, access_opset=None):
        self.pg = pg
        self.access_opset = access_opset
        self.context_obj = context_obj

    def __getitem__(self, key):
        rtup, vtup = utilities.ref_value_split(key)
        obj = self.context_obj
        for ref in rtup:
            obj = self.pg.object_references[obj][ref]
        if self.access_opset is not None:
            self.access_opset.add((obj, vtup))
        return self.pg._resolve_parameter(obj, vtup)

    math = np
    const = constants.constants_floats


class ParameterGraphFunctionView(object):
    """
    View object for ParameterGraph calls from parameter resolving functions
    """

    def __init__(self, pg, context_obj, active_assign_dict):
        self.pg = pg
        self.active_assign_dict = active_assign_dict
        self.access_dict = dict()
        self.assign_dict = dict()
        self.context_obj = context_obj

    def __getitem__(self, key):
        rtup, vtup = utilities.ref_value_split(key)
        obj = self.context_obj
        for ref in rtup:
            obj = self.pg.object_references[obj][ref]
        val = self.pg._resolve_parameter(obj, vtup)
        self.access_dict[(obj, vtup)] = val
        return val

    def __setitem__(self, key, val):
        rtup, vtup = utilities.ref_value_split(key)
        try:
            rseq = self.active_assign_dict[(rtup, vtup)]
            self.assign_dict[rseq[-1], vtup] = val
        except KeyError:
            print("Not active (dbg)", self.context_obj, key)
        # obj = self.context_obj
        # rseq = [obj]
        # for ref in rtup:
        #    obj = self.object_references[obj][ref]
        #    rseq.append(obj)
        # rseq = tuple(rseq)
        return
