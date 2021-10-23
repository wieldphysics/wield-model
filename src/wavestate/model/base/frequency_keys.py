#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""
try:
    from collections.abc import Mapping as MappingABC
except ImportError:
    from collections import Mapping as MappingABC


class FrequencyKey(object):
    __slots__ = ("F_dict", "_cache_hash", "_cache_tup")

    def __init__(self, F_dict):
        self.F_dict = F_dict

    def DC_is(self):
        return not self.F_dict

    def ord(self):
        return sum(abs(v) for v in self.F_dict.values())

    def __hash__(self):
        try:
            return self._cache_hash
        except AttributeError:
            pass
        self._cache_hash = hash(self.hash_tuple())
        return self._cache_hash

    def hash_tuple(self):
        try:
            return self._cache_tup
        except AttributeError:
            pass
        self._cache_tup = tuple(
            sorted((f, n) for f, n in list(self.F_dict.items()) if n != 0)
        )
        return self._cache_tup

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.F_dict == other.F_dict

    def __getitem__(self, Fname):
        return self.F_dict[Fname]

    # def frequency(self):
    #    F_sum = 0
    #    for Fname, n in list(self.F_dict.items()):
    #        if n != 0:
    #            #TODO maybe move the .val somewhere else
    #            F_sum += Fname.F_Hz.val * n
    #    return F_sum

    def __repr__(self):
        l = tuple(sorted(((Fname, n) for Fname, n in list(self.F_dict.items()))))
        flist = []
        for Fname, n in l:
            if n == 1:
                flist.append("+" + Fname)
            elif n == -1:
                flist.append("-" + Fname)
            elif n > 1:
                flist.append("+" + str(n) + Fname)
            elif n < -1:
                flist.append(str(n) + Fname)
        if not flist:
            flist.append("0")
        return "".join(flist)

    def __add__(self, other):
        F_dict = dict(self.F_dict)
        for Fname, n in list(other.F_dict.items()):
            current_idx = self.F_dict.get(Fname, 0)
            new_idx = current_idx + n
            if new_idx == 0 and Fname in F_dict:
                del F_dict[Fname]
            else:
                F_dict[Fname] = new_idx
        return self.__class__(F_dict)

    def __sub__(self, other):
        F_dict = dict(self.F_dict)
        for Fname, n in list(other.F_dict.items()):
            current_idx = self.F_dict.get(Fname, 0)
            new_idx = current_idx - n
            if new_idx == 0 and Fname in F_dict:
                del F_dict[Fname]
            else:
                F_dict[Fname] = new_idx
        return self.__class__(F_dict)

    def __rmul__(self, other):
        F_dict = dict()
        for Fname, n in list(self.F_dict.items()):
            F_dict[Fname] = other * n
        return self.__class__(F_dict)

    def __neg__(self):
        F_dict = dict()
        for Fname, n in list(self.F_dict.items()):
            current_idx = self.F_dict[Fname]
            F_dict[Fname] = -current_idx
        return self.__class__(F_dict)
