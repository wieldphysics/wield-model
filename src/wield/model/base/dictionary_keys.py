#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""
try:
    from collections.abc import Mapping as MappingABC
except ImportError:
    from collections import Mapping as MappingABC
from wield import declarative


class DictKey(MappingABC):
    __slots__ = ("_dict", "_cache_hash")

    def __init__(self, __dict=None, **kwargs):
        self._dict = kwargs
        if __dict:
            self._dict.update(__dict)

    def copy_update(self, **kwargs):
        newdict = dict(self._dict)
        newdict.update(**kwargs)
        return self.__class__(**newdict)

    def __hash__(self):
        try:
            return self._cache_hash
        except AttributeError:
            pass
        l = tuple(sorted(self._dict.items()))
        self._cache_hash = hash(l)
        return self._cache_hash

    def __iter__(self):
        return iter(self._dict)

    def __getitem__(self, key):
        return self._dict[key]

    def get(self, key, default=declarative.NOARG):
        if default is declarative.NOARG:
            return self._dict[key]
        return self._dict.get(key, default)

    def __len__(self):
        return len(self._dict)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self._dict == other._dict

    def __lt__(self, other):
        if not isinstance(other, self.__class__):
            return False
        # TODO this is probably megaslow. Should likely use id or hash
        l1 = tuple(sorted(self._dict.items()))
        o1 = tuple(sorted(self._dict.items()))
        return l1 < o1

    # @repr_compat
    def __repr__(self):
        l = tuple(sorted(self._dict.items()))
        # print(unicode(repr(l), 'utf-8'))
        l2 = ["{0}:{1}".format(i, j) for i, j in l]
        return "DK{{{0}}}".format("|".join(l2))

    def __or__(self, other):
        cp = dict(self._dict)
        cp.update(other._dict)
        return self.__class__(cp)

    def iteritems(self):
        return list(self._dict.items())

    def kv_contains(self, k, v):
        v2 = self._dict.get(k, declarative.NOARG)
        return v == v2

    def __and__(self, other):
        if len(self) > len(other):
            larger = self
            smaller = other
        else:
            larger = other
            smaller = self
        newdict = dict()
        for k, v in list(smaller.items()):
            if larger.kv_contains(k, v):
                newdict[k] = v
        return self.__class__(**newdict)

    def contains(self, other):
        for k, v in list(other.items()):
            if not self.kv_contains(k, v):
                return False
        return True

    def without_keys(self, *keys):
        cp = dict(self._dict)
        for key in keys:
            del cp[key]
        return self.__class__(**cp)

    def purge_keys(self, *keys):
        cp = dict(self._dict)
        for key in keys:
            try:
                del cp[key]
            except KeyError:
                pass
        return self.__class__(**cp)

    def replace_keys(self, key_dict, *more_key_dicts):
        cp = dict(self._dict)
        for key, val in list(key_dict.items()):
            cp[key] = val
        if more_key_dicts:
            for key_dict in more_key_dicts:
                for key, val in list(key_dict.items()):
                    cp[key] = val
        return self.__class__(**cp)

    def subkey_has(self, other):
        try:
            for k, v in list(other.items()):
                v2 = self._dict[k]
                if v != v2:
                    return False
        except KeyError:
            return False
        return True

    def __sub__(self, other):
        cp = dict(self._dict)
        for k, v in list(other._dict.items()):
            assert cp[k] == v
            del cp[k]
        return self.__class__(cp)

    def __deepcopy__(self, memo):
        # by the immutibility of this object and anything it stores, it is OK to return the same thing
        return self

    def __copy__(self):
        # by the immutibility of this object and anything it stores, it is OK to return the same thing
        return self
