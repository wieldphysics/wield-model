# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals, absolute_import


class LinkageAlgorithm(object):
    def __init__(self, bg_algo):
        self.bg_algo = bg_algo
        self.pbg = self.bg_algo.pbg
        return

class LinkageAlgorithmView(BondGraphAlgorithmView):
    def __init__(self, **kw):
        super(LinkageAlgorithmView, self).__init__(**kw)

class LinkageAlgorithmManipulator(object):

    def link_transform # internal

    #is a dictionary of known basis elements as of this linkage
    def link_basis('A!i')

    def links_used()
    #returns the links which must be supplied, will be a mixture of types. May
    #make some way to query this if that becomes useful

    #stores cachable values for future algorithms
    def cache[]
