# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals

from ..pgraph import (
    ParameterObject,
    #port_types,
)


class Frequency(ParameterObject):
    def __init__(self):
        super(Frequency, self).__init__()

        with self._internal():
            self['frequency[Hz]'] = None
            self['order_mechanical'] = 0
            self['order_optical'] = 2

            @self.deco_one_one('order_optical')
            def order_signal(order_optical):
                return 2 * order_optical

            @self.deco_one_one('order_signal')
            def order_electrical(order_signal):
                return order_signal

            #set the default span to be the value
            self.set_assign(
                kto = 'frequency_span[Hz]',
                kfrom = 'frequency[Hz]',
            )

        self.values_settable([
            'frequency[Hz]',
            'frequency_span[Hz]',
            'order_optical',
            'order_signal',
            'order_electrical',
            'order_mechanical',
        ])


class FrequencySpan(ParameterObject):
    def __init__(self):
        super(FrequencySpan, self).__init__()

        #should be a list or iterable of the frequency names
        self['frequencies']   = None
        self['mechanical'] = True
        self['electrical'] = True
        self['optical']    = True
        self['signal']     = True

        @self.deco_one_one('frequencies')
        def order(frequencies):
            return len(frequencies) + 1


class FrequencySuppress(ParameterObject):
    def __init__(self):
        super(FrequencySuppress, self).__init__()

        #should be a dictionary of the orders, or lists of the orders
        #if lists, all of the lengths must match
        self['suppress']   = None

        #indicates the the values should be plus-minus'ed
        self['balanced']   = True

        self['mechanical'] = True
        self['optical']    = True
        self['signal']     = True

class FrequencyAliases(ParameterObject):
    pass

