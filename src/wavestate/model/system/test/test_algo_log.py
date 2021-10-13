# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals, absolute_import
from os import path

from wavestate.model.system import algo_log


def test_log():
    def investigate_f(log):
        print("HMMM")

    print()
    log = algo_log.LoggingAlgorithm(
        hints = dict(
            log_level = 8,
            log_folder = path.split(__file__)[0]
        ),
        filters = {
            r'tester' : dict(investigate = False),
            r'warn2' : dict(investigate = True),
        }
    )
    log.warn("test")
    log.investigate(
        'warn2',
        'checks if investigations are working',
        func = investigate_f
    )
    log.warn(5, "test")
    with log.heading('warn2'):
        log.warn(5, "test")
        log.warn(5, "test")
        with log.reference('test/'):
            log.warn(5, "test")
        log.investigate(
            'check',
            'checks if investigations are working',
            func = investigate_f
        )
        log.investigate(
            'tester',
            'checks if investigations are working',
            func = investigate_f
        )

