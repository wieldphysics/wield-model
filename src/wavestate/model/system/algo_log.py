#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import sys
import time
import logging
import re

import contextlib
from wavestate import declarative

from wavestate.utilities.strings import padding_remove

class LoggingAlgorithm(object):
    def __init__(
        self,
        hints = None,
        filters = dict(),
        **kwargs
    ):
        self.mtime_start = time.time()

        if hints is None:
            hints = kwargs
        else:
            hints = dict(hints)
            hints.update(kwargs)

        if isinstance(hints, (list, tuple)):
            usehints = dict()
            for hdict in hints:
                usehints.update(hdict)
        else:
            usehints = dict(hints)

        self.hints = usehints
        self.filter_hints = dict()

        compiled_filters = dict()
        for re_tup, f_hints in filters.items():
            if isinstance(re_tup, (tuple, list)):
                (header_re, ref_re) = re_tup
                header_re = re.compile(header_re)
                ref_re = re.compile(ref_re)
            else:
                header_re = re.compile(re_tup)
                ref_re = None
            compiled_filters[header_re, ref_re] = f_hints

        self.filters = compiled_filters

        #holds a heading for the logging, as well as sets tabbing
        self.log_header_stack = ""
        self.log_header_short_stack = ()
        self.log_header_stack_len = 0
        #indicates how far into the header has been printed yet.
        #for the live log
        self.log_header_printed = 0
        self.log_ref = ''

        #log entries
        self._logs = []
        self.log_number = 0

        #investigations to view
        self._investigations = dict()
        return

    def hint_has(self, hname):
        return hname in self.hints

    def hint_setdefault(self, hname, hval):
        self.hints.setdefault(hname, hval)
        return

    def hint(self, *args, **kwargs):
        superarg = []
        for arg in args:
            if isinstance(arg, (list, tuple)):
                superarg.extend(arg)
            else:
                superarg.append(arg)

        for key in superarg:
            key = key.format(**kwargs)

            ret = self.filter_hints.get(key, declarative.NOARG)
            if ret is not declarative.NOARG:
                return ret

            ret = self.hints.get(key, declarative.NOARG)
            if ret is not declarative.NOARG:
                return ret
        return kwargs['default']

    def log(
        self,
        *args,
        **kwargs
    ):
        """
        First argument is the level, should include a log group, which must be one of
        ['info', 'debug', 'warn', 'progress', 'investigate']
        """
        level = args[0]
        if isinstance(level, int):
            level = args[0]
            args = args[1:]
            group = kwargs.setdefault('group', 'info')
        else:
            level = -1
            group = kwargs.setdefault('group', 'debug')
            #TODO print line and file upon hint request
            #args = args

        header = self.log_header_stack
        kwargs['header'] = header
        kwargs['ref'] = self.log_ref
        kwargs['time'] = time.time()
        kwargs['time_start'] = self.mtime_start

        if self.hint('log_off', default = False):
            return

        kwargs['args'] = args
        #TODO, merge if consecutive with the same parameters
        if args:
            save_logs = self.hint(['save_logs'], default = False)
            if save_logs:
                self._logs.append(
                    kwargs
                )

        self.log_number += 1

        #FOR LIVE PRINTING
        if group == 'info':
            log_mod_level = logging.INFO
            group_character = 'I'
            level_limit = self.hint([
                'log_level_info',
                'log_level',
            ], default = 0)
        elif group == 'debug':
            log_mod_level = logging.DEBUG
            group_character = 'D'
            level_limit = self.hint([
                'log_level_debug',
                'log_level',
            ], default = 0)
        elif group == 'warn':
            log_mod_level = logging.WARNING
            group_character = 'W'
            level_limit = self.hint([
                'log_level_warn',
                'log_level',
            ], default = 0)
        elif group == 'investigate':
            log_mod_level = logging.INFO
            group_character = 'I'
            level_limit = self.hint([
                'log_level_investigate',
                'log_level',
            ], default = 0)
        elif group == 'progress':
            log_mod_level = logging.INFO
            group_character = 'P'
            level_limit = self.hint([
                'log_level_progress',
                'log_level',
            ], default = 0)
        else:
            raise RuntimeError("Unrecognized log grouping")

        if self.hint('log_print', default = True) and level <= level_limit:
            hint_log_stdout = self.hint('log_stdout', default = True)
            if hint_log_stdout not in [None, True, False]:
                lfile = hint_log_stdout
            else:
                lfile = sys.stdout

            header = self.log_header_stack
            header_len = self.log_header_stack_len

            prefix = "{}{} {: >6.2f} {}".format(
                level if level >= 0 else '-',
                group_character,
                kwargs['time'] - kwargs['time_start'],
                '  ' * header_len
            )

            #TODO, make these take a header argument
            if not self.hint('logging_use', default = False):
                def pfunc(*args, **kwargs):
                    print(*args, **kwargs)
            else:
                def pfunc(*args, **kwargs):
                    kwargs.pop('file', None)
                    logging.log(log_mod_level + 9 - level, *args, **kwargs)

            if header_len > self.log_header_printed:
                pfunc(
                    "{}{}".format(
                        '-' * (len(prefix)),
                        header
                    ),
                    file = lfile
                )
                self.log_header_printed = header_len
                #tag that the header has been printed

            hint_log_stderr = self.hint('log_stderr', default = True)
            if hint_log_stderr and group == 'warn':
                if hint_log_stderr not in [None, True, False]:
                    lfile = hint_log_stderr
                else:
                    lfile = sys.stderr
            else:
                lfile = sys.stdout

            if not args:
                return

            arg_lines = [[]]
            for arg in args:
                if isinstance(arg, str):
                    if '\n' in arg:
                        arg = padding_remove(arg)
                    arg_spl = arg.split('\n')
                    arg_lines[-1].append(arg_spl[0])
                    for subline in arg_spl[1:]:
                        arg_lines.append([subline])
                else:
                    arg_lines[-1].append(arg)

            #TODO, have pfunc do this splitting
            pfunc(
                prefix, *arg_lines[0],
                file = lfile
            )
            for argsl in arg_lines[1:]:
                pfunc(
                    ' ' * len(prefix), *argsl,
                    file = lfile
                )
        return

    def investigate(self, name, desc, func, offline = False):
        with self.heading(name):
            kwargs = dict()
            kwargs['group'] = 'investigate'
            if offline:
                offline_notify = " (available offline)"
            else:
                offline_notify = ""
            self.log(5, desc + offline_notify, **kwargs)
            if offline:
                self.investigations[self.log_header_stack] = func
            if self.hint('investigate', default = False):
                self.info(
                    5,
                    "investigating: {}".format(desc),
                    **kwargs
                )
                func(self)

    def debug(self, *args, **kwargs):
        kwargs['group'] = 'debug'
        self.log(*args, **kwargs)

    def warn(self, *args, **kwargs):
        kwargs['group'] = 'warn'
        self.log(*args, **kwargs)

    def info(self, *args, **kwargs):
        kwargs['group'] = 'info'
        self.log(*args, **kwargs)

    def progress(self, *args, **kwargs):
        kwargs['group'] = 'progress'
        self.log(*args, **kwargs)

    @contextlib.contextmanager
    def heading(self, header, header_short = None):
        save_stack = self.log_header_stack
        self.log_header_stack = save_stack + header + ":"
        save_stack_short = self.log_header_short_stack
        if header_short is None:
            header_short = header
        self.log_header_short_stack = save_stack_short + (header_short,)
        self.log_header_stack_len += 1

        #TODO, could defer this until logging actually occurs
        filter_hints = dict(self.filter_hints)
        filters_save = self.filters
        filters_new = None
        for (header_re, ref_re), hints in self.filters.items():
            if header_re.search(self.log_header_stack):
                if ref_re is None or ref_re.search(self.log_ref):
                    filter_hints.update(hints)
                    if filters_new is None:
                        filters_new = dict(self.filters)
                    filters_new.pop((header_re, ref_re))
        if filters_new is not None:
            self.filters = filters_new
        filter_hints_save = self.filter_hints
        self.filter_hints = filter_hints

        if self.hint('log_heading', default = False):
            #print a null info, which will print the heading, but nothing under it
            self.info()
        #TODO, auto print header on command?
        yield
        self.filters = filters_save
        self.filter_hints = filter_hints_save
        self.log_header_short_stack = save_stack_short
        self.log_header_stack = save_stack
        self.log_header_stack_len -= 1
        if self.log_header_printed > self.log_header_stack_len:
            self.log_header_printed = self.log_header_stack_len

    @contextlib.contextmanager
    def reference(self, ref):
        save_ref = self.log_ref
        self.log_ref = ref

        #TODO, could defer this until logging actually occurs
        filter_hints = dict(self.filter_hints)
        for (header_re, ref_re), hints in self.filters.items():
            if header_re.match(self.log_header_stack):
                if ref_re is None or ref_re.match(self.log_ref):
                    filter_hints.update(hints)
        filter_hints_save = self.filter_hints
        self.filter_hints = filter_hints

        #TODO, could defer this until logging actually occurs
        filter_hints = dict(self.filter_hints)
        filters_save = self.filters
        filters_new = None
        for (header_re, ref_re), hints in self.filters.items():
            if header_re.search(self.log_header_stack):
                if ref_re is None or ref_re.search(self.log_ref):
                    filter_hints.update(hints)
                    if filters_new is None:
                        filters_new = dict(self.filters)
                    filters_new.pop((header_re, ref_re))
        if filters_new is not None:
            self.filters = filters_new
        filter_hints_save = self.filter_hints
        self.filter_hints = filter_hints

        #TODO, auto print ref on command?
        yield
        self.filters = filters_save
        self.filter_hints = filter_hints_save
        self.log_ref = save_ref

