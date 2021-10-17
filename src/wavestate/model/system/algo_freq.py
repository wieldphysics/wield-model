# -*- coding: utf-8 -*-
"""
"""

import itertools
import numpy as np
try:
    from collections import abc
except ImportError:
    import collections as abc

from .. import base
from .. import _base
from .. import optics


class FrequenciesAlgorithm(object):
    def __init__(
            self,
            pbg,
            log = _base.log_default,
    ):
        self.pbg = pbg
        self.log = log
        self.frequencies_root = self.pbg.referred_path(self.pbg.root, 'frequencies/')

        self._assemble_frequencies()
        self._map_frequencies()

        #TODO, allow the AC key to be missing, but disable AC analysis if it is
        #AC key is removed
        #self.AC_root  = self.pbg.referred_path(self.pbg.root, 'AC/')
        #self.AC_fkey  = base.FrequencyKey({'AC': 1})
        #self.AC_fspan = self.pbg.resolve_parameter(self.AC_root, 'frequency_span[Hz]')
        return

    def _assemble_frequencies(self):
        frequencies = dict()
        freq_basis = dict()

        wavelengths = dict()
        wlen_basis = dict()
        spans = dict()
        suppressors = dict()
        wavelength_aliases_obj = dict()
        frequency_aliases_obj = dict()
        for fname, F_obj in self.pbg.object_references[self.frequencies_root].items():
            #print(fname, F_obj)
            if isinstance(F_obj, base.Frequency):
                frequencies[fname] = F_obj
                freq_basis[fname] = base.FrequencyKey({fname: 1})
            elif isinstance(F_obj, base.FrequencySuppress):
                suppressors[fname] = F_obj
            elif isinstance(F_obj, base.FrequencySpan):
                spans[fname] = F_obj
            elif isinstance(F_obj, optics.OpticalFrequency):
                wavelengths[fname] = F_obj
                wlen_basis[fname] = base.FrequencyKey({fname: 1})
            elif isinstance(F_obj, optics.OpticalFrequencyAliases):
                wavelength_aliases_obj[fname] = F_obj
            elif isinstance(F_obj, base.FrequencyAliases):
                frequency_aliases_obj[fname] = F_obj

        self.frequencies = frequencies
        self.freq_basis = freq_basis
        self.wavelengths = wavelengths
        self.wlen_basis = wlen_basis

        frequency_values = dict()
        frequency_spans  = dict()
        for fname, F_obj in self.frequencies.items():
            fval = self.pbg.resolve_parameter(F_obj, 'frequency[Hz]')
            fspan = self.pbg.resolve_parameter(F_obj, 'frequency_span[Hz]')
            if fval is None:
                raise RuntimeError("Frequency must be set")
            if fspan is None:
                fspan = fval

            frequency_values[fname] = fval
            frequency_spans[fname] = fspan
        self.frequency_values = frequency_values
        self.frequency_spans = frequency_spans

        wavenumber_values = dict()
        for wname, W_obj in self.wavelengths.items():
            wnval = self.pbg.resolve_parameter(W_obj, 'wavenumber[1_m]')
            wavenumber_values[wname] = wnval
            if wnval is None:
                raise RuntimeError("wavelength must be set")
        self.wavenumber_values = wavenumber_values

        sc_dict_optical = dict()
        sc_dict_signal = dict()
        sc_dict_electrical = dict()
        sc_dict_mechanical = dict()
        #TODO, allow wavelength spans
        sc_dict_wavenumbers = dict()
        for sname, span_obj in spans.items():
            freqs = self.pbg.resolve_parameter(span_obj, 'frequencies')
            fnames = frozenset(freqs)
            order = self.pbg.resolve_parameter(span_obj, 'order')

            if self.pbg.resolve_parameter(span_obj, 'optical'):
                assert(fnames not in sc_dict_optical)
                sc_dict_optical[fnames] = order

            if self.pbg.resolve_parameter(span_obj, 'signal'):
                assert(fnames not in sc_dict_signal)
                sc_dict_signal[fnames] = order

            if self.pbg.resolve_parameter(span_obj, 'electrical'):
                assert(fnames not in sc_dict_electrical)
                sc_dict_electrical[fnames] = order

            if self.pbg.resolve_parameter(span_obj, 'mechanical'):
                assert(fnames not in sc_dict_mechanical)
                sc_dict_mechanical[fnames] = order

        fkey_suppress_optical = set()
        fkey_suppress_signal = set()
        fkey_suppress_mechanical = set()
        #TODO, allow wavelength spans
        fkey_suppress_wavelength = set()

        for sname, suppress_obj in suppressors.items():
            balanced = self.pbg.resolve_parameter(suppress_obj, 'balanced')
            suppress = self.pbg.resolve_parameter(suppress_obj, 'suppress')
            suppress = dict(suppress)

            test_v = next(iter(suppress.values()))

            #TODO improve error handling
            def add_fkey(fkey):
                if self.pbg.resolve_parameter(suppress_obj, 'optical'):
                    fkey_suppress_optical.add(fkey)
                    if balanced:
                        fkey_suppress_optical.add(-fkey)

                if self.pbg.resolve_parameter(suppress_obj, 'signal'):
                    fkey_suppress_signal.add(fkey)

                if self.pbg.resolve_parameter(suppress_obj, 'mechanical'):
                    fkey_suppress_mechanical.add(fkey)

            if isinstance(test_v, (list, tuple)):
                klist = list(suppress.keys())
                vlistlist = [suppress[k] for k in klist]
                for vlist in zip(*vlistlist):
                    fkey = base.FrequencyKey(dict(zip(klist, vlist)))
                    add_fkey(fkey)
            else:
                fkey = base.FrequencyKey(suppress)
                add_fkey(fkey)

        #set some default spans
        for fname, F_obj in frequencies.items():
            sckey = frozenset([fname])
            if sckey not in sc_dict_optical:
                sc_dict_optical.setdefault(sckey, self.pbg.resolve_parameter(F_obj, 'order_optical'))

        self.freq_set_optical = self._complete_frequencies(
            self.freq_basis,
            sc_dict = sc_dict_optical,
            suppress = fkey_suppress_optical,
            balanced = True,
        )

        #set some default spans
        for fname, F_obj in frequencies.items():
            sckey = frozenset([fname])
            if sckey not in sc_dict_signal:
                sc_dict_signal.setdefault(sckey, self.pbg.resolve_parameter(F_obj, 'order_signal'))

        self.freq_set_signal = self._complete_frequencies(
            self.freq_basis,
            sc_dict = sc_dict_signal,
            suppress = fkey_suppress_signal,
            balanced = False,
        )

        #set some default spans
        for fname, F_obj in frequencies.items():
            sckey = frozenset([fname])
            if sckey not in sc_dict_electrical:
                sc_dict_electrical.setdefault(sckey, self.pbg.resolve_parameter(F_obj, 'order_electrical'))

        self.freq_set_electrical = self._complete_frequencies(
            self.freq_basis,
            sc_dict = sc_dict_electrical,
            suppress = fkey_suppress_signal,
            balanced = False,
        )

        #set some default spans
        for fname, F_obj in frequencies.items():
            sckey = frozenset([fname])
            if sckey not in sc_dict_mechanical:
                sc_dict_mechanical.setdefault(sckey, self.pbg.resolve_parameter(F_obj, 'order_mechanical'))

        self.freq_set_mechanical = self._complete_frequencies(
            self.freq_basis,
            sc_dict = sc_dict_mechanical,
            suppress = fkey_suppress_mechanical,
            balanced = False,
        )

        #set some default spans
        for fname, F_obj in wavelengths.items():
            sckey = frozenset([fname])
            if sckey not in sc_dict_wavenumbers:
                sc_dict_wavenumbers.setdefault(sckey, self.pbg.resolve_parameter(F_obj, 'order'))
        #print('wavenum: ', sc_dict_wavenumbers)

        freq_set_wavenumbers = self._complete_frequencies(
            self.wlen_basis,
            sc_dict_wavenumbers,
            default_restrict = 0,
            balanced = False,
        )
        freq_set_wavenumbers.remove(base.FrequencyKey({}))
        self.freq_set_wavenumbers = freq_set_wavenumbers

        #TODO, check for double mappings and warn appropriately
        wavelength_aliases = {}
        for aname, wAO in wavelength_aliases_obj.items():
            wkeydict = self.pbg.resolve_parameter(wAO, 'to')
            wk = base.FrequencyKey(wkeydict)
            for name in self.pbg.resolve_parameter(wAO, 'names'):
                wavelength_aliases[name] = wk
        self.wavelength_aliases = wavelength_aliases

        frequency_aliases = {}
        for aname, fAO in frequency_aliases_obj.items():
            fkeydict = self.pbg.resolve_parameter(fAO, 'to')
            fk = base.FrequencyKey(fkeydict)
            for name in self.pbg.resolve_parameter(fAO, 'names'):
                frequency_aliases[name] = fk
        self.frequency_aliases = frequency_aliases
        return

    def default_optical_wavekey(self):
        #TODO, make a better defaulting system
        assert(len(self._pa_algo.fs.freq_set_wavenumbers) == 1)
        return next(iter(self._pa_algo.fs.freq_set_wavenumbers))

    def parameter_to_fk(self, fparam):
        #check for mappings first since they are unhashable and fail for the alias check
        if isinstance(fparam, abc.Mapping):
            for fname, m in fparam.items():
                assert(fname in self.frequencies)
            return base.FrequencyKey(fparam)

        #first check if there is a direct alias name
        mapping = self.frequency_aliases.get(fparam, None)
        if mapping is not None:
            return mapping

        if isinstance(fparam, str):
            assert(fparam in self.frequencies)
            return base.FrequencyKey({fparam : 1})
        raise RuntimeError("Frequency parameter {} not recognized".format(fparam))

    def parameter_to_wk(self, wparam):
        if isinstance(wparam, base.FrequencyKey):
            assert(wparam in self.freq_set_wavenumbers)
            return wparam
        elif isinstance(wparam, abc.Mapping):
            for fname, m in wparam.items():
                assert(fname in self.frequencies)
            return base.FrequencyKey(wparam)

        #first check if there is a direct alias name
        mapping = self.wavelength_aliases.get(wparam, None)
        if mapping is not None:
            return mapping

        if isinstance(wparam, str):
            assert(wparam in self.frequencies)
            return base.FrequencyKey({wparam : 1})
        raise RuntimeError("Wavelength parameter {} not recognized".format(wparam))

    def _complete_frequencies(
            self,
            freq_basis,
            sc_dict,
            suppress = set(),
            default_restrict = None,
            balanced = False,
    ):
        fnames_all = list(freq_basis.keys())
        F_keys = set([base.FrequencyKey({})])

        for N_combinations in range(1, 1+len(fnames_all)):
            for fnames in itertools.combinations(fnames_all, N_combinations):
                fnames = frozenset(fnames)
                restrict = sc_dict.get(fnames, default_restrict)

                def add_one_each(Fkey, skip_check = None):
                    if Fkey in F_keys:
                        return

                    if Fkey in suppress:
                        return

                    for fname in fnames:
                        if fname is skip_check:
                            continue
                        fd = dict(Fkey.F_dict)
                        del fd[fname]
                        if base.FrequencyKey(fd) not in F_keys:
                            #all subkeys must be allowed
                            return

                    if restrict is not None and Fkey.ord() > restrict:
                        return

                    #print("Adding:", Fkey)
                    F_keys.add(Fkey)

                    for fname in fnames:
                        Ford = Fkey.F_dict.get(fname, 0)
                        if Ford > 0:
                            add_one_each(Fkey + freq_basis[fname], skip_check = fname)
                        elif Ford < 0:
                            add_one_each(Fkey - freq_basis[fname], skip_check = fname)
                        else:
                            raise RuntimeError("Cannot Happen (BUG!)")

                if balanced:
                    for Fvect in itertools.product(*[(freq_basis[fname], -freq_basis[fname]) for fname in fnames]):
                        Fsum = np.sum(Fvect)
                        add_one_each(Fsum)
                else:
                    Fsum = np.sum([freq_basis[fname] for fname in fnames])
                    add_one_each(Fsum)
                #print("BUILDING: ", fnames, restrict)
                ##add_one_each(base.FrequencyKey({}))
                ##TODO, recursive algorithm to recurse down all paths, check if
                ##lower groups as disallowed
        return F_keys

    def print_frequencies(self):
        print("OPTICAL: ", self.freq_set_optical)
        print("SIGNAL: ", self.freq_set_signal)
        print("MECH: ", self.freq_set_mechanical)
        print("WLEN: ", self.freq_set_wavenumbers)

    def _map_frequencies(self):
        """
        Generates the mapping from FrequencyKeys to the numeric frequencies
        returned
        """
        freq_set_all = set()
        freq_set_all.update(self.freq_set_optical)
        freq_set_all.update(self.freq_set_signal)
        freq_set_all.update(self.freq_set_mechanical)

        freq_map = {}
        span_map = {}
        for fd in freq_set_all:
            fval = 0
            fspan = 0
            for fname, ford in fd.F_dict.items():
                fval = fval + self.frequency_values[fname] * ford
                fspan = fspan + self.frequency_spans[fname] * ford
            freq_map[fd] = fval
            span_map[fd] = fspan
        self.frequency_map = freq_map
        self.frequency_span_map = span_map

        wnum_map = {}
        wlen_map = {}
        for wd in self.freq_set_wavenumbers:
            wnval = 0
            for fname, ford in wd.F_dict.items():
                wnval = wnval + self.wavenumber_values[fname] * ford
            wnum_map[wd] = wnval
            wlen_map[wd] = (2 * np.pi)/wnval
        self.wavenumber_map = wnum_map
        self.wavelength_map = wlen_map



        
