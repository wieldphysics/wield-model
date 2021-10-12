# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals, absolute_import
import declarative

from transient.numrep import constants

from .. import base
from .. import _base

from . import (
    algo_freq,
    algo_bg,
)


class PhysicsAlgorithm(object):
    check_build = True

    def __init__(
        self,
        sys = None,
        pbg = None,
        log = _base.log_default,
        previous = None,
    ):
        if pbg is None:
            pbg = base.ParameterGraph(sys)
        self.log = log
        self.pbg = pbg
        self.previous = previous
        return

    def regenerate(self, pbg):
        return self.__class__(
            pbg = pbg,
            log = self.log,
            previous = self,
        )

    _fs = None
    @property
    def fs(self):
        if self._fs is not None:
            return self._fs
        else:
            self._fs = algo_freq.FrequenciesAlgorithm(self.pbg, log = self.log)
        return self._fs

    _bg = None
    @property
    def bg(self):
        if self._bg is not None:
            return self._bg
        else:
            self._bg = algo_bg.BondGraphAlgorithm(self.pbg, log = self.log)
        return self._bg

    _mm = None
    @property
    def mm(self):
        if self._mm is not None:
            return self._mm
        else:
            from . import algo_alm
            if self.previous is not None:
                previous = self.previous.mm
            else:
                previous = None
            self._mm = algo_alm.ModeMatchingAlgorithm(self, previous = previous)
        return self._mm

    _dc = None
    @property
    def dc(self):
        if self._dc is not None:
            return self._dc
        else:
            self._setup_optical()
            self._setup_mechanical()
            self._setup_electrical()
            self._setup_signal()

            from . import algo_DC
            self._dc = algo_DC.PhysicsDCAlgorithm(self)
        return self._dc

    _ac = None
    @property
    def ac(self):
        if self._ac is not None:
            return self._ac
        else:
            from . import algo_AC
            self._ac = algo_AC.PhysicsACAlgorithm(self, dc = self.dc)
        return self._dc

    _ga = None
    @property
    def ga(self):
        if self._ga is not None:
            return self._ga
        from . import algo_graphs
        self._ga = algo_graphs.GraphPlottingAlgorithm(self)
        return self._ga

    def _setup_optical(self):
        optical_settings = wavestate.bunch.Bunch(
            polarization    = True,
            quantum         = 'pm',
            transverse      = None,
            transverse_mech = 'optickle',
        )

        kg_frequency_DC_optical = base.KeyGroup(
            'frequency',
            sorted(self.fs.freq_set_optical,
                   key = lambda fd : self.fs.frequency_map[fd])
        )
        kg_wavenumber = base.KeyGroup(
            'wavenumber',
            sorted(self.fs.freq_set_wavenumbers,
                   key = lambda wd : -self.fs.wavenumber_map[wd])
        )
        optical_basis_DC = dict(
            quantum    = base.kg_twophoton_pm,
            frequency  = kg_frequency_DC_optical,
            wavenumber = kg_wavenumber,
        )

        if optical_settings.polarization:
            optical_basis_DC['polarization'] = base.kg_polarization

        if optical_settings.transverse is None:
            pass
        elif optical_settings.transverse == 'HG':
            optical_basis_DC['HGX'] = base.KeyGroup('HGX', [0])
            optical_basis_DC['HGY'] = base.KeyGroup('HGY', [0])
        elif optical_settings.transverse == 'LG':
            optical_basis_DC['LG'] = base.KeyGroup('LG', [(0, 0)])

        self.optical_settings = optical_settings
        self.optical_basis_DC = optical_basis_DC

        fvals = []
        for fd in kg_frequency_DC_optical.enumerated:
            fvals.append(self.fs.frequency_span_map[fd])
        self._basis_frequency_values_DC_optical = fvals

        wvals = []
        for wd in kg_wavenumber.enumerated:
            wvals.append(self.fs.wavenumber_map[wd])
        self.basis_wavenumber_values = wvals

        optical_basis_AC = dict(optical_basis_DC)
        optical_basis_AC['quantum'] = base.kg_twophoton_pmAC
        self.optical_basis_AC = optical_basis_AC
        return

    def _setup_mechanical(self):
        #HMM
        self.kg_frequency_DC_mechanical = base.KeyGroup('frequency', sorted(
            self.fs.freq_set_mechanical,
            key = lambda fd : self.fs.frequency_map[fd]
        ))
        self.mechanical_basis_DC = dict(
            quantum    = base.kg_twophoton_pm,
            frequency  = self.kg_frequency_DC_mechanical,
        )
        fvals_DC = []
        for fd in self.kg_frequency_DC_mechanical.enumerated:
            fvals_DC.append(self.fs.frequency_span_map[fd])
        self._basis_frequency_values_DC_mech = fvals_DC
        mechanical_basis_AC = dict(self.mechanical_basis_DC)
        mechanical_basis_AC['quantum'] = base.kg_twophoton_pmAC
        self.mechanical_basis_AC = mechanical_basis_AC

    def _setup_electrical(self):
        #HMM
        self.kg_frequency_DC_electrical = base.KeyGroup('frequency', sorted(
            self.fs.freq_set_electrical,
            key = lambda fd : self.fs.frequency_map[fd]
        ))
        self.electrical_basis_DC = dict(
            quantum    = base.kg_twophoton_pm,
            frequency  = self.kg_frequency_DC_electrical,
        )
        fvals_DC = []
        for fd in self.kg_frequency_DC_electrical.enumerated:
            fvals_DC.append(self.fs.frequency_span_map[fd])
        self._basis_frequency_values_DC_electrical = fvals_DC

        electrical_basis_AC = dict(self.electrical_basis_DC)
        electrical_basis_AC['quantum'] = base.kg_twophoton_pmAC
        self.electrical_basis_AC = electrical_basis_AC

    def _setup_signal(self):
        #HMM
        self.kg_frequency_DC_signal = base.KeyGroup('frequency', sorted(
            self.fs.freq_set_signal,
            key = lambda fd : self.fs.frequency_map[fd]
        ))
        self.signal_basis_DC = dict(
            quantum    = base.kg_twophoton_pm,
            frequency  = self.kg_frequency_DC_signal,
        )
        fvals_DC = []
        for fd in self.kg_frequency_DC_signal.enumerated:
            fvals_DC.append(self.fs.frequency_span_map[fd])
        self._basis_frequency_values_DC_signal = fvals_DC

        signal_basis_AC = dict(self.signal_basis_DC)
        signal_basis_AC['quantum'] = base.kg_twophoton_pmAC
        self.signal_basis_AC = signal_basis_AC


class PhysicsAlgorithmView(algo_bg.LinkageAlgorithmView):
    _pa_algo = None

    constants = constants.constants_floats

    def __init__(self, pa_algo, **kw):
        super(PhysicsAlgorithmView, self).__init__(**kw)
        self._pa_algo = pa_algo
        self.check_build = pa_algo.check_build

    def optical_settings(self):
        return self._pa_algo.optical_settings

    def default_optical_freqkey(self):
        return base.FrequencyKey({})

    def default_optical_wavekey(self):
        return self._pa_algo.fs.default_optical_wavekey()

    @property
    def log(self):
        self._pa_algo.log

    def configure_optical_wavenumber(self, wdict):
        #TODO, check that the wavenumber exists
        #assert(len(self._pa_algo.fs.freq_set_wavenumbers) == 1)
        return base.FrequencyKey(wdict)

    def basis_wavenumbers(self, with_keys = False):
        #print("FVs: ", self._pa_algo._basis_frequency_values_DC_optical)
        if with_keys:
            return zip(
                self._pa_algo.optical_basis_DC['wavenumber'].enumerated,
                self._pa_algo.basis_wavenumber_values
            )
        else:
            return self._pa_algo.basis_wavenumber_values

    def parameter_to_fk(self, fparam):
        return self._pa_algo.fs.parameter_to_fk(fparam)

    def parameter_to_wk(self, wparam):
        return self._pa_algo.fs.parameter_to_wk(wparam)

    def optical_frequency_allowed(self, fk):
        return fk in self._pa_algo.fs.freq_set_optical
