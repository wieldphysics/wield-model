"""
"""
from __future__ import print_function, division, unicode_literals, absolute_import
from wavestate.utilities.np.semidense import KeyGroup


kg_polarization = KeyGroup(
    'polarization',
    ('S', 'P'),
)

kg_twophoton_qp = KeyGroup(
    'quantum',
    ('q', 'p'),
    identifier = 'qp',
)

kg_twophoton_pm = KeyGroup(
    'quantum',
    ('+', '-'),
    identifier = 'pm',
)

kg_twophoton_qpAC = KeyGroup(
    'quantum',
    ('qAC', 'pAC'),
    identifier = 'qpAC',
)

kg_twophoton_pmAC = KeyGroup(
    'quantum',
    ('+AC', '-AC'),
    identifier = 'pmAC',
)

