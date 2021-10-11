# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals, absolute_import
import yaml
import re


# HACK: fix loading number in scientific notation
#from PyGWINC
#
# https://stackoverflow.com/questions/30458977/yaml-loads-5e-6-as-string-and-not-a-number
#
# An apparent bug in python-yaml prevents it from regognizing
# scientific notation as a float.  The following is a modified version
# of the parser that recognize scientific notation appropriately.
yaml_loader = yaml.SafeLoader
yaml_loader.add_implicit_resolver(
    'tag:yaml.org,2002:float',
    re.compile('''^(?:
     [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
    |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
    |\\.[0-9_]+(?:[eE][-+][0-9]+)?
    |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
    |[-+]?\\.(?:inf|Inf|INF)
    |\\.(?:nan|NaN|NAN))$''', re.X),
    list('-+0123456789.'))


def yaml_load(fname):
    d = yaml.load(fname, Loader=yaml_loader)
    return d


