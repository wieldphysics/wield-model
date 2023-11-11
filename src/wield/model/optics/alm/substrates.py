#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

from ...base import FrequencyKey


substrates = dict(
    fused_silica={
        FrequencyKey({"Nd1064": 1}): 1.4496,
        FrequencyKey({"Nd1064": 2}): 1.4607,
        FrequencyKey({"1550": 1}): 1.440,
        FrequencyKey({"1550": 2}): 1.4538,
        # https://refractiveindex.info/?shelf=glass&book=fused_silica&page=Malitson
        # n^2-1=\frac{0.6961663λ^2}{λ^2-0.0684043^2}+\frac{0.4079426λ^2}{λ^2-0.1162414^2}+\frac{0.8974794λ^2}{λ^2-9.896161^2}
        "Sellmeier": None,
    },
    silicon={
        FrequencyKey({"1550": 1}): 3.4850,
        FrequencyKey({"1550": 2}): 3.6950,
    },
    BK7={
        FrequencyKey({"Nd1064": 1}): 1.5066,
        FrequencyKey({"Nd1064": 2}): 1.5195,
    },
    vacuum={
        FrequencyKey({"Nd1064": 1}): 1.0,
        FrequencyKey({"Nd1064": 2}): 1.0,
        FrequencyKey({"1550": 1}): 1.0,
        FrequencyKey({"1550": 2}): 1.0,
    },
    nitrogen={
        FrequencyKey({"Nd1064": 1}): 1.0002952,
        FrequencyKey({"Nd1064": 2}): 1.0002994,
    },
    PPKTP={
        FrequencyKey({"Nd1064": 1}): 1.8302,
        FrequencyKey({"Nd1064": 2}): 1.7779,
        FrequencyKey({"1550": 1}): 1.7282,
        FrequencyKey({"1550": 2}): 1.7497,
    },
)
