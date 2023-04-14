#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""
import pytest
from wield.model.system import algo_log as ws_algo_log


@pytest.fixture
def algo_log(tpath):
    return ws_algo_log.LoggingAlgorithm(
        log_level=9,
        log_folder=tpath,
    )
