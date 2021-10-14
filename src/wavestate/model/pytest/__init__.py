# -*- coding: utf-8 -*-
"""
"""
import pytest
from wavestate.model.system import algo_log as ws_algo_log

@pytest.fixture
def algo_log(tpath):
    return ws_algo_log.LoggingAlgorithm(
        log_level = 9,
        log_folder = tpath,
    )
