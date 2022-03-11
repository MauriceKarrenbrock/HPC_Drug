# -*- coding: utf-8 -*-
# pylint: disable=missing-docstring
# pylint: disable=redefined-outer-name
# pylint: disable=wildcard-import
# pylint: disable=unused-wildcard-import
# pylint: disable=no-self-use
#############################################################
# Copyright (c) 2021-2021 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# AGPL v3.0 License                                         #
#############################################################
"""pytest fixtures definitions and other general test configurations
"""
try:  #python>=3.7
    import importlib.resources as importlib_resources
except ImportError:  # python3.6
    import importlib_resources

import pytest


@pytest.fixture(scope='session')
def get_data_dir():
    """Will return the path of a file in integration_tests/data
    """
    with importlib_resources.path('tests.files4tests',
                                  'dummy.txt') as f:
        p = f.resolve()
    return p.parent
