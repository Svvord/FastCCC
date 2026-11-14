import os

import pytest


def pytest_addoption(parser):
    parser.addoption('--integration', action='store_true', default=False,
                     help='run integration tests (download public data, run full FastCCC workflows)')


def pytest_collection_modifyitems(config, items):
    enabled = config.getoption('--integration') or os.environ.get('FASTCCC_INTEGRATION') == '1'
    if enabled:
        return
    skip = pytest.mark.skip(reason='integration test: pass --integration or set FASTCCC_INTEGRATION=1')
    for item in items:
        if 'integration' in item.keywords:
            item.add_marker(skip)
