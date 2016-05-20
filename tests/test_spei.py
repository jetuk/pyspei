from __future__ import unicode_literals
from distutils import dir_util
import pytest
import os
import numpy as np


@pytest.fixture
def datadir(tmpdir, request):
    '''
    Fixture responsible for searching a folder with the same name of test
    module and, if available, moving all contents to a temporary directory so
    tests can use them freely.
    '''
    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)

    if os.path.isdir(test_dir):
        dir_util.copy_tree(test_dir, bytes(tmpdir))
    return tmpdir


def test_tampa(datadir):
    from pyspei.cli import main

    output_file = str(datadir.join('test_out.txt'))
    # Calculate SPEI using CLI method
    main(12, str(datadir.join('tampa.txt')), output_file)

    # Read correct values
    with open(str(datadir.join('tampa_spei_12.txt'))) as fh_test, open(output_file) as fh_calc:
        # Check name
        assert fh_test.readline() == fh_calc.readline()
        np.testing.assert_allclose(float(fh_test.readline()), float(fh_calc.readline()))
        assert fh_test.readline() == fh_calc.readline()
        assert fh_test.readline() == fh_calc.readline()
        # data computed from C version
        test_data = np.array([float(v) for v in fh_test.readlines()])
        spei_data = np.array([float(v) for v in fh_calc.readlines()])

        np.testing.assert_allclose(test_data, spei_data, rtol=1e-4, atol=1e-4)
