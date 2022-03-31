import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import unittest
import pyee
import numpy as np


class PyeeBasicTestSuite(unittest.TestCase):
    """Basic Pyee test cases."""

    def setUp(self):
        self.data = {'general': {'simdir': '.'},
                     'simulation': {'modes': [0]},
                     'geometry': {'axi': 0, 'Lz': 2.0, 'Lx': 2.0, 'nz': 2, 'nx': 2},
                     'conditions': {'jz': ['example', 1], 'jx': ['example', 1], 'jy': ['example', 1]}}

    def test_prepareData_default(self):
        """ Default data has the basic fields """
        data = pyee.pre.prepareData()
        keys = data.keys()
        self.assertIn('general', keys)
        self.assertIn('simulation', keys)
        self.assertIn('geometry', keys)

    def test_constructProblem(self):
        """ Matrix and vector construction """
        A, b = pyee.core.constructProblem(self.data, 0)

    def test_getField(self):
        """ get Field """
        A, b = pyee.core.constructProblem(self.data, 0)
        solution = np.linalg.solve(A, b)
        Z, R, Ez = pyee.core.getField(self.data, solution, 'Ez')
        print(Ez)


if __name__ == '__main__':
    unittest.main()
