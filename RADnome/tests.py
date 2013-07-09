#!/usr/bin/env python
# encoding: utf-8

import os
import sys

sys.path.insert(0, os.path.abspath('..'))  # Seriously?! This is fucking ugly.

import glob
import unittest
from RADnome import GenerateRADnome

class TestPipeline(unittest.TestCase):

    def setUp(self):
        module_dir = os.path.dirname(GenerateRADnome.__file__)
        base_dir = os.path.split(module_dir)[0]
        self.test_data = os.path.join(base_dir, "tests/data")
        pass

    def test_rainbow(self):

        P = GenerateRADnome.RunPipeline()

        fq1 = os.path.join(self.test_data, "fq1.10k.fq")
        fq2 = os.path.join(self.test_data, "fq2.10k.fq")

        a = P.run_cluster_cmd(fq1)
        b = P.run_cluster_cmd(fq2)

        c = P.run_div_cmd(fq1)
        d = P.run_div_cmd(fq2)

        e = P.run_merge_cmd(fq1)
        f = P.run_merge_cmd(fq2)

        results = [a, b, c, d, e, f]
        self.assertEqual(results, [1, 1, 1, 1, 1, 1])

    def test_READnomes(self):

        P = GenerateRADnome.RunPipeline()

        fq1 = os.path.join(self.test_data, "fq1.10k.fq")
        fq2 = os.path.join(self.test_data, "fq2.10k.fq")

        P.make_READnome(fq1, "test.R1")
        P.make_READnome(fq2, "test.R2")

        # results = [a, b, c, d, e, f]
        # self.assertEqual(results, [1, 1, 1, 1, 1, 1])


    # def tearDown(self):
    #     outs = os.path.join(self.test_data, "*.out")
    #     [os.remove(p) for p in glob.glob(outs)]



if __name__ == '__main__':
    unittest.main()