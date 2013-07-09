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
        self.run_ID = 'test'
        pass

    def test_pipeline(self):
        """Run comlete pipline with test data."""

        P = GenerateRADnome.RunPipeline()

        #  SETUP PATHS TO OUTPUT FILES
        fq1 = os.path.join(self.test_data, "fq1.10k.fq")
        fq2 = os.path.join(self.test_data, "fq2.10k.fq")

        R1 = os.path.join(self.test_data, "test.R1")
        R2 = os.path.join(self.test_data, "test.R2")

        run_ID = self.run_ID

        #  RUN STEPS OF PIPLINE
        a = P.run_cluster_cmd(fq1)
        b = P.run_cluster_cmd(fq2)

        c = P.run_div_cmd(fq1, out_path=self.test_data)
        d = P.run_div_cmd(fq2, out_path=self.test_data)

        e = P.run_merge_cmd(fq1, out_path=self.test_data)
        f = P.run_merge_cmd(fq2, out_path=self.test_data)

        g = P.make_READnome(fq1, R1, out_path=self.test_data)
        h = P.make_READnome(fq2, R2, out_path=self.test_data)

        i = P.run_bowtie2(fq1, cores=1, out_path=self.test_data)
        j = P.run_bowtie2(fq2, cores=1, out_path=self.test_data)

        l = P.sam_to_sorted_sam(fq1, out_path=self.test_data)
        m = P.sam_to_sorted_sam(fq2, out_path=self.test_data)

        n = P.ascContigs(fq1, fq2, run_ID, min_mapq=3,
                         out_path=self.test_data)

        o = P.make_RADnome(fq1, fq2, run_ID, out_path=self.test_data)

        results =                 [a, b, c, d, e, f, g, h, i, j, l, m, n, o]
        self.assertEqual(results, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])


    def tearDown(self):
        """Remove intermediate files from test/data directory."""

        files = os.path.join(self.test_data, "*.txt")
        [os.remove(p) for p in glob.glob(files)]

        files = os.path.join(self.test_data, "fq*.10k.fq.*")
        [os.remove(p) for p in glob.glob(files)]

        files = os.path.join(self.test_data, "{}*".format(self.run_ID))
        [os.remove(p) for p in glob.glob(files)]







if __name__ == '__main__':
    unittest.main()