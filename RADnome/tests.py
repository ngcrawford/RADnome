#!/usr/bin/env python
# encoding: utf-8

import os
import sys

sys.path.insert(0, os.path.abspath('..'))  # Seriously Python?! This is fucking ugly.

import glob
import shlex
import shutil

import unittest
from subprocess import Popen, PIPE
import pipeline
from RADnome import GenerateRADnome



class TestCmds(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        """Define shared variables, classes, and decompress tar.gz"""

        # Setup Variables
        self.run_ID = 'test_run'
        self.min_mapq = 3
        self.min_depth = 1

        # Setup Classes
        self.P = pipeline.RunPipeline()
        self.C = GenerateRADnome.ContigAssembler()
        self.M = GenerateRADnome.MergeAssemblies()
        self.G = GenerateRADnome.GenerateRADnome()

        # Decompress tar.gz (using TarFile module caused issues)
        module_dir = os.path.dirname(GenerateRADnome.__file__)
        self.base_dir = os.path.split(module_dir)[0]
        self.test_data = os.path.join(self.base_dir, "tests/data")
        self.data_archive = os.path.join(self.test_data, "data.tar.gz")

        os.chdir(self.test_data)
        cli = "tar -xzf {0}".format("data.tar.gz")
        cli_list = shlex.split(cli)
        line, err = Popen(cli_list,
                          stdin=PIPE,
                          stdout=PIPE,
                          stderr=PIPE).communicate()

        # Make tmp dir (TO DO redo with 'real' tempdir)
        os.mkdir(os.path.join(self.base_dir, "tests/data/tmp"))

        # Reset working directory
        os.chdir(self.base_dir)

    @classmethod
    def tearDownClass(self):
        """Clean up the mess (e.g., data and tmp dir)"""

        test_data = os.path.join(self.base_dir, "tests/data/data")
        shutil.rmtree(test_data)
        tmp = os.path.join(self.base_dir, "tests/data/tmp")
        shutil.rmtree(tmp)

    def test_cluster_cmd(self):
        """Test rainbow clustering wrapper."""

        fq1 = os.path.join(self.test_data, "data/fq1.10k.fq")
        z = self.P.run_cluster_cmd(fq1)
        self.assertTrue(z)

    def test_div_cmd(self):
        """Test rainbow dividing wrapper."""

        fq1 = os.path.join(self.test_data, "data/rainbow/fq1.10k.fq")
        z = self.P.run_div_cmd(fq1, out_path="data/tmp")
        self.assertTrue(z)

    def test_merge_cmd(self):
        """Test rainbow merging wrapper."""

        fq1 = os.path.join(self.test_data, "data/rainbow/fq1.10k.fq")
        z = self.P.run_merge_cmd(fq1, out_path="data/tmp")
        self.assertTrue(z)

    def test_READnome_cmd(self):
        """Test generating READnome method."""

        asm = os.path.join(self.test_data, "data/rainbow/fq1.10k.fq.asm.out")
        fa = os.path.join(self.test_data, "tmp/test.fa")

        buff = 500
        run_ID = self.run_ID

        z = self.G.make_READnome(asm, fa, buff, run_ID)
        self.assertTrue(z)

    def test_ascContigs(self):
        """Test associating contigs method."""

        bam1 = os.path.join(self.test_data, "data/alignments/fq1.10k.fq.sorted.bam")
        sam2 = os.path.join(self.test_data, "data/alignments/fq2.10k.fq.sorted.sam")
        contig_positions = os.path.join(self.test_data, "data/fastas/test_run.R1_to_R2_contig_associations.pkl")

        run_ID = self.run_ID
        min_mapq = self.min_mapq
        min_depth = self.min_depth

        z = self.M.associate_contigs(bam1, sam2, min_mapq, min_depth, run_ID)
        self.assertTrue(z)

    def test_contigs_2_RADnome(self):
        """Test generating RADnome"""

        def seq_len(fq):
            with open(fq,'rU') as fin:
                fin.readline()
                return len(fin.readline().strip())

        rad1 = os.path.join(self.test_data, "data/fastas/fq1.10k.fq.asm.fa")
        rad2 = os.path.join(self.test_data, "data/fastas/fq2.10k.fq.asm.fa")

        r1_contig_len = seq_len(os.path.join(self.test_data, "data/fq1.10k.fq"))
        r2_contig_len = seq_len(os.path.join(self.test_data, "data/fq2.10k.fq"))
        contig_2_contig_dict = os.path.join(self.test_data, "data/fastas/test_run.R1_to_R2_contig_associations.pkl")

        R2_starts = os.path.join(self.test_data, 'data/alignments/{}.R2.contig_start_pos.txt')
        remaining_R1s = os.path.join(self.test_data, 'data/alignments/{}.R1.contig_start_pos.no_pass.txt')
        run_name = self.run_ID
        N_padding = 500
        insert_size = 50
        proportion = 0.8
        overlap = 10

        out_path = os.path.join(self.test_data, "tmp")

        z = self.G.contigs_2_RADnome(rad1,
                                     rad2,
                                     r1_contig_len,
                                     r2_contig_len,
                                     contig_2_contig_dict,
                                     run_name,
                                     N_padding,
                                     insert_size,
                                     proportion,
                                     overlap,
                                     out_path,
                                     R2_starts,
                                     remaining_R1s)
        self.assertTrue(z)


class TestPipeline(unittest.TestCase):

    def setUp(self):
        module_dir = os.path.dirname(GenerateRADnome.__file__)
        base_dir = os.path.split(module_dir)[0]
        self.test_data = os.path.join(base_dir, "tests/data")
        self.run_ID = 'test'
        pass

    def test_pipeline(self):
        """Run complete pipline with test data."""

        P = pipeline.RunPipeline()

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

        n = P.ascContigs(fq1, fq2,  min_mapq=3, min_depth=10, 
                         run_ID=run_ID, out_path=self.test_data)


        #def make_RADnome(self, fq1, fq2, run_ID, N_padding, insert_size,
        #             proportion, overlap, out_path=None)

        o = P.make_RADnome(fq1, fq2, run_ID, N_padding=500,
                           insert_size=50, proportion=0.8,
                           overlap=10, out_path=self.test_data)

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

