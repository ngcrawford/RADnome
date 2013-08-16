#!/usr/bin/env python
# encoding: utf-8

"""
File: pipeline.py
Author: Nicholas Crawford

Created by Nicholas Crawford on Mon Jun 17 15:35:51 PDT 2013
Copyright (c) 2012 Nicholas G. Crawford All rights reserved.

Description:

Runs the entire RADnome generation code as a pipeline.
"""

import os
import sys
import glob
import shlex
import shutil
import datetime
from subprocess import Popen, PIPE
from RADnome.GenerateRADnome import *


class RunPipeline(object):
    """docstring for RunPipeline"""
    def __init__(self):
        super(RunPipeline, self).__init__()


    def run_cluster_cmd(self, fq_id):

        cli = "rainbow cluster -1 {0}".format(fq_id)
        cli_list = shlex.split(cli)
        line, err = Popen(cli_list,
                          stdin=PIPE,
                          stdout=open("{}.cluster.out".format(fq_id), 'w'),
                          stderr=PIPE).communicate()
        return 1


    def run_div_cmd(self, fq_id, out_path=None):

        if out_path != None:
            fout_name = os.path.split(fq_id)[-1]
            fout = os.path.join(out_path, fout_name)
        else:
            fout = os.path.split(fq_id)[-1]

        cli = "rainbow div -i {0}.cluster.out -o {1}.div.out".format(fq_id, fout)
        cli_list = shlex.split(cli)
        line, err = Popen(cli_list,
                          stdin=PIPE,
                          stdout=PIPE,
                          stderr=PIPE).communicate()
        return 1


    def run_merge_cmd(self, fq_id, out_path=None):

        if out_path != None:
            fout_name = os.path.split(fq_id)[-1]
            fout = os.path.join(out_path, fout_name)
        else:
            fout = os.path.split(fq_id)[-1]

        cli = "rainbow merge -a -i {0}.div.out -o {1}.asm.out".format(fq_id, fout)
        cli_list = shlex.split(cli)
        line, err = Popen(cli_list,
                          stdin=PIPE,
                          stdout=PIPE,
                          stderr=PIPE).communicate()
        return 1

    def run_gzip(self, fq_id):

        cli = "gzip {0}".format(fq_id)
        cli_list = shlex.split(cli)
        line, err = Popen(cli_list,
                          stdin=PIPE,
                          stdout=PIPE,
                          stderr=PIPE).communicate()
        return 1

    def make_READnome(self, fq_id, run_ID, buff=500, out_path=None):

        if out_path != None:
            fout_name = os.path.split(fq_id)[-1]
            fout = os.path.join(out_path, fout_name)
        else:
            fout = os.path.split(fq_id)[-1]

        asm = "{}.asm.out".format(fout)
        fa = "{}.asm.fa".format(fout)

        G = GenerateRADnome()
        G.make_READnome(asm, fa, buff, run_ID, out_path)
        return 1

    def run_bowtie2(self, fq_id, cores, out_path=None):

        if out_path != None:
            fout_name = os.path.split(fq_id)[-1]
            fout = os.path.join(out_path, fout_name)
        else:
            fout = os.path.split(fq_id)[-1]

        cli = "bowtie2-build -f {0}.asm.fa {0}.asm".format(fout)
        cli_list = shlex.split(cli)
        line, err = Popen(cli_list,
                          stdin=PIPE,
                          stdout=PIPE,
                          stderr=PIPE).communicate()

        cli = "bowtie2 \
               --very-sensitive \
               --end-to-end \
               -p {1} \
               -x {0}.asm \
               -U {2}".format(fout, cores, fq_id)

        cli_list = shlex.split(cli)
        line, err = Popen(cli_list,
                          stdin=PIPE,
                          stdout=open("{}.sam".format(fout), 'w'),
                          stderr=PIPE).communicate()

        err_out = open('{}.bwt2.error.log'.format(fq_id), 'w')
        [err_out.write(line) for line in err]
        return 1

    def create_faidx(self, fa, out_path=None):

        if out_path != None:
            fout_name = os.path.split(fa)[-1]
            fout = os.path.join(out_path, fout_name)
        else:
            fout = os.path.split(fa)[-1]

        cli = "samtools faidx {}".format(fa)
        cli_list = shlex.split(cli)
        line, err = Popen(cli_list,
                          stdin=PIPE,
                          stdout=PIPE,
                          stderr=PIPE).communicate()


    def sam_to_sorted_sam(self, fq_id, out_path=None):

        if out_path != None:
            fout_name = os.path.split(fq_id)[-1]
            fout = os.path.join(out_path, fout_name)
        else:
            fout = os.path.split(fq_id)[-1]

        cli = "samtools view -bS {}.sam".format(fout)
        cli_list = shlex.split(cli)
        line, err = Popen(cli_list,
                          stdin=PIPE,
                          stdout=open("{}.bam".format(fout), 'wb'),
                          stderr=PIPE).communicate()

        cli = "samtools sort {0}.bam {0}.sorted".format(fout)
        cli_list = shlex.split(cli)
        line, err = Popen(cli_list,
                          stdin=PIPE,
                          stdout=PIPE,
                          stderr=PIPE).communicate()

        cli = "samtools index {0}.sorted.bam".format(fout)
        cli_list = shlex.split(cli)
        line, err = Popen(cli_list,
                          stdin=PIPE,
                          stdout=PIPE,
                          stderr=PIPE).communicate()

        cli = "samtools view -h {0}.sorted.bam".format(fout)
        cli_list = shlex.split(cli)
        line, err = Popen(cli_list,
                          stdin=PIPE,
                          stdout=open("{0}.sorted.sam".format(fout), 'w'),
                          stderr=PIPE).communicate()
        return 1

    def ascContigs(self, fq1, fq2, run_ID, min_mapq, min_depth, out_path=None):

        if out_path != None:
            fout1 = os.path.split(fq1)[-1]
            fout2 = os.path.split(fq2)[-1]

            fout1 = os.path.join(out_path, fout1)
            fout2 = os.path.join(out_path, fout2)

            run_ID = os.path.join(out_path, run_ID)
        else:
            fout1 = os.path.split(fq1)[-1]
            fout2 = os.path.split(fq2)[-1]

        bam1 = "{}.sorted.bam".format(fout1)
        sam2 = "{}.sorted.sam".format(fout2)

        contig_positions = "{}.R1.contig_start_pos.txt".format(run_ID)

        M = MergeAssemblies()
        log = M.associate_contigs(bam1, sam2, min_mapq, min_depth, run_ID)
        return 1

    def make_RADnome(self, fq1, fq2, run_ID, N_padding, insert_size,
                     proportion, overlap, out_path=None):

        def seq_len(fq):
            with open(fq,'rU') as fin:
                fin.readline()
                return len(fin.readline().strip())

        if out_path != None:
            fout1 = os.path.split(fq1)[-1]
            fout2 = os.path.split(fq2)[-1]

            fout1 = os.path.join(out_path, fout1)
            fout2 = os.path.join(out_path, fout2)

            run_ID = os.path.join(out_path, run_ID)
        else:
            fout1 = os.path.split(fq1)[-1]
            fout2 = os.path.split(fq2)[-1]

        rad1 = "{}.asm.fa".format(fout1)
        rad2 = "{}.asm.fa".format(fout2)

        r1_contig_len = seq_len(fq1)
        r2_contig_len = seq_len(fq2)

        contig_2_contig_dict = "{}.R1_to_R2_contig_associations.pkl".format(run_ID)
        run_name = run_ID
        # N_padding = 500
        # insert_size = 50
        # proportion = 0.8
        # overlap = 10

        self.create_faidx(rad1)
        self.create_faidx(rad2)

        G = GenerateRADnome()
        G.contigs_2_RADnome(rad1, rad2, r1_contig_len,
                  r2_contig_len, contig_2_contig_dict,
                  run_name, N_padding, insert_size,
                  proportion, overlap)
        return 1

    def merge_fastas(self, run_ID):
        Paired_nome_out = run_ID + ".paired_contigs.fa"
        Singtons_R1_nome_out = run_ID + ".singleton_R1_contigs.fa"
        Singtons_R2_nome_out = run_ID + ".singleton_R2_contigs.fa"
        Contig_nome_out = run_ID + ".assembled_contigs.fa"

        cli = "cat {0} {1} {2} {3}".format(Paired_nome_out, Contig_nome_out,
                                       Singtons_R1_nome_out, Singtons_R2_nome_out)
        cli_list = shlex.split(cli)
        line, err = Popen(cli_list,
                          stdin=PIPE,
                          stdout=open("{0}.RADnome.fa".format(run_ID), 'w'),
                          stderr=PIPE).communicate()
        return 1

    def tidy_dir(self, fq1, run_ID):

        if os.path.exists('bowtie2/') is False:
            os.mkdir("bowtie2")

        a = [shutil.move(f, 'bowtie2/') for f in glob.glob('*.bt2')]

        if os.path.exists('alignments/') is False:
            os.mkdir("alignments")

        a = [shutil.move(f, 'alignments/') for f in glob.glob('*.sam')]
        a = [shutil.move(f, 'alignments/') for f in glob.glob('*.bam')]
        a = [shutil.move(f, 'alignments/') for f in glob.glob('*.bai')]
        a = [shutil.move(f, 'alignments/') for f in glob.glob('*.fidx')]


        if os.path.exists('fastas/') is False:
            os.mkdir("fastas")

        fastas = glob.glob('*.fa')
        fastas.remove("{}.RADnome.fa".format(run_ID))

        a = [shutil.move(f, 'fastas/') for f in fastas]
        a = [shutil.move(f, 'fastas/') for f in glob.glob('*.fai')]
        a = [shutil.move(f, 'fastas/') for f in glob.glob('*.contig_start*')]
        a = [shutil.move(f, 'fastas/') for f in glob.glob('*.pkl')]

        if os.path.exists('rainbow/') is False:
            os.mkdir("rainbow")

        a = [shutil.move(f, 'rainbow/') for f in glob.glob('*.out')]

        if os.path.exists('logs/') is False:
            os.mkdir("logs")

        a = [shutil.move(f, 'logs/') for f in glob.glob('*.log')]
        return 1

    def pipeline(self, fq1, fq2, run_ID, N_padding, insert_size,
                     proportion, overlap, cores):

        # -----------
        # Run Rainbow
        # -----------

        sys.stdout.write("""Generating RADnome: '{0}'' on {1}\n\n""".format(run_ID, datetime.date.today()))

        sys.stdout.write("Step 1: Running Rainbow Cluster ...\n")
        self.run_cluster_cmd(fq1)
        self.run_cluster_cmd(fq2)

        sys.stdout.write("Step 2: Running Rainbow Div ...\n")
        self.run_div_cmd(fq1)
        self.run_div_cmd(fq2)

        sys.stdout.write("Step 3: Running Rainbow Merge ...\n")
        self.run_merge_cmd(fq1)
        self.run_merge_cmd(fq2)

        # -------------------------------
        # Generate READnomes and RADnomes
        # -------------------------------

        sys.stdout.write("Step 4: Creating R1 and R2 READnomes ...\n")
        self.make_READnome(fq1, "{}.R1".format(run_ID))
        self.make_READnome(fq2, "{}.R2".format(run_ID))


        sys.stdout.write("Step 5: Running Bowtie2 ...\n")

        sys.stdout.write("  .. aligning {}\n".format(fq1))
        self.run_bowtie2(fq1, cores)

        sys.stdout.write("  .. aligning {}\n".format(fq2))
        self.run_bowtie2(fq2, cores)


        sys.stdout.write("Step 6: Creating sorted BAMs ...\n")
        self.sam_to_sorted_sam(fq1)
        self.sam_to_sorted_sam(fq2)

        sys.stdout.write("Step 7: Associating contigs ...\n")
        self.ascContigs(fq1, fq2, run_ID, min_depth=1, min_mapq=3)

        sys.stdout.write("Step 8: Create RADnome ...\n")
        self.make_RADnome(fq1, fq2, run_ID, N_padding,
                          insert_size, proportion, overlap,)
        self.merge_fastas(run_ID)

        sys.stdout.write("Step 9: Organize directory ...\n")
        self.tidy_dir(fq1, run_ID)