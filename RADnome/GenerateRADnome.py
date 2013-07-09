#!/usr/bin/env python
# encoding: utf-8

"""
File: GenerateRADnome.py
Author: Nicholas Crawford

Created by Nicholas Crawford on Mon Jun 17 15:35:51 PDT 2013
Copyright (c) 2012 Nicholas G. Crawford All rights reserved.

Description...

Add later.
"""


import os
import sys
import gzip
import math
import shlex
import glob
import shutil
import numpy
import pysam
import argparse
import datetime
import textwrap
import collections
import cPickle as pickle
from fileindex import FileIndex
from subprocess import Popen, PIPE
from collections import defaultdict, Counter


class SamLine(object):
    def __init__(self, fh):

        line = fh.readline().split("\t") or [None]

        if line[0].startswith("@"):
            self.qname = "@"

        elif line == ['']:
            self.qname = None

        else:
            self.qname = line[0]
            self.flag = int(line[1])
            self.rname = line[2]
            self.pos = int(line[3])
            self.mapq = int(line[4])
            self.ciagr = line[5]
            self.mrnm = line[6]
            self.mpos = int(line[7])
            self.tlen = int(line[8])
            self.seq = line[9]
            self.qual = line[10]


class Logging(object):
    """Class to hold logging methods."""
    def __init__(self):
        super(Logging, self).__init__()

        self.today = datetime.date.today()

    def __associate_contigs_log__(self, results_dict):

        results_dict["prop_asc"] =  float(results_dict["passed_filter"]) / results_dict["total_queries"]

        txt = """\
        Ran 'ascContigs' on {0[date]}

          Query file: {0[query_sam]}
          Searched file: {0[searched_sam]}
          Contig Positions file: {0[contig_positions]}
          Output Pickle file: {0[pickle]}

          Minimum MAPQ value: {0[min_mapq]}

        Results:

          Total reads processed: {0[total_queries]}
          Paired reads passing filter: {0[passed_filter]} ( = {0[prop_asc]:%} )

        """.format(results_dict)

        txt = textwrap.dedent(txt)

        outfile = self.today.strftime('ascRADnom.%m%d%y.log')
        outfile = open(outfile, 'w')
        outfile.write(txt)

        pass


    def __contigs_2_RADnome__(self, results_dict):

        results_dict["mean_unique_asc"] = float(results_dict["unique_associations"]) / results_dict['rad_frag_count']
        results_dict["mode_diff"] = float(results_dict["mode_diff"]) / results_dict['potential_rad_frags']

        txt = """\
        Ran 'makeRADnome' as {0[run_name]} on {0[date]}

          R1 file: {0[R1]}
            contig length: {0[R1_len]} bp

          R2 file: {0[R2]}
            contig length: {0[R2_len]} bp

          Output File: {0[fout]}

          Minimum proportion of reads contributing to mode: {0[proportion]}

          Padding between RADfrags: {0[N_padding]} bp
          Padding between R1/R2 contigs: {0[insert_size]} bp

        Basic Stats:

          Tested RADfrags: {0[potential_rad_frags]:g}
          Total Associated RADfrags: {0[rad_frag_count]}
          Mean unique associations per RADfrag: {0[mean_unique_asc]:g}
           (e.g., The mean number of R2 contigs associated with
            a particular R1 contig.)
          Mean proportion of reads contributing to mode: {0[mode_diff]:g}

        """.format(results_dict)

        txt = textwrap.dedent(txt)

        outfile = self.today.strftime('makeRADnome.%m%d%y.log')
        outfile = open(outfile, 'w')
        outfile.write(txt)
        pass


class GenerateRADnome(Logging):
    """Docstring for GenerateRADGenome"""
    def __init__(self):
        super(GenerateRADnome, self).__init__()

    run_info = {"Total clusters": 0,
                'seq_lengths': []
                }

    def __split_len__(self, seq, length):
        return [seq[i:i+length] for i in range(0, len(seq), length)]

    def __open_files__(self, fin, how):
        """Tests if gzip and opens appropriately."""

        if type(fin) == file or type(fin) == gzip.GzipFile:
            return fin

        if fin.endswith('.gz') is True:
            return gzip.open(fin, how)

        else:
            return open(fin, how)

    def _append_to_pseudo_genome_(self, seq, fout, pos, span=80):

        # PROCESS SHORT CHUNKS
        if (pos + len(seq)) <= span:
            fout.write(seq)
            return pos + len(seq)

        else:
            first_chunk_size = span - pos                               # deterimine size of first chunk
            first_chunk = seq[:first_chunk_size]                        # slice out chunk from seq
            chunks = [first_chunk]                                      # convert to list
            chunks += self.__split_len__(seq[first_chunk_size:], span)  # append additional chunks
            [fout.write(c+"\n") for c in chunks[:-1]]                   # write chunks
            fout.write(chunks[-1])                                      # write last chunks without '\n'

            return len(chunks[-1])

    # def get_seq(self, ID, index_path):
    #     """Select sequence from indexed fastq file."""

    #     cli = "cdbyank {0} -a {1}".format(index_path, ID)
    #     cli_parts = shlex.split(cli)
    #     line, err = Popen(cli_parts, stdin=PIPE, stdout=PIPE, stderr=PIPE).communicate()
    #     ID, seq, plus, quals = line.strip().split("\n")
    #     return (ID, seq, plus, quals)


    def make_pseudo_genome(self, rainbow_clustered_fin, pseudo_genome_fout, Ns, run_name):
        """Provide input and output file names as well as the number of Ns to insert
           between contigs.
        """

        rainbow_clustered_fin = self.__open_files__(rainbow_clustered_fin, 'rb')

        with rainbow_clustered_fin as fin:

            # SETUP COUNTERS, ETC.
            id_count = 0
            current_cluster = collections.defaultdict(list)
            previous_pos = 0
            seq_start = 0

            # PREP OUTPUT FILE
            fout = self.__open_files__(pseudo_genome_fout,'w')
            fout.write('>{0}\n'.format(run_name))

            # PREP CONTIG POSITIONS OUTPUT FILE
            contig_starts_log = open('{}.contig_start_pos.txt'.format(run_name), 'w')

            # ITERATE OVER RAINBOW ASSEMBLY
            for count, line in enumerate(fin):

                if line.startswith('E') is True:
                    self.run_info["Total clusters"] += 1

                    if id_count != 0:

                        seq = current_cluster["S"][0]
                        self.run_info['seq_lengths'].append(len(seq))

                        previous_pos = self._append_to_pseudo_genome_(seq, fout, previous_pos)

                        contig_starts_log.write("{}\n".format(seq_start))
                        seq_start += len(seq) + Ns

                        Ns2add = "N" * Ns
                        previous_pos = self._append_to_pseudo_genome_(Ns2add, fout, previous_pos)

                    current_cluster = collections.defaultdict(list)
                    id_count += 1

                line_parts = line.strip().split(" ")
                if len(line_parts) == 2:
                    key, value = line_parts
                    current_cluster[key].append(value)

            fout.write("\n")
            fout.close()
            contig_starts_log.close()

    def filter_contigs(self, c, p):
        """"""

        mx = float(c.most_common()[0][1])
        rmdr = sum([i for k, i in c.most_common()[1:]])
        diff = 1 - (rmdr / mx)
        return diff

    def __get_fasta_header__(self, fasta_path):
        header_line = open(fasta_path, 'rU').readline()
        return header_line.strip(">").strip()

    def contigs_2_RADnome(self, rad1, rad2, r1_contig_len,
                          r2_contig_len, contig_2_contig_dict,
                          run_name, N_padding, insert_size, proportion):

        run_results = {
            "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
            "pickle": contig_2_contig_dict,
            'R1': rad1,
            'R1_len': r1_contig_len,
            'R2': rad2,
            'R2_len': r2_contig_len,
            'fout': None,
            'run_name':run_name,
            'N_padding': N_padding,
            "insert_size": insert_size,
            "proportion": proportion,
            "rad_frag_count": 0,
            "potential_rad_frags": 0,
            "unique_associations": 0,
            "mode_diff": 0,
            }

        r1_name = self.__get_fasta_header__(rad1)
        r2_name = self.__get_fasta_header__(rad2)

        r1 = pysam.Fastafile(rad1)
        r2 = pysam.Fastafile(rad2)

        # today = datetime.date.today()
        fout = run_name + ".RADnome.fa"
        run_results['fout'] = fout
        fout = open(fout, 'w')

        c2c = pickle.load(open(contig_2_contig_dict, 'rb'))

        #fout = open('test.PE_RADnome.fa','w')
        count = 0
        rad_frag_count = 0
        previous_pos = 0
        for key, value in c2c.iteritems():
            run_results["potential_rad_frags"] += 1
            if key == 0:           # skip zero key as this may be where 'unalignable fragments' get put.
                continue

            # CALCULATE MODE
            counts = collections.Counter(value)
            mode = counts.most_common()[0]
            run_results['unique_associations'] += len(counts)

            # WRITE FASTA HEADER
            if count == 0:
                fout.write('>{}\n'.format(run_name))

            if self.filter_contigs(counts, proportion) >= proportion:

                run_results["rad_frag_count"] += 1
                run_results["mode_diff"] += self.filter_contigs(counts, proportion)

                # SETUP DNA FRAGMENTS
                bigNs = "N" * N_padding
                smallNs = "N" * insert_size
                contig1 = r1.fetch(r1_name, key, key+r1_contig_len)
                contig2 = r2.fetch(r2_name, mode[0], mode[0]+r2_contig_len)

                # ADD FRAGMENTS TO RADNOME
                previous_pos = self._append_to_pseudo_genome_(bigNs, fout, previous_pos, span=80)
                previous_pos = self._append_to_pseudo_genome_(contig1, fout, previous_pos, span=80)
                previous_pos = self._append_to_pseudo_genome_(smallNs, fout, previous_pos, span=80)
                previous_pos = self._append_to_pseudo_genome_(contig2, fout, previous_pos, span=80)

            count += 1

        fout.close()

        # Send data for logging
        self.__contigs_2_RADnome__(run_results)


class MergeAssemblies(Logging):
    """docstring for MergeAssemblies"""
    def __init__(self):
        super(MergeAssemblies, self).__init__()

    @staticmethod
    def get_hit_pos(s):

        try:
            return int(s.pos)

        except AttributeError:
            return 'unpaired'

    @staticmethod
    def round_bp_pos(value, nearest=5):
        if value != "unpaired":
            return int(round(float(value) / nearest) * nearest)
        else:
            return 'unpaired'


    def associate_contigs(self, sam1, sam2, contig_positions, min_mapq, run_ID):

        run_results = {
            "total_queries": 0,
            "passed_filter": 0,
            "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
            "min_mapq": min_mapq,
            "query_sam": sam1,
            "searched_sam": sam2,
            "contig_positions": contig_positions,
            "pickle": None,
            }

        contig_2_contig_dict = defaultdict(list)

        #  CREATE SAM INDEX IF IT DOESN'T EXIST
        if not os.path.exists(sam2 + FileIndex.ext):
            FileIndex.create(sam2, lambda fh: SamLine(fh).qname, allow_multiple=True)

        fi = FileIndex(sam2, SamLine, allow_multiple=True)
        ma = MergeAssemblies()

        contig_starts = [int(l.strip()) for l in open(contig_positions, 'rU')]

        with open(sam1, 'rU') as qs:

            for count, q in enumerate(qs):

                if q.startswith("@"):                      # skip header lines
                    continue

                # SPLIT LINE AND IDENTIFY QUERY POSITION.
                q = q.strip().split("\t")
                query_pos = int(q[3])
                query_pos = ma.round_bp_pos(query_pos)     # this could be more sophisticated.

                q_id = q[0][:-1] + "2"                     # create query id (e.g., ends with "2")
                q_mapq = int(q[4])                         # and, mapping quality.


                # SEARCH FOR QUERY AND PARSE RESULTS
                for s in fi[q_id]:

                    if (q_mapq > min_mapq) and (s.mapq > min_mapq):

                        run_results["passed_filter"] += 1

                        hit_pos = ma.get_hit_pos(s)
                        hit_pos = ma.round_bp_pos(hit_pos)
                        contig_2_contig_dict[query_pos].append(hit_pos)


        #CREATE PICKLE OUTPUT FILE
        today = datetime.date.today()

        path = os.path.split(sam1)[0]
        pkl_output_file_name = os.path.join(path, '{}.R1_to_R2_contig_associations.pkl'.format(run_ID))
        pkl_out = open(pkl_output_file_name, 'wb')
        pickle.dump(contig_2_contig_dict, pkl_out)

        # UPDATE RUN RESULTS AND GENERATE LOGFILE
        run_results["total_queries"] = count
        run_results['pickle'] = pkl_output_file_name
        self.__associate_contigs_log__(run_results)        # format logging results



    def __generate_test_data__(self, sam1, sam2, contig_positions, min_mapq):


        run_results = {
            "total_queries": 0,
            "passed_filter": 0,
            "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
            "min_mapq": min_mapq,
            "query_sam": sam1,
            "searched_sam": sam2,
            "contig_positions": contig_positions,
            "pickle": None,
            }

        contig_2_contig_dict = defaultdict(list)

        #  CREATE SAM INDEX IF IT DOESN'T EXIST
        if not os.path.exists(sam2 + FileIndex.ext):
            FileIndex.create(sam2, lambda fh: SamLine(fh).qname, allow_multiple=True)

        fi = FileIndex(sam2, SamLine, allow_multiple=True)
        ma = MergeAssemblies()

        contig_starts = [int(l.strip()) for l in open(contig_positions, 'rU')]

        fq1_out = open("test_fq1.10k.fq",'w')
        fq2_out = open("test_fq2.10k.fq",'w')

        with open(sam1, 'rU') as qs:

            for count, q in enumerate(qs):

                if q.startswith("@"):                      # skip header lines
                    continue

                # SPLIT LINE AND IDENTIFY QUERY POSITION.
                q = q.strip().split("\t")
                query_pos = int(q[3])
                query_pos = ma.round_bp_pos(query_pos)     # this could be more sophisticated.

                q_id = q[0][:-1] + "2"                     # create query id (e.g., ends with "2")
                q_mapq = int(q[4])                         # and, mapping quality.

                if count > 10000: break

                # SEARCH FOR QUERY AND PARSE RESULTS
                for s in fi[q_id]:

                    fq1 = (q[0], q[9], q[10])
                    fq2 = (s.qname, s.seq, s.qual)

                    fq1_line = "@{}\n{}\n+\n{}\n".format(*fq1)
                    fq2_line = "@{}\n{}\n+\n{}\n".format(*fq2)

                    fq1_out.write(fq1_line)
                    fq2_out.write(fq2_line)


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
        G.make_pseudo_genome(asm, fa, buff, run_ID)
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
        return 1

    def create_faidx(self, fa, out_path=None):

        # if out_path != None:
        #     fout_name = os.path.split(fa)[-1]
        #     fout = os.path.join(out_path, fout_name)
        # else:
        #     fout = os.path.split(fa)[-1]

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

        cli = "samtools view -h {0}.sorted.bam".format(fout)
        cli_list = shlex.split(cli)
        line, err = Popen(cli_list,
                          stdin=PIPE,
                          stdout=open("{0}.sorted.sam".format(fout), 'w'),
                          stderr=PIPE).communicate()
        return 1

    def ascContigs(self, fq1, fq2, run_ID, min_mapq, out_path=None):

        if out_path != None:
            fout1 = os.path.split(fq1)[-1]
            fout2 = os.path.split(fq2)[-1]

            fout1 = os.path.join(out_path, fout1)
            fout2 = os.path.join(out_path, fout2)

            run_ID = os.path.join(out_path, run_ID)
        else:
            fout1 = os.path.split(fq1)[-1]
            fout2 = os.path.split(fq2)[-1]

        sam1 = "{}.sorted.sam".format(fout1)
        sam2 = "{}.sorted.sam".format(fout2)

        contig_positions = "{}.R1.contig_start_pos.txt".format(run_ID)

        M = MergeAssemblies()
        log = M.associate_contigs(sam1, sam2, contig_positions, min_mapq, run_ID)
        return 1

    def make_RADnome(self, fq1, fq2, run_ID, out_path=None):

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
        N_padding = 500
        insert_size = 50
        proportion = 0.8

        print 'here', rad1

        self.create_faidx(rad1)
        self.create_faidx(rad2)

        G = GenerateRADnome()
        G.contigs_2_RADnome(rad1, rad2, r1_contig_len,
                  r2_contig_len, contig_2_contig_dict,
                  run_name, N_padding, insert_size, proportion)
        return 1

    def tidy_dir(self, fq1):

        if os.path.exists('bowtie2/') is False:
            os.mkdir("bowtie2")

        a = [shutil.move(f, 'bowtie2/') for f in glob.glob('*.bt2')]

        if os.path.exists('alignments/') is False:
            os.mkdir("alignments")

        a = [shutil.move(f, 'alignments/') for f in glob.glob('*.sam')]
        a = [shutil.move(f, 'alignments/') for f in glob.glob('*.bam')]
        a = [shutil.move(f, 'alignments/') for f in glob.glob('*.fidx')]


        if os.path.exists('fastas/') is False:
            os.mkdir("fastas")

        a = [shutil.move(f, 'fastas/') for f in glob.glob('*.fa')]
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

    def pipeline(self, fq1, fq2, run_ID, cores=3):

        # -----------
        # Run Rainbow
        # -----------

        sys.stdout.write("""Generating {0} on {1}\n\n""".format(run_ID, datetime.date.today()))

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
        self.ascContigs(fq1, fq2, run_ID, min_mapq=3)

        sys.stdout.write("Step 8: Create RADnome ...\n")
        self.make_RADnome(fq1, fq2, run_ID)

        sys.stdout.write("Step 9: Organize directory ...\n")
        self.tidy_dir(fq1)

if __name__ == '__main__':
    pass

    # fq1 = "/home/ngcrawford/Data/Nearctic_Turtles/test_fq1.fq"
    # fq2 = "/home/ngcrawford/Data/Nearctic_Turtles/test_fq2.fq"

    # fq1 = "/Users/ngcrawford/Desktop/RADs/CalAcd/test_dir/test_fq1.10k.fq"
    # fq2 = "/Users/ngcrawford/Desktop/RADs/CalAcd/test_dir/test_fq2.10k.fq"

    # R = RunRainbow()
    # R.run_rainbow(fq1, fq2, 'test_run')

    # G = MergeAssemblies()
    # sam1 = "/home/ngcrawford/Data/Nearctic_Turtles/test_data/neartic.trachemys.1.50k.sam"
    # sam2 = "/home/ngcrawford/Data/Nearctic_Turtles/READ_aligments/sams/neartic.trachemys.2.sam"
    # cntg_pos = "/home/ngcrawford/Data/Nearctic_Turtles/RADnome/neartic.trachemys.1.contig_start_pos.txt"
    # G.__generate_test_data__(sam1, sam2, cntg_pos, min_mapq=3)




    # ---------------------------------
    #  Make test data
    # ---------------------------------

