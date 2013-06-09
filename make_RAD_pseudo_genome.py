#!/usr/bin/env python
# encoding: utf-8

import sys
import gzip
import shlex
import numpy
import pysam
import argparse
import datetime
import collections
from subprocess import Popen, PIPE

"""
File: make_RAD_pseudo_genome.py
Author: Nicholas Crawford

Created by Nicholas Crawford on Thu May  2 14:52:44 EDT 2013
Copyright (c) 2012 Nicholas G. Crawford All rights reserved.

Description...

Add later.

"""


def get_args():
    """Parse sys.argv"""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-n', '--N-padding',
        type=int,
        default=500,
        help="Number of N's to insert between each cluster/contig.")

    parser.add_argument(
        '-r', '--run-name',
        type=str,
        default='RADnome',
        help="Name of run. Becomes fasta name.")

    parser.add_argument(
        'infile',
        nargs='?',
        default=sys.stdin,
        help='Path to input. (default is STDIN)')

    parser.add_argument(
        'outfile',
        nargs='?',
        default=sys.stdout,
        help='Path to output. (default is STOUT)')

    args = parser.parse_args()
    return args


class GenerateRADGenome(object):
    """Docstring for GenerateRADGenome"""
    def __init__(self):
        super(GenerateRADGenome, self).__init__()

    run_info = {"Total clusters": 0,
                'seq_lengths': []
                }

    def __split_len__(self, seq, length):
        return [seq[i:i+length] for i in range(0, len(seq), length)]

    def __open_files__(self, fin, how):
        """Tests if gzip."""

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
            [fout.write(c+"\n") for c in chunks[:-1]]             # write chunks
            fout.write(chunks[-1])                                # write last chunks without '\n'

            return len(chunks[-1])

    def get_seq(self, ID, index_path):
        """Select sequence from indexed fastq file."""

        cli = "cdbyank {0} -a {1}".format(index_path, ID)
        cli_parts = shlex.split(cli)
        line, err = Popen(cli_parts, stdin=PIPE, stdout=PIPE, stderr=PIPE).communicate()
        ID, seq, plus, quals = line.strip().split("\n")
        return (ID, seq, plus, quals)


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

            # PREP OUTPUT FIEL
            fout = self.__open_files__(pseudo_genome_fout,'w')
            fout.write('>{0}\n'.format(run_name))

            for count, line in enumerate(fin):

                if line.startswith('E') is True:
                    self.run_info["Total clusters"] += 1

                    if id_count != 0:

                        seq = current_cluster["S"][0]
                        self.run_info['seq_lengths'].append(len(seq))

                        previous_pos = self._append_to_pseudo_genome_(seq, fout, previous_pos)

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

    def make_log_file(self):

        self.run_info['Input'] = args.infile
        self.run_info['Output'] = args.outfile
        self.run_info['Run Name'] = args.run_name
        self.run_info['Padding (bases)'] = args.N_padding
        self.run_info['Date'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

        with open("{}.log".format(self.run_info['Run Name']), 'w') as log_file:

            log_file.write("Ran make_RAD_pseudo_genome.py with the following parameters:\n\n")


            for key in ['Date', 'Run Name', 'Input', 'Output', 'Padding (bases)', 'Total clusters', 'seq_lengths']:
                value = self.run_info[key]
                if key is 'seq_lengths':
                    values = numpy.array(value)
                    log_file.write("  Mean Seq Length: {0} +/- {1} stdev\n".format(values.mean(), values.std()))

                else:
                    log_file.write("  {0}: {1}\n".format(key, value))


    def contigs_2_rad_frags(bam1, bam2):


        
        pass





if __name__ == '__main__':

    args = get_args()
    G = GenerateRADGenome()
    #G.make_pseudo_genome(args.infile, args.outfile, args.N_padding, args.run_name)
    #G.make_log_file()

    G.contigs_2_rad_frags(bam1, bam2)

