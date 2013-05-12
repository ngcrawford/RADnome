#!/usr/bin/env python
# encoding: utf-8

import sys
import shlex
import argparse
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
        'infile',
        type=argparse.FileType('r'),
        nargs='?',
        default=sys.stdin,
        help='Path to input. (default is STDIN)')

    parser.add_argument(
        'outfile',
        type=argparse.FileType('w'),
        nargs='?',
        default=sys.stdout,
        help='Path to output. (default is STOUT)')

    args = parser.parse_args()
    return args


class GenerateRADGenome(object):
    """docstring for GenerateRADGenome"""
    def __init__(self):
        super(GenerateRADGenome, self).__init__()

    def get_seq(self, ID, index_path):
        """Select sequence from indexed fastq file."""

        cli = "cdbyank {0} -a {1}".format(index_path, ID)
        cli_parts = shlex.split(cli)
        line, err = Popen(cli_parts, stdin=PIPE, stdout=PIPE, stderr=PIPE).communicate()
        ID, seq, plus, quals = line.strip().split("\n")
        return (ID, seq, plus, quals)

    def __split_len__(self, seq, length):
        return [seq[i:i+length] for i in range(0, len(seq), length)]

    def _append_to_pseudo_genome_(self, seq, fout, pos, span=80):

        if pos + len(seq) > span:
            size = span - pos

            chunks = [seq[:size]]
            chunks += self.__split_len__(seq[size:], span)

            [sys.stdout.write(c+"\n") for c in chunks[:-1]]
            sys.stdout.write(chunks[-1])

        return len(chunks[-1])

    def make_pseudo_genome(self, rainbow_clustered_fin, pseudo_genome_fout, Ns):

        with rainbow_clustered_fin as fin:

            id_count = 0
            current_cluster = collections.defaultdict(list)
            previous_pos = 0

            fout = pseudo_genome_fout
            fout.write('>rainbow_clustered\n'.format(fin))
            for count, line in enumerate(fin):

                if line.startswith('E') is True:

                    if id_count != 0:

                        seq = current_cluster["S"][0]
                        previous_pos = self._append_to_pseudo_genome_(seq, fout, previous_pos)

                        Ns2add = "N"*Ns
                        previous_pos = self._append_to_pseudo_genome_(Ns2add, fout, previous_pos)

                    current_cluster = collections.defaultdict(list)
                    id_count += 1

                line_parts = line.strip().split(" ")
                if len(line_parts) == 2:
                    key, value = line_parts
                    current_cluster[key].append(value)

            fout.close()


if __name__ == '__main__':

    args = get_args()
    G = GenerateRADGenome()
    G.make_pseudo_genome(args.infile, args.outfile, args.N_padding)

