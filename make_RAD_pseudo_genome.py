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
import cPickle as pickle
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

    parser = argparse.ArgumentParser(prog='RADnome')

    parser.add_argument(
        '-n', '--N-padding',
        type=int,
        default=500,
        help="Number of N's to insert between each cluster/contig.")




    parser_




    parser.add_argument(
        '-i', '--insert-size',
        type=int,
        default=50,
        help="Size of insert when libraries were generated. Used to deterimine\
              how many N's to insert between associated contigs.")

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
            [fout.write(c+"\n") for c in chunks[:-1]]                   # write chunks
            fout.write(chunks[-1])                                      # write last chunks without '\n'

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


    def filter_contigs(self, mode, total):
        """Mode must match total."""

        if float(mode[-1]) / len(total) == 1.0:
            return True
        else:
            return False

    def contigs_2_RADnome(self, rad1, rad2, fout, contig_2_contig_dict, args):

        r1 = pysam.Fastafile(rad1)
        r2 = pysam.Fastafile(rad2)
        fout = open(fout,'w')

        c2c = contig_2_contig_dict

        #fout = open('test.PE_RADnome.fa','w')
        count = 0
        previous_pos = 0
        for key, value in c2c.iteritems():
            # if count > 10: break
            if key == 0:           # skip zero key as I think this may be where 'unalignable fragments' get put.
                continue
            
            # CALCULATE MODE    
            counts = collections.Counter(value)
            mode = counts.most_common()[0]

            # WRITE FASTA HEADER
            if count == 0:
                fout.write('>{}\n'.format(args.run_name))
            
            if self.filter_contigs(mode, value) is True:

                # SETUP DNA FRAGMENTS
                bigNs = "N" * args.N_padding
                smallNs = "N" * args.insert_size
                contig1 = r1.fetch("r1.asm.RADnome", key, key+95)
                contig2 = r2.fetch("r2.asm", mode[0], mode[0]+100)

                # ADD FRAGMENTS TO RADNOME
                previous_pos = self._append_to_pseudo_genome_(bigNs, fout, previous_pos, span=80)
                previous_pos = self._append_to_pseudo_genome_(contig1, fout, previous_pos, span=80)
                previous_pos = self._append_to_pseudo_genome_(smallNs, fout, previous_pos, span=80)
                previous_pos = self._append_to_pseudo_genome_(contig2, fout, previous_pos, span=80)
      
            count += 1

        fout.close()

if __name__ == '__main__':

    args = get_args()
    G = GenerateRADGenome()
    # G.make_pseudo_genome(args.infile, args.outfile, args.N_padding, args.run_name)
    # G.make_log_file()

    rad1 = "/Users/ngcrawford/Dropbox/JMP_CAS/Nick/RADnome/r1.asm.fa"
    rad2 = "/Users/ngcrawford/Dropbox/JMP_CAS/Nick/RADnome/r2.asm.fa"
    c2c = pickle.load(open('contig_2_contig_associations_dict.pkl', 'rb'))
    fout = 'test.RADnome.fa'

    G.contigs_2_RADnome(rad1, rad2, fout, c2c, args)



