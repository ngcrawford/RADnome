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
import shlex
import numpy
import pysam
import argparse
import datetime
import collections
import cPickle as pickle
from fileindex import FileIndex
from subprocess import Popen, PIPE
from collections import defaultdict, Counter

class GenerateRADnome(object):
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

    def contigs_2_RADnome(self, rad1, rad2, r1_contig_len, r2_contig_len, fout, contig_2_contig_dict, run_name, N_padding, insert_size):

        r1 = pysam.Fastafile(rad1)
        r2 = pysam.Fastafile(rad2)
        fout = open(fout, 'w')

        c2c = pickle.load(open(contig_2_contig_dict, 'rb'))

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
                fout.write('>{}\n'.format(run_name))
            
            if self.filter_contigs(mode, value) is True:

                # SETUP DNA FRAGMENTS
                bigNs = "N" * N_padding
                smallNs = "N" * insert_size
                contig1 = r1.fetch("r1.asm.RADnome", key, key+r1_contig_len)
                contig2 = r2.fetch("r2.asm", mode[0], mode[0]+r2_contig_len)

                # ADD FRAGMENTS TO RADNOME
                previous_pos = self._append_to_pseudo_genome_(bigNs, fout, previous_pos, span=80)
                previous_pos = self._append_to_pseudo_genome_(contig1, fout, previous_pos, span=80)
                previous_pos = self._append_to_pseudo_genome_(smallNs, fout, previous_pos, span=80)
                previous_pos = self._append_to_pseudo_genome_(contig2, fout, previous_pos, span=80)
      
            count += 1

        fout.close()
        self.make_log_file()


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


class MergeAssemblies(object):
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

    @staticmethod
    def associate_contigs(sam1, sam2, contig_positions):

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
                
                if count % 50000 == 0:
                    print count
                
                # SEARCH FOR QUERY AND PARSE RESULTS
                for s in fi[q_id]:
                    hit_pos = ma.get_hit_pos(s)
                    hit_pos = ma.round_bp_pos(hit_pos)
                    contig_2_contig_dict[query_pos].append(hit_pos)

        print 'dumping assoications'
        associations = open('contig_2_contig_associations_dict.pkl', 'wb')
        pickle.dump(contig_2_contig_dict, associations)


if __name__ == '__main__':



    queries = 'r1.sorted.sam'
    sam = 'r2.sorted.sam'
    contig_positions = "r1.asm.contig_start_pos.txt"




