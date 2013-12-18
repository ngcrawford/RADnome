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
import pysam
import bisect
import datetime
import textwrap
import collections
# import numpy as np
import cPickle as pickle
from fileindex import FileIndex
from collections import defaultdict


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


    def __make_READnome_log__(self, results_dict, path):

        results_dict['contig_count'] = len(results_dict["seq_lengths"])
        results_dict['mean_contig_len'] = sum(results_dict["seq_lengths"]) / float(len(results_dict["seq_lengths"]))

        txt = """\
        Ran 'makeREADnome' as {0[run_name]} on {0[date]}

          Rainbow Assembly: {0[fin]}

          Pseudo Genome: {0[fout]}

          Length of N padding: {0[Ns]:g}

        Basic Stats:

          Total Contigs: {0[contig_count]:g}
          Mean Contig length: {0[mean_contig_len]}

        """.format(results_dict)

        txt = textwrap.dedent(txt)
        fin = os.path.split(results_dict["fin"])[-1]
        outfile = self.today.strftime('makeREADnome.{}.%m%d%y.log'.format(fin))
        outfile = os.path.join(path, outfile)
        outfile = open(outfile, 'w')
        outfile.write(txt)
        outfile.close()
        pass

    def __associate_contigs_log__(self, results_dict, path):

        results_dict["prop_asc"] = float(results_dict["passed_filter"]) / results_dict["total_queries"]

        txt = """\
        Ran 'ascContigs' on {0[date]}

          Query file: {0[query_sam]}
          Searched file: {0[searched_sam]}
          Output Pickle file: {0[pickle]}

          Minimum MAPQ value: {0[min_mapq]}

        Results:

          Total reads processed: {0[total_queries]}
          Paired reads passing filter: {0[passed_filter]}
          Unpaired R1s: {0[failed_filter]}

        """.format(results_dict)

        txt = textwrap.dedent(txt)

        outfile = self.today.strftime('ascRADnom.%m%d%y.log')

        outfile = os.path.join(path, outfile)
        outfile = open(outfile, 'w')
        outfile.write(txt)
        outfile.close()
        pass


    def __contigs_2_RADnome__(self, results_dict, path):

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

        Counts:

            Total associated contigs: {0[associated_contigs]}
            Total assembled contigs: {0[assembled_contigs]}
            Total R1 singleton contigs: {0[R1_singletons]}
            Total R2 singleton contigs: {0[R2_singletons]}

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
        outfile = os.path.join(path, outfile)
        outfile = open(outfile, 'w')
        outfile.write(txt)
        outfile.close()
        pass


class GenerateRADnome(Logging):
    """Docstring for GenerateRADGenome"""
    def __init__(self):
        super(GenerateRADnome, self).__init__()

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


    def make_READnome(self, rainbow_clustered_fin, pseudo_genome_fout, Ns, run_name, out_path=None):
        """Provide input and output file names as well as the number of Ns to insert
           between contigs.
        """

        results_dict = {"date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
                        "fin": rainbow_clustered_fin,
                        "fout": pseudo_genome_fout,
                        "Ns": Ns,
                        "run_name": run_name,
                        "seq_lengths": []
                        }

        rainbow_clustered_fin = self.__open_files__(rainbow_clustered_fin, 'rb')

        with rainbow_clustered_fin as fin:

            # SETUP COUNTERS, ETC.
            id_count = 0
            current_cluster = collections.defaultdict(list)
            previous_pos = 0
            seq_start = 0

            # PREP OUTPUT FILE
            fout = self.__open_files__(pseudo_genome_fout, 'w')
            fout.write('>{0}\n'.format(run_name))

            # PREP CONTIG POSITIONS OUTPUT FILE
            path = os.path.split(pseudo_genome_fout)[0]
            contig_starts_log = os.path.join(path, '{}.contig_start_pos.txt'.format(run_name))
            contig_starts_log = open(contig_starts_log, 'w')

            # ITERATE OVER RAINBOW ASSEMBLY
            for count, line in enumerate(fin):

                if line.startswith('E') is True:

                    if id_count != 0:

                        seq = current_cluster["S"][0]

                        if seq == 'N': # skip empty sequences.
                            continue

                        results_dict["seq_lengths"].append(len(seq))

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

            self.__make_READnome_log__(results_dict, path)
            return 1

    def make_READnome_from_Stacks_tsv(self, tsv_fin, pseudo_genome_fout, Ns, run_name, out_path=None):
        """Provide input and output file names as well as the number of Ns to insert
           between contigs.
        """

        results_dict = {"date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
                        "fin": tsv_fin,
                        "fout": pseudo_genome_fout,
                        "Ns": Ns,
                        "run_name": run_name,
                        "seq_lengths": []
                        }

        tsv_fin = self.__open_files__(tsv_fin, 'rb')

        with tsv_fin as fin:

            # SETUP COUNTERS, ETC.
            id_count = 0
            current_cluster = collections.defaultdict(list)
            previous_pos = 0
            seq_start = 0

            # PREP OUTPUT FILE
            fout = self.__open_files__(pseudo_genome_fout, 'w')
            fout.write('>{0}\n'.format(run_name))

            # PREP CONTIG POSITIONS OUTPUT FILE
            path = os.path.split(pseudo_genome_fout)[0]
            contig_starts_log = os.path.join(path, '{}.contig_start_pos.txt'.format(run_name))
            contig_starts_log = open(contig_starts_log, 'w')

            # ITERATE OVER TSV CONTIGS FILE
            for count, line in enumerate(fin):

                line_parts = line.strip().split()
                seq = line_parts[8]

                results_dict["seq_lengths"].append(len(seq))

                previous_pos = self._append_to_pseudo_genome_(seq, fout, previous_pos)

                contig_starts_log.write("{}\n".format(seq_start))
                seq_start += len(seq) + Ns

                Ns2add = "N" * Ns
                previous_pos = self._append_to_pseudo_genome_(Ns2add, fout, previous_pos)

            fout.write("\n")
            fout.close()

            self.__make_READnome_log__(results_dict, path)
            return 1

    def filter_contigs(self, c, p):
        """"""

        mx = float(c.most_common()[0][1])
        rmdr = sum([i for k, i in c.most_common()[1:]])
        diff = 1 - (rmdr / mx)
        return diff

    def __get_fasta_header__(self, fasta_path):
        header_line = open(fasta_path, 'rU').readline()
        return header_line.strip(">").strip()

    def __make_consensus__(self, seq):
        cons = ""
        for s in seq:
            if s == 'A':
                cons += 'T'
            elif s == 'T':
                cons += 'A'
            elif s == 'C':
                cons += 'G'
            elif s == 'G':
                cons += 'C'
            else:
                cons += s

        return cons

    # def test(self):
    #     unique_items = set([item for sublist in c2c.values() for item in sublist])
    #     min_i = min(z)
    #     max_i = max(z)
    #     full_set = set(xrange(min_i, max_i + 1, 600))
    #     print sorted(list(full_set - unique_items))

    def write_contig_start_stop_info(self, outfile, run_name, current_position, start_stop_info, pseudo_genome_type):

        if pseudo_genome_type == "RADnome":
            bigNs, contig1, smallNs, contig2 = start_stop_info
            current_position += len(bigNs) + len(contig1) + len(smallNs) + len(contig2)
            start = (current_position - (len(contig1) + len(smallNs) + len(contig2))) + 1
            stop = current_position
            line = '{0}_NonOverlapping {1} {2}\n'.format(run_name, start, stop)
            outfile.write(line)

        elif pseudo_genome_type == 'Contignome':
            bigNs, assembled_contig = start_stop_info
            current_position += len(bigNs) + len(assembled_contig)
            start = (current_position - len(assembled_contig)) + 1
            stop = current_position
            line = '{0}_Assembled {1} {2}\n'.format(run_name, start, stop)
            outfile.write(line)

        elif pseudo_genome_type == 'R1nome':
            bigNs, contig1 = start_stop_info
            current_position += len(bigNs) + len(contig1)
            start = (current_position - len(contig1)) + 1
            stop = current_position
            line = '{0}_R1_Singletons {1} {2}\n'.format(run_name, start, stop)
            outfile.write(line)

        elif pseudo_genome_type == 'R2nome':
            bigNs, contig1 = start_stop_info
            current_position += len(bigNs) + len(contig1)
            start = (current_position - len(contig1)) + 1
            stop = current_position
            line = '{0}_R2_Singletons {1} {2}\n'.format(run_name, start, stop)
            outfile.write(line)

        return current_position

    def __find_nearest__(self, array, value):
        n = [abs(i-value) for i in array]
        idx = n.index(min(n))
        return array[idx]

    def contigs_2_RADnome(self, rad1, rad2, r1_contig_len,
                          r2_contig_len, contig_2_contig_dict,
                          run_name, N_padding, insert_size,
                          proportion, overlap, out_path=None,
                          R2_starts='{}.R2.contig_start_pos.txt',
                          remaining_R1s="{}.R1.contig_start_pos.no_pass.txt"):

        C = ContigAssembler()
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
            "associated_contigs": 0,
            'assembled_contigs': 0,
            'R1_singletons': 0,
            'R2_singletons': 0,
            }

        r1_name = self.__get_fasta_header__(rad1)
        r2_name = self.__get_fasta_header__(rad2)

        r1 = pysam.Fastafile(rad1)
        r2 = pysam.Fastafile(rad2)

        # MAKE OUTPUT FILENAMES
        RADnome_out = run_name + ".paired_contigs.fa"
        Singltons_R1_nome_out = run_name + ".singleton_R1_contigs.fa"
        Singltons_R2_nome_out = run_name + ".singleton_R2_contigs.fa"
        Contig_nome_out = run_name + ".assembled_contigs.fa"

        contig_positions_out = run_name + '.RADnome.contig_positions.txt'

        run_results['RADnome_out'] = RADnome_out

        if out_path is not None:

            RADnome_out = os.path.join(out_path, RADnome_out)
            RADnome_out = open(RADnome_out, 'w')

            Singltons_R1_nome_out = os.path.join(out_path, Singltons_R1_nome_out)
            Singltons_R1_nome_out = open(Singltons_R1_nome_out, 'w')

            Singltons_R2_nome_out = os.path.join(out_path, Singltons_R2_nome_out)
            Singltons_R2_nome_out = open(Singltons_R2_nome_out, 'w')

            Contig_nome_out = os.path.join(out_path, Contig_nome_out)
            Contig_nome_out = open(Contig_nome_out, 'w')

            contig_positions_out = os.path.join(out_path, contig_positions_out)
            contig_positions_out = open(contig_positions_out, 'w')

        else:

            RADnome_out = open(RADnome_out, 'w')
            Singltons_R1_nome_out = open(Singltons_R1_nome_out, 'w')
            Singltons_R2_nome_out = open(Singltons_R2_nome_out, 'w')
            Contig_nome_out = open(Contig_nome_out, 'w')
            contig_positions_out = open(contig_positions_out, 'w')

        c2c = pickle.load(open(contig_2_contig_dict, 'rb'))

        count = 0
        rad_frag_count = 0
        RADnome_pos = 0
        singletons_pos_R1 = 0
        contigs_pos = 0
        seq_start = 1
        r2_unassociated_contigs = []
        r2_assoicated_contigs = []

        contignome_current_position = 0
        radnome_current_position = 0
        r1_current_position = 0
        r2_current_position = 0

        R2_starts = R2_starts.format(run_name)
        R2_starts = tuple(int(i.strip()) for i in open(R2_starts,'rU'))


        # ITERATE OVER R1 CONTIG KEYS IN CONTIG-2-CONTIG DICTIONARY
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

                RADnome_out.write('>{}_NonOverlapping\n'.format(run_name))
                Contig_nome_out.write('>{}_Assembled\n'.format(run_name))
                Singltons_R1_nome_out.write('>{}_R1_Singletons\n'.format(run_name))
                Singltons_R2_nome_out.write('>{}_R2_Singletons\n'.format(run_name))

            # ASSOCIATE CONTIG PAIRS
            if self.filter_contigs(counts, proportion) >= proportion:

                run_results["rad_frag_count"] += 1
                run_results["mode_diff"] += self.filter_contigs(counts, proportion)

                # SETUP DNA FRAGMENTS
                bigNs = "N" * N_padding
                smallNs = "N" * insert_size
                contig1 = r1.fetch(r1_name, key, key+r1_contig_len)

                r2_start = self.__find_nearest__(R2_starts, mode[0])
                contig2 = r2.fetch(r2_name, r2_start, r2_start+r2_contig_len)

                r2_assoicated_contigs.append(contig2)

                if len(counts.keys()) > 1:

                    try:
                        non_mode_contigs = counts.keys()
                        non_mode_contigs.remove(mode[0])
                        non_mode_contigs = [self.__find_nearest__(R2_starts, s) for s in non_mode_contigs]
                        r2_unassociated_contigs.extend(non_mode_contigs)

                    except:
                        pass

                # RAINBOW orientation is reversed for the R1 reads
                # so it is necessary to flip contig1's oriention
                # ToDo: Make this generalized based on the SAM alignments
                contig1 = contig1[::-1]
                contig1 = self.__make_consensus__(contig1)

                assembly_results = C.assemble_contigs(contig1[::-1], contig2[::-1], overlap)

                if assembly_results is not None:

                    assembled_contig = assembly_results

                    # ADD FRAGMENTS TO CONTIG NOME
                    pos1 = self._append_to_pseudo_genome_(bigNs, Contig_nome_out, contigs_pos, span=80)
                    pos2 = self._append_to_pseudo_genome_(assembled_contig, Contig_nome_out, pos1, span=80)
                    contigs_pos = pos2

                    contignome_current_position = self.write_contig_start_stop_info(contig_positions_out, run_name, contignome_current_position, \
                                                                                (bigNs, assembled_contig), 'Contignome')
                    run_results["assembled_contigs"] += 1

                else:
                    # ADD FRAGMENTS TO RADNOME
                    r_pos1 = self._append_to_pseudo_genome_(bigNs, RADnome_out, RADnome_pos, span=80)
                    r_pos2 = self._append_to_pseudo_genome_(contig1, RADnome_out, r_pos1, span=80)
                    r_pos3 = self._append_to_pseudo_genome_(smallNs, RADnome_out, r_pos2, span=80)
                    r_pos4 = self._append_to_pseudo_genome_(contig2, RADnome_out, r_pos3, span=80)
                    RADnome_pos = r_pos4

                    radnome_current_position = self.write_contig_start_stop_info(contig_positions_out, run_name, radnome_current_position,
                                                                                (bigNs, contig1, smallNs, contig2), 'RADnome')

                    seq_start += r_pos4
                    run_results["associated_contigs"] += 1

            # PROCESS R1s
            else:

                r2_unassociated_contigs.extend(counts.keys()) # update R2 ids

                bigNs = "N" * N_padding
                contig1 = r1.fetch(r1_name, key, key+r1_contig_len)

                pos1 = self._append_to_pseudo_genome_(bigNs, Singltons_R1_nome_out, singletons_pos_R1, span=80)
                pos2 = self._append_to_pseudo_genome_(contig1, Singltons_R1_nome_out, pos1, span=80)

                r1_current_position = self.write_contig_start_stop_info(contig_positions_out, run_name, r1_current_position, (bigNs, contig1), 'R1nome')

                singletons_pos_R1 = pos2
                run_results["R1_singletons"] += 1

            count += 1

        # print 'un-associated R2s', len(r2_unassociated_contigs), len(set(r2_unassociated_contigs))
        # print 'associated R2s', len(r2_assoicated_contigs), len(set(r2_assoicated_contigs))

        # print 'set dif', len(list(set(r2_unassociated_contigs) - set(r2_assoicated_contigs)))

        # PROCESS REMAINING R1s
        #remaining_R1s = "{}.R1.contig_start_pos.no_pass.txt"
        for start in open(remaining_R1s.format(run_name),'r'):

            start = int(start)
            contig1 = r1.fetch(r1_name, start, start+r1_contig_len)
            pos1 = self._append_to_pseudo_genome_(bigNs, Singltons_R1_nome_out, singletons_pos_R1, span=80)
            pos2 = self._append_to_pseudo_genome_(contig1, Singltons_R1_nome_out, pos1, span=80)

            r1_current_position = self.write_contig_start_stop_info(contig_positions_out, run_name, r1_current_position, (bigNs, contig1), 'R1nome')

            singletons_pos_R1 = pos2
            run_results["R1_singletons"] += 1

        # PROCESS R2s
        singletons_pos_R2 = 0
        r2_unassociated_contigs = set(r2_unassociated_contigs)

        for c in r2_unassociated_contigs:
            pos1 = self._append_to_pseudo_genome_(bigNs, Singltons_R2_nome_out, singletons_pos_R2, span=80)
            pos2 = self._append_to_pseudo_genome_(contig1, Singltons_R2_nome_out, pos1, span=80)

            r2_current_position = self.write_contig_start_stop_info(contig_positions_out, run_name, r2_current_position, (bigNs, contig1), 'R2nome')

            singletons_pos_R2 = pos2

        run_results["R2_singletons"] = len(r2_unassociated_contigs)

        Singltons_R1_nome_out.write("\n")
        Singltons_R2_nome_out.write("\n")
        Contig_nome_out.write("\n")
        RADnome_out.write("\n")

        Singltons_R1_nome_out.close()
        Singltons_R2_nome_out.close()
        Contig_nome_out.close()
        RADnome_out.close()

        # Send data for logging
        if out_path == None:
            out_path = os.getcwd()

        self.__contigs_2_RADnome__(run_results, out_path)
        return 1


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

    @staticmethod
    def closest_starting_pos(positions, value, array_size):
        
        insert_point = bisect.bisect(positions, value)

        if insert_point >= array_size:
            return positions[insert_point-1]

        else:

            left, right = positions[insert_point-1], positions[insert_point]
            
            if abs(left - insert_point) > abs(right - insert_point):
                return right
            else:
                return left

    def associate_contigs(self, bam1, sam2, min_mapq, run_ID, R1_starts=None, R2_starts=None):
        # print bam1, sam2, min_mapq, run_ID
        # sys.exit()

        run_results = {
            "total_queries": 0,
            "passed_filter": 0,
            "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
            "min_mapq": min_mapq,
            "query_sam": bam1,
            "searched_sam": sam2,
            # "contig_positions": contig_positions,
            "pickle": None,
            "failed_filter":0,
            }

        contig_2_contig_dict = defaultdict(list)

        #  CREATE SAM INDEX IF IT DOESN'T EXIST
        if not os.path.exists(sam2 + FileIndex.ext):
            FileIndex.create(sam2, lambda fh: SamLine(fh).qname, allow_multiple=True)

        fi = FileIndex(sam2, SamLine, allow_multiple=True)
        ma = MergeAssemblies()

        R1_sam = pysam.Samfile(bam1,'rb')
        path = os.path.split(bam1)[0]

        if R1_starts == None:
            R1_starts = os.path.join(path, "{}.R1.contig_start_pos.txt".format(run_ID))

        if R2_starts == None:
            R2_starts = os.path.join(path, "{}.R2.contig_start_pos.txt".format(run_ID))
        
        R2_starts = [int(l.strip()) for l in open(R2_starts,'rU')]

        R2_starts_len = len(set(R2_starts))

        low_depth_R1_starts = os.path.join(path, "{}.R1.contig_start_pos.no_pass.txt".format(run_ID))
        low_depth_R1_starts = open(low_depth_R1_starts, 'w')

        all_R1s = set()

        # ITERATE OVER START POSITIONS
        for count, start in enumerate(open(R1_starts,'rU')):
            query_pos = int(start.strip())
            all_R1s.update([query_pos])

            if query_pos is 0:
                continue

            #print "{}.R1".format(run_ID), query_pos-5, query_pos+5
            # samtools view alignments/Trachs_Merged.1.fq.sorted.bam Trachs_Merged_stacks.R1:590-600
            reads = R1_sam.fetch("{}.R1".format(run_ID), query_pos-5, query_pos+5)
            #reads = [r for r in reads]                  # unpack interator

            for q in reads:
                q_id = q.qname[:-1] + "2"               # create query id (e.g., ends with "2")
                q_mapq = q.mapq                         # and, mapping quality.

                # SEARCH FOR QUERY AND PARSE RESULTS
                for s in fi[q_id]:

                    if (q_mapq > min_mapq) and (s.mapq > min_mapq):

                        hit_pos = ma.get_hit_pos(s)
                        #print hit_pos, ma.closest_starting_pos(R2_starts, hit_pos, R2_starts_len)
                        hit_pos = ma.closest_starting_pos(R2_starts, hit_pos, R2_starts_len)
                        #hit_pos = ma.round_bp_pos(hit_pos)
                        contig_2_contig_dict[query_pos].append(hit_pos)


        # print 'R1s total count', len(all_R1s)
        # print 'R1s associated:', len(contig_2_contig_dict.keys())
        # print 'R1s un-associated:', len(all_R1s.difference(set(contig_2_contig_dict.keys())))
        
        # print 'R2s total count', R2_starts_len 
        # print 'R2s associated:', len(set([item for sublist in contig_2_contig_dict.values() for item in sublist]))


        run_results["passed_filter"] = len(contig_2_contig_dict.keys())

        unassociated_R1s = all_R1s.difference(set(contig_2_contig_dict.keys()))
        run_results["failed_filter"] = len(unassociated_R1s)

        [low_depth_R1_starts.write("{}\n".format(p)) for p in unassociated_R1s]

        #CREATE PICKLE OUTPUT FILE
        today = datetime.date.today()

        path = os.path.split(bam1)[0]
        pkl_output_file_name = os.path.join(path, '{}.R1_to_R2_contig_associations.pkl'.format(run_ID))
        pkl_out = open(pkl_output_file_name, 'wb')
        pickle.dump(contig_2_contig_dict, pkl_out)

        low_depth_R1_starts.close()

        # UPDATE RUN RESULTS AND GENERATE LOGFILE
        run_results["total_queries"] = count
        run_results['pickle'] = pkl_output_file_name
        self.__associate_contigs_log__(run_results, path)        # format logging results
        return 1

    def __generate_test_data__(self, bam1, sam2, contig_positions, min_mapq):

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

        fq1_out = open("test_fq1.10k.fq", 'w')
        fq2_out = open("test_fq2.10k.fq", 'w')

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


class ContigAssembler(object):
    """Code for assembling contigs."""
    def __init__(self):
        super(ContigAssembler, self).__init__()

    def make_pairs(self, a, b):
        """Make Pairs"""
        a = [i for i in a]
        b = [i for i in b]
        return zip(a, b)

    def is_match(self, pair):
        """Test if pair is match."""

        if pair[0] == pair[1]:
            return 1
        else:
            return 0

    def get_max_overlap(self, mismatch_list):
        """Identify the maximum overlap that has zero mismatches."""

        pos = None # Ensures that 'None' is returned if there is no valid overlap
        for count, i in enumerate(mismatch_list):
            if i is 0:
                pos = count
        return pos

    def assemble_contigs(self, z, y, overlap=None):
        """ Do the assembly. Only merge contigs if the overlap produces zero mismatches."""

        mismatch_list = []
        chars = len(z) + 1

        for i in range(1, chars):
            a, b = z[(-1 * i):], y[:i]

            algn = map(self.is_match, self.make_pairs(a, b))

            padding = " "*((chars - 1) - i)
            #print algn.count(0),"{}{}".format(padding, y) # print statements show what is going on.

            match_str = ''.join(map(str, algn))
            #print algn.count(0), "{}{}".format(padding, match_str)

            mismatch_list.append(algn.count(0))

        # identify position of max overlap
        mpos = self.get_max_overlap(mismatch_list)

        # merge reads.
        if mpos is not None:

            assembled_contigs = z + y[mpos+1:]
            overlapped = mpos

            if overlapped >= overlap:
                return assembled_contigs
            else:
                return None
        else:
            return mpos


if __name__ == '__main__':
    G = GenerateRADnome()

    tsv_fin ="/Users/ngcrawford/Dropbox/JMP_CAS/Nick/Stacks_output/batch_1.catalog.tags_R2.tsv"
    pseudo_genome_fout = 'stacksnome.fa'
    Ns = 500
    run_name = 'test_run'

    G.make_READnome_from_Stacks_tsv(tsv_fin, pseudo_genome_fout, Ns, run_name)


    pass

