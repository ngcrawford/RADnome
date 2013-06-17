#!/usr/bin/env python
# encoding: utf-8


"""
File: RADnome.py
Author: Nicholas Crawford

Created by Nicholas Crawford on Mon Jun 17 15:35:51 PDT 2013
Copyright (c) 2012 Nicholas G. Crawford All rights reserved.

Description...

Add later.

"""

import sys
import argparse
import GenerateRADnome


def generate_READnome(args):
    G = GenerateRADnome.GenerateRADnome()
    G.make_pseudo_genome(args.infile, args.outfile, args.N_padding, args.run_name)
    pass


def generate_RADnome(args):
    G = GenerateRADnome.GenerateRADnome()
    G.contigs_2_RADnome(args.readnomes[0], args.readnomes[1], args.contig_lengths[0],
                        args.contig_lengths[1], args.outfile, args.pickle_dict,
                        args.run_name, args.N_padding, args.insert_size)
    pass

def assoicate_contigs(args):
    M = GenerateRADnome.MergeAssemblies()
    M.associate_contigs(args.sams[0], args.sams[1], args.pos_file)


def get_args():
    """Parse sys.argv"""

    parser = argparse.ArgumentParser(prog='RADnome')

    # SETUP SUB-PARSERS
    subparsers = parser.add_subparsers(title='subcommands',
                                       help='sub-command help')

    parser_readnome = subparsers.add_parser(
        'READnome',
        help='Create pseudo-genome from single end reads assembled \
              with Rainbow (Chong 2012).')

    parser_ascContigs = subparsers.add_parser(
        'ascContigs',
        help='Match up contigs based on read mapping.')

    parser_RADnome = subparsers.add_parser(
        'RADnome',
        help='Create pseudo-genome for associated contigs.')


    # SETUP READnome SUB-COMMAND
    parser_readnome.add_argument(
        '-n', '--N-padding',
        type=int,
        default=500,
        help="Number of N's to insert between each cluster/contig.")

    parser_readnome.add_argument(
        '-r', '--run-name',
        type=str,
        default='RADnome',
        help="Name of run. Becomes fasta name.")

    parser_readnome.add_argument(
        'infile',
        nargs='?',
        default=sys.stdin,
        help='Path to input. (default is STDIN)')

    parser_readnome.add_argument(
        'outfile',
        nargs='?',
        default=sys.stdout,
        help='Path to output. (default is STOUT)')

    # SETUP RADnome SUB-COMMAND
    parser_RADnome.add_argument(
        '-n', '--N-padding',
        type=int,
        default=500,
        help="Number of N's to insert between each cluster/contig.")

    parser_RADnome.add_argument(
        '-i', '--insert-size',
        type=int,
        default=50,
        help="Size of insert when libraries were generated. Used to deterimine\
              how many N's to insert between associated contigs.")

    parser_RADnome.add_argument(
        '-r', '--run-name',
        type=str,
        default='RADnome',
        help="Name of run. Becomes fasta name.")

    parser_RADnome.add_argument(
        '--readnomes',
        type=str,
        nargs=2,
        help="Forward and Reverse READnome files. Note that order matters.")

    parser_RADnome.add_argument(
        '--pickle-dict',
        type=str,
        default=None,
        help="'Pickled' dictionary from sam_index run.")

    parser_RADnome.add_argument(
        '-cl','--contig-lengths',
        type=int,
        nargs=2,
        help="Length of contigs in forward and reverse files. Note that order matters.")

    parser_RADnome.add_argument(
        'outfile',
        nargs='?',
        default=sys.stdout,
        help='Path to output. (default is STOUT)')

    # SETUP ascContigs SUB-COMMAND
    parser_ascContigs.add_argument(
        '-s', '--sams',
        type=str,
        nargs=2,
        default='Path to file containing start positions.',
        help="Forward and Reverse SAM files. Note that order matters.")

    parser_ascContigs.add_argument(
        '-p', '--pos-file',
        type=str,
        help="Forward and Reverse SAM files. Note that order matters.")

    # 
    parser_readnome.set_defaults(func=generate_READnome)
    parser_ascContigs.set_defaults(func=assoicate_contigs)
    parser_RADnome.set_defaults(func=generate_RADnome)

    args = parser.parse_args()
    args.func(args)
    return args

print get_args()
