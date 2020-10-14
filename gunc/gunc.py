#!/usr/bin/env python3
import os
import sys
import argparse
from . import external_tools
from .get_scores import chim_score
from ._version import get_versions


def parse_args(args):
    """Parse Arguments

    Arguments:
        args {List} -- List of args supplied to script.

    Returns:
        {Namespace} -- assigned args

    """
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)

    parser.add_argument('-d', '--database_file',
                        help='Diamond database reference file.',
                        required=True,
                        metavar='')
    group.add_argument('-i', '--input_file',
                       help='Input file in FASTA fna format..',
                       metavar='')
    group.add_argument('-g', '--gene_calls',
                       help='Input genecalls FASTA faa format..',
                       metavar='')
    parser.add_argument('-p', '--threads',
                        help='number of CPU threads.',
                        default='4',
                        metavar='')
    parser.add_argument('-t', '--temp_dir',
                        help='Directory to store temporary files.',
                        default=os.getcwd(),
                        metavar='')
    parser.add_argument('-o', '--out_dir',
                        help='Output dir.',
                        default=os.getcwd(),
                        metavar='')
    parser.add_argument('-s', '--sensitive',
                        help='Run with high sensitivity',
                        action='store_true',
                        default=False)
    parser.add_argument('-v', '--version',
                        help='Print version number and exit.',
                        action='version',
                        version=get_versions()['version'])
    args = parser.parse_args()
    return args


def create_dir(path):
    """Create a directory

    Will create a directory if it doesnt already exist.

    Arguments:
        path {str} -- directory path
    """
    if not os.path.exists(path):
        os.makedirs(path)


def main():
    args = parse_args(sys.argv[1:])

    if args.input_file:
        input_basename = os.path.basename(args.input_file)
        prodigal_outfile = os.path.join(args.out_dir,
                                        f'{input_basename}.genecalls.faa')
        external_tools.prodigal(args.input_file, prodigal_outfile)
        diamond_inputfile = prodigal_outfile
    else:
        input_basename = os.path.basename(args.gene_calls)
        diamond_inputfile = args.gene_calls

    diamond_outfile = os.path.join(args.out_dir,
                                   f'{input_basename}.diamond.out')
    external_tools.diamond(diamond_inputfile,
                           args.threads,
                           args.temp_dir,
                           args.database_file,
                           diamond_outfile)
    gene_count = external_tools.get_record_count_in_fasta(diamond_inputfile)

    print('[INFO] Calculating scores for each tax-level:')
    df = chim_score(diamond_outfile, gene_count, args.sensitive)
    df.to_csv(f'{diamond_outfile}.chimerism_scores', index=False, sep='\t')


if __name__ == "__main__":
    main()
