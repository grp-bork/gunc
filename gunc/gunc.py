#!/usr/bin/env python3
import os
import sys
import argparse
from . import gunc_database
from . import external_tools
from . import visualisation as vis
from .get_scores import chim_score
from ._version import get_versions


def parse_args(args):
    """Parse Arguments

    Arguments:
        args {List} -- List of args supplied to script.

    Returns:
        {Namespace} -- assigned args

    """
    description = ('Tool for detection of chimerism and '
                   'contamination in prokaryotic genomes.\n')
    parser = argparse.ArgumentParser(description=description,)
    subparsers = parser.add_subparsers(title='GUNC subcommands',
                                       metavar='',
                                       dest='cmd')
    run = subparsers.add_parser('run',
                                help='Run chimerism detection.')
    run_group = run.add_mutually_exclusive_group(required=True)
    download_db = subparsers.add_parser('download_db',
                                        help='Download GUNC db.')
    vis = subparsers.add_parser('plot',
                                help='Create interactive visualisation.',
                                formatter_class=lambda prog:
                                    argparse.ArgumentDefaultsHelpFormatter(prog,
                                                                           max_help_position=100))

    run.add_argument('-d', '--db_file',
                     help='Diamond database reference file.',
                     default=os.environ.get('GUNC_DB'),
                     metavar='')
    run_group.add_argument('-i', '--input_file',
                           help='Input file in FASTA fna format.',
                           metavar='')
    run_group.add_argument('-g', '--gene_calls',
                           help='Input genecalls FASTA faa format.',
                           metavar='')
    run.add_argument('-p', '--threads',
                     help='number of CPU threads. Default: 4',
                     default='4',
                     metavar='')
    run.add_argument('-t', '--temp_dir',
                     help='Directory to store temp files. Default: cwd',
                     default=os.getcwd(),
                     metavar='')
    run.add_argument('-o', '--out_dir',
                     help='Output dir.  Default: cwd',
                     default=os.getcwd(),
                     metavar='')
    run.add_argument('-s', '--sensitive',
                     help='Run with high sensitivity',
                     action='store_true',
                     default=False)
    vis.add_argument('-d', '--diamond_file',
                     help='GUNC diamond outputfile.',
                     required=True,
                     metavar='')
    vis.add_argument('-g', '--gunc_file',
                     help='GUNC scores outputfile.',
                     metavar='')
    vis.add_argument('-o', '--out_dir',
                     help='Output directory.',
                     default=os.getcwd(),
                     metavar='')
    vis.add_argument('-t', '--tax_levels',
                     help='Tax levels to display (comma-seperated).',
                     default='kingdom,phylum,family,genus,contig',
                     metavar='')
    vis.add_argument('-r', '--remove_minor_clade_level',
                     help='Tax levels at which to remove minor clades.',
                     default='kingdom',
                     metavar='')
    vis.add_argument('-c', '--contig_display_num',
                     help='Numper of contigs to visualise.',
                     default=1000,
                     type=int,
                     metavar='')
    download_db.add_argument('path',
                             help='Download database to given direcory.',
                             metavar='dest_path')
    parser.add_argument('-v', '--version',
                        help='Print version number and exit.',
                        action='version',
                        version=get_versions()['version'])
    if not args:
        parser.print_help(sys.stderr)
        sys.exit()
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


def start_checks():
    """Checks if tool dependencies are available."""
    if not external_tools.check_if_tool_exists('diamond'):
        sys.exit('[ERROR] Diamond not found..')
    else:
        diamond_ver = external_tools.check_diamond_version()
        if diamond_ver != '2.0.4':
            sys.exit(f'[ERROR] Diamond version is {diamond_ver}, not 2.0.4')
    if not external_tools.check_if_tool_exists('prodigal'):
        sys.exit('[ERROR] Prodigal not found..')
    if not external_tools.check_if_tool_exists('zgrep'):
        sys.exit('[ERROR] zgrep not found..')


def get_genecount_from_gunc_output(gunc_file):
    """Extract gene count from GUNC output file."""
    with open(gunc_file, 'r') as f:
        next(f)
        return f.readline().split('\t')[1]


def main():
    args = parse_args(sys.argv[1:])
    if args.cmd == 'download_db':
        gunc_database.get_db(args.path)

    if args.cmd == 'run':
        start_checks()

        if not args.db_file:
            sys.exit('[WARNING] database_file argument missing.')

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
                               args.db_file,
                               diamond_outfile)
        gene_count = external_tools.get_record_count_in_fasta(diamond_inputfile)

        print('[INFO] Calculating scores for each tax-level:')
        df = chim_score(diamond_outfile, gene_count, args.sensitive)
        df.to_csv(f'{diamond_outfile}.chimerism_scores', index=False, sep='\t')

    if args.cmd == 'plot':
        if not args.gunc_file:
            gunc_file = f'{os.path.abspath(args.diamond_file)}.chimerism_scores'
            if not os.path.isfile(gunc_file):
                sys.exit('[ERROR] GUNC file not found!')
        else:
            gunc_file = args.gunc_file
        gene_count = get_genecount_from_gunc_output(gunc_file)
        genome_name = os.path.basename(args.diamond_file).replace('.diamond.out', '')
        viz_html = vis.create_viz_from_diamond_file(args.diamond_file,
                                                    gene_count,
                                                    args.tax_levels,
                                                    args.contig_display_num,
                                                    args.remove_minor_clade_level)
        viz_html_path = os.path.join(args.out_dir, f'{genome_name}.viz.html')
        with open(viz_html_path, 'w') as f:
            f.write(viz_html)


if __name__ == "__main__":
    main()
