#!/usr/bin/env python3
import os
import sys
import glob
import json
import argparse
from . import gunc_database
from . import external_tools
from . import visualisation as vis
from .get_scores import chim_score
from ._version import get_versions
from .external_tools import get_record_count_in_fasta as record_count


def parse_args(args):
    """Parse Arguments

    Arguments:
        args {List} -- List of args supplied to script.

    Returns:
        {Namespace} -- assigned args

    """
    description = ('Tool for detection of chimerism and '
                   'contamination in prokaryotic genomes.\n')
    parser = argparse.ArgumentParser(description=description)
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
                     help='DiamondDB reference file. Default: GUNC_DB envvar',
                     default=os.environ.get('GUNC_DB'),
                     metavar='')
    run_group.add_argument('-i', '--input_file',
                           help='Input file in FASTA fna format.',
                           metavar='')
    run_group.add_argument('-f', '--input_dir',
                           help='Input dir with files in FASTA fna format.',
                           metavar='')
    run.add_argument('-e', '--file_suffix',
                     help='Suffix of files in input_dir. Default: .fa',
                     default='.fa',
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
                     help='Run with high sensitivity. Default: False',
                     action='store_true',
                     default=False)
    vis.add_argument('-d', '--diamond_file',
                     help='GUNC diamond outputfile.',
                     required=True,
                     metavar='')
    vis.add_argument('-g', '--gunc_gene_count_file',
                     help='GUNC gene_counts.json file.',
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
                     help='Tax level at which to remove minor clades.',
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
    args = parser.parse_args(args)
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
        sys.exit('[ERROR] Diamond 2.0.4 not found..')
    else:
        diamond_ver = external_tools.check_diamond_version()
        if diamond_ver != '2.0.4':
            sys.exit(f'[ERROR] Diamond version is {diamond_ver}, not 2.0.4')
    if not external_tools.check_if_tool_exists('prodigal'):
        sys.exit('[ERROR] Prodigal not found..')
    if not external_tools.check_if_tool_exists('zgrep'):
        sys.exit('[ERROR] zgrep not found..')


def get_genecount_from_gunc_output(gene_counts_file, basename):
    """Extract gene count from GUNC gene_counts.json file."""
    with open(gene_counts_file) as f:
        data = json.load(f)
    return int(data[basename])


def get_fna_files_in_dir(directory, suffix):
    """Get files in directory that end in suffix."""
    return glob.glob(os.path.join(directory, f'*{suffix}'))


def merge_genecalls(genecall_files, out_dir):
    merged_outfile = os.path.join(out_dir, 'merged.genecalls.faa')
    with open(merged_outfile, 'w') as ofile:
        for file in genecall_files:
            with open(file, 'r') as infile:
                genome_name = os.path.basename(file).replace('.genecalls.faa',
                                                             '')
                for line in infile:
                    if line.startswith('>'):
                        contig_name = line.split(' ')[0]
                        line = f'{contig_name}_-_{genome_name}\n'
                    ofile.write(line)
    return merged_outfile


def split_diamond_output(diamond_outfile, out_dir):
    outfiles = []
    output = {}
    with open(diamond_outfile, 'r') as f:
        for line in f:
            genome_name = line.split('\t')[0].split('_-_')[1]
            line = line.replace(f'_-_{genome_name}', '')
            output[genome_name] = output.get(genome_name, '') + line
    for genome_name in output:
        outfile = os.path.join(out_dir, f'{genome_name}.diamond.out')
        outfiles.append(outfile)
        with open(outfile, 'w') as ofile:
            ofile.write(output[genome_name])
    return outfiles


def write_json(data, filename):
    with open(filename, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=4)


def main():
    args = parse_args(sys.argv[1:])
    if args.cmd == 'download_db':
        gunc_database.get_db(args.path)

    if args.cmd == 'run':
        start_checks()

        if not args.db_file:
            sys.exit('[WARNING] database_file argument missing.')

        genes_called = {}

        if args.input_dir:
            fnas = get_fna_files_in_dir(args.input_dir, args.file_suffix)
            if len(fnas) == 0:
                sys.exit(f'[ERROR] No fastas found in {args.input_dir} '
                         f'with ending {args.file_suffix}')
        elif args.input_file:
            fnas = [args.input_file]

        if args.gene_calls:
            input_basename = os.path.basename(args.gene_calls)
            diamond_inputfile = args.gene_calls
            genes_called[input_basename] = record_count(diamond_inputfile)
            diamond_outfile = os.path.join(args.out_dir,
                                           f'{input_basename}.diamond.out')
        else:
            gene_calls_out_dir = os.path.join(args.out_dir, 'gene_calls')
            create_dir(gene_calls_out_dir)
            genecall_files = []
            for i, fna in enumerate(fnas):
                print(f'[INFO] Running Prodigal {i}/{len(fnas)}', flush=True)
                input_basename = os.path.basename(fna)
                prodigal_outfile = os.path.join(gene_calls_out_dir,
                                                f'{input_basename}.genecalls.faa')
                external_tools.prodigal(fna, prodigal_outfile)
                genes_called[input_basename] = record_count(prodigal_outfile)
                genecall_files.append(prodigal_outfile)
            diamond_inputfile = merge_genecalls(genecall_files,
                                                gene_calls_out_dir)
            diamond_outfile = os.path.join(args.out_dir,
                                           f'merged.diamond.out')
            genes_called_outfile = os.path.join(gene_calls_out_dir,
                                                f'gene_counts.json')
            write_json(genes_called, genes_called_outfile)

        external_tools.diamond(diamond_inputfile,
                               args.threads,
                               args.temp_dir,
                               args.db_file,
                               diamond_outfile)

        diamond_out_dir = os.path.join(args.out_dir, 'diamond_output')
        create_dir(diamond_out_dir)
        if args.input_dir:
            diamond_outfiles = split_diamond_output(diamond_outfile,
                                                    diamond_out_dir)
            os.remove(diamond_outfile)
            os.remove(diamond_inputfile)
        else:
            diamond_outfiles = [diamond_outfile]

        for diamond_outfile in diamond_outfiles:
            basename = os.path.basename(diamond_outfile).replace('.diamond.out',
                                                                 '')
            print(f'[INFO] Calculating GUNC scores for {basename}:')
            gunc_out_dir = os.path.join(args.out_dir,
                                        'gunc_output')
            gunc_out_file = os.path.join(gunc_out_dir,
                                         f'{basename}.chimerism_scores')
            create_dir(gunc_out_dir)
            gene_call_count = genes_called[basename]
            df = chim_score(diamond_outfile, gene_call_count, args.sensitive)
            df.to_csv(gunc_out_file,
                      index=False,
                      sep='\t')

    if args.cmd == 'plot':
        if not args.gunc_gene_count_file:
            diamond_file_path = os.path.abspath(args.diamond_file)
            gunc_dir = os.path.dirname(os.path.dirname(diamond_file_path))
            gene_counts_file = os.path.join(gunc_dir,
                                            "gene_calls/gene_counts.json")
            if not os.path.isfile(gene_counts_file):
                sys.exit('[ERROR] GUNC gene_counts.json file not found!')
        else:
            gene_counts_file = args.gunc_gene_count_file
        basename = os.path.basename(args.diamond_file).replace('.diamond.out',
                                                               '')
        genes_called = get_genecount_from_gunc_output(gene_counts_file,
                                                      basename)
        print(genes_called)
        viz_html = vis.create_viz_from_diamond_file(args.diamond_file,
                                                    genes_called,
                                                    args.tax_levels,
                                                    args.contig_display_num,
                                                    args.remove_minor_clade_level)
        viz_html_path = os.path.join(args.out_dir, f'{basename}.viz.html')
        with open(viz_html_path, 'w') as f:
            f.write(viz_html)


if __name__ == "__main__":
    main()
