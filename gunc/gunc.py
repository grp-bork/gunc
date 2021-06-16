#!/usr/bin/env python3
import os
import sys
import gzip
import glob
import json
import argparse
import pandas as pd
import multiprocessing
from . import checkm_merge
from . import gunc_database
from . import external_tools
from datetime import datetime
from . import visualisation as vis
from .get_scores import chim_score
from ._version import get_versions
from .external_tools import get_record_count_in_fasta as record_count


def parse_args(args):
    """Parse Arguments

    Arguments:
        args (List): List of args supplied to script.

    Returns:
        Namespace: assigned args

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
    merge_checkm = subparsers.add_parser('merge_checkm',
                                         help='Merge GUNC and CheckM outputs.')
    vis = subparsers.add_parser('plot',
                                help='Create interactive visualisation.')

    run.add_argument('-r', '--db_file',
                     help='DiamondDB reference file. Default: GUNC_DB envvar',
                     default=os.environ.get('GUNC_DB'),
                     metavar='')
    run_group.add_argument('-i', '--input_fasta',
                           help='Input file in FASTA format.',
                           metavar='')
    run_group.add_argument('-f', '--input_file',
                           help='File with paths to FASTA format files.',
                           metavar='')
    run_group.add_argument('-d', '--input_dir',
                           help='Input dir with files in FASTA format.',
                           metavar='')
    run.add_argument('-e', '--file_suffix',
                     help='Suffix of files in input_dir. Default: .fa',
                     default='.fa',
                     metavar='')
    run.add_argument('-g', '--gene_calls',
                     help='Input files are FASTA faa format. Default: False',
                     action='store_true',
                     default=False)
    run.add_argument('-t', '--threads',
                     help='number of CPU threads. Default: 4',
                     default='4',
                     metavar='')
    run.add_argument('-o', '--out_dir',
                     help='Output dir.  Default: cwd',
                     default=os.getcwd(),
                     metavar='')
    run.add_argument('--temp_dir',
                     help='Directory to store temp files. Default: cwd',
                     default=os.getcwd(),
                     metavar='')
    run.add_argument('--sensitive',
                     help='Run with high sensitivity. Default: False',
                     action='store_true',
                     default=False)
    run.add_argument('--detailed_output',
                     help='Output scores for every taxlevel. Default: False',
                     action='store_true',
                     default=False)
    run.add_argument('--use_species_level',
                     help=('Allow species level to be picked as maxCSS. '
                           'Default: False'),
                     action='store_true',
                     default=False)
    run.add_argument('--min_mapped_genes',
                     help=('Dont calculate GUNC score if number of mapped '
                           'genes is below this value. Default: 11'),
                     type=int,
                     default=11,
                     metavar='')
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
                     help='Number of contigs to visualise.',
                     default=1000,
                     type=int,
                     metavar='')
    download_db.add_argument('path',
                             help='Download database to given direcory.',
                             metavar='dest_path')
    download_db.add_argument('-db', '--database',
                             help='Which db to download. progenomes,gtdb Default: progenomes',
                             default='progenomes',
                             metavar='')
    merge_checkm.add_argument('-g', '--gunc_file',
                              help='Path of GUNC.maxCSS_level.tsv file.',
                              required=True,
                              metavar='')
    merge_checkm.add_argument('-c', '--checkm_file',
                              help='CheckM qa output file',
                              required=True,
                              metavar='')
    merge_checkm.add_argument('-o', '--out_dir',
                              help='Output directory for merged file',
                              default=os.getcwd(),
                              metavar='')
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
        path (str): directory path
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


def get_files_in_dir_with_suffix(directory, suffix):
    """Get files in directory that end in suffix."""
    files = glob.glob(os.path.join(directory, f'*{suffix}'))
    if len(files) == 0:
        sys.exit(f'[ERROR] No files found in {directory} with ending {suffix}')
    return files


def openfile(filename):
    """Open file whether gzipped or not."""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'r')
    else:
        return open(filename, 'r')


def merge_genecalls(genecall_files, out_dir, file_suffix):
    """Merge genecall files.

    Merges fastas together to run diamond more efficiently.
    Adds the name of the file to each record (delimiter '_-_')
    so they can be separated after diamond mapping.

    Arguments:
        genecall_files (list): Paths of genecall fastas to merge
        out_dir (str): Directory to put the merged file
        file_suffix (str): Suffix of input files

    Returns:
        str: path of the merged file
    """
    merged_outfile = os.path.join(out_dir, 'merged.genecalls.faa')
    with open(merged_outfile, 'w') as ofile:
        for file in genecall_files:
            if os.path.isfile(file):
                with openfile(file) as infile:
                    genome_name = os.path.basename(file).replace('.genecalls.faa',
                                                                 '')
                    if genome_name.endswith(file_suffix):
                        genome_name = genome_name[:-len(file_suffix)]
                    for line in infile:
                        if line.startswith('>'):
                            contig_name = line.split(' ')[0]
                            line = f'{contig_name}_-_{genome_name}\n'
                        ofile.write(f'{line.strip()}\n')
    return merged_outfile


def split_diamond_output(diamond_outfile, out_dir, db):
    """Split diamond output into per-sample files.

    Separate diamond output file into the constituent sample files.
    This uses the identifiers that were added by :func:`merge_genecalls`

    Arguments:
        diamond_outfile (str): path to the diamond file to be split
        out_dir (str): Directory in which to put the split files.
        db (str): Which db is used: progenomes or gtdb

    Returns:
        list: Of the split file paths
    """
    outfiles = []
    output = {}
    with open(diamond_outfile, 'r') as f:
        for line in f:
            genome_name = line.split('\t')[0].split('_-_')[1]
            line = line.replace(f'_-_{genome_name}', '')
            output[genome_name] = output.get(genome_name, '') + line
    for genome_name in output:
        outfile = os.path.join(out_dir, f'{genome_name}.diamond.{db}.out')
        outfiles.append(outfile)
        with open(outfile, 'w') as ofile:
            ofile.write(output[genome_name])
    return outfiles


def write_json(data, filename):
    """Write json data to filename."""
    with open(filename, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=4)


def get_paths_from_file(input_file):
    """Extract paths from a text file."""
    if input_file.endswith('.gz'):
        sys.exit('[ERROR] Input should be list of filepaths: use -i instead?.')
    with open(input_file, 'r') as f:
        paths = f.readlines()
    if paths[0].startswith('>'):
        sys.exit('[ERROR] Input should be list of filepaths: use -i instead?.')
    return [path.strip() for path in paths]


def run_from_gene_calls(faas, out_dir, file_suffix):
    """Prepare genecalls for running diamond."""
    print(f'[START] {datetime.now().strftime("%H:%M:%S")} Merging genecalls..',
          flush=True)
    genes_called = {}
    for faa in faas:
        basename = os.path.basename(faa).split(file_suffix)[0]
        if basename.endswith('.genecalls'):
            basename = basename.split('.genecalls')[0]
        genes_called[basename] = record_count(faa)
    diamond_inputfile = merge_genecalls(faas, out_dir, file_suffix)
    print(f'[END]   {datetime.now().strftime("%H:%M:%S")} Finished Merging..',
          flush=True)
    return genes_called, diamond_inputfile


def run_from_fnas(fnas, out_dir, file_suffix, threads):
    """Call genes and prepare diamond input.

    Args:
        fnas (list): Of fasta files
        out_dir (str): Directory in which to put genecalls
        file_suffix (str): Suffix of input files
        threads (int): Numper of parallel prodigal instances to run

    Returns:
        tuple:
            - genes_called (dict): {input_file: gene count}
            - diamond_inputfile (str): merged input file for diamond

    """
    print(f'[START] {datetime.now().strftime("%H:%M:%S")} Running Prodigal..',
          flush=True)
    create_dir(out_dir)
    genecall_files = []
    genes_called = {}
    prodigal_info = []
    for fna in fnas:
        basename = os.path.basename(fna).split(file_suffix)[0]
        prodigal_outfile = os.path.join(out_dir, f'{basename}.genecalls.faa')
        prodigal_info.append((fna, prodigal_outfile))
    p = multiprocessing.Pool(int(threads))
    p.map(run_prodigal, prodigal_info)
    for fna, prodigal_outfile in prodigal_info:
        basename = os.path.basename(fna).split(file_suffix)[0]
        genes_called[basename] = record_count(prodigal_outfile)
        genecall_files.append(prodigal_outfile)
    diamond_inputfile = merge_genecalls(genecall_files, out_dir, file_suffix)
    print(f'[END]   {datetime.now().strftime("%H:%M:%S")} Finished Prodigal..',
          flush=True)
    return genes_called, diamond_inputfile


def run_prodigal(prodigal_info):
    fna, prodigal_outfile = prodigal_info
    external_tools.prodigal(fna, prodigal_outfile)


def run_diamond(infile, threads, temp_dir, db_file, out_dir, db):
    """Run diamond and split ouput.

    Runs diamond on infile and if needed splits the constitiuent samples.

    Args:
        infile (str): Path to gene calls fasta to run diamond on.
        threads (int): Number of threads to use for diamond.
        temp_dir (str): Path of tempdir for diamond running.
        db_file (str): Path to diamond database file (GUNC_DB).
        out_dir (str): Path of directory in which to put the output files.
        db (str): Which db is used: progenomes or gtdb

    Returns:
        list: Of diamond output files.
    """
    print(f'[START] {datetime.now().strftime("%H:%M:%S")} Running Diamond..',
          flush=True)
    outfile = os.path.join(out_dir, f'{os.path.basename(infile)}.diamond.{db}.out')
    external_tools.diamond(infile, threads, temp_dir, db_file, outfile)

    if infile.endswith('merged.genecalls.faa'):
        out_dir = os.path.join(out_dir, 'diamond_output')
        create_dir(out_dir)
        diamond_outfiles = split_diamond_output(outfile, out_dir, db)
        os.remove(outfile)
        os.remove(infile)
    else:
        diamond_outfiles = [outfile]
    print(f'[END]   {datetime.now().strftime("%H:%M:%S")} Finished Diamond..',
          flush=True)
    return diamond_outfiles


def run_gunc(diamond_outfiles, genes_called, out_dir, sensitive,
             detailed_output, db, min_mapped_genes, use_species_level):
    """Call GUNC scores on diamond output files.

    Outputs a pandas.DataFrame with one line per inputfile
    (taxlevel with highest CSS score).
    If detailed_output = True, files with all taxlevels are written.

    Args:
        diamond_outfiles (list): Paths of diamond output files.
        genes_called (dict): filename: genecount.
        out_dir (str): Path to output directory.
        sensitive (bool): Run with high sensitivity.
        detailed_output (bool): Output scores for every taxlevel.
        db (str): Which db to use: progenomes or gtdb
        min_mapped_genes (int): Minimum number of mapped genes
                                at which to calculate scores
        use_species_level (bool): Allow species level to be picked as maxCSS

    Returns:
        pandas.DataFrame: One line per inputfile Gunc scores
    """
    print(f'[START] {datetime.now().strftime("%H:%M:%S")} Running scoring..',
          flush=True)
    gunc_output = []
    for diamond_file in diamond_outfiles:
        basename = os.path.basename(diamond_file).split('.diamond.')[0]
        gene_call_count = genes_called[basename]
        detailed, single = chim_score(diamond_file, gene_call_count,
                                      sensitive, min_mapped_genes,
                                      use_species_level, db)
        if detailed_output:
            detailed_gunc_out_dir = os.path.join(out_dir, 'gunc_output')
            detailed_gunc_out_file = os.path.join(detailed_gunc_out_dir,
                                                  f'{basename}.{db}.all_levels.tsv')
            create_dir(detailed_gunc_out_dir)
            detailed.to_csv(detailed_gunc_out_file,
                            index=False,
                            sep='\t',
                            na_rep='nan')
        gunc_output.append(single)
    print(f'[END]   {datetime.now().strftime("%H:%M:%S")} Finished scoring..',
          flush=True)
    return pd.concat(gunc_output).sort_values('genome')


def check_for_duplicate_filenames(fnas, file_suffix):
    """Check if there are duplicate filenames in fna list.

    As the output files are based on the names of the input files, if
    there are duplicate filenames the output files wil be overwritten.

    Args:
        fnas (list): fna file paths
        file_suffix (str): File suffix
    """
    names = [os.path.basename(fna).split(file_suffix)[0] for fna in fnas]
    duplicated_items = set([x for x in names if names.count(x) > 1])
    if len(duplicated_items) > 0:
        sys.exit(f'Filenames appear more than once: {duplicated_items}')


def remove_missing_fnas(fnas):
    existing_fnas = []
    for fna in fnas:
        if os.path.isfile(fna):
            existing_fnas.append(fna)
        else:
            print(f'[WARNING] {fna} not found')
    return existing_fnas


def run(args):
    """Run entire GUNC workflow."""
    if not os.path.isdir(args.out_dir):
        sys.exit(f'[ERROR] Output Directory {args.out_dir} doesnt exist.')

    if args.input_dir:
        fastas = get_files_in_dir_with_suffix(args.input_dir, args.file_suffix)
    elif args.input_file:
        fastas = get_paths_from_file(args.input_file)
        fastas = remove_missing_fnas(fastas)
    elif args.input_fasta:
        fastas = [args.input_fasta]

    check_for_duplicate_filenames(fastas, args.file_suffix)

    if args.gene_calls:
        gene_calls_out_dir = args.out_dir
        genes_called, diamond_input = run_from_gene_calls(fastas,
                                                          gene_calls_out_dir,
                                                          args.file_suffix)
    else:
        gene_calls_out_dir = os.path.join(args.out_dir, 'gene_calls')
        genes_called, diamond_input = run_from_fnas(fastas,
                                                    gene_calls_out_dir,
                                                    args.file_suffix,
                                                    args.threads)

    genes_called_outfile = os.path.join(gene_calls_out_dir, 'gene_counts.json')
    write_json(genes_called, genes_called_outfile)

    if 'gtdb' in os.path.basename(args.db_file):
        db = 'gtdb_95'
    else:
        db = 'progenomes_2.1'

    diamond_outfiles = run_diamond(diamond_input, args.threads, args.temp_dir,
                                   args.db_file, args.out_dir, db)
    if not args.gene_calls and len(diamond_outfiles) != len(fastas):
        diamond_outfiles = add_empty_diamond_output(args.out_dir,
                                                    fastas,
                                                    args.file_suffix)
    gunc_output = run_gunc(diamond_outfiles, genes_called, args.out_dir,
                           args.sensitive, args.detailed_output, db,
                           args.min_mapped_genes, args.use_species_level)
    gunc_out_file = os.path.join(args.out_dir, f'GUNC.{db}.maxCSS_level.tsv')
    gunc_output.to_csv(gunc_out_file, index=False, sep='\t', na_rep='nan')


def add_empty_diamond_output(diamond_outdir, fnas, file_suffix):
    """Create placeholder file for missing diamond output.

    If an input doesnt map to reference it is missing from the final ouput.
    Creating a placeholder file means it is not silently lost.
    TODO: There are better ways to do this

    Args:
        diamond_outdir (str): Path to dir containing diamond output.
        fnas (list): All input file paths
        file_suffix (str): suffix of the input files

    Returns:
        list: Paths of diamond output files incl. empty placeholders
    """
    expected_diamond_outfiles = []
    for fna in fnas:
        basename = os.path.basename(fna).split(file_suffix)[0]
        diamond_outfile = os.path.join(diamond_outdir,
                                       'diamond_output',
                                       f'{basename}.diamond.out')
        if not os.path.isfile(diamond_outfile):
            print(f'[WARNING] no genes mapped to reference: {basename}')
            open(diamond_outfile, 'a').close()
        expected_diamond_outfiles.append(diamond_outfile)
    return expected_diamond_outfiles


def get_gene_count_file(args):
    """Search for gunc gencount file."""
    if not args.gunc_gene_count_file:
        diamond_file_path = os.path.abspath(args.diamond_file)
        gunc_dir = os.path.dirname(os.path.dirname(diamond_file_path))
        gene_counts_file = os.path.join(gunc_dir,
                                        "gene_calls/gene_counts.json")
        if not os.path.isfile(gene_counts_file):
            sys.exit('[ERROR] GUNC gene_counts.json file not found!')
    else:
        gene_counts_file = args.gunc_gene_count_file
    return gene_counts_file


def get_genecount_from_gunc_output(gene_counts_file, basename):
    """Extract gene count from GUNC gene_counts.json file."""
    with open(gene_counts_file) as f:
        data = json.load(f)
    return int(data[basename])


def plot(args):
    """Run visualisation function."""
    basename = os.path.basename(args.diamond_file).split('.diamond.')[0]
    genes_called = get_genecount_from_gunc_output(get_gene_count_file(args),
                                                  basename)
    viz_html = vis.create_viz_from_diamond_file(args.diamond_file,
                                                genes_called,
                                                args.tax_levels,
                                                args.contig_display_num,
                                                args.remove_minor_clade_level)
    viz_html_path = os.path.join(args.out_dir, f'{basename}.viz.html')
    with open(viz_html_path, 'w') as f:
        f.write(viz_html)


def merge_checkm(args):
    """Merge gunc output with checkm output."""
    merged = checkm_merge.merge_checkm_gunc(args.checkm_file, args.gunc_file)
    outfile = os.path.join(args.out_dir, 'GUNC_checkM.merged.tsv')
    merged.to_csv(outfile, sep='\t', index=False)


def main():
    args = parse_args(sys.argv[1:])
    start_time = datetime.now()
    print(f'[START] {start_time.strftime("%H:%M:%S %Y-%m-%d")}')
    if args.cmd == 'download_db':
        gunc_database.get_db(args.path, args.database)
    if args.cmd == 'run':
        start_checks()
        if not args.db_file:
            sys.exit('[WARNING] database_file argument missing.')
        else:
            if not os.path.isfile(args.db_file):
                sys.exit('[WARNING] database_file not found.')

        run(args)
    if args.cmd == 'plot':
        plot(args)
    if args.cmd == 'merge_checkm':
        merge_checkm(args)
    end_time = datetime.now()
    run_time = str(end_time - start_time).split('.')[0]
    print(f'[INFO]  {end_time.strftime("%H:%M:%S")} Runtime: {run_time}')
    print(f'[END]   {end_time.strftime("%H:%M:%S %Y-%m-%d")}')


if __name__ == "__main__":
    main()
