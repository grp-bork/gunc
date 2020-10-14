#!/usr/bin/env python3
import sys
import argparse
import subprocess


def parse_args(args):
    """Parse Arguments

    Arguments:
        args {list} -- List of args supplied to script.

    Returns:
        {Namespace} -- assigned args

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--database_file',
                        help='Diamond database reference file.',
                        required=True,
                        metavar='')
    parser.add_argument('-i', '--input_file',
                        help='Input file in FASTA format..',
                        required=True,
                        metavar='')
    parser.add_argument('-p', '--diamond_threads',
                        help='number of CPU threads used for diamond search.',
                        default='10',
                        metavar='')
    parser.add_argument('-t', '--temp_dir',
                        help='Directory to store temporary files.',
                        required=True,
                        metavar='')
    parser.add_argument('-o', '--out_file',
                        help='Output file.',
                        required=True,
                        metavar='')
    args = parser.parse_args()
    return args


def get_record_count_in_fasta(fasta_file):
    """Count number of records in a fasta file.

    Arguments:
        fasta_file {str} -- full path of fasta file

    Returns:
        str -- count of records
    """
    print('[INFO] Counting fasta records..')
    return subprocess.check_output(f'zgrep -c ">" {fasta_file}',
                                   shell=True,
                                   universal_newlines=True).strip()


def prodigal(input_file, out_file):
    """Run prodigal

    Arguments:
        input_file {str} -- full path of input fasta
        out_file {str} -- fullpath of output file
    """
    try:
        print('[INFO] Running Prodigal..', flush=True)
        subprocess.check_output(['prodigal',
                                 '-i', input_file,
                                 '-a', out_file,
                                 '-o', '/dev/null',
                                 '-p', 'meta',
                                 '-q'],
                                universal_newlines=True)
    except subprocess.CalledProcessError:
        sys.exit('[ERROR] Failed to run Prodigal')


def diamond(input_file, threads, temp_dir, database_file, out_file):
    """Run diamond.

    Arguments:
        input_file {str} -- full path of gene calls
        threads {int} -- number of threads to use
        temp_dir {str} -- path of directory to use for tmp files
        database_file {str} -- full path fof diamond db file
        out_file {str} -- full path of output file
    """
    try:
        print('[INFO] Running Diamond..', flush=True)
        subprocess.check_output(['diamond', 'blastp',
                                 '--query', input_file,
                                 '--threads', threads,
                                 '--max-target-seqs', '1',
                                 '--masking', '0',
                                 '--evalue', '1',
                                 '--tmpdir', temp_dir,
                                 '--db', database_file,
                                 '--out', out_file],
                                universal_newlines=True)
    except subprocess.CalledProcessError:
        sys.exit('[ERROR] Failed to run Diamond')


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    diamond(args.input_file,
            args.threads,
            args.temp_dir,
            args.database_file,
            args.out_file)
