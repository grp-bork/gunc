#!/usr/bin/env python3
import os
import sys
import gzip
import glob
import json
import logging
import argparse
import multiprocessing
from . import get_scores
from datetime import datetime
from . import external_tools
from . import __version__
from .external_tools import get_record_count_in_fasta as record_count

logger = logging.getLogger("gunc.py")


def parse_args(args):
    """Parse Arguments

    Arguments:
        args (List): List of args supplied to script.

    Returns:
        Namespace: assigned args

    """
    description = (
        "Tool for detection of chimerism and contamination in prokaryotic genomes.\n"
    )

    parser = argparse.ArgumentParser(description=description)
    subparsers = parser.add_subparsers(title="GUNC subcommands", metavar="", dest="cmd")
    run = subparsers.add_parser("run", help="Run chimerism detection.")
    run_group = run.add_mutually_exclusive_group(required=True)
    download_db = subparsers.add_parser("download_db", help="Download GUNC db.")
    merge_checkm = subparsers.add_parser(
        "merge_checkm", help="Merge GUNC and CheckM outputs."
    )
    vis = subparsers.add_parser("plot", help="Create interactive visualisation.")
    summarise = subparsers.add_parser(
        "summarise",
        aliases=["rescore"],
        help="Re-score genomes using a different contamination cutoff (alias: rescore).",
    )

    run.add_argument(
        "-r",
        "--db_file",
        help="DiamondDB reference file. Default: GUNC_DB envvar",
        default=os.environ.get("GUNC_DB"),
        metavar="",
    )
    run.add_argument(
        "--custom_genome2taxonomy",
        help="If using a custom DB, use this to provide taxonomy of genomes in the custom DB.",
        default=None,
        metavar=""
    )
    run_group.add_argument(
        "-i", "--input_fasta", help="Input file in FASTA format.", metavar=""
    )
    run_group.add_argument(
        "-f", "--input_file", help="File with paths to FASTA format files.", metavar=""
    )
    run_group.add_argument(
        "-d", "--input_dir", help="Input dir with files in FASTA format.", metavar=""
    )
    run.add_argument(
        "-e",
        "--file_suffix",
        help="Suffix of files in input_dir. Default: .fa",
        default=".fa",
        metavar="",
    )
    run.add_argument(
        "-g",
        "--gene_calls",
        help="Input files are FASTA faa format. Default: False",
        action="store_true",
        default=False,
    )
    run.add_argument(
        "-t",
        "--threads",
        help="number of CPU threads. Default: 4",
        default="4",
        metavar="",
    )
    run.add_argument(
        "-o",
        "--out_dir",
        help="Output dir.  Default: cwd",
        default=os.getcwd(),
        metavar="",
    )
    run.add_argument(
        "--temp_dir",
        help="Directory to store temp files. Default: cwd",
        default=os.getcwd(),
        metavar="",
    )
    run.add_argument(
        "--sensitive",
        help="Run with high sensitivity. Default: False",
        action="store_true",
        default=False,
    )
    run.add_argument(
        "--detailed_output",
        help="Output scores for every taxlevel. Default: False",
        action="store_true",
        default=False,
    )
    run.add_argument(
        "--contig_taxonomy_output",
        help="Output assignments for each contig. Default: False",
        action="store_true",
        default=False,
    )
    run.add_argument(
        "--use_species_level",
        help="Allow species level to be picked as maxCSS. Default: False",
        action="store_true",
        default=False,
    )
    run.add_argument(
        "--min_mapped_genes",
        help=(
            "Don't calculate GUNC score if number of mapped "
            "genes is below this value. Default: 11"
        ),
        type=int,
        default=11,
        metavar="",
    )
    run.add_argument(
        "-v",
        "--verbose",
        help="Verbose output for debugging",
        action="store_true",
        default=False,
    )
    vis.add_argument(
        "-d",
        "--diamond_file",
        help="GUNC diamond outputfile.",
        required=True,
        metavar="",
    )
    vis.add_argument(
        "-g", "--gunc_gene_count_file", help="GUNC gene_counts.json file.", metavar=""
    )
    vis.add_argument(
        "-o", "--out_dir", help="Output directory.", default=os.getcwd(), metavar=""
    )
    vis.add_argument(
        "-t",
        "--tax_levels",
        help="Tax levels to display (comma-separated).",
        default="kingdom,phylum,family,genus,contig",
        metavar="",
    )
    vis.add_argument(
        "-r",
        "--remove_minor_clade_level",
        help="Tax level at which to remove minor clades.",
        default="kingdom",
        metavar="",
    )
    vis.add_argument(
        "-c",
        "--contig_display_num",
        help="Number of contigs to visualise. [0 plots all contigs]",
        default=1000,
        type=int,
        metavar="",
    )
    vis.add_argument(
        "-l",
        "--contig_display_list",
        help="Comma-separated list of contig names to plot",
        default=None,
        metavar="",
    )
    vis.add_argument(
        "-v",
        "--verbose",
        help="Verbose output for debugging",
        action="store_true",
        default=False,
    )
    download_db.add_argument(
        "path", help="Download database to given directory.", metavar="dest_path"
    )
    download_db.add_argument(
        "-db",
        "--database",
        help="Which db to download. progenomes_2.1, progenomes_3, gtdb_95, gtdb_214, test_data. Default: progenomes_2.1",
        default="progenomes_2.1",
        metavar="",
    )
    download_db.add_argument(
        "-v",
        "--verbose",
        help="Verbose output for debugging",
        action="store_true",
        default=False,
    )
    merge_checkm.add_argument(
        "-g",
        "--gunc_file",
        help="Path of GUNC.maxCSS_level.tsv file.",
        required=True,
        metavar="",
    )
    merge_checkm.add_argument(
        "-c", "--checkm_file", help="CheckM qa output file", required=True, metavar=""
    )
    merge_checkm.add_argument(
        "-o",
        "--out_dir",
        help="Output directory for merged file",
        default=os.getcwd(),
        metavar="",
    )
    merge_checkm.add_argument(
        "-v",
        "--verbose",
        help="Verbose output for debugging",
        action="store_true",
        default=False,
    )
    summarise.add_argument(
        "-m",
        "--max_csslevel_file",
        help="MaxCSS output file from GUNC (e.g. GUNC.progenomes_2.1.maxCSS_level.tsv).",
        required=True,
        metavar="FILE",
    )
    summarise.add_argument(
        "-d",
        "--gunc_detailed_output_dir",
        help="GUNC detailed output dir (e.g. gunc_output).",
        required=True,
        metavar="DIR",
    )
    summarise.add_argument(
        "-c",
        "--contamination_cutoff",
        help="Alternative cutoff to use.",
        default=0.05,
        type=float,
        metavar="FLOAT",
    )
    summarise.add_argument(
        "-o",
        "--output_file",
        help="File in which to write output with rescored pass.GUNC column.",
        required=True,
        metavar="FILE",
    )
    summarise.add_argument(
        "-v",
        "--verbose",
        help="Verbose output for debugging",
        action="store_true",
        default=False,
    )
    check = subparsers.add_parser(
        "check", help="Validate environment and input files before running."
    )
    check.add_argument(
        "-r",
        "--db_file",
        help="DiamondDB reference file to validate. Default: GUNC_DB envvar",
        default=os.environ.get("GUNC_DB"),
        metavar="FILE",
    )
    check.add_argument(
        "--custom_genome2taxonomy",
        help="Custom genome-to-taxonomy TSV file to validate.",
        default=None,
        metavar="FILE",
    )
    check.add_argument(
        "-o",
        "--out_dir",
        help="Output directory to check for write access.",
        default=None,
        metavar="DIR",
    )
    check.add_argument(
        "-v",
        "--verbose",
        help="Verbose output for debugging",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-v",
        "--version",
        help="Print version number and exit.",
        action="version",
        version=__version__,
    )
    if not args:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args(args)
    return args


def create_dir(path):
    """Create a directory

    Will create a directory if it doesn't already exist.

    Arguments:
        path (str): directory path
    """
    if not os.path.exists(path):
        os.makedirs(path)


_GENOME2TAX_COLS = ["genome", "kingdom", "phylum", "class", "order", "family", "genus", "species"]


def _check_result(label, ok, detail=""):
    """Print a single check result line."""
    status = "PASS" if ok else "FAIL"
    line = f"[{status}] {label}"
    if detail:
        line += f": {detail}"
    print(line)
    return ok


def validate_custom_genome2taxonomy(path):
    """Validate a custom genome-to-taxonomy TSV file.

    Arguments:
        path (str): Path to the TSV file.

    Returns:
        list: List of (label, ok, detail) tuples for each check performed.
    """
    import pandas as pd

    results = []

    if not os.path.isfile(path):
        results.append(("custom_genome2taxonomy exists", False, path))
        return results
    results.append(("custom_genome2taxonomy exists", True, path))

    if not os.access(path, os.R_OK):
        results.append(("custom_genome2taxonomy readable", False, ""))
        return results
    results.append(("custom_genome2taxonomy readable", True, ""))

    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as e:
        results.append(("custom_genome2taxonomy parseable as TSV", False, str(e)))
        return results
    results.append(("custom_genome2taxonomy parseable as TSV", True, ""))

    missing = [c for c in _GENOME2TAX_COLS if c not in df.columns]
    if missing:
        results.append((
            "custom_genome2taxonomy has required columns",
            False,
            f"missing: {', '.join(missing)}; expected: {', '.join(_GENOME2TAX_COLS)}",
        ))
    else:
        results.append(("custom_genome2taxonomy has required columns", True, ", ".join(_GENOME2TAX_COLS)))

    n_rows = len(df)
    if n_rows == 0:
        results.append(("custom_genome2taxonomy has data rows", False, "file is empty"))
    else:
        results.append(("custom_genome2taxonomy has data rows", True, f"{n_rows} genomes"))

    if n_rows > 0 and "genome" in df.columns:
        n_empty = df["genome"].isna().sum() + (df["genome"] == "").sum()
        if n_empty > 0:
            results.append(("custom_genome2taxonomy genome column has no empty values", False, f"{n_empty} empty values"))
        else:
            results.append(("custom_genome2taxonomy genome column has no empty values", True, ""))

    return results


def run_check(args):
    """Validate environment and input files without running the pipeline."""
    all_ok = True

    # Tool dependencies
    for tool in ("diamond", "prodigal", "grep", "zcat"):
        ok = external_tools.check_if_tool_exists(tool)
        all_ok = _check_result(f"{tool} on PATH", ok) and all_ok

    # Diamond version
    skip_ver = bool(os.environ.get("GUNC_SKIP_DIAMOND_VERSION_CHECK"))
    if skip_ver:
        _check_result(
            "diamond version",
            True,
            "skipped (GUNC_SKIP_DIAMOND_VERSION_CHECK set)",
        )
    else:
        required = external_tools.REQUIRED_DIAMOND_VERSION
        ok = external_tools.check_diamond_version_correct()
        detail = (
            external_tools.check_diamond_version() if external_tools.check_if_tool_exists("diamond") else "diamond not found"
        )
        all_ok = _check_result(
            f"diamond version == {required}", ok, detail
        ) and all_ok

    # Database file
    if args.db_file:
        ok = os.path.isfile(args.db_file)
        all_ok = _check_result("db_file exists", ok, args.db_file) and all_ok
        if ok:
            ok = os.access(args.db_file, os.R_OK)
            all_ok = _check_result("db_file readable", ok) and all_ok
            if ok:
                size_mb = os.path.getsize(args.db_file) / 1024 / 1024
                ok = size_mb > 0.001
                all_ok = _check_result("db_file non-empty", ok, f"{size_mb:.1f} MB") and all_ok
    else:
        all_ok = _check_result("db_file", False, "not provided and GUNC_DB env var not set") and all_ok

    # Custom genome2taxonomy
    if args.custom_genome2taxonomy:
        for label, ok, detail in validate_custom_genome2taxonomy(args.custom_genome2taxonomy):
            all_ok = _check_result(label, ok, detail) and all_ok

    # Output directory
    if args.out_dir:
        ok = os.path.isdir(args.out_dir)
        all_ok = _check_result("out_dir exists", ok, args.out_dir) and all_ok
        if ok:
            ok = os.access(args.out_dir, os.W_OK)
            all_ok = _check_result("out_dir writable", ok) and all_ok

    print()
    if all_ok:
        print("All checks passed.")
    else:
        print("One or more checks failed.")
        sys.exit(1)


def start_checks():
    """Checks if tool dependencies are available."""
    if not external_tools.check_if_tool_exists("diamond"):
        logger.error("Diamond not found.")
        sys.exit(1)
    if not os.environ.get("GUNC_SKIP_DIAMOND_VERSION_CHECK"):
        if not external_tools.check_diamond_version_correct():
            actual = external_tools.check_diamond_version()
            required = external_tools.REQUIRED_DIAMOND_VERSION
            logger.error(
                f"Diamond version {actual} found but {required} is required. "
                "Set GUNC_SKIP_DIAMOND_VERSION_CHECK=1 to bypass this check."
            )
            sys.exit(1)
    if not external_tools.check_if_tool_exists("prodigal"):
        logger.error("Prodigal not found.")
        sys.exit(1)
    if not external_tools.check_if_tool_exists("grep"):
        logger.error("grep not found.")
        sys.exit(1)
    if not external_tools.check_if_tool_exists("zcat"):
        logger.error("zcat not found.")
        sys.exit(1)


def get_files_in_dir_with_suffix(directory, suffix):
    """Get files in directory that end in suffix."""
    files = glob.glob(os.path.join(directory, f"*{suffix}"))
    if len(files) == 0:
        logger.error(
            f"No files found in {directory} with ending {suffix}. "
            f"Use --file_suffix to specify the correct extension "
            f"(e.g. --file_suffix .fna or --file_suffix .fasta)."
        )
        sys.exit(1)
    return files


def openfile(filename):
    """Open file whether gzipped or not."""
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")


def merge_genecalls(genecall_files, out_dir, file_suffix):
    """Merge genecall files.

    Merges fastas together to run diamond more efficiently.
    Adds the name of the file to each record (delimiter '/')
    so they can be separated after diamond mapping.

    Arguments:
        genecall_files (list): Paths of genecall fastas to merge
        out_dir (str): Directory to put the merged file
        file_suffix (str): Suffix of input files

    Returns:
        str: path of the merged file
    """
    merged_outfile = os.path.join(out_dir, "genecalls.merged.faa")
    with open(merged_outfile, "w") as ofile:
        for file in genecall_files:
            if os.path.isfile(file):
                with openfile(file) as infile:
                    genome_name = os.path.basename(file).replace(".genecalls.faa", "")
                    if genome_name.endswith(file_suffix):
                        genome_name = genome_name[: -len(file_suffix)]
                    for line in infile:
                        if line.startswith(">"):
                            contig_name = line.split(" ")[0]
                            line = f"{contig_name}/{genome_name}\n"
                        ofile.write(f"{line.strip()}\n")
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
    with open(diamond_outfile, "r") as f:
        for line in f:
            query = line.split("\t")[0]
            contig_name, genome_name = query.rsplit("/", 1)
            line = line.replace(query, contig_name, 1)
            output[genome_name] = output.get(genome_name, "") + line
    for genome_name in output:
        outfile = os.path.join(out_dir, f"{genome_name}.diamond.{db}.out")
        outfiles.append(outfile)
        with open(outfile, "w") as ofile:
            ofile.write(output[genome_name])
    return outfiles


def write_json(data, filename):
    """Write json data to filename."""
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4)


def get_paths_from_file(input_file):
    """Extract paths from a text file."""
    if input_file.endswith(".gz"):
        logger.error("Input looks like a FASTA file. Use -i/--input_fasta for a single genome or -d/--input_dir for a directory.")
        sys.exit(1)
    with open(input_file, "r") as f:
        paths = f.readlines()
    if not paths:
        logger.error(f"Input file is empty: {input_file}")
        sys.exit(1)
    if paths[0].startswith(">"):
        logger.error("Input looks like a FASTA file. Use -i/--input_fasta for a single genome or -d/--input_dir for a directory.")
        sys.exit(1)
    return [path.strip() for path in paths if path.strip()]


def run_from_gene_calls(faas, out_dir, file_suffix):
    """Prepare genecalls for running diamond."""
    logger.info("START Merging genecalls")
    genes_called = {}
    for faa in faas:
        basename = os.path.basename(faa).split(file_suffix)[0]
        if basename.endswith(".genecalls"):
            basename = basename.split(".genecalls")[0]
        genes_called[basename] = record_count(faa)
    diamond_inputfile = merge_genecalls(faas, out_dir, file_suffix)
    logger.info("END  Merging genecalls")
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
    logger.info("START Prodigal")
    create_dir(out_dir)
    genecall_files = []
    genes_called = {}
    prodigal_info = []
    for fna in fnas:
        basename = os.path.basename(fna).split(file_suffix)[0]
        prodigal_outfile = os.path.join(out_dir, f"{basename}.genecalls.faa")
        prodigal_info.append((fna, prodigal_outfile))
    with multiprocessing.Pool(int(threads)) as p:
        p.map(run_prodigal, prodigal_info)
    for fna, prodigal_outfile in prodigal_info:
        basename = os.path.basename(fna).split(file_suffix)[0]
        genes_called[basename] = record_count(prodigal_outfile)
        if os.path.isfile(prodigal_outfile) and os.path.getsize(prodigal_outfile) > 0:
            genecall_files.append(prodigal_outfile)
        else:
            logger.warning(f"Prodigal failed for {fna}")
    if genecall_files:
        diamond_inputfile = merge_genecalls(genecall_files, out_dir, file_suffix)
        logger.info("END   Prodigal")
        return genes_called, diamond_inputfile
    else:
        logger.error("No genecalls to run.")
        sys.exit(1)


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
    logger.info("START Diamond")
    outfile = os.path.join(out_dir, f"{os.path.basename(infile)}.diamond.{db}.out")
    external_tools.diamond(infile, threads, temp_dir, db_file, outfile)

    if infile.endswith("genecalls.merged.faa"):
        out_dir = os.path.join(out_dir, "diamond_output")
        create_dir(out_dir)
        diamond_outfiles = split_diamond_output(outfile, out_dir, db)
        os.remove(outfile)
        os.remove(infile)
    else:
        diamond_outfiles = [outfile]
    logger.info("END   Diamond")
    return diamond_outfiles


def run_gunc(
    diamond_outfiles,
    genes_called,
    out_dir,
    sensitive,
    detailed_output,
    db,
    min_mapped_genes,
    use_species_level,
    contig_taxonomy_output,
    custom_genome2taxonomy,
):
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
        contig_taxonomy_output (bool): Output contig assignments file
        custom_genome2taxonomy (str): genome2taxonomy tsv if using custom db

    Returns:
        pandas.DataFrame: One line per inputfile Gunc scores
    """
    import pandas as pd

    logger.info("START Scoring")
    gunc_output = []
    for diamond_file in diamond_outfiles:
        basename = os.path.basename(diamond_file).split(".diamond.")[0]
        gene_call_count = genes_called[basename]
        detailed, single = get_scores.chim_score(
            diamond_file,
            gene_call_count,
            sensitive,
            min_mapped_genes,
            use_species_level,
            db,
            custom_genome2taxonomy,
        )
        if detailed_output or contig_taxonomy_output:
            detailed_gunc_out_dir = os.path.join(out_dir, "gunc_output")
            create_dir(detailed_gunc_out_dir)
        if detailed_output:
            detailed_gunc_out_file = os.path.join(
                detailed_gunc_out_dir, f"{basename}.{db}.all_levels.tsv"
            )
            detailed.to_csv(detailed_gunc_out_file, index=False, sep="\t", na_rep="nan")
        if contig_taxonomy_output:
            contig_assignments_out_file = os.path.join(
                detailed_gunc_out_dir, f"{basename}.contig_assignments.tsv"
            )
            contig_assignments = create_contig_assignments(
                diamond_file, gene_call_count, custom_genome2taxonomy
            )
            if contig_assignments is not None:
                contig_assignments.to_csv(
                    contig_assignments_out_file, index=False, sep="\t"
                )
            else:
                logger.warning(f"Failed to parse {diamond_file}")
        gunc_output.append(single)
    logger.info("END   Scoring")
    if not gunc_output:
        logger.error("No scores produced. Check that input genomes mapped to the reference database.")
        sys.exit(1)
    result = pd.concat(gunc_output).sort_values("genome")
    rrs = result["reference_representation_score"].dropna()
    n_low_rrs = (rrs < 0.3).sum()
    if n_low_rrs > 0:
        logger.warning(
            f"{n_low_rrs} out of {len(rrs)} genomes were not well represented in the "
            f"GUNC reference (reference_representation_score < 0.3). "
            f"Please interpret GUNC clade separation scores for these with caution."
        )
    return result


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
        logger.error(f"Filenames appear more than once: {duplicated_items}")
        sys.exit(1)


def remove_missing_fnas(fnas):
    existing_fnas = []
    for fna in fnas:
        if os.path.isfile(fna):
            existing_fnas.append(fna)
        else:
            logger.warning(f"{fna} not found")
    return existing_fnas


def create_contig_assignments(diamond_file, gene_count, custom_genome2taxonomy=None):
    """Create contig gene assignment dataframe.

    Arguments:
        diamond_file (str): GUNC diamond output file path
        gene_count (int): Count of genes in original fasta
        custom_genome2taxonomy (str): genome2taxonomy tsv if using custom db

    Returns:
        pandas.DataFrame: with contig,
                               tax_level,
                               assignment,
                               count_of_genes_assigned columns
    """
    import pandas as pd

    diamond_basename = os.path.basename(diamond_file)
    db = get_scores.detect_db_from_filename(diamond_basename)
    tax_data, genome_name, cutoff = get_scores.get_base_data_for_plotting(
        diamond_file, gene_count, db=db, custom_genome2taxonomy=custom_genome2taxonomy
    )
    assignments = []
    if len(tax_data) == 0:
        return None
    for contig in tax_data["contig"].unique():
        contig_data = tax_data[tax_data["contig"] == contig]
        for tax_level in get_scores.TAX_LEVELS:
            counts = contig_data[tax_level].value_counts().to_dict()
            for assignment in counts:
                assignments.append(
                    {
                        "contig": contig,
                        "tax_level": tax_level,
                        "assignment": assignment,
                        "count_of_genes_assigned": counts[assignment],
                    }
                )
    return pd.DataFrame(assignments)


def run(args):
    """Run entire GUNC workflow."""
    if not os.path.isdir(args.out_dir):
        logger.error(f"Output directory {args.out_dir} doesn't exist.")
        sys.exit(1)
    if not os.path.isdir(args.temp_dir):
        logger.error(f"Temporary directory {args.temp_dir} doesn't exist.")
        sys.exit(1)
    if args.custom_genome2taxonomy:
        if not os.path.isfile(args.custom_genome2taxonomy):
            logger.error(f"--custom_genome2taxonomy file not found: {args.custom_genome2taxonomy}")
            sys.exit(1)

    if args.input_dir:
        fastas = get_files_in_dir_with_suffix(args.input_dir, args.file_suffix)
    elif args.input_file:
        fastas = get_paths_from_file(args.input_file)
    elif args.input_fasta:
        fastas = [args.input_fasta]

    fastas = remove_missing_fnas(fastas)
    if not fastas:
        logger.error("No input files found.")
        sys.exit(1)
    check_for_duplicate_filenames(fastas, args.file_suffix)

    if args.gene_calls:
        gene_calls_out_dir = args.out_dir
        genes_called, diamond_input = run_from_gene_calls(
            fastas, gene_calls_out_dir, args.file_suffix
        )
    else:
        gene_calls_out_dir = os.path.join(args.out_dir, "gene_calls")
        genes_called, diamond_input = run_from_fnas(
            fastas, gene_calls_out_dir, args.file_suffix, args.threads
        )

    genes_called_outfile = os.path.join(gene_calls_out_dir, "gene_counts.json")
    write_json(genes_called, genes_called_outfile)

    db_file = os.path.basename(args.db_file)
    if args.custom_genome2taxonomy:
        db = os.path.splitext(db_file)[0]
    else:
        db = get_scores.detect_db_from_filename(db_file)

    diamond_outfiles = run_diamond(
        diamond_input, args.threads, args.temp_dir, args.db_file, args.out_dir, db
    )
    if not args.gene_calls and len(diamond_outfiles) != len(fastas):
        diamond_outfiles = add_empty_diamond_output(
            args.out_dir, fastas, args.file_suffix, db
        )
    gunc_output = run_gunc(
        diamond_outfiles,
        genes_called,
        args.out_dir,
        args.sensitive,
        args.detailed_output,
        db,
        args.min_mapped_genes,
        args.use_species_level,
        args.contig_taxonomy_output,
        args.custom_genome2taxonomy,
    )
    gunc_out_file = os.path.join(args.out_dir, f"GUNC.{db}.maxCSS_level.tsv")
    gunc_output.to_csv(gunc_out_file, index=False, sep="\t", na_rep="nan")
    logger.info(f"Output written to: {gunc_out_file}")


def add_empty_diamond_output(diamond_outdir, fnas, file_suffix, db):
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
    count_existing_diamond_files = 0
    expected_diamond_outfiles = []
    for fna in fnas:
        basename = os.path.basename(fna).split(file_suffix)[0]
        diamond_outfile = os.path.join(
            diamond_outdir, "diamond_output", f"{basename}.diamond.{db}.out"
        )
        if not os.path.isfile(diamond_outfile):
            logger.warning(f"No genes mapped to reference: {basename}")
            open(diamond_outfile, "a").close()
        else:
            count_existing_diamond_files += 1
        expected_diamond_outfiles.append(diamond_outfile)
    if count_existing_diamond_files == 0:
        logger.error("No diamond output files.")
        sys.exit(1)
    else:
        logger.info(f"{count_existing_diamond_files}/{len(fnas)} run successfully with diamond")
    return expected_diamond_outfiles


def get_gene_count_file(args):
    """Search for gunc gencount file."""
    if not args.gunc_gene_count_file:
        diamond_file_path = os.path.abspath(args.diamond_file)
        gunc_dir = os.path.dirname(os.path.dirname(diamond_file_path))
        gene_counts_file = os.path.join(gunc_dir, "gene_calls/gene_counts.json")
        if not os.path.isfile(gene_counts_file):
            logger.error("GUNC gene_counts.json file not found!")
            sys.exit(1)
    else:
        gene_counts_file = args.gunc_gene_count_file
    return gene_counts_file


def get_genecount_from_gunc_output(gene_counts_file, basename):
    """Extract gene count from GUNC gene_counts.json file."""
    try:
        with open(gene_counts_file) as f:
            data = json.load(f)
        return int(data[basename])
    except (KeyError, ValueError) as e:
        logger.error(f"Could not read gene count for {basename} from {gene_counts_file}: {e}")
        sys.exit(1)


def plot(args):
    """Run visualisation function."""
    from . import visualisation as vis

    basename = os.path.basename(args.diamond_file).split(".diamond.")[0]
    genes_called = get_genecount_from_gunc_output(get_gene_count_file(args), basename)
    viz_html = vis.create_viz_from_diamond_file(
        args.diamond_file,
        genes_called,
        args.tax_levels,
        args.contig_display_num,
        args.contig_display_list,
        args.remove_minor_clade_level,
    )

    viz_html_path = os.path.join(args.out_dir, f"{basename}.viz.html")
    with open(viz_html_path, "w") as f:
        f.write(viz_html)
    logger.info(f"Output written to: {viz_html_path}")


def merge_checkm(args):
    """Merge gunc output with checkm output."""
    from . import checkm_merge

    merged = checkm_merge.merge_checkm_gunc(args.checkm_file, args.gunc_file)
    outfile = os.path.join(args.out_dir, "GUNC_checkM.merged.tsv")
    merged.to_csv(outfile, sep="\t", index=False)
    logger.info(f"Output written to: {outfile}")


def get_scores_using_supplied_cont_cutoff(detail_file, cutoff=0.05):
    import pandas as pd

    df = pd.read_csv(detail_file, sep="\t", header=0)
    df = df.drop(["pass.GUNC"], axis=1)
    max_CSS = df.iloc[[0]].to_dict("records")[0]
    df = df[df["taxonomic_level"] != "species"]
    df = df[df["contamination_portion"] > cutoff]
    if len(df) > 0:
        max_CSSidx = df["clade_separation_score"].idxmax()
        if not pd.isna(max_CSSidx):
            max_CSS = df.loc[[max_CSSidx]].to_dict("records")[0]
            if max_CSS["clade_separation_score"] > get_scores.CSS_CHIMERIC_THRESHOLD:
                max_CSS[f"pass.GUNC_{cutoff}"] = False
                return max_CSS
    max_CSS[f"pass.GUNC_{cutoff}"] = True
    return max_CSS


def summarise(args):
    import pandas as pd

    max_csslevel_file = pd.read_csv(args.max_csslevel_file, sep="\t", header=0)
    max_csslevel_file = max_csslevel_file.to_dict("index")
    db = os.path.basename(args.max_csslevel_file).replace("GUNC.", "").replace(".maxCSS_level.tsv", "")
    logger.debug(max_csslevel_file)
    for row in max_csslevel_file:
        pass_gunc = max_csslevel_file[row]["pass.GUNC"]
        logger.debug(pass_gunc)
        if pass_gunc is True or pass_gunc == "True" or pd.isna(pass_gunc):
            max_csslevel_file[row][f"pass.GUNC_{args.contamination_cutoff}"] = True
        elif pass_gunc is False or pass_gunc == "False":
            bin_name = max_csslevel_file[row]["genome"]
            detail_file = os.path.join(
                args.gunc_detailed_output_dir, bin_name + f".{db}.all_levels.tsv"
            )
            if not os.path.isfile(detail_file):
                logger.error(f"Detail file not found: {detail_file}")
                sys.exit(1)
            new_row = get_scores_using_supplied_cont_cutoff(
                detail_file, args.contamination_cutoff
            )
            max_csslevel_file[row] = new_row
    out_df = pd.DataFrame.from_dict(max_csslevel_file, orient="index")
    logger.debug(out_df)
    out_df.to_csv(args.output_file, index=False, sep="\t", na_rep="nan")
    logger.info(f"Output written to: {args.output_file}")


def main():
    args = parse_args(sys.argv[1:])
    if args.verbose:
        logging.basicConfig(
            format="%(asctime)s : [%(levelname)7s] : %(name)s:%(lineno)s %(funcName)20s() : %(message)s",
            datefmt="%H:%M:%S",
            level=logging.DEBUG,
        )
    else:
        logging.basicConfig(
            format="%(asctime)s : %(message)s",
            datefmt="%H:%M:%S",
            level=logging.INFO,
        )
    logger = logging.getLogger("gunc.py")
    logger.debug(args)
    start_time = datetime.now()
    logger.info(f'START {start_time.strftime("%Y-%m-%d")}')
    if args.cmd == "check":
        run_check(args)
    elif args.cmd == "download_db":
        from . import gunc_database

        gunc_database.get_db(args.path, args.database)
    elif args.cmd == "run":
        start_checks()
        if not args.db_file:
            logger.error("Database file (-r/--db_file) is required. Set it or export the GUNC_DB environment variable.")
            sys.exit(1)
        else:
            if not os.path.isfile(args.db_file):
                logger.error(f"Database file not found: {args.db_file}")
                sys.exit(1)

        run(args)
    elif args.cmd == "plot":
        plot(args)
    elif args.cmd == "merge_checkm":
        merge_checkm(args)
    elif args.cmd in ("summarise", "rescore"):
        summarise(args)
    end_time = datetime.now()
    run_time = str(end_time - start_time).split(".")[0]
    logger.info(f"END   Runtime: {run_time}")


if __name__ == "__main__":
    main()
