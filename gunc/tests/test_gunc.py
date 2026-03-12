import argparse
import pytest
import pandas as pd
from ..gunc import (
    get_genecount_from_gunc_output, start_checks, parse_args,
    get_files_in_dir_with_suffix, get_scores_using_supplied_cont_cutoff,
    summarise, split_diamond_output,
)
from ..get_scores import CSS_CHIMERIC_THRESHOLD, detect_db_from_filename
from importlib.resources import files as _pkg_files


def resource_filename(package, resource):
    return str(_pkg_files(package).joinpath(resource))


gunc_gene_counts = resource_filename("gunc.tests", "test_data/tiny_test.gene_counts.json")
test_data_dir = resource_filename("gunc.tests", "test_data/")


def test_get_genecount_from_gunc_output():
    assert get_genecount_from_gunc_output(gunc_gene_counts, "tiny_test.faa") == 35


def test_start_checks():
    with pytest.raises(SystemExit):
        start_checks()


def test_parse_args():
    with pytest.raises(SystemExit):
        parse_args(["-h"])
    with pytest.raises(SystemExit):
        parse_args(["-f"])
    parser = parse_args(["run", "-f", "test_path"])
    assert parser.input_file == "test_path"
    assert parser.sensitive is False


def test_get_files_in_dir_with_suffix():
    assert len(get_files_in_dir_with_suffix(test_data_dir, ".json")) == 1


# --- helpers for summarise / get_scores_using_supplied_cont_cutoff tests ---

_DETAIL_COLS = [
    "genome", "n_genes_called", "n_genes_mapped", "n_contigs",
    "taxonomic_level", "proportion_genes_retained_in_major_clades",
    "genes_retained_index", "clade_separation_score", "contamination_portion",
    "n_effective_surplus_clades", "mean_hit_identity",
    "reference_representation_score", "pass.GUNC",
]


def _detail_row(genome, tax_level, css, cont_portion, pass_gunc=True):
    return {
        "genome": genome, "n_genes_called": 100, "n_genes_mapped": 90,
        "n_contigs": 10, "taxonomic_level": tax_level,
        "proportion_genes_retained_in_major_clades": 0.9,
        "genes_retained_index": 0.9, "clade_separation_score": css,
        "contamination_portion": cont_portion, "n_effective_surplus_clades": 1.5,
        "mean_hit_identity": 90.0, "reference_representation_score": 0.8,
        "pass.GUNC": pass_gunc,
    }


def _write_detail_tsv(path, rows):
    pd.DataFrame(rows, columns=_DETAIL_COLS).to_csv(path, sep="\t", index=False, na_rep="nan")


def _write_maxcss_tsv(path, rows):
    pd.DataFrame(rows, columns=_DETAIL_COLS).to_csv(path, sep="\t", index=False, na_rep="nan")


def _make_summarise_args(tmp_path, maxcss_file, cutoff=0.05):
    return argparse.Namespace(
        max_csslevel_file=str(maxcss_file),
        contamination_cutoff=cutoff,
        gunc_detailed_output_dir=str(tmp_path),
        output_file=str(tmp_path / "out.tsv"),
    )


# --- get_scores_using_supplied_cont_cutoff ---

def test_get_scores_using_supplied_cont_cutoff_chimeric(tmp_path):
    """When a tax level has CSS > threshold and cont_portion > cutoff, genome is chimeric."""
    detail_file = tmp_path / "genome_a.progenomes_2.1.all_levels.tsv"
    _write_detail_tsv(detail_file, [
        _detail_row("genome_a", "kingdom", CSS_CHIMERIC_THRESHOLD + 0.1, 0.3),
        _detail_row("genome_a", "phylum",  0.2, 0.1),
        _detail_row("genome_a", "class",   0.1, 0.05),
        _detail_row("genome_a", "order",   0.1, 0.05),
        _detail_row("genome_a", "family",  0.1, 0.05),
        _detail_row("genome_a", "genus",   0.1, 0.05),
        _detail_row("genome_a", "species", 0.1, 0.05),
    ])
    result = get_scores_using_supplied_cont_cutoff(str(detail_file), cutoff=0.05)
    assert result["pass.GUNC_0.05"] is False


def test_get_scores_using_supplied_cont_cutoff_clean(tmp_path):
    """When no tax level exceeds the contamination cutoff, genome passes."""
    detail_file = tmp_path / "genome_b.progenomes_2.1.all_levels.tsv"
    _write_detail_tsv(detail_file, [
        _detail_row("genome_b", "kingdom", 0.1, 0.01),
        _detail_row("genome_b", "phylum",  0.1, 0.01),
        _detail_row("genome_b", "class",   0.1, 0.01),
        _detail_row("genome_b", "order",   0.1, 0.01),
        _detail_row("genome_b", "family",  0.1, 0.01),
        _detail_row("genome_b", "genus",   0.1, 0.01),
        _detail_row("genome_b", "species", 0.1, 0.01),
    ])
    result = get_scores_using_supplied_cont_cutoff(str(detail_file), cutoff=0.05)
    assert result["pass.GUNC_0.05"] is True


def test_get_scores_using_supplied_cont_cutoff_below_threshold(tmp_path):
    """When cont_portion passes the cutoff but CSS <= threshold, genome still passes."""
    detail_file = tmp_path / "genome_c.progenomes_2.1.all_levels.tsv"
    _write_detail_tsv(detail_file, [
        _detail_row("genome_c", "kingdom", CSS_CHIMERIC_THRESHOLD - 0.01, 0.3),
        _detail_row("genome_c", "phylum",  0.1, 0.3),
        _detail_row("genome_c", "class",   0.1, 0.3),
        _detail_row("genome_c", "order",   0.1, 0.3),
        _detail_row("genome_c", "family",  0.1, 0.3),
        _detail_row("genome_c", "genus",   0.1, 0.3),
        _detail_row("genome_c", "species", 0.1, 0.3),
    ])
    result = get_scores_using_supplied_cont_cutoff(str(detail_file), cutoff=0.05)
    assert result["pass.GUNC_0.05"] is True


# --- summarise ---

def test_summarise_pass_gunc_true(tmp_path):
    """Genomes with pass.GUNC=True are carried through as passing."""
    maxcss = tmp_path / "GUNC.progenomes_2.1.maxCSS_level.tsv"
    _write_maxcss_tsv(maxcss, [_detail_row("genome_a", "genus", 0.1, 0.05, pass_gunc=True)])
    args = _make_summarise_args(tmp_path, maxcss)
    summarise(args)
    out = pd.read_csv(args.output_file, sep="\t")
    assert bool(out.loc[0, "pass.GUNC_0.05"]) is True


def test_summarise_pass_gunc_nan(tmp_path):
    """Genomes with pass.GUNC=NaN (could not be scored) are treated as passing."""
    maxcss = tmp_path / "GUNC.progenomes_2.1.maxCSS_level.tsv"
    _write_maxcss_tsv(maxcss, [_detail_row("genome_b", "genus", float("nan"), float("nan"), pass_gunc=float("nan"))])
    args = _make_summarise_args(tmp_path, maxcss)
    summarise(args)
    out = pd.read_csv(args.output_file, sep="\t")
    assert bool(out.loc[0, "pass.GUNC_0.05"]) is True


def test_summarise_pass_gunc_false_rescores(tmp_path):
    """Genomes with pass.GUNC=False are rescored against the detail file."""
    db = "progenomes_2.1"
    maxcss = tmp_path / f"GUNC.{db}.maxCSS_level.tsv"
    _write_maxcss_tsv(maxcss, [_detail_row("genome_c", "genus", CSS_CHIMERIC_THRESHOLD + 0.1, 0.3, pass_gunc=False)])

    # Write a detail file where the genome is still chimeric at the given cutoff
    detail_file = tmp_path / f"genome_c.{db}.all_levels.tsv"
    _write_detail_tsv(detail_file, [
        _detail_row("genome_c", "kingdom", CSS_CHIMERIC_THRESHOLD + 0.1, 0.3),
        _detail_row("genome_c", "phylum",  0.2, 0.3),
        _detail_row("genome_c", "class",   0.1, 0.3),
        _detail_row("genome_c", "order",   0.1, 0.3),
        _detail_row("genome_c", "family",  0.1, 0.3),
        _detail_row("genome_c", "genus",   0.1, 0.3),
        _detail_row("genome_c", "species", 0.1, 0.3),
    ])
    args = _make_summarise_args(tmp_path, maxcss)
    summarise(args)
    out = pd.read_csv(args.output_file, sep="\t")
    assert bool(out.loc[0, "pass.GUNC_0.05"]) is False


def test_summarise_pass_gunc_false_passes_at_cutoff(tmp_path):
    """A genome that failed at default cutoff can pass at a stricter contamination cutoff."""
    db = "progenomes_2.1"
    maxcss = tmp_path / f"GUNC.{db}.maxCSS_level.tsv"
    _write_maxcss_tsv(maxcss, [_detail_row("genome_d", "genus", CSS_CHIMERIC_THRESHOLD + 0.1, 0.3, pass_gunc=False)])

    # Write a detail file where contamination_portion is below the cutoff at all levels
    detail_file = tmp_path / f"genome_d.{db}.all_levels.tsv"
    _write_detail_tsv(detail_file, [
        _detail_row("genome_d", "kingdom", CSS_CHIMERIC_THRESHOLD + 0.1, 0.01),
        _detail_row("genome_d", "phylum",  0.2, 0.01),
        _detail_row("genome_d", "class",   0.1, 0.01),
        _detail_row("genome_d", "order",   0.1, 0.01),
        _detail_row("genome_d", "family",  0.1, 0.01),
        _detail_row("genome_d", "genus",   0.1, 0.01),
        _detail_row("genome_d", "species", 0.1, 0.01),
    ])
    args = _make_summarise_args(tmp_path, maxcss)
    summarise(args)
    out = pd.read_csv(args.output_file, sep="\t")
    assert bool(out.loc[0, "pass.GUNC_0.05"]) is True


# --- detect_db_from_filename (T4) ---

def test_detect_db_from_filename_gtdb214():
    assert detect_db_from_filename("genome.diamond.gtdb_214.out") == "gtdb_214"
    assert detect_db_from_filename("genome.diamond.gtdb214.out") == "gtdb_214"


def test_detect_db_from_filename_gtdb95():
    assert detect_db_from_filename("genome.diamond.gtdb_95.out") == "gtdb_95"
    assert detect_db_from_filename("genome.diamond.gtdb95.out") == "gtdb_95"


def test_detect_db_from_filename_progenomes3():
    assert detect_db_from_filename("genome.diamond.progenomes_3.out") == "progenomes_3"
    assert detect_db_from_filename("genome.diamond.progenomes3.out") == "progenomes_3"


def test_detect_db_from_filename_progenomes2_default():
    assert detect_db_from_filename("genome.diamond.progenomes_2.1.out") == "progenomes_2.1"
    assert detect_db_from_filename("genome.diamond.out") == "progenomes_2.1"


# --- split_diamond_output round-trip (T3) ---

def _make_diamond_line(query, target="1747.SAMN03982942"):
    return f"{query}\t{target}\t98.6\t145\t2\t0\t1\t145\t139\t283\t9.1e-78\t296.6\n"


def test_split_diamond_output_round_trip(tmp_path):
    """merge_genecalls encodes genome/contig; split_diamond_output decodes back."""
    merged = tmp_path / "merged.diamond.out"
    with open(merged, "w") as f:
        f.write(_make_diamond_line(">contig_1/genome_a"))
        f.write(_make_diamond_line(">contig_2/genome_a"))
        f.write(_make_diamond_line(">contig_3/genome_b"))

    out_dir = tmp_path / "split"
    out_dir.mkdir()
    outfiles = split_diamond_output(str(merged), str(out_dir), "progenomes_2.1")

    assert len(outfiles) == 2
    genome_a_file = out_dir / "genome_a.diamond.progenomes_2.1.out"
    genome_b_file = out_dir / "genome_b.diamond.progenomes_2.1.out"
    assert genome_a_file.exists()
    assert genome_b_file.exists()

    lines_a = genome_a_file.read_text().splitlines()
    assert len(lines_a) == 2
    assert lines_a[0].startswith(">contig_1\t")

    lines_b = genome_b_file.read_text().splitlines()
    assert len(lines_b) == 1
    assert lines_b[0].startswith(">contig_3\t")


def test_split_diamond_output_slash_in_contig_name(tmp_path):
    """Contig names containing '/' are handled correctly via rsplit."""
    merged = tmp_path / "merged.diamond.out"
    with open(merged, "w") as f:
        f.write(_make_diamond_line("path/to/contig_1/genome_x"))

    out_dir = tmp_path / "split"
    out_dir.mkdir()
    outfiles = split_diamond_output(str(merged), str(out_dir), "progenomes_2.1")

    assert len(outfiles) == 1
    genome_x_file = out_dir / "genome_x.diamond.progenomes_2.1.out"
    assert genome_x_file.exists()
    lines = genome_x_file.read_text().splitlines()
    assert lines[0].startswith("path/to/contig_1\t")
