import pytest
import numpy as np
import pandas as pd
from ..get_scores import (
    read_diamond_output, get_n_effective_surplus_clades, calc_conditional_entropy,
    create_base_data, get_stats, get_abundant_lineages_cutoff, calc_contamination_portion,
    mean, calc_mean_hit_identity, is_chimeric, get_scores_for_taxlevel, chim_score,
    read_genome2taxonomy_reference,
)
from importlib.resources import files as _pkg_files

def resource_filename(package, resource):
    return str(_pkg_files(package).joinpath(resource))

diamond_output = resource_filename("gunc.tests", "test_data/tiny_test.diamond.out")
diamond_df = read_diamond_output(diamond_output)
ref_base_data_path = resource_filename("gunc.tests", "test_data/tiny_test.base_data")
ref_base_data = pd.read_csv(ref_base_data_path)


def test_read_diamond_output():
    empty_file = resource_filename("gunc.tests", "__init__.py")
    assert len(read_diamond_output(empty_file)) == 0
    diamond_df = read_diamond_output(diamond_output)
    assert len(diamond_df) == 17
    assert len(diamond_df.columns) == 4


def test_get_n_effective_surplus_clades():
    assert get_n_effective_surplus_clades([1]) == 0
    assert get_n_effective_surplus_clades([2, 8, 1, 1, 3]) == pytest.approx(
        1.84, rel=1e-2
    )


def test_calc_expected_conditional_entropy():
    s = pd.Series
    assert calc_conditional_entropy(s([0, 0, 1, 1]), s([1, 1, 0, 0])) == 0
    assert calc_conditional_entropy(s([0, 0, 1, 1]), s([0, 1, 0, 1])) > 0.1


def test_create_base_data():
    new_base_data = create_base_data(diamond_df, "progenomes_2.1", None)
    pd.testing.assert_frame_equal(
        ref_base_data.sort_values("query").reset_index(drop=True),
        new_base_data.sort_values("query").reset_index(drop=True),
    )


def test_get_stats():
    assert get_stats(diamond_df) == (17, 15)


def test_get_abundant_lineages_cutoff():
    assert get_abundant_lineages_cutoff(True, 3) == 10
    assert get_abundant_lineages_cutoff(False, 10) == 0.2


def test_calc_contamination_portion():
    counts = pd.Series({"a": 1, "b": 2, "c": 3, "d": 4})
    assert calc_contamination_portion(counts) == 0.6


def test_mean():
    assert mean([1, 2, 3]) == 2
    assert mean([0]) == 0
    assert mean([1, 1, 1]) == 1


def test_calc_mean_hit_identity():
    assert calc_mean_hit_identity([1, 2, 3]) == 0.02


def test_is_chimeric():
    assert is_chimeric(0) is False
    assert is_chimeric(1) is True
    assert is_chimeric(0.44) is False
    assert is_chimeric(0.45) is False
    assert is_chimeric(0.46) is True


def test_get_scores_for_taxlevel():
    # Test kingdom level
    data = get_scores_for_taxlevel(
        ref_base_data, "kingdom", 0.34, "test", 35, 17, 15, 11
    )
    assert data["genome"] == "test"
    assert data["n_contigs"] == 15
    assert data["n_genes_called"] == 35
    assert data["n_genes_mapped"] == 17
    assert data["taxonomic_level"] == "kingdom"
    assert data["clade_separation_score"] == 0
    assert data["contamination_portion"] == 0
    assert data["n_effective_surplus_clades"] == 0
    assert data["proportion_genes_retained_in_major_clades"] == 1
    assert round(data["mean_hit_identity"], 2) == 0.92
    assert round(data["genes_retained_index"], 2) == 0.49
    assert round(data["reference_representation_score"], 2) == 0.45
    assert data["pass.GUNC"]

    # Test species level
    data = get_scores_for_taxlevel(
        ref_base_data, "species", 0.34, "test", 35, 17, 15, 11
    )
    assert data["genome"] == "test"
    assert data["n_contigs"] == 15
    assert data["n_genes_called"] == 35
    assert data["n_genes_mapped"] == 17
    assert data["taxonomic_level"] == "species"
    assert round(data["clade_separation_score"], 2) == 1
    assert round(data["contamination_portion"], 2) == 0.35
    assert round(data["n_effective_surplus_clades"], 2) == 1.28
    assert data["proportion_genes_retained_in_major_clades"] == 1
    assert round(data["mean_hit_identity"], 2) == 0.92
    assert round(data["genes_retained_index"], 2) == 0.49
    assert round(data["reference_representation_score"], 2) == 0.45
    assert data["pass.GUNC"] is False

    # Test case where nothing is left after minor clade filtering
    data = get_scores_for_taxlevel(ref_base_data, "class", 1000, "test", 35, 17, 15, 11)
    assert data["genome"] == "test"
    assert data["n_contigs"] == 15
    assert data["n_genes_called"] == 35
    assert data["n_genes_mapped"] == 17
    assert data["taxonomic_level"] == "class"
    np.testing.assert_equal(data["clade_separation_score"], np.nan)
    np.testing.assert_equal(data["contamination_portion"], np.nan)
    np.testing.assert_equal(data["n_effective_surplus_clades"], np.nan)
    np.testing.assert_equal(data["proportion_genes_retained_in_major_clades"], np.nan)
    np.testing.assert_equal(data["mean_hit_identity"], np.nan)
    np.testing.assert_equal(data["genes_retained_index"], np.nan)
    np.testing.assert_equal(data["reference_representation_score"], np.nan)
    np.testing.assert_equal(data["pass.GUNC"], np.nan)


_GENOME2TAX_COLS = ["genome", "kingdom", "phylum", "class", "order", "family", "genus", "species"]


def test_read_genome2taxonomy_reference_progenomes2():
    df = read_genome2taxonomy_reference("progenomes_2.1", None)
    assert len(df) > 0
    assert list(df.columns) == _GENOME2TAX_COLS


def test_read_genome2taxonomy_reference_progenomes3():
    df = read_genome2taxonomy_reference("progenomes_3", None)
    assert len(df) > 0
    assert list(df.columns) == _GENOME2TAX_COLS


def test_read_genome2taxonomy_reference_gtdb95():
    df = read_genome2taxonomy_reference("gtdb_95", None)
    assert len(df) > 0
    assert list(df.columns) == _GENOME2TAX_COLS


def test_read_genome2taxonomy_reference_gtdb214():
    df = read_genome2taxonomy_reference("gtdb_214", None)
    assert len(df) > 0
    assert list(df.columns) == _GENOME2TAX_COLS


def test_read_genome2taxonomy_reference_custom(tmp_path):
    custom_tsv = tmp_path / "custom.tsv"
    import pandas as pd
    pd.DataFrame([
        {"genome": "g1", "kingdom": "Bacteria", "phylum": "p", "class": "c",
         "order": "o", "family": "f", "genus": "g", "species": "s"},
    ]).to_csv(custom_tsv, sep="\t", index=False)
    df = read_genome2taxonomy_reference("progenomes_2.1", str(custom_tsv))
    assert len(df) == 1
    assert df.iloc[0]["genome"] == "g1"


def test_read_genome2taxonomy_reference_unknown_db():
    with pytest.raises(SystemExit):
        read_genome2taxonomy_reference("unknown_db", None)


def test_chim_score():
    diamond_file_path = resource_filename(
        "gunc.tests", "test_data/test_genome.fa.diamond.out"
    )
    data, _ = chim_score(diamond_file_path, genes_called=1832)

    expected_data_path = resource_filename(
        "gunc.tests", "test_data/test_genome.fa.diamond.out.chimerism_scores"
    )
    expected_data = pd.read_csv(expected_data_path, sep="\t")
    pd.testing.assert_frame_equal(data, expected_data, check_dtype=False)
