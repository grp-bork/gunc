import pytest
import pandas as pd
from ..get_scores import *
from pkg_resources import resource_filename

diamond_output = resource_filename(__name__,
                                   'test_data/tiny_test.diamond.out')
diamond_df = read_diamond_output(diamond_output)
ref_base_data_path = resource_filename(__name__,
                                       'test_data/tiny_test.base_data')
ref_base_data = pd.read_csv(ref_base_data_path)


def test_get_chimerism_score():
    assert get_chimerism_score([1]) == 0
    assert get_chimerism_score([2, 8, 1, 1, 3]) == pytest.approx(1.84,
                                                                 rel=1e-2)


def test_calc_completeness_score():
    assert calc_completeness_score([0, 0, 1, 1], [1, 1, 0, 0]) == 1
    assert calc_completeness_score([0, 0, 1, 1], [0, 1, 0, 1]) == 0


# def test_calc_mean_clade_separation_score():
#     # This is awkward as its supposed to be a bit random
#     mean = calc_mean_clade_separation_score([0, 0, 1, 1], [1, 1, 0, 0])
#     print(mean)
#     assert mean > .25
#     assert mean < .75


def test_create_base_data():
    new_base_data = create_base_data(diamond_df)
    pd.testing.assert_frame_equal(ref_base_data, new_base_data)


def test_get_stats():
    assert get_stats(diamond_df) == (17, 15)


def test_get_abundant_lineages_cutoff():
    assert get_abundant_lineages_cutoff(True, 3) == 10
    assert get_abundant_lineages_cutoff(False, 10) == 0.2


def test_calc_contamination_portion():
    counts = pd.Series({'a': 1, 'b': 2, 'c': 3, 'd': 4})
    assert calc_contamination_portion(counts) == 0.6


def test_mean():
    assert mean([1, 2, 3]) == 2
    assert mean([0]) == 0
    assert mean([1, 1, 1]) == 1


def test_calc_mean_hit_identity():
    assert calc_mean_hit_identity([1, 2, 3]) == 0.02


def test_column_to_list():
    test_list = ['k115_6143_1',
                 'k115_4699_1',
                 'k115_7226_1',
                 'k115_6505_1',
                 'k115_6866_1',
                 'k115_1087_1',
                 'k115_8310_1',
                 'k115_7589_1',
                 'k115_2172_1',
                 'k115_3977_1',
                 'k115_726_1',
                 'k115_1450_1',
                 'k115_2_1',
                 'k115_2_2',
                 'k115_364_1',
                 'k115_364_2',
                 'k115_6868_1']
    print(diamond_df)
    assert column_to_list(diamond_df, 'query') == test_list


def test_calc_clade_separation_score():
    assert calc_clade_separation_score(0, 10) == 0
    assert calc_clade_separation_score(1, 10) == 10


def test_determine_adjustment():
    assert determine_adjustment(1, 2, 3) == 0
    assert determine_adjustment(2, 1, 3) == 1
    assert determine_adjustment(2, 1, 0) == 0


def test_is_chimeric():
    assert is_chimeric(0) is False
    assert is_chimeric(1) is True
    assert is_chimeric(0.39) is False
    assert is_chimeric(0.4) is False
    assert is_chimeric(0.41) is True


def test_get_scores_for_taxlevel():
    data = get_scores_for_taxlevel(ref_base_data,
                                   'kingdom',
                                   0.34,
                                   'test',
                                   35,
                                   17,
                                   15)
    assert data['genome'] == 'test'
    assert data['n_contigs'] == 15
    assert data['n_genes_called'] == 35
    assert data['n_genes_mapped'] == 17
    assert data['taxonomic_level'] == 'kingdom'
    assert data['clade_separation_score'] == 0
    assert data['contamination_portion'] == 0
    assert data['n_effective_clades'] == 0
    assert data['genes_retained_in_major_clades'] == 1
    assert round(data['mean_hit_identity'], 2) == 0.92
    assert round(data['mean_random_clade_separation_score'], 2) == pytest.approx(1, rel=1e-1)
    assert round(data['genes_retained_index'], 2) == 2.06
    assert round(data['out_of_reference_score'], 2) == 1.89
    assert data['adjustment'] == 0
    assert data['clade_separation_score_adjusted'] == 0
    assert data['chimeric'] is False

    data = get_scores_for_taxlevel(ref_base_data,
                                   'specI',
                                   0.34,
                                   'test',
                                   35,
                                   17,
                                   15)
    assert data['genome'] == 'test'
    assert data['n_contigs'] == 15
    assert data['n_genes_called'] == 35
    assert data['n_genes_mapped'] == 17
    assert data['taxonomic_level'] == 'specI'
    assert round(data['clade_separation_score'], 2) == 1
    assert round(data['contamination_portion'], 2) == 0.35
    assert round(data['n_effective_clades'], 2) == 1.28
    assert data['genes_retained_in_major_clades'] == 1
    assert round(data['mean_hit_identity'], 2) == 0.92
    assert round(data['mean_random_clade_separation_score'], 2) == pytest.approx(0.92, rel=1e-1)
    assert round(data['genes_retained_index'], 2) == 2.06
    assert round(data['out_of_reference_score'], 2) == 1.89
    assert data['adjustment'] == 0
    assert data['clade_separation_score_adjusted'] == 0
    assert data['chimeric'] is False


def test_chim_score():
    diamond_file_path = resource_filename(__name__,
                                          'test_data/test_genome.fa.diamond.out')
    data = chim_score(diamond_file_path, gene_count=1832)

    expected_data_path = resource_filename(__name__,
                                           'test_data/test_genome.fa.diamond.out.chimerism_scores')
    expected_data = pd.read_csv(expected_data_path, sep='\t')
    mrcss = 'mean_random_clade_separation_score'
    pd.testing.assert_frame_equal(data.drop(mrcss, axis=1),
                                  expected_data.drop(mrcss, axis=1))
    assert data[mrcss].tolist() == pytest.approx(expected_data[mrcss].tolist(), rel=1e-1)
