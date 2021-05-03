#!/usr/bin/env python3
import os
import sys
import scipy
import scipy.stats
import numpy as np
import pandas as pd
from collections import OrderedDict
from pkg_resources import resource_filename


def read_diamond_output(file_path):
    """Read in Diamond output.

    Args:
        file_path (str): Full path to diamond output.

    Returns:
        pandas.DataFrame: diamond output
    """
    df = pd.read_csv(file_path,
                     sep='\t',
                     header=None,
                     usecols=[0, 1, 2],
                     names=['query', 'genome', 'id'])
    if len(df) == 0:
        return(df)
    else:
        df['contig'] = df['query'].str.rpartition('_')[[0]]
        return df


def get_n_effective_surplus_clades(counts):
    """Calculate Inverse Simpson Index.

    Inverse Simpson Index of all clade labels distribution
    -1 (as 1 genome is expected)

    Arguments:
        counts (pandas.Series): Counts of taxons in a taxlevel.

    Returns:
        float: Score describing the extent of chimerism,
        i.e. the effective number of surplus clades represented at a taxlevel.
    """
    if len(counts) == 0:
        return 0
    else:
        denom = sum(counts) ** 2
        return 1 / sum([x ** 2 / denom for x in counts]) - 1


def expected_entropy_estimate(probabilities, sample_count):
    """Compute the expected entropy estimate sampling N elements underlying
    probabilities `p`


    Arguments:
        probabilities (numpy.ndarray): probabilities (p)
                                       (assumed to sum to 1.0)
        sample_count (int): number of samples (N)

    Returns:
        entropy (float): expected entropy
    """
    entropy = 0.0
    for probability in probabilities:
        sample_array = np.arange(1, sample_count)
        entropy += np.sum(scipy.special.comb(sample_count, sample_array)
                          * np.power(probability, sample_array)
                          * np.power(1 - probability,
                                     sample_count - sample_array)
                          # should be sample_array/sample_count
                          # but we moved the 1/sample_count out:
                          * sample_array
                          * np.log(sample_array / sample_count))
    return -entropy / sample_count


def calc_expected_conditional_entropy(contigs, taxons):
    """Compute the expected measured conditional entropy under the null
    hypothesis that there is no relationship between contig membership and
    taxonomy

    When the bucket is large enough, the estimates are expected to be close
    enough to the global estimate that we will no longer compute the estimate
    via the more costly `expected_entropy_estimate` function

    Arguments:
        contigs (pandas.Series): contig names in data
        taxons (pandas.Series): taxons in data

    Returns:
        float: expected measured conditional entropy
    """

    MAX_BUCKET_SIZE = 500
    counts = taxons.value_counts()
    taxon_entropy = scipy.stats.entropy(counts.values)
    if taxon_entropy == 0:
        return 0.0

    bucket_sizes = contigs.value_counts().value_counts()
    nr_elements = (bucket_sizes * bucket_sizes.index)
    contribution = nr_elements / nr_elements.sum()

    taxon_probability = counts / counts.sum()
    taxon_probability = taxon_probability.values

    total_entropy = 0.0
    for bucket_size, contig_count in contribution.iteritems():
        if bucket_size <= MAX_BUCKET_SIZE:
            total_entropy += (contig_count
                              * expected_entropy_estimate(taxon_probability,
                                                          bucket_size))
        else:
            total_entropy += contig_count * taxon_entropy
    return total_entropy


def read_genome2taxonomy_reference(db):
    """Read in genome2taxonomy reference data.

    This file translates the genomes in the diamond reference file to
    taxonomy with seven taxonomic levels from kingdom to species.

    Arguments:
        db (str): Which db to use: progenomes or gtdb

    Returns:
        pandas.DataFrame: genome2taxonomy reference
    """
    if db == 'progenomes_2.1':
        genome2taxonomy = resource_filename(__name__,
                                            'data/genome2taxonomy_ref.tsv')
    elif db == 'gtdb_95':
        genome2taxonomy = resource_filename(__name__,
                                            'data/genome2taxonomy_gtdbref.tsv')
    else:
        sys.exit(f'[ERROR] {db} unknown. Allowed: progenomes_2.1, gtdb_95')
    return pd.read_csv(genome2taxonomy, sep='\t')


def create_base_data(diamond_df, db):
    """Assign taxonomies to diamond output.

    Arguments:
        diamond_df (pandas.DataFrame): diamond output
        db (str): Which db to use: progenomes or gtdb

    Returns:
        pandas.DataFrame: merged dataframe
    """
    genome2taxonomy_df = read_genome2taxonomy_reference(db)
    return pd.merge(diamond_df, genome2taxonomy_df, on='genome', how="inner")


def get_stats(diamond_df):
    """Get summary statistics of diamond mapping output.

    Arguments:
        diamond_df (pandas.DataFrame): diamond output

    Returns:
        tuple:
        (number of genes mapped to GUNC reference DB,
        number of contigs containing mapped genes)
    """
    if len(diamond_df) == 0:
        return (0, 0)
    genes_mapped = len(diamond_df)
    contig_count = diamond_df['contig'].nunique()
    return (genes_mapped, contig_count)


def get_abundant_lineages_cutoff(sensitive, genes_mapped):
    """Determine cutoff for abundant lineages.

    Removal of all genes coming from clades consisting of <2% of all mapped
    genes is intended to reduce noise introduced by genes mapping to a wide
    range of clades due to their poor representation in the reference.
    In sensitive mode that value is reduced to just 10 genes.

    Arguments:
        sensitive (bool): sets cutoff to 10 if true
        genes_mapped (int): total number of genes mapped by diamond to GUNC DB

    Returns:
        float: cutoff used to determine an abundant lineage
    """
    if sensitive:
        return 10
    else:
        return 0.02 * genes_mapped


def calc_contamination_portion(counts):
    """Calculate contamination portion

    Arguments:
        counts (pandas.Series): Counts of taxons in a taxlevel.

    Returns:
        float:
        portion of genes assigning to all clades
        except the one with most genes.
    """
    return 1 - counts.max() / counts.sum()


def mean(list_of_numbers):
    """Get mean of a list of numbers.

    Arguments:
        list_of_numbers (list): numbers to get mean of

    Returns:
        float: mean of numbers
    """
    return sum(list_of_numbers) / len(list_of_numbers)


def calc_mean_hit_identity(identity_scores):
    """Calculate mean hit identity score.

    Calculates the mean identity with which genes in abundant lineages (>2%)
    hit genes in the reference.

    Arguments:
        identity_scores (list): list of identity scores

    Returns:
        float: the mean of identity scores
    """
    if len(identity_scores) == 0:
        return 0
    return mean(identity_scores) / 100


def calc_conditional_entropy(contigs, taxons):
    """Compute conditional entropy

    Arguments:
        contigs (pandas.Series): IDs of contigs
        taxons (pandas.Series): IDs of taxonomic clade assignments

    Returns:
        float: measured conditional entropy
    """
    cross = pd.crosstab(contigs, taxons)
    MI = scipy.stats.entropy(cross.values.ravel(),
                             np.outer(cross.sum(1), cross.sum()).ravel())
    H_t = scipy.stats.entropy(taxons.value_counts())
    return H_t - MI


def calc_clade_separation_score(contamination_portion,
                                conditional_entropy,
                                expected_conditional_entropy):
    """Get clade separation score (CSS).

    CSS = 0, if contamination_portion = 0 or H(T|C) <= H(T|R)
    CSS = 1 - H(T|C)/H(T|R), else

    Arguments:
        contamination_portion (float): GUNC contamination portion,
            when equal to 0 means all genes map to the same taxonomic clade.
        conditional_entropy (float): H(T|C),
                                     the entropy of taxonomic clade
                                     labels given their contig assignment.
        expected_conditional_entropy (float): H(T|R),
            the expected value of H(T|C) given identical contig size
            distribution and given there is no relationship between taxonomic
            clade and contig labels.

    Returns:
        float: GUNC CSS
    """
    if contamination_portion == 0:
        return 0
    elif contamination_portion == np.nan:
        return np.nan
    elif expected_conditional_entropy == 0:
        return np.nan
    else:
        return (1 - conditional_entropy / expected_conditional_entropy
                if conditional_entropy <= expected_conditional_entropy
                else 0)


def determine_adjustment(genes_retained_index):
    """Determine if adjustment is necessary.

    Adjustment of GUNC CSS score is done by setting it to 0 when there are
    <40% of all called genes retained in abundant lineages/clades representing
    >2% of all mapped genes.

    Arguments:
        genes_retained_index (float): proportion of all called genes retained
                                      in abundant lineages (>2%).

    Returns:
        int: 1 or 0. If 1, keeps CSS as it was. If 0, sets CSS to zero.
    """
    if genes_retained_index > 0.4:
        return 1
    else:
        return 0


def is_chimeric(clade_separation_score_adjusted):
    """Determine if chimeric.

    The cutoff of 0.45 was identified using benchmarks and is used to call
    a genome chimeric/contaminated if CSS is higher than this cutoff.

    Arguments:
        clade_separation_score_adjusted (float): score

    Returns:
        bool: is genome chimeric
    """
    if clade_separation_score_adjusted > 0.45:
        return True
    else:
        return False


def get_scores_for_taxlevel(base_data, tax_level, abundant_lineages_cutoff,
                            genome_name, genes_called, genes_mapped,
                            contig_count, min_mapped_genes):
    """Run chimerism check.

    Calculates the various scores needed to determine if genome ic chimeric.

    Arguments:
        base_data (pandas.DataFrame): Diamond output merged with taxonomy table
        tax_level (str): tax level to run
        abundant_lineages_cutoff (float): Cutoff value for abundant lineages
        genome_name (str): Name of input genome
        genes_called (int): Number of genes called by prodigal and
                            used by diamond for mapping to GUNC DB
        genes_mapped (int): Number of genes mapped to GUNC DB by diamond
        contig_count (int): Count of contigs
        min_mapped_genes (int): Minimum number of mapped genes
                                at which to calculate scores

    Returns:
        OrderedDict: scores for chosen taxlevel
    """
    tax_data = base_data.groupby(tax_level).filter(
        lambda x: len(x) > abundant_lineages_cutoff)
    if len(tax_data) < min_mapped_genes:
        return OrderedDict({'genome': genome_name,
                            'n_genes_called': genes_called,
                            'n_genes_mapped': genes_mapped,
                            'n_contigs': contig_count,
                            'taxonomic_level': tax_level,
                            'proportion_genes_retained_in_major_clades': np.nan,
                            'genes_retained_index': np.nan,
                            'clade_separation_score': np.nan,
                            'contamination_portion': np.nan,
                            'n_effective_surplus_clades': np.nan,
                            'mean_hit_identity': np.nan,
                            'reference_representation_score': np.nan,
                            'pass.GUNC': np.nan})
    counts = tax_data[tax_level].value_counts()
    contigs = tax_data['contig']
    taxons = tax_data[tax_level]

    n_effective_surplus_clades = get_n_effective_surplus_clades(counts)
    contamination_portion = calc_contamination_portion(counts)
    mean_hit_identity = calc_mean_hit_identity(tax_data['id'].tolist())
    expected_conditional_entropy = calc_expected_conditional_entropy(contigs,
                                                                     taxons)
    conditional_entropy = calc_conditional_entropy(contigs, taxons)
    clade_separation_score = calc_clade_separation_score(contamination_portion,
                                                         conditional_entropy,
                                                         expected_conditional_entropy)

    portion_genes_retained = len(tax_data) / genes_mapped
    genes_retained_index = genes_mapped / genes_called * portion_genes_retained
    reference_representation_score = genes_retained_index * mean_hit_identity
    adjustment = determine_adjustment(genes_retained_index)
    clade_separation_score_adjusted = clade_separation_score * adjustment
    chimeric = is_chimeric(clade_separation_score_adjusted)
    passGUNC = not chimeric
    return OrderedDict({'genome':
                        genome_name,
                        'n_genes_called':
                            genes_called,
                        'n_genes_mapped':
                            genes_mapped,
                        'n_contigs':
                            contig_count,
                        'taxonomic_level':
                            tax_level,
                        'proportion_genes_retained_in_major_clades':
                            portion_genes_retained,
                        'genes_retained_index':
                            genes_retained_index,
                        'clade_separation_score':
                            clade_separation_score_adjusted,
                        'contamination_portion':
                            contamination_portion,
                        'n_effective_surplus_clades':
                            n_effective_surplus_clades,
                        'mean_hit_identity':
                            mean_hit_identity,
                        'reference_representation_score':
                            reference_representation_score,
                        'pass.GUNC':
                            passGUNC})


def chim_score(diamond_file_path, genes_called=0, sensitive=False,
               min_mapped_genes=11, use_species_level=False, db='progenomes_2.1',
               plot=False):
    """Get chimerism scores for a genome.

    Arguments:
        diamond_file_path (str): Full path to diamond output

    Keyword Arguments:
        genes_called (int):
        sensitive (bool): Run in sensitive mode (default: (False))
        min_mapped_genes (int): Minimum number of mapped genes
                                at which to calculate scores (default: (11)
        use_species_level (bool): Allow species level to be selected for maxCSS
                                  (default: (False))
        plot (bool): Return data needed for plotting (default: (False))
        db (str): Which db to use: progenomes or gtdb (default: (progenomes)

    Returns:
        pandas.DataFrame: GUNC scores
    """
    diamond_df = read_diamond_output(diamond_file_path)
    base_data = create_base_data(diamond_df, db)
    genes_mapped, contig_count = get_stats(diamond_df)

    genome_name = os.path.basename(diamond_file_path).split('.diamond.')[0]
    abundant_lineages_cutoff = get_abundant_lineages_cutoff(sensitive,
                                                            genes_mapped)
    if plot:
        return base_data, genome_name, abundant_lineages_cutoff

    scores = []
    for tax_level in ['kingdom', 'phylum', 'class',
                      'order', 'family', 'genus', 'species']:
        scores.append(get_scores_for_taxlevel(base_data,
                                              tax_level,
                                              abundant_lineages_cutoff,
                                              genome_name,
                                              genes_called,
                                              genes_mapped,
                                              contig_count,
                                              min_mapped_genes))
    df = pd.DataFrame(scores).round(2)
    df['pass.GUNC'] = df['pass.GUNC'].astype(str)
    if use_species_level:
        max_CSSidx = df['clade_separation_score'].idxmax()
    else:
        max_CSSidx = df[:-1]['clade_separation_score'].idxmax()
    max_CSS = df.iloc[[0] if pd.isna(max_CSSidx) else [max_CSSidx]]
    return df, max_CSS
