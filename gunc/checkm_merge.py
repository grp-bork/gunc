#!/usr/bin/env python3
import sys
import pandas as pd


def read_tsv(tsv_file):
    """Read tsv to pandas.DataFrame

    Args:
        tsv_file (str): Pathof tsv file to read.

    Returns:
        pandas.DataFrame: of tsv file
    """
    return pd.read_csv(tsv_file, sep='\t')


def merge_checkm_gunc(checkm_file, gunc_file):
    """Merge checkM and gunc outputs.

    Args:
        checkm_file (str): Path of qa.tsv file form checkm
        gunc_file (str): Path of gunc_scores.tsv file

    Returns:
        pandas.DataFrame: Merged scores
    """
    checkm = read_tsv(checkm_file)
    gunc = read_tsv(gunc_file)
    output = []
    for guncdata in gunc.itertuples():
        try:
            samplename = guncdata.genome
        except AttributeError:
            sys.exit(f'[ERROR] Invalid input file: {gunc_file}')
        checkmdata = checkm[checkm['Bin Id'] == samplename]
        if len(checkmdata) == 1:
            checkmdata = checkmdata.to_dict(orient='records')[0]
        elif len(checkmdata) == 0:
            print(f'[WARNING] {samplename} not found in {checkm_file}.')
            continue
        else:
            sys.exit(f'[ERROR] {samplename} appears more '
                     f'than once in {checkm_file}')
        MIMAG_medium = (checkmdata['Completeness'] >= 50
                        and checkmdata['Contamination'] < 10)
        MIMAG_high = (checkmdata['Completeness'] >= 90
                      and checkmdata['Contamination'] < 5)
        passGUNC = guncdata.clade_separation_score < 0.45
        line = {'genome': samplename,
                'GUNC.n_contigs': guncdata.n_contigs,
                'GUNC.n_genes_called': guncdata.n_genes_called,
                'GUNC.n_genes_mapped': guncdata.n_genes_mapped,
                'GUNC.divergence_level': guncdata.taxonomic_level,
                'GUNC.contamination_portion': guncdata.contamination_portion,
                'GUNC.n_effective_surplus_clades': guncdata.n_effective_surplus_clades,
                'GUNC.RRS': guncdata.reference_representation_score,
                'GUNC.CSS': guncdata.clade_separation_score,
                'checkM.lineage': checkmdata['Marker lineage'],
                'checkM.genome_size': checkmdata.get('Genome size (bp)', ''),
                'checkM.GC': checkmdata.get('GC', ''),
                'checkM.coding_density': checkmdata.get('Coding density', ''),
                'checkM.N50_contigs': checkmdata.get('N50 (contigs)', ''),
                'checkM.completeness': checkmdata['Completeness'],
                'checkM.contamination': checkmdata['Contamination'],
                'checkM.strain_heterogeneity': checkmdata['Strain heterogeneity'],
                'pass.MIMAG_medium': MIMAG_medium,
                'pass.MIMAG_high': MIMAG_high,
                'pass.GUNC': passGUNC}
        output.append(line)
    return(pd.DataFrame(output))
