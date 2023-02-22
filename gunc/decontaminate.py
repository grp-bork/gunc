import gzip
import pandas as pd


def choose_decontamination_level(
    df,
    decont_species_level=False,
):
    """Decide to decontaminate or not and if yes, choose level.

    This function takes "detailed" GUNC scores per level and returns the
    taxonomic level with the highest product of the two scores (CSS & RRS)
    which also doesn't pass GUNC quality control (GUNC.pass == False).

    Args:
        df (pandas.DataFrame): detailed GUNC scores per level
        decont_species_level (bool): A boolean indicating whether to consider
                                     the species level as a decontamination level.

    Returns:
        decontamination level
    """
    if not decont_species_level:
        df = df[:-1]
    df = df[df["pass.GUNC"] == "False"]
    if df.shape[0] > 0:
        df["CSS_RRS_product"] = (
            df["clade_separation_score"] * df["reference_representation_score"]
        )
        return [True, df.loc[df["CSS_RRS_product"].idxmax(), "taxonomic_level"]]
    else:
        return [False, "kingdom"]


def get_contig_gene_counts(base_data, decontamination_level, min_contig_gene_percent=0):
    output = []
    for name, group in base_data.groupby("contig"):
        genes_in_contig = len(group)
        clades_in_contig = set(list(group[decontamination_level]))
        if len(clades_in_contig) == 1:
            output.append([name, genes_in_contig, clades_in_contig.pop()])
        else:
            counts = group[decontamination_level].value_counts()
            total_gene_count = counts.sum()
            clade_with_max_genecount = counts.index[0]
            most_gene_count = counts.iloc[0]
            second_most_geneCount = counts.iloc[1]
            if most_gene_count / total_gene_count > min_contig_gene_percent:
                if most_gene_count != second_most_geneCount:
                    output.append([name, genes_in_contig, clade_with_max_genecount])
                else:
                    group = group[
                        group[decontamination_level].isin(
                            counts[counts == most_gene_count].index.to_list()
                        )
                    ]
                    mean_id = group.groupby(decontamination_level).mean()
                    max_mean_id = mean_id.iloc[0][0]
                    second_max_id = mean_id.iloc[1][0]
                    if max_mean_id != second_max_id:
                        output.append([name, genes_in_contig, mean_id.idxmax()[0]])
                    else:
                        group = group[
                            group[decontamination_level].isin(
                                mean_id[mean_id["id"] == max_mean_id].index.to_list()
                            )
                        ]
                        output.append(
                            [
                                name,
                                genes_in_contig,
                                mean_id[mean_id == max_mean_id]
                                .sort_index()
                                .idxmax()[0],
                            ]
                        )
            else:
                print(
                    f"[WARNING] {name} with {genes_in_contig} genes is being removed @ {decontamination_level}"
                )
    return pd.DataFrame(output, columns=["contig", "n_genes", "clade"]).set_index(
        "contig"
    )


def transform_base_data(
    base_data, decontamination_level, min_genes_for_decontamination
):
    contig_gene_counts = get_contig_gene_counts(base_data, decontamination_level)

    clade_gene_sums = (
        contig_gene_counts.groupby("clade")["n_genes"]
        .sum()
        .reset_index()
        .sort_values("n_genes", ascending=False)
    )

    clades_to_decontaminate = clade_gene_sums[
        clade_gene_sums["n_genes"] > min_genes_for_decontamination
    ]["clade"].to_list()

    contigs_in_each_clade = (
        contig_gene_counts.reset_index().groupby("clade").agg({"contig": list})
    )

    deconted_clade_contigs = {
        clade: contigs_in_each_clade.loc[clade, "contig"]
        for clade in clades_to_decontaminate
    }

    return deconted_clade_contigs


def get_gene_count(gene_calls, decontaminated_clade_contigs):
    gene_ids = []
    with open(gene_calls, "r") as f:
        for line in f:
            if line.startswith(">"):
                gene_ids.append(line.strip()[1:].split(" ")[0].rsplit("_", 1)[0])
    return len([x for x in gene_ids if x in decontaminated_clade_contigs])


def openfile(filename):
    """Open file whether gzipped or not."""
    if filename.endswith(".gz"):
        return gzip.open(filename, "r")
    else:
        return open(filename, "r")


def create_decontaminated_bin_fasta(fna, out_file, decontaminated_clade_contigs):
    with open(out_file, "w") as of:
        with openfile(fna) as f:
            for line in f:
                if line.startswith(">"):
                    writeline = False
                    contig_id = line.strip()[1:]
                    if contig_id in decontaminated_clade_contigs:
                        writeline = True
                if writeline:
                    of.write(line)
