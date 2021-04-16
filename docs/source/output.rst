===========
GUNC Output
===========

Output files
------------

Normal Output
   *GUNC.{DB version}maxCSS_level.tsv* - output file with scores for a taxonomic level with the highest CSS score (or the level closest to kingdom if multiple maxima).

Detailed Output
   *{PREFIX}.all_levels.tsv* - output file with results for each taxonomic level.

Output Columns
--------------

genome
   name of input genome
n_genes_called
   number of genes called by prodigal or directly provided by the user.
n_genes_mapped
   number of genes mapped by diamond into GUNC refDB.
n_contigs
   number of contigs containing mapped genes.
taxonomic_level
   taxonomic clade labels at this taxonomic level were used to calculate values in all following columns. For each genome, all scores at six levels (species level can be added using a command-line option) are calculated.
proportion_genes_retained_in_major_clades
   only major clades that have >2% of all mapped genes assigned to them are retained to calculate other scores. Value of this column is n_genes_retained/n_genes_mapped.
genes_retained_index
   ``n_genes_mapped/n_genes_called * proportion_genes_retained_in_major_clades``, i.e. a portion of all called genes retained in major clades.
clade_separation_score
   a result of applying a formula explained in GUNC paper to taxonomy and contig labels of genes retained in major clades. Ranges from 0 to 1 and is set to 0 when genes_retained index is <0.4 because that is too few genes left.
contamination_portion
   a portion of genes retained in major clades assigned to all clades except the one clade with the highest proportion of genes assigned to it.
n_effective_surplus_clades
   an Inverse Simpson Index of fractions of all clades - 1 (as 1 genome is expected). It is a score describing the extent of chimerism, i.e. the effective number of surplus clades represented at a taxlevel.
mean_hit_identity
   the mean identity with which genes in abundant lineages (>2%) hit genes in the reference.
reference_representation_score
   ``genes_retained_index * mean_hit_identity``
   Estimates how well a genome is represented in the GUNC DB.
pass.GUNC
   If a genome passes GUNC analysis it means it is likely to not be chimeric (or that chimerism cannot be detected especially when its reference representation (RRS) is low).
   A genome passes if clade_separation_score <= 0.45, a cutoff benchmarked using simulated genomes.


.. note::
   Please note that most of genomes having reference_representation_score >0.5 (roughly) are labelled as passing GUNC filters not necessarily because they are non-chimeric but rather because they are so poorly represented in the reference that it is much more difficult to judge.

.. note::
   MIMAG_medium and MIMAG_high filters are incomplete; the MIMAG standard additionally requires data on rRNA and tRNA counts
