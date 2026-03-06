=====
Usage
=====

Quickstart
----------

.. code-block:: shell

    gunc run -i genome.fa -r /path/to/db

This will run gunc on genome.fa with outputs going to the current working directory.

Main Commands
-------------

 * :code:`gunc run` The main functionality of GUNC, runs chimerism detection.
 * :code:`gunc check` Validate your environment and input files before running a job.
 * :code:`gunc plot` Produce an interactive plot using the output from :code:`gunc run`
 * :code:`gunc merge_checkm` Produce a merged file combining the outputs of GUNC and checkM
 * :code:`gunc rescore` (alias: :code:`gunc summarise`) Re-score genomes from a previous :code:`gunc run` using a different contamination cutoff
 * :code:`gunc download_db` Download the GUNC database (required to run :code:`gunc run`)

Any of the above commands can be run with :code:`-h` to get function specific information.

------------

GUNC accepts a progenomes or GTDB based reference database via the :code:`--db_file` option. Four databases are available: ``progenomes_2.1`` (default), ``progenomes_3``, ``gtdb_95``, and ``gtdb_214``. All can be downloaded using the :code:`gunc download_db` command (see below). Note that using GTDB will lead to higher resource requirements and longer run times; in accuracy benchmarks, the performance of GTDB and the default proGenomes-derived GUNC database performed very similarly.


GUNC RUN
--------

Run chimerism detection.

Required Flags
^^^^^^^^^^^^^^

 * :code:`--db_file` Path to the GUNC database file. Can be set as environment variable GUNC_DB.

One of the following is required. If flag :code:`--gene_calls` is not set gene calling will be done using prodigal with option "-p meta".

 * :code:`--input_dir` Input dir with files in FASTA format.
 * :code:`--input_file` Input file containing paths to FASTA format files.
 * :code:`--input_fasta` Input file in FASTA format.

Optional Flags
^^^^^^^^^^^^^^

 * :code:`--file_suffix` Suffix of input files. Default: ``.fa``. Use this if your files end in ``.fna`` or ``.fasta``.
 * :code:`--gene_calls` Input is FASTA faa format genecalls.
 * :code:`--use_species_level` Allow species level to be picked as maxCSS. Default: False
 * :code:`--min_mapped_genes` Don't calculate GUNC score if number of mapped genes is below this value. Default: 11
 * :code:`--threads` Number of CPU threads.
 * :code:`--temp_dir` Directory to store temporary files. Default: Current working directory.
 * :code:`--out_dir` Directory in which to put output. Default: Current working directory.
 * :code:`--sensitive` Run with high sensitivity. (Uses a different cutoff to determine an abundant lineage)
 * :code:`--detailed_output` Output scores for every tax_level.
 * :code:`--contig_taxonomy_output` Output taxonomic assignment for each contig.
 * :code:`--custom_genome2taxonomy` Path to a custom genome-to-taxonomy TSV file when using a custom database (see below).

Custom Database Format
^^^^^^^^^^^^^^^^^^^^^^

When using :code:`--custom_genome2taxonomy`, the file must be a tab-separated TSV with the following columns in this order:

.. code-block:: text

    genome  kingdom  phylum  class  order  family  genus  species

The ``genome`` column must match the sequence identifiers used in your diamond database. No header transformation is needed — the column names must appear exactly as shown.

------------

GUNC CHECK
----------

Validate your environment and input files before submitting a long job. Useful on clusters where tool availability varies between login and compute nodes.

.. code-block:: shell

    gunc check -r /path/to/gunc.dmnd --out_dir /path/to/output

Each check prints ``[PASS]`` or ``[FAIL]`` and exits with code 1 if any check fails.

Optional Flags
^^^^^^^^^^^^^^

 * :code:`-r / --db_file` Diamond database file to validate. Default: ``GUNC_DB`` env var.
 * :code:`--custom_genome2taxonomy` Custom genome-to-taxonomy TSV file to validate (checks existence, readability, correct columns, and non-empty ``genome`` column).
 * :code:`-o / --out_dir` Output directory to check for write access.

------------


GUNC PLOT
---------

Create interactive plot to visualise chimerism.

Required Flags
^^^^^^^^^^^^^^

 * :code:`--diamond_file` GUNC diamond outputfile. (one of the output files in :code:`diamond_output` produced by :code:`gunc run`)

Optional Flags
^^^^^^^^^^^^^^

 * :code:`--gunc_gene_count_file` GUNC gene_counts.json file. (Not needed if :code:`--diamond` file is in the file structure made by :code:`gunc run`)
 * :code:`--out_dir` Output directory.  Default: Current working directory.
 * :code:`--tax_levels` Tax levels to display (comma-separated). (default: kingdom,phylum,family,genus,contig)
 * :code:`--remove_minor_clade_level` Tax level at which to remove minor clades. (default: kingdom)
 * :code:`--contig_display_num` Number of contigs to visualise. (default: 1000, 0 plots all contigs)
 * :code:`--contig_display_list` Comma-separated list of contig names to plot.

------------


GUNC MERGE_CHECKM
-----------------

Merge outputs of GUNC and checkM. Both should have been run on the same input files. CheckM qa should be run with :code:`-f qa.tsv -o 2 --tab_table` parameters. If run without :code:`-o 2` the extra columns will be empty.

Required Flags
^^^^^^^^^^^^^^

 * :code:`--gunc_file` Path of GUNC maxCSS output file (``GUNC.{db}.maxCSS_level.tsv``).
 * :code:`--checkm_file` CheckM output (qa.tsv) file (run :code:`checkm qa` with :code:`-o 2 --tab_table` parameters).

Optional Flags
^^^^^^^^^^^^^^

 * :code:`--out_dir` Output directory.  Default: Current working directory.

------------


GUNC RESCORE (summarise)
------------------------

Re-score genomes from a previous :code:`gunc run` using a different contamination cutoff. Both :code:`gunc rescore` and :code:`gunc summarise` are accepted; ``summarise`` is kept for backward compatibility. This is useful when you want to apply a stricter or more lenient threshold without re-running the full pipeline.

Example
^^^^^^^

.. code-block:: shell

    gunc summarise \
        -m GUNC.progenomes_2.1.maxCSS_level.tsv \
        -d gunc_output/ \
        -c 0.10 \
        -o GUNC.rescore_0.10.tsv

This re-scores all genomes from the previous run, marking those where no taxonomic level exceeds a contamination portion of 0.10 (at CSS > 0.45) as passing.

Required Flags
^^^^^^^^^^^^^^

 * :code:`-m / --max_csslevel_file` MaxCSS output file from a previous :code:`gunc run` (``GUNC.{db}.maxCSS_level.tsv``).
 * :code:`-d / --gunc_detailed_output_dir` Directory containing the per-genome all-levels detail files (the ``gunc_output`` directory from :code:`gunc run`).
 * :code:`-o / --output_file` Output file path for the rescored results.

Optional Flags
^^^^^^^^^^^^^^

 * :code:`-c / --contamination_cutoff` Contamination portion threshold to use for pass/fail. Default: 0.05.

.. note::

   ``gunc rescore`` uses a simplified re-evaluation: it selects the representative
   taxonomic level by filtering to levels where ``contamination_portion > cutoff`` first,
   then picking the highest CSS among those. This differs subtly from a full
   :code:`gunc run`, where the level with the globally highest CSS is always selected
   regardless of contamination portion. In edge cases, ``rescore`` can be more stringent
   than a full re-run. For definitive results at a different threshold, re-running the
   full pipeline is recommended.

------------


GUNC DOWNLOAD_DB
----------------

Required Flags
^^^^^^^^^^^^^^

 * :code:`positional argument` Download database to given directory.

Optional Flags
^^^^^^^^^^^^^^

 * :code:`--db` Which database to download. Options: ``progenomes_2.1`` (default), ``progenomes_3``, ``gtdb_95``, ``gtdb_214``.

------------


Special Flags
-------------

 * :code:`--version` Print version number and exit.
 * :code:`--help` Print help message and exit.

------------


Output Column Definitions
-------------------------

The main output file (``GUNC.{db}.maxCSS_level.tsv``) contains one row per genome at the taxonomic level with the highest clade separation score. Columns:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Column
     - Description
   * - ``genome``
     - Input genome filename (without extension).
   * - ``n_genes_called``
     - Number of genes called by Prodigal.
   * - ``n_genes_mapped``
     - Number of genes that mapped to the reference database.
   * - ``n_contigs``
     - Number of contigs in the genome.
   * - ``taxonomic_level``
     - The taxonomic level at which the maxCSS score was obtained.
   * - ``proportion_genes_retained_in_major_clades``
     - Fraction of mapped genes assigned to clades above the abundance cutoff.
   * - ``genes_retained_index``
     - Proportion of called genes that were retained (mapped / called).
   * - ``clade_separation_score``
     - Core GUNC metric (CSS). Measures how unevenly genes are distributed across clades at this taxonomic level. Range 0–1; values above 0.45 indicate chimerism.
   * - ``contamination_portion``
     - Fraction of mapped genes assigned to non-dominant clades (i.e. estimated contamination).
   * - ``n_effective_surplus_clades``
     - Effective number of extra clades beyond the dominant one (entropy-based).
   * - ``mean_hit_identity``
     - Mean amino acid identity of diamond hits to the reference.
   * - ``reference_representation_score``
     - How well the reference database covers this genome. Low values (< 0.3) mean the genome is poorly represented in the reference and the GUNC score is unreliable — treat such genomes as unscored regardless of ``pass.GUNC``.
   * - ``pass.GUNC``
     - ``True`` if the genome passes (CSS ≤ 0.45 or contamination portion ≤ 0.05); ``False`` if chimeric; ``nan`` if too few genes were mapped to score.

------------

Caveats and Limitations
-----------------------

**Horizontal gene transfer (HGT)**

GUNC cannot distinguish HGT islands from genuine contamination. Both result in genes
mapping to a secondary taxonomic clade, which raises the CSS and contamination portion.
Organisms with known extensive HGT — including members of Firmicutes, Thermotogae, and
Halobacteria — may generate false positive chimeric calls. If your genome belongs to one
of these groups, interpret a failing GUNC score alongside other quality metrics (e.g.
CheckM completeness/contamination, coverage uniformity) before discarding it.

**Low reference representation**

If ``reference_representation_score`` is below 0.3, the reference database contains few
close relatives of this genome. GUNC will still report a score, but it is unreliable —
the genome is effectively unscored. A warning is printed at the end of ``gunc run`` if
any genomes fall below this threshold.
