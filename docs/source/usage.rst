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
 * :code:`gunc plot` Produce an interactive plot using the output from :code:`gunc run`
 * :code:`gunc merge_checkm` Produce a merged file combining the outputs of GUNC and checkM
 * :code:`gunc download_db` Download the GUNC database (required to run :code:`gunc run`)

Any of the above commands can be run with :code:`-h` to get function specific information.

------------

GUNC accepts either a progenomes or GTDB based reference database via the :code:`--db_file` option. Both can be downladed using the :code:`gunc download_db` command (see below). Note that using GTDB will lead to higher resource requirements and longer run times; in accuracy benchmarks, the performance of GTDB and the default proGenome-derived GUNC database performed very similarly.


GUNC RUN
--------

Run chimerism detection.

Required Flags
^^^^^^^^^^^^^^

 * :code:`--db_file` Path to the GUNC database file. Can be set as environment variable GUNC_DB.

One of the following is required. If flad :code:`--gene_calls` is not set gene calling will be done using prodigal with option "-p meta".

 * :code:`--input_dir` Input dir with files in FASTA format.
 * :code:`--input_file` Input file containing paths to FASTA format files.
 * :code:`--file_suffix` Only needed if suffix of input files is not the default .fa.
 * :code:`--input_fasta` Input file in FASTA format.

Optional Flags
^^^^^^^^^^^^^^

 * :code:`--gene_calls` Input is FASTA faa format genecalls.
 * :code:`--use_species_level` Allow species level to be picked as maxCSS. Default: False
 * :code:`--min_mapped_genes` Dont calculate GUNC score if number of mapped genes is below this value. Default: 11
 * :code:`--threads` Number of CPU threads.
 * :code:`--temp_dir` Directory to store temporary files. Default: Current working directory.
 * :code:`--out_dir` Directory in which to put output. Default: Current working directory.
 * :code:`--sensitive` Run with high sensitivity. (Uses a different cutoff to determine an abundant lineage)
 * :code:`--detailed_output` Output scores for every tax_level.

------------

GUNC PLOT
---------

Create interactive plot to visualise chimerism.

Required Flags
^^^^^^^^^^^^^^

 * :code:`--diamond_file` GUNC diamond outputfile. (one of the output files in `diamond_output` produced by :code:`gunc run`)

Optional Flags
^^^^^^^^^^^^^^

 * :code:`--gunc_gene_count_file` GUNC gene_counts.json file. (Not needed if `--diamond` file is in the file structure made by :code:`gunc run`)
 * :code:`--out_dir` Output directory.  Default: Current working directory.
 * :code:`--tax_levels` Tax levels to display (comma-seperated). (default: kingdom,phylum,family,genus,contig)
 * :code:`--remove_minor_clade_level` Tax level at which to remove minor clades. (default: kingdom)
 * :code:`--contig_display_num` Number of contigs to visualise. (default: 1000)

------------


GUNC MERGE_CHECKM
-----------------

Merge outputs of GUNC and checkM. Both should have been run on the same input files. CheckM qa should be run with `-f qa.tsv -o 2 --tab_table` parameters. If run without `-o 2` the extra columns will be empty.

Required Flags
^^^^^^^^^^^^^^

 * :code:`--gunc_file` Path of gunc_scores.tsv file.
 * :code:`--checkm_file` CheckM output (qa.tsv) file (run `checkm qa` with `-o 2 --tab_table` parameters).

Optional Flags
^^^^^^^^^^^^^^

 * :code:`--out_dir` Output directory.  Default: Current working directory.

------------


GUNC DOWNLOAD_DB
----------------

Required Flags
^^^^^^^^^^^^^^

 * :code:`positional argument` Download database to given directory.

Optional Flags
^^^^^^^^^^^^^^

 * :code:`--db` Which db to download (progenomes or gtdb). Default: progenomes

------------


Special Flags
-------------

 * :code:`--version` Print version number and exit.
 * :code:`--help` Print help message and exit.

