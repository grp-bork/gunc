=====
Usage
=====

Quickstart
----------

.. code-block:: shell

    gunc -i genome.fa -d /path/to/db

Required Flags
--------------

One of the following it required. If contigs are supplied the gene calls will be done using prodigal.

 * :code:`--input_file` Input file in FASTA fna format.
 * :code:`--gene_calls` Input genecalls FASTA faa format.

Optional Flags
--------------

 * :code:`--database_file` Diamond database reference file. Optionally can be set using an env var: GUNC_DB
 * :code:`--threads` Number of CPU threads.
 * :code:`--temp_dir` Directory to store temporary files. Default: Current working directory.
 * :code:`--out_dir` Directory in which to put output. Default: Current working directory.
 * :code:`--sensitive` Run with high sensitivity. (Uses a different cutoff to determine an abundant lineage)

Special Flags
-------------

 * :code:`--config` Config file path. All options can optionally be provided using a config file. The config file can be located at :code:`~/.gunc` or provided via command line with the :code:`--config` option

 * :code:`--download_db` Download database to given direcory.
 * :code:`--version` Print version number and exit.
 * :code:`--help` Print hekp message and exit.

