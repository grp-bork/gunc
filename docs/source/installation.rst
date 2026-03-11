============
Installation
============

The recommended method of installation is through `bioconda <https://anaconda.org/bioconda/gunc>`_ ::

    $ conda install -c bioconda gunc

Download GUNC DB
----------------
Before use the GUNC DB needs to be downloaded. The default progenomes2.1 DB (~13G) can be downloaded with::

    $ gunc download_db /path/to/output/dir/

Alternative databases (progenomes3, gtdb95, gtdb214) can be downloaded with the :code:`--db` flag::

    $ gunc download_db /path/to/output/dir/ --db gtdb214

The path to the GUNC_DB can be supplied as a command line argument
to :code:`gunc run` or set as an environment variable: :code:`GUNC_DB`

Verifying your installation
----------------------------

A minimal test dataset (a tiny diamond database and two small test genomes)
can be downloaded to verify that your GUNC installation is working correctly::

    $ gunc download_db /path/to/test_dir/ --db test_data

This downloads four files:

- :code:`ci_test.dmnd` — minimal diamond database
- :code:`chimeric.faa` — gene calls from a chimeric genome (expected: ``pass.GUNC=False``)
- :code:`clean.faa` — gene calls from a clean genome (expected: ``pass.GUNC=True``)
- :code:`genome2taxonomy.tsv` — taxonomy reference for the test database

Run both test genomes together with::

    $ TEST_DIR=/path/to/test_dir
    $ gunc run --gene_calls \
        --input_dir ${TEST_DIR} \
        --file_suffix .faa \
        --db_file ${TEST_DIR}/ci_test.dmnd \
        --custom_genome2taxonomy ${TEST_DIR}/genome2taxonomy.tsv \
        --out_dir ./gunc_test_out

.. note::
   Use absolute paths (or a variable as shown above) rather than relative paths
   to avoid path resolution issues.

:code:`--input_dir` with :code:`--file_suffix .faa` picks up only ``chimeric.faa``
and ``clean.faa``, ignoring the ``.dmnd`` and ``.tsv`` files in the same directory.

In the output TSV check that ``chimeric`` has :code:`pass.GUNC=False` and
``clean`` has :code:`pass.GUNC=True`.

Manual Installation
-------------------

Alternatively gunc can be install using pip::

    $ pip install gunc

GUNC DB should be downloaded by following the above instructions.
Separately prodigal and diamond 2.1.24 need to be installed.
The required diamond version is enforced at startup; set the environment variable
``GUNC_SKIP_DIAMOND_VERSION_CHECK=1`` to bypass this check if needed.


Containers
----------

Both Singularity and docker containers are available by following instructions here: `LINK <https://biocontainers.pro/#/tools/gunc>`_
