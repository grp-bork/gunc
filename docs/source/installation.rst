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

Manual Installation
-------------------

Alternatively gunc can be install using pip::

    $ pip install gunc

GUNC DB should be downloaded by following the above instructions.
Separately prodigal and diamond 2.1.4 need to be installed


Containers
----------

Both Singularity and docker containers are available by following instructions here: `LINK <https://biocontainers.pro/#/tools/gunc>`_
