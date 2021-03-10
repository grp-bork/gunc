============
Installation
============

The recommended method of installation is through `bioconda <https://anaconda.org/bioconda/gunc>`_ ::

    $ conda install -c bioconda gunc

Download GUNC DB
----------------
Before use the GUNC DB needs to be downloaded (13G)::

    $ gunc download_db /path/to/output/dir/

The path to the GUNC_DB can be supplied as a command line argument
to :code:`gunc run` or set as an environment variable: :code:`GUNC_DB`

Manual Installation
-------------------

Alternatively gunc can be install using pip::

    $ pip install gunc

GUNC DB should be downloaded by following the above instructions.
Separately prodigal and diamond 2.0.4 need to be installed


Containers
----------

Both Singularity and docker containers are available by following instructions here: `LINK <https://biocontainers.pro/#/tools/gunc>`_
