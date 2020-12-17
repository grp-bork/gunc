
.. image:: GUNC_LOGO.svg
    :width: 400px
    :align: left
    :alt: GUNC

|
|
|
|
|

GUNC Documentation
==================

.. toctree::
   :maxdepth: 1

   installation
   usage
   output
   functions
   changelog
   contributing
   datasets
   authors
   citations

.. image:: https://img.shields.io/pypi/v/gunc.svg
        :target: https://pypi.python.org/pypi/gunc
.. image:: https://anaconda.org/bioconda/gunc/badges/version.svg
        :target: https://anaconda.org/bioconda/gunc
.. image:: https://anaconda.org/bioconda/gunc/badges/downloads.svg
        :target: https://anaconda.org/bioconda/gunc
.. image:: https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400
        :target: https://choosealicense.com/licenses/gpl-3.0/
.. image:: https://anaconda.org/bioconda/gunc/badges/installer/conda.svg
        :target: https://conda.anaconda.org/bioconda



* Free software: GNU General Public License v3 or later

Genome UNClutterer (GUNC) is a tool for detection of chimerism and contamination in prokaryotic genomes resulting from mis-binning of genomic contigs from unrelated lineages. It does so by applying an entropy based score on taxonomic assignment and contig location of all genes in a genome.

.. figure:: /_static/GUNC_PLOT_example.png
  :width: 100%
  :alt: Example of GUNC visualisation function.
  :align: center

  Example of GUNC visualisation.
