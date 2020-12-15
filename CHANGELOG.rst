=========
Changelog
=========

v0.1.3
------

Features
^^^^^^^^
 - Added option to download the GUNC_DB file
 - Added option to merge GUNC output with checkM output
 - Added option to create interactive HTML based visualisation
 - Added option to run all fastas in a directory
 - Added option to provide input filepaths in a file
 - Added min_mapped_genes option so scores are not calculated when there are not enough genes
 - Added use_species_level option for determining tax_level with maxCSS score
 - Can now accept gzipped fna files (with .gz ending)
 - Allow GUNC_DB to be supplied using an env var
 - Updated arguments to a subcommand structure
 - Complete rewrite of how scores are calculated
 - Gene calling is now done in parallel

Bugfixes
^^^^^^^^
 - genome2taxonomy was not included in pip package
 - GUNC failed if nothing left after minor clade filtering
 - If duplicate filenames were in input, output files were overwritten

Other
^^^^^
 - Documentation updates
 - sklearn dependency removed
 - Added the bioconda recipe to repo
 - Added check for zgrep, prodigal and diamond
 - Changed output names to match those in paper
 - Fixed diamond version to 2.0.4 (needs to be compatable with GUNC_DB)
 - Better quality LOGOs
 - Diamond logs are silenced
 - Timestamps added to log output


Initial Release v0.1.2 (2020-10-14)
-----------------------------------

 - First release
