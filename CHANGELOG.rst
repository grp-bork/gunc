=========
Changelog
=========

v1.0.5
------

Features
^^^^^^^^
 - Added option to plot specific contigs
 - Added option to plot ALL contigs

Bugfixes
^^^^^^^^
 - Handle cases in which no genes are called with a more useful error message.
 - Added check if temp_dir exists before starting analysis.

Other
^^^^^
 - Removed dependency on `zgrep` (for compatability with nf-core tests)

v1.0.4
------

Features
^^^^^^^^
 - Added `contig_taxonomy_output` option to output detailed taxonomy assignment count per contig.

Bugfixes
^^^^^^^^
 - Fix version of dependancy in conda recipe: requests>=2.22.0

v1.0.3
------


Bugfixes
^^^^^^^^
 - Running with genecalls as iput failed
 - GUNC plot contig_display_num displayed a defined number of genes instead of contigs


v1.0.2
------


Features
^^^^^^^^
 - GUNC can now be run with GTDB database
 - Added option to download GTDB_GUNC database
 - Input file options can be gene_calls (faa) instead of fna if `--gene_calls` flag is set
 - Input genecalls can be gzipped
 - Output maxCSS file is now sorted


Bugfixes
^^^^^^^^
 - Fix version of dependancy: requests>=2.22.0 (older versions not compatable)
 - Better error message if gunc_db does not exist
 - `checkm_merge` didnt work with unless checkm qa was run with `-o 2`

Other
^^^^^
 - Documentation updates
    - Links to synthetic datasets added
    - Citations for diamond and prodigal added
    - Clarified how checkM should be run for `checkm_merge`
    - Corrected command for `download_db`
 - Check if fasta is given with `-f` option instead of list of filepaths


v1.0.1
------

Bugfixes
^^^^^^^^
 - Running from genecounts failed
 - Fixed case where pass.GUNC output was converted to ints
 - Fixed silently ignoring input samples that did not map to reference

Other
^^^^^
 - Better error message if ouput_dir doesnt exist
 - Documentation updates


v1.0.0
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
 - Inputs that dont map any genes to GUNC_DB were silently missing in output

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
