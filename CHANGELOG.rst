=========
Changelog
=========

v1.0.7
------

Summary
^^^^^^^

This release adds support for two new reference databases (ProGenomes 3, GTDB r214) and
a custom database option. A new ``gunc check`` subcommand validates your environment before
submitting a long job, and ``gunc rescore`` is introduced as a clearer alias for 
``gunc summarise``. A warning is now emitted when genomes have low reference representation scores.
Packaging has been modernised to ``pyproject.toml`` and the CI pipeline updated.

Features
^^^^^^^^
 - Added support for progenomes_3 and gtdb_214 reference databases.
 - Added ``--custom_genome2taxonomy`` option to allow use of a custom reference database.
 - Diamond version requirement removed.
 - Added ``gunc rescore`` as the preferred name for the ``summarise`` subcommand; ``gunc summarise`` remains as a backward-compatible alias.
 - Added ``gunc check`` subcommand to validate environment (tool dependencies, database file, custom genome-to-taxonomy TSV format, output directory write access) without running the pipeline.
 - All subcommands (``run``, ``plot``, ``merge_checkm``, ``summarise``) now log the output file path on completion.
 - ``--file_suffix`` error message now suggests the correct flag usage when no files are found.
 - Fixed ``metavar="\\b"`` hack in ``summarise`` argparse definitions; replaced with meaningful placeholders (``FILE``, ``DIR``, ``FLOAT``).
 - Documentation: added ``gunc summarise`` section with worked example; fixed ``--file_suffix`` incorrectly listed as required; fixed ``--gunc_file`` help referencing ``gunc_scores.tsv`` (actual filename is ``GUNC.{db}.maxCSS_level.tsv``); added ``--custom_genome2taxonomy`` file format spec; added output column definitions table; updated DB names to underscore convention throughout.

Bugfixes
^^^^^^^^
 - Fixed ``summarise`` subcommand incorrectly marking all genomes as passing GUNC.
 - Fixed ``pass.GUNC`` column being silently converted to strings in output TSV; ``summarise`` now uses proper NaN detection instead of string comparison.
 - Fixed ``summarise`` not rescoring genomes with boolean ``False`` in ``pass.GUNC``; previously only the string ``"False"`` was matched, so boolean values (the normal case) were silently skipped.
 - Fixed genome identity corruption in ``split_diamond_output`` when contig names contain ``/``; now uses ``rsplit`` to always extract the genome name from the last path segment.
 - Fixed DB detection logic duplicated across three code paths with subtly different ordering; extracted into single ``detect_db_from_filename()`` function.
 - Fixed ``prodigal()`` leaving partial output files on disk when gene calling fails; partial files are now removed so the caller's size check correctly excludes failed genomes.
 - Fixed ``extract_node_data()`` in visualisation missing colour entries for ``class`` and ``order`` tax levels, causing ``KeyError`` when non-default ``--tax_levels`` are used.
 - Extracted ``plot=True`` path from ``chim_score()`` into dedicated ``get_base_data_for_plotting()`` function; ``chim_score()`` now has a single consistent return type.
 - Fixed empty diamond output files not being named correctly when a genome fails to map ( thanks to @pamelaferretti ).
 - Fixed edge case where contamination score was incorrectly calculated when contamination portion was NaN.
 - Fixed crash when no genes were called or mapped to the reference database.
 - Fixed shell injection risk in ``get_record_count_in_fasta``.

Other
^^^^^
 - Removed versioneer; version is now statically set.
 - Fixed 8 flake8 errors: import ordering in ``get_scores.py`` and ``visualisation.py``, trailing whitespace in ``gunc.py``, spurious ``f``-string prefixes in ``gunc_database.py``.
 - Extracted ``CSS_CHIMERIC_THRESHOLD = 0.45`` and ``TAX_LEVELS`` as named constants in ``get_scores.py``; replaced all three scattered hardcoded copies of the threshold and tax level list across ``gunc.py``, ``checkm_merge.py``, and ``visualisation.py``.
 - Fixed all ``sys.exit(string)`` calls in ``visualisation.py`` and ``get_scores.py`` to use ``logger.error()`` + ``sys.exit(1)`` consistently with the rest of the codebase; added module-level logger to ``get_scores.py``.
 - Fixed ``add_empty_diamond_output()`` using ``print()`` for progress output; now uses ``logger.info()``.
 - Fixed ``check_diamond_version()`` using ``shell=True``; now uses list-form subprocess call.
 - Added guard against empty ``gunc_output`` list before ``pd.concat()`` in ``run_gunc()`` to give a clear error instead of a cryptic ``ValueError``.
 - Reference data files renamed to reflect database version (e.g. ``genome2taxonomy_pg2.1ref.tsv``).
 - Documentation updated: diamond version, all four database options, ``--custom_genome2taxonomy`` flag.
 - Migrated packaging from ``setup.py`` + ``setup.cfg`` + ``MANIFEST.in`` + ``requirements.txt`` to a single ``pyproject.toml`` (PEP 621); fixed ``package_data`` paths, license field (GPLv3), dropped ``universal=1``, and added minimum version pins for numpy (>=1.20), scipy (>=1.7), and plotly (>=5.0).
 - Replaced all ``from module import *`` in test files with explicit named imports; marked network-dependent tests in ``test_gunc_database.py`` with ``@pytest.mark.integration``; added ``conftest.py`` registering the ``integration`` marker.
 - Added tests for ``summarise()``, ``get_scores_using_supplied_cont_cutoff()``, ``read_genome2taxonomy_reference()`` (all 4 DBs + custom + unknown), ``split_diamond_output()`` round-trip, and ``detect_db_from_filename()``.


v1.0.6
------

Features
^^^^^^^^
 - GUNC will now report the number of failed genomes if some fail when running in batch mode.

Bugfixes
^^^^^^^^
 - Fixed cases where number of mapped genes might be different when running smaller genomes individually vs grouped.
 - Fixed cases where gunc would fail if _-_ was in a sequence identifier.
 - Fixed cases where contig colour in visualisation was incorrect.
 - Fixed incorrect ordering of taxonomic levels in plot if user supplies them in an incorrect order.
 - Allow GUNC to continue if some genomes fail diamond mapping.
 - Fixed crash if diamond fails to map anything to GUNC reference.
 - Fixed bug where GUNC goes into infinite loop if an input filename is merged.fa

Other
^^^^^
 - class and order are now included in contig assignment output file (@aaronmussig)
 - Fix pandas deprecation waring (@tmaklin)
 - Added fixes for pandas FutureWarning (changing behaviour of Series.idxmax)


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
