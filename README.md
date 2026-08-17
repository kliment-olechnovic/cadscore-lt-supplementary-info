# Supplementary materials for the CAD-score-LT publication

This repository contains materials supporting the CAD-score-LT publication: software source codes, as well as the benchmarking and example workflows, and the results of those workflows.

## Software source codes

The following subdirectories contain the source codes of the benchmarked software to be published:

* [software_source_codes/cadscore-lt_v0.9.196/](software_source_codes/cadscore-lt_v0.9.196/) contains CAD-score-LT version 0.9.196

* [software_source_codes/cadscorelt-extra-data-utilities-0.9.1/](software_source_codes/cadscorelt-extra-data-utilities-0.9.1/) contains extra data-preparation utilities for CAD-score-LT version 0.9.1

For future releases of CAD-score-LT the main access point is at [https://www.voronota.com/expansion_lt_cadscore/](https://www.voronota.com/expansion_lt_cadscore/).

## Benchmarking and example workflows

The following workflow directories cover the bencharking and examples presented in the CAD-score-LT publication:

* [supplementary_workflows/comparison_with_different_scores](supplementary_workflows/comparison_with_different_scores) covers
    - benchmarking CAD-score-LT against previous CAD-score implementations
    - benchmarking CAD-score-LT integrated chain remapping
    - benchmarking CAD-score-LT against other scores

* [supplementary_workflows/massivefold_clustering_example](supplementary_workflows/massivefold_clustering_example) covers
    - example application to a MassiveFold data set

* [supplementary_workflows/large_complex_clustering_example](supplementary_workflows/large_complex_clustering_example) covers
    - example application to a set of large models with complex stoichiometry

* [supplementary_workflows/mixed_complex_analysis](supplementary_workflows/mixed_complex_analysis) covers
    - example application to predictions of a protein-nucleic acid complex

* [supplementary_workflows/large_complex_clustering_example](supplementary_workflows/large_complex_clustering_example) covers
    - example application to protein-small molecule ligand docking models

Every workflow directory contains a ``workflow.bash`` script that defines environmental variables and runs the workflow scripts.

The workflow scripts are intended to work on Linux or in a Linux-like environment.

The workflow scripts automatically download and prepare most of the software, but
some software packages need to be installed manually or available beforehand:

* R (any version from last 20 years) required for all the workflows

* OpenStructure (recommended version 2.11.1) required for [supplementary_workflows/comparison_with_different_scores](supplementary_workflows/comparison_with_different_scores)

