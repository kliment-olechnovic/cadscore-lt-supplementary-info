# Supplementary materials for the CAD-score-LT publication

This repository contains materials supporting the CAD-score-LT publication:

* software source code
* benchmarking and example workflows, together with the results produced by those workflows

## Software source code

The following subdirectories contain the source code of the software versions benchmarked in the publication:

* [software_source_codes/cadscore-lt_v0.9.196/](software_source_codes/cadscore-lt_v0.9.196/) contains CAD-score-LT version 0.9.196

* [software_source_codes/cadscorelt-extra-data-utilities-0.9.1/](software_source_codes/cadscorelt-extra-data-utilities-0.9.1/) contains version 0.9.1 of the extra data-preparation utilities for CAD-score-LT

The main access point for future releases of CAD-score-LT is https://www.voronota.com/expansion_lt_cadscore/.

## Benchmarking and example workflows

The following workflow directories cover the benchmarking and examples presented in the CAD-score-LT publication:

* [supplementary_workflows/comparison_with_different_scores](supplementary_workflows/comparison_with_different_scores) covers:
    - benchmarking CAD-score-LT against previous CAD-score implementations
    - benchmarking CAD-score-LT integrated chain remapping
    - benchmarking CAD-score-LT against other scores

* [supplementary_workflows/massivefold_clustering_example](supplementary_workflows/massivefold_clustering_example) covers:
    - an example application to a MassiveFold data set

* [supplementary_workflows/large_complex_clustering_example](supplementary_workflows/large_complex_clustering_example) covers:
    - an example application to a set of large models with complex stoichiometry

* [supplementary_workflows/mixed_complex_analysis](supplementary_workflows/mixed_complex_analysis) covers:
    - an example application to predictions of a protein--nucleic acid complex

* [supplementary_workflows/protein_ligand_docking_analysis](supplementary_workflows/protein_ligand_docking_analysis) covers:
    - an example application to protein--small-molecule ligand docking models

Each workflow directory contains a `workflow.bash` script that defines environment variables and runs the workflow scripts.

The workflow scripts are intended to run on Linux or in a Linux-like environment.

The workflow scripts automatically download and prepare most of the required software. However, some software packages need to be installed manually or made available beforehand:

* R (any version released within the last 20 years) is required for all workflows

* OpenStructure (recommended version 2.11.1) is required for [supplementary_workflows/comparison_with_different_scores](supplementary_workflows/comparison_with_different_scores)
