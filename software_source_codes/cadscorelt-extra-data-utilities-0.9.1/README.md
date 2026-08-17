# Extra data-preparation utilities for CAD-score-LT

This repository contains tools that can be utilized in preparing protein-ligand complex data for comparison with [CAD-score-LT](https://www.voronota.com/expansion_lt_cadscore/index.html).

The available tools are:

* [prepare-canonical-receptor-ligand-mmcif.py](./prepare-canonical-receptor-ligand-mmcif/README.md) is a command-line utility for preparing molecular structure files with consistent ligand atom naming.

* [group-equivalent-ligand-atoms.py](./group-equivalent-ligand-atoms/README.md) is a command-line utility that reads a molecular structure from file, identifies unique ligands, computes RDKit atom symmetry classes, and outputs standardized equivalence-based atom names.

