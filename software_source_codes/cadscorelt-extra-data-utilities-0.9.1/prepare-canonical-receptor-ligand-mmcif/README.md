# prepare-canonical-receptor-ligand-mmcif.py

`prepare-canonical-receptor-ligand-mmcif.py` is a command-line utility for preparing macromolecular structure files with consistent ligand atom naming.

The script:

* Reads a PDB or mmCIF structure.
* Identifies ligands (non-polymer, non-water residues).
* Removes all hydrogen atoms from the structure.
* Assigns deterministic canonical atom names to ligand atoms.
* Optionally merges ligand structures from SDF or MOL or MOL2 files.
* Writes the resulting structure in mmCIF format.
* If the ligand file contains multiple molecules, multiple output mmCIF files are produced.

The tool is designed to produce stable, reproducible ligand atom names that are independent of the input atom naming.

# Quick Reference

Convert PDB to canonicalized mmCIF:

```
python prepare-canonical-receptor-ligand-mmcif.py input_complex.pdb output_complex.cif
```

Add ligand:

```
python prepare-canonical-receptor-ligand-mmcif.py input_receptor.cif output_complex.cif --ligand-file input_ligand.sdf
```

Process docking poses:

```
python prepare-canonical-receptor-ligand-mmcif.py input_receptor.cif output_complex.cif --ligand-file input_poses.sdf
```

Verbose mode:

```
python prepare-canonical-receptor-ligand-mmcif.py input.cif output.cif --verbose
```

# Key Features

### Canonical ligand atom naming

Ligand atoms are renamed deterministically using RDKit canonical atom ranking.

Example:

```
C1
C2
C3
N1
O1
```

### Automatic ligand detection

Ligands are defined as residues that are:

* not protein
* not nucleic acid
* not water

Ions and small molecules are included.

### Hydrogen removal

All hydrogens are removed from:

* the receptor
* existing ligands
* externally supplied ligands

### Multi-model ligand input

If an SDF or MOL2 file contains multiple molecules, the script generates:

```
output_pose1.cif
output_pose2.cif
output_pose3.cif
```

### PDB numbering preservation

When converting from PDB input:

* original chain IDs are preserved
* original residue numbers are preserved
* mmCIF `label_*` fields are rewritten to match `auth_*` fields

This ensures compatibility with software expecting PDB-style numbering.

# Installation

## Python requirements

Python 3.8 pr later recommended.

Required packages:

```
gemmi
rdkit
```

### Install with pip

```
pip install gemmi rdkit
```

# Command-line Usage

Basic syntax:

```
python prepare-canonical-receptor-ligand-mmcif.py INPUT_STRUCTURE OUTPUT_CIF [options]
```

# Required Arguments

## Input structure

Supported formats:

* `.pdb`
* `.ent`
* `.cif`
* `.mmcif`

Example:

```
receptor.pdb
structure.cif
```

## Output path

Target mmCIF file path.

Example:

```
complex.cif
```

# Optional Arguments

## External ligand file

```
--ligand-file FILE
```

Supported formats:

```
.sdf
.mol
.mol2
```

Example:

```
--ligand-file ligand.sdf
```

## Ligand residue name

```
--ligand-resname NAME
```

Default:

```
LIG
```

Example:

```
--ligand-resname DRG
```

## Ligand chain identifier

```
--ligand-chain CHAIN
```

Default:

```
L
```

Example:

```
--ligand-chain X
```

## Ligand residue number

```
--ligand-seq NUMBER
```

Default:

```
1
```

Example:

```
--ligand-seq 501
```

## Model index

```
--model-index INDEX
```

Used when the structure contains multiple models.

Default:

```
0
```

Example:

```
--model-index 1
```

## Verbose output

```
--verbose
```

Prints detailed information about:

* residues processed
* hydrogen removal
* atom renaming

Example:

```
--verbose
```

# Basic Examples

## Example 1 — Canonicalize ligands in an existing structure

```
python prepare-canonical-receptor-ligand-mmcif.py input.cif output.cif
```

Result:

* ligand atoms renamed
* hydrogens removed
* output written as mmCIF


## Example 2 — Convert PDB to mmCIF

```
python prepare-canonical-receptor-ligand-mmcif.py receptor.pdb receptor.cif
```

Features:

* preserves chain IDs
* preserves residue numbers
* removes hydrogens

## Example 3 — Add a ligand from SDF

```
python prepare-canonical-receptor-ligand-mmcif.py receptor.cif complex.cif \
  --ligand-file ligand.sdf
```

Result:

* ligand inserted as residue `LIG`
* atom names canonicalized

## Example 4 — Use custom ligand identifiers

```
python prepare-canonical-receptor-ligand-mmcif.py receptor.cif complex.cif \
  --ligand-file ligand.sdf \
  --ligand-resname INH \
  --ligand-chain X \
  --ligand-seq 100
```

Resulting ligand:

```
Chain X
Residue INH 100
```

## Example 5 — Multi-model SDF input

```
python prepare-canonical-receptor-ligand-mmcif.py receptor.cif complex.cif \
  --ligand-file poses.sdf
```

If the SDF contains:

```
pose1
pose2
pose3
```

Outputs:

```
complex_pose1.cif
complex_pose2.cif
complex_pose3.cif
```

# Ligand Detection Rules

Residues are considered ligands if they are:

* not amino acid
* not nucleic acid
* not water

# Hydrogen Handling

Hydrogens are always removed from everywhere.

# Atom Naming Strategy

Ligand atoms are renamed using RDKit canonical atom ranking.

Procedure:

1. Construct RDKit molecule from atomic coordinates
2. Infer bonds
3. Compute canonical atom order
4. Assign element-based names

Example output:

```
C1
C2
C3
N1
O1
S1
```

Naming is deterministic but may not match PDB CCD atom names.

# Multi-model Ligand Input

### SDF

Each molecule in the file becomes a separate output structure.

### MOL2

Multiple molecules are supported if separated by:

```
@<TRIPOS>MOLECULE
```

### MOL

Single molecule only.

# Output File Naming

For multi-model ligand files:

```
output.cif
```

becomes

```
output_model1.cif
output_model2.cif
```

or

```
output_pose1.cif
output_pose2.cif
```

depending on molecule names.

Invalid characters in names are sanitized.

# Logging Output

Typical summary:

```
Input structure: receptor.cif
Input ligand file: poses.sdf
Ligand model: pose1
Output: complex_pose1.cif
Ligands processed: 3
Atoms renamed: 27
Hydrogens removed: 1482
```

With `--verbose`, additional residue-level messages appear.

# Special Cases

## No ligands present

The structure is written unchanged (except hydrogens removed).

## External ligand only

If the receptor has no ligand but a ligand file is supplied:

* the ligand is inserted
* canonical naming applied

## Multi-model SDF

Each ligand model produces a separate mmCIF.

## Invalid ligand file

The script stops with an error if:

* no molecules can be read
* ligand coordinates are missing

# Limitations

### Bond inference

Bond detection uses coordinate-based heuristics and may fail for:

* unusual metal complexes
* distorted geometries

Fallback naming ensures the script continues.

### Naming is not CCD-compatible

Atom names are deterministic but may not match:

```
PDB Chemical Component Dictionary
```

### No docking or alignment

External ligand coordinates are inserted directly without alignment.

# Exit Behavior

Normal completion:

```
exit code 0
```

Errors:

* unreadable ligand file
* invalid structure
* missing coordinates

