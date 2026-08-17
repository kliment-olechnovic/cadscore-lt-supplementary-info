#!/bin/bash

set -u

cd "$WORKROOTDIR"

###########################################################################

PDBID="$1"

TARGET_RECEPTOR_INFILE="${CASF2016_INPUT_DIR}/coreset/${PDBID}/${PDBID}_protein.pdb"
TARGET_LIGAND_INFILE="${CASF2016_INPUT_DIR}/coreset/${PDBID}/${PDBID}_ligand.mol2"
MODEL_LIGAND_INFILE="${CASF2016_INPUT_DIR}/decoys_docking/${PDBID}_decoys.mol2"
RMSD_INFILE="${CASF2016_INPUT_DIR}/decoys_docking/${PDBID}_rmsd.dat"

if [ ! -s "$TARGET_RECEPTOR_INFILE" ] || [ ! -s "$TARGET_LIGAND_INFILE" ] || [ ! -s "$MODEL_LIGAND_INFILE" ] || [ ! -s "$RMSD_INFILE" ] 
then
	echo "No files for ${PDBID}."
	exit 1
fi

ROOT_OUTDIR="./output/prepared_input_files/${PDBID}"

if [ -d "$ROOT_OUTDIR" ]
then
	echo "Skipping already prepared ${PDBID}."
	exit 0
fi

TARGET_OUTDIR="${ROOT_OUTDIR}/target"
mkdir -p "$TARGET_OUTDIR"

MODELS_OUTDIR="${ROOT_OUTDIR}/models"
mkdir -p "$MODELS_OUTDIR"

source ./workvenv/bin/activate

python3 "${CADSCORELT_EXTRAUTILS_DIR}/prepare-canonical-receptor-ligand-mmcif/prepare-canonical-receptor-ligand-mmcif.py" "$TARGET_RECEPTOR_INFILE" "${TARGET_OUTDIR}/target.cif" \
  --remove-ligands-from-receptor --ligand-file "$TARGET_LIGAND_INFILE" --ligand-resname LIG --ligand-chain x --ligand-seq 1

python3 "${CADSCORELT_EXTRAUTILS_DIR}/prepare-canonical-receptor-ligand-mmcif/prepare-canonical-receptor-ligand-mmcif.py" "$TARGET_RECEPTOR_INFILE" "${MODELS_OUTDIR}/model.cif" \
  --remove-ligands-from-receptor --ligand-file "$MODEL_LIGAND_INFILE" --ligand-resname LIG --ligand-chain x --ligand-seq 1

TARGET_OUTFILE="$(find ${TARGET_OUTDIR} -type f)"

python3 "${CADSCORELT_EXTRAUTILS_DIR}/group-equivalent-ligand-atoms/group-equivalent-ligand-atoms.py" "$TARGET_OUTFILE" > "${ROOT_OUTDIR}/equivalent_atoms.txt"

{
echo "model rmsd"
cat "$RMSD_INFILE" | tail -n +2 | awk '{print "model_" $1 ".cif " $2}'
} \
| sed 's/\s\+/\t/g' \
> "${ROOT_OUTDIR}/rmsd.tsv"

###########################################################################

exit 0
