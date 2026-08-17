#!/bin/bash

set -u

################################################################################

mkdir -p "$INPUT_DATA_FILES_DIR"
cd "$INPUT_DATA_FILES_DIR"

cat "${WORKROOTDIR}/input_data_lists/urls_of_model_archives.txt" \
| while read -r MAURL
do
	ARCHIVEFILE="$(basename ${MAURL})"
	TNAME="$(basename ${ARCHIVEFILE} .tar.gz)"
	SEQFILE="sequences/${TNAME}.fasta"
	
	if [ ! -d "$TNAME" ]
	then
		[ -s "$ARCHIVEFILE" ] || wget -q "$MAURL"
		[ -s "$ARCHIVEFILE" ] || { echo >&2 "Error: failed to download ${MAURL}"; exit 1; }
		tar -xf "$ARCHIVEFILE"
		[ -d "$TNAME" ] || { echo >&2 "Error: failed to extract ${ARCHIVEFILE}"; exit 1; }
		rm -f "$ARCHIVEFILE"
	fi
	
	if [ ! -s "$SEQFILE" ]
	then
		mkdir -p "$(dirname ${SEQFILE})"
		SEQURL="$(echo "https://predictioncenter.org/casp16/target.cgi?target=${TNAME}&view=sequence" | sed 's|o\&view|\&view|')"
		curl -s "$SEQURL" > "$SEQFILE"
	fi
done

################################################################################

exit 0
