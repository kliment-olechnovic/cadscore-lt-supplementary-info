#!/bin/bash

set -u

################################################################################

mkdir -p "${WORKROOTDIR}/input_model_files/"
cd "${WORKROOTDIR}/input_model_files/"

TNAME="H2424"
MAURL="https://predictioncenter.org/download_area/CASP17/predictions/oligo/${TNAME}.tar.gz"

ARCHIVEFILE="$(basename ${MAURL})"

if [ ! -d "$TNAME" ]
then
	[ -s "$ARCHIVEFILE" ] || wget -q "$MAURL"
	[ -s "$ARCHIVEFILE" ] || { echo >&2 "Error: failed to download ${MAURL}"; exit 1; }
	tar -xf "$ARCHIVEFILE"
	[ -d "$TNAME" ] || { echo >&2 "Error: failed to extract ${ARCHIVEFILE}"; exit 1; }
	rm -f "$ARCHIVEFILE"
fi

################################################################################

exit 0
