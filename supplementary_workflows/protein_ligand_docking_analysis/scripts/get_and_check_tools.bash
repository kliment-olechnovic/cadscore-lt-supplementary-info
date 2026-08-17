#!/bin/bash

set -u

################################################################################

mkdir -p "${WORKROOTDIR}/tools"
cd "${WORKROOTDIR}/tools"

################################################################################

CADSCORELT_TGZFILE="cadscore-lt_v${CADSCORELT_VERSION}.tar.gz"

[ -s "$CADSCORELT_TGZFILE" ] || wget "https://github.com/kliment-olechnovic/voronota/releases/download/v${VORONOTA_VERSION}/${CADSCORELT_TGZFILE}"

[ -d "$CADSCORELT_DIR" ] || tar -xf "$CADSCORELT_TGZFILE"

[ -d "$CADSCORELT_DIR" ] || { echo >&2 "Error: failed to prepare ${CADSCORELT_DIR}"; exit 1; }

################################################################################


if [ ! -s "${CADSCORELT_DIR}/cadscore-lt" ]
then
	cd "$CADSCORELT_DIR"
	g++ -std=c++17 -Ofast -march=native -fopenmp -I ./src -o ./cadscore-lt ./src/cadscore_lt.cpp
	cd - >/dev/null
fi

[ -s "${CADSCORELT_DIR}/cadscore-lt" ] || { echo >&2 "Error: failed to build ${CADSCORELT_DIR}/cadscore-lt"; exit 1; }

################################################################################

mkdir -p "${WORKROOTDIR}/tools/extras"
cd "${WORKROOTDIR}/tools/extras"

UTILS_ARCH_TGZFILE="v${CADSCORELT_EXTRAUTILS_VERSION}.tar.gz"

[ -s "$UTILS_ARCH_TGZFILE" ] || wget "https://github.com/kliment-olechnovic/cadscorelt-extra-data-utilities/archive/refs/tags/v${CADSCORELT_EXTRAUTILS_VERSION}.tar.gz"

[ -d "$CADSCORELT_EXTRAUTILS_DIR" ] || tar -xf "$UTILS_ARCH_TGZFILE"

[ -d "$CADSCORELT_EXTRAUTILS_DIR" ] || { echo >&2 "Error: failed to prepare extra utilities"; exit 1; }

################################################################################

exit 0

