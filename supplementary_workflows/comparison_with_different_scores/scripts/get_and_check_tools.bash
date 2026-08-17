#!/bin/bash

set -u

################################################################################

mkdir -p "${WORKROOTDIR}/tools"
cd "${WORKROOTDIR}/tools"

OLDCADSCORE_TGZFILE="cadscore_${OLDCADSCORE_VERSION}.tar.gz"
VORONOTA_TGZFILE="voronota_${VORONOTA_VERSION}.tar.gz"
CADSCORELT_TGZFILE="cadscore-lt_v${CADSCORELT_VERSION}.tar.gz"

[ -s "$OLDCADSCORE_TGZFILE" ] || wget "https://github.com/kliment-olechnovic/old_cadscore/releases/download/${OLDCADSCORE_VERSION}/${OLDCADSCORE_TGZFILE}"
[ -s "$VORONOTA_TGZFILE" ] || wget "https://github.com/kliment-olechnovic/voronota/releases/download/v${VORONOTA_VERSION}/${VORONOTA_TGZFILE}"
[ -s "$CADSCORELT_TGZFILE" ] || wget "https://github.com/kliment-olechnovic/voronota/releases/download/v${VORONOTA_VERSION}/${CADSCORELT_TGZFILE}"

[ -d "$OLDCADSCORE_DIR" ] || tar -xf "$OLDCADSCORE_TGZFILE"
[ -d "$VORONOTA_DIR" ] || tar -xf "$VORONOTA_TGZFILE"
[ -d "$CADSCORELT_DIR" ] || tar -xf "$CADSCORELT_TGZFILE"

[ -d "$OLDCADSCORE_DIR" ] || { echo >&2 "Error: failed to prepare ${OLDCADSCORE_DIR}"; exit 1; }
[ -d "$VORONOTA_DIR" ] || { echo >&2 "Error: failed to prepare ${VORONOTA_DIR}"; exit 1; }
[ -d "$CADSCORELT_DIR" ] || { echo >&2 "Error: failed to prepare ${CADSCORELT_DIR}"; exit 1; }

################################################################################

if [ ! -s "${OLDCADSCORE_DIR}/bin/voroprot2" ]
then
	cd "$OLDCADSCORE_DIR"
	g++ -O3 -o bin/voroprot2 src/*.cpp
	cd - >/dev/null
fi

if [ ! -s "${VORONOTA_DIR}/voronota" ]
then
	cd "$VORONOTA_DIR"
	g++ -O3 -std=c++11 -o voronota $(find ./src/ -name '*.cpp')
	cd - >/dev/null
fi

if [ ! -s "${VORONOTA_DIR}/expansion_js/voronota-js" ]
then
	cd "${VORONOTA_DIR}/expansion_js"
	g++ -std=c++14 -I"./src/dependencies" -O3 -o "./voronota-js" $(find ./src/ -name '*.cpp')
	cd - >/dev/null
fi

if [ ! -s "${CADSCORELT_DIR}/cadscore-lt" ]
then
	cd "$CADSCORELT_DIR"
	g++ -std=c++17 -Ofast -march=native -fopenmp -I ./src -o ./cadscore-lt ./src/cadscore_lt.cpp
	cd - >/dev/null
fi

[ -s "${OLDCADSCORE_DIR}/bin/voroprot2" ] || { echo >&2 "Error: failed to build ${OLDCADSCORE_DIR}/bin/voroprot2"; exit 1; }
[ -s "${VORONOTA_DIR}/voronota" ] || { echo >&2 "Error: failed to build ${VORONOTA_DIR}/voronota"; exit 1; }
[ -s "${VORONOTA_DIR}/expansion_js/voronota-js" ] || { echo >&2 "Error: failed to build ${VORONOTA_DIR}/expansion_js/voronota-js"; exit 1; }
[ -s "${CADSCORELT_DIR}/cadscore-lt" ] || { echo >&2 "Error: failed to build ${CADSCORELT_DIR}/cadscore-lt"; exit 1; }

################################################################################

if [[ ! -s "$OPENSTRUCTURE_APP" ]] || [[ "$(${OPENSTRUCTURE_APP} --version)" != "OpenStructure"* ]]
then
	echo >&2 "Error: OpenStructure application not available at ${OPENSTRUCTURE_APP}."
	exit 1
fi

################################################################################

exit 0

