#!/usr/bin/env bash

set -e

TOOLS_DIR=${1:-tools}

mkdir -p "${TOOLS_DIR}"
cd "${TOOLS_DIR}"

URL_PREFIX="https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64"
wget -q "${URL_PREFIX}/bigWigAverageOverBed"
wget -q "${URL_PREFIX}/wigToBigWig"
chmod +x *
