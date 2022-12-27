#!/bin/bash

set -e

source /broad/software/scripts/useuse
use Samtools
use .igvtools-2.4.16
use .r-3.6.0-bioconductor

cd ~/ProductionCode/

./scripts/binning/bin.sh $@
