#!/bin/bash

set -e

# Parse inputs and store parmeters in appropriate variables
while [ $# -gt 0 ]
do
	case $1 in
		--bam_filename)
		  shift
		  BAM_FILENAME=$1
		  shift
		  ;;
		--binned_bed_filename)
		  shift
		  BINNED_BED_FILENAME=$1
		  shift
		  ;;
		*)
			echo "ERROR: Unknown flag"
			exit 1
      ;;
	esac
done

SCRIPT_DIR=$(dirname "$0")
"${SCRIPT_DIR}/SubmitRescaleBinnedFiles.R" --bam_filename $BAM_FILENAME --binned_bed_filename $BINNED_BED_FILENAME
