#!/bin/bash

set -e

# use these commands for Broad cluster only
#source /broad/software/scripts/useuse
#use Samtools

# possible that will have multiple inputs so using a case structure
# Parse inputs and store parmeters in appropriate variables
while [ $# -gt 0 ]
do
	case $1 in
		-f)
			shift
			BAM=$1
			shift
			;;
		*)
			echo "ERROR: Unknown flag"
			exit 1
      ;;
	esac
done

# Check if BAM exists
if [ ! -f $BAM ]; then
	echo "ERROR: BAM not found"
	exit 2
fi

pairs=$(cat <(samtools view -H $BAM) <(samtools view $BAM | head -n 1) | samtools view -c -f 1)
echo $pairs

seq=$(samtools view $BAM | head -n 100 | cut -f 10)
LEN=${#seq}
echo $((LEN/100))

