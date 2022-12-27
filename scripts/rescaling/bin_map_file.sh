#!/bin/bash
# make sure to be running this from ProductionCode folder
export PATH=${PATH}:./tools
set -e

MAPURL=F
# Parse inputs and store parmeters in appropriate variables
while [ $# -gt 0 ]
do
	case $1 in
		--mapURL)
		  shift
		  MAPURL=$1
		  shift
		  ;;
		--mapFile)
		  shift
		  MAPFILE=$1
		  shift
		  ;;
		--binSize)
		  shift
		  BINSIZE=$1
		  shift
		  ;;
		--referenceFile)
		  shift
		  REF=$1
		  shift
		  ;;
		--chromSizes)
		  shift
		  CHROMSIZES=$1
		  shift
		  ;;
		--saveFilename)
		  shift
		  SAVENAME=$1
		  shift
		  ;;
		*)
			echo "ERROR: Unknown flag"
			exit 1
      ;;
	esac
done

TMP_DIR=$(mktemp -dp .)

if [ "$MAPURL" != "F" ]; then
  # Download map file to TMP_DIR
  echo 'Dowloading mappability track'
  wget -nv -P "${TMP_DIR}" ${MAPURL}
  # update location of MAPFILE
  MAPFILE="${TMP_DIR}/${MAPFILE}"
fi

# does it have a gz in the filename?
if [[ $MAPFILE == *.gz ]]; then
  gunzip $MAPFILE
  MAPFILE=${MAPFILE//.gz/}
fi

# convert bedgraph to bed
if [[ $MAPFILE == *.bedgraph ]]; then
  NEWFILE=${MAPFILE//graph/}
  mv $MAPFILE $NEWFILE
  MAPFILE=$NEWFILE
fi

# convert bedgraph to bigwig first
if [[ $MAPFILE == *.bed ]]; then
  echo "this is a bedgraph file"
  WIGFILE="${TMP_DIR}/map_count.wig"
  BIGWIGFILE="${TMP_DIR}/map_out_temp.bw"
  igvtools count -w ${BINSIZE} "${MAPFILE}" "${WIGFILE}" "${CHROMSIZES}"
  wigToBigWig "${WIGFILE}" "${CHROMSIZES}" "${BIGWIGFILE}"
  MAPFILE="${BIGWIGFILE}"
fi

#then make map file
paste $REF <(bigWigAverageOverBed $MAPFILE $REF stdout | cut -f 5) | \
awk 'BEGIN{FS="\t"; OFS="\t"}; {print $1, $2, $3, $7}' >$SAVENAME

rm -rf "${TMP_DIR}"
