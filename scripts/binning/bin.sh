#!/usr/bin/env bash
# bin.sh bins the reads of a BAM file
# It uses the actual fragment length for paired reads (estimates 200 bp for single-end)
# Automatically detects BAM read length and if single/paired (assumes homogeneity of BAM)
# Currently supported genomes: mm9, mm10, hg18, hg19, hg38
# Converts BAM to WIG to BigWig using igvtools count
# BigWig is averaged over a reference file to create a BED
# Expects the full path of to the BAM file (-f)
# User specifies bin size (-w [default=5000])
# User specifies genome (-g [default=hg19])
# User specifies modal length (-r [default=200])
# User specifies output filename (-n [default=output])
# User can specify suffix (-s [default=binned])

usage(){
	echo "usage: $(basename $0) [-f infile] [-g genome] [-w binSize]"
	echo -e "Bins the reads of a BAM file\n\
It uses the actual fragment length for paired reads (estimates 200 bp for single-end)\n\
Automatically detects BAM read length and if single/paired (assumes homogeneity of BAM)\n\
Currently supported genomes: mm9, mm10, hg18, hg19, hg38\n\
Converts BAM to WIG to BigWig using igvtools count which implicitly filters these flags: 4,512,1024 \n\
BigWig is averaged over a reference file to create a bedgraph"
	echo -e "-f, -file \e[4minfile\e[0m\n\t Path to input BAM"
	echo -e "-g, -genome \e[4mgenome\e[0m\n\t Reference genome (default: hg19)"
	echo -e "-w, -window-size \e[4mbinSize\e[0m\n\t Size of bins (default: 5000)"
	echo -e "-r, -modal-length \e[4mmodalLen\e[0m\n\t Estimated length of fragment for single-end read (default: 200)"
	echo -e "-n, -name \e[4mprefix\e[0m\n\t Optional prefix for output file name (default: output)"
	echo -e "-s, -suffix \e[4msuffix\e[0m\n\t Optional suffix for output file name (default: binned)"
}

# set -e # Exit on error
set -u # Exit if undeclared variable is used
# set -x # Trace executions
set -o pipefail # Catch pipe fails

# Defaults
bin=5000
gen="hg19"
modalLen=200
prefix="output"
suffix="binned"

# Parse inputs and store parmeters in appropriate variables
while [ $# -gt 0 ]
do
	case $1 in
		-h | --help)
			usage
			exit
			;;
		-f | -file)
			shift
			bam="$1"
			;;
		-w | -window-size)
			shift
    		bin="$1"
			;;
		-g | -genome)
			shift
  			gen="$1"
			;;
		-n | -name)
			shift
			prefix="$1"
			;;
		-r | -modal-length)
			shift
			modalLen="$1"
			;;
		-s | -suffix)
			shift
			suffix="$1"
			;;
		*)
			echo "ERROR: Unknown flag"
			echo "Run $(basename $0) --help for usage instructions"
			exit 1
			;;
	esac
	shift
done

SCRIPT_DIR=$(dirname "$0")
REFS_DIR="${SCRIPT_DIR}/../../references"

# Check if BAM exists
if [ ! -f "$bam" ]; then
	echo "ERROR: BAM not found"
	exit 2
fi

# Check if reference directory exists
if [ ! -d "$REFS_DIR" ]; then
	echo "ERROR: Reference directory not found"
	exit 2
fi

# Check if valid genome
AVAIL_GEN=(mm9 mm10 hg18 hg19 hg38)

if [[ ! "${AVAIL_GEN[@]}" =~ "$gen" ]]; then
	echo "ERROR: Genome assembly unavailable"
	exit 2
fi

# Check if chrom sizes file exists
if [ ! -f "${REFS_DIR}/${gen}.chrom.sizes" ]; then
	echo "ERROR: Chrom sizes not found"
	exit 2
fi

# Make reference if not available
ref="${REFS_DIR}/${gen}_${bin}_tiles.bed"
if [ ! -f $ref ]; then
	"${SCRIPT_DIR}/makeBins.R" $gen $bin
fi

# Check first read
pairs=$(cat <(samtools view -H $bam) <(samtools view $bam | head -n 1) | samtools view -c -f 1)
paired=$([ "$pairs" == 1 ] && echo "true" || echo "false")

if [[ $paired == true ]]; then
	igv_args="--pairs"
else
	seq=$(samtools view $bam | head -n 1 | cut -f 10)
	LEN=${#seq}

	igv_args="-e $(($modalLen - $LEN))"
fi

# IGVTools count
echo -e "Running command:\n\
    [igvtools count -w $bin --minMapQuality 1 $igv_args $bam ${prefix}.wig $gen]"

igvtools count -w $bin --minMapQuality 1 $igv_args $bam ${prefix}.wig $gen

# WigToBigWig
echo -e "Running command:\n\
    [wigToBigWig ${prefix}.wig "${REFS_DIR}/${gen}.chrom.sizes" ${prefix}.bw]"

wigToBigWig ${prefix}.wig "${REFS_DIR}/${gen}.chrom.sizes" ${prefix}.bw

# BigWigAverageOverBed
echo -e "Running command:\n\
    [paste $ref <(bigWigAverageOverBed ${prefix}.bw $ref stdout | cut -f 5) | awk 'BEGIN{FS="\t"; OFS="\t"}; {print \$1, \$2, \$3, \$7}' > ${prefix}_${suffix}.bed]"
paste $ref <(bigWigAverageOverBed ${prefix}.bw $ref stdout | cut -f 5) | awk 'BEGIN{FS="\t"; OFS="\t"}; {print $1, $2, $3, $7}' > ${prefix}_${suffix}.bed

# Clean up intermediate files
rm ${prefix}.wig ${prefix}.bw
