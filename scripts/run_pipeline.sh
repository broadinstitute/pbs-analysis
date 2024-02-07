#!/bin/bash

set -e

# Required inputs:
# -f, --is_input_control, --params_output, --pbs_output
# -f is the filename of the BAM, IS_INPUT_CONTROL indicates whether the BAM is a WCE/input control
# PARAMS_OUTPUT and PBS_OUTPUT are the filenames for the final outputs of the pipeline


# DEFAULTS:
BIN_SIZE=5000
GENOME="hg19"
BIN_OUTPUT_PREFIX="output"
BIN_OUTPUT_SUFFIX="binned"
RESCALED_OUTPUT_SUFFIX="map_scaled"
SNR_OUTPUT_SUFFIX="snr.txt"
CNV_OUTPUT_SUFFIX="binned_final"
CNV_FLAG_FILENAME_SUFFIX="cnv_flag.txt"
CNV_RATIOS_FILENAME=NULL

SCRIPTS_DIR=$(dirname "$0")
REFS_DIR="${SCRIPTS_DIR}/../references"


while [ $# -gt 0 ]
do
    case $1 in
        -f)
            shift
            BAM=$1
            shift
            ;;
        -w)
            shift
            BIN_SIZE=$1
            shift
            ;;
        -g)
            shift
            GENOME=$1
            shift
            ;;
        -n)
            shift
            BIN_OUTPUT_PREFIX=$1
            shift
            ;;
        -s)
            shift
            BIN_OUTPUT_SUFFIX=$1
            shift
            ;;
        --rescaled_output_suffix)
            shift
            RESCALED_OUTPUT_SUFFIX=$1
            shift
            ;;
        --snr_output_suffix)
            shift
            SNR_OUTPUT_SUFFIX=$1
            shift
            ;;
        --is_input_control)
            shift
            IS_INPUT_CONTROL=$1
            shift
            ;;
        --cnv_ratios_filename)
            shift
            CNV_RATIOS_FILENAME=$1
            shift
            ;;
        --gc_content_filename)
            shift
            GC_CONTENT_FILENAME=$1
            shift
            ;;
        --cnv_output_suffix)
            shift
            CNV_OUTPUT_SUFFIX=$1
            shift
            ;;
        --cnv_flag_filename_suffix)
            shift
            CNV_FLAG_FILENAME_SUFFIX=$1
            shift
            ;;
        --cnv_rescale_success_output)
            shift
            CNV_RESCALE_SUCCESS_OUTPUT=$1
            shift
            ;;
        --bypass_cnv_rescaling_step)
            shift
            BYPASS_CNV_RESCALING_STEP=$1
            shift
            ;;
        --params_output)
            shift
            PARAMS_OUTPUT=$1
            shift
            ;;
        --pbs_output)
            shift
            PBS_OUTPUT=$1
            shift
            ;;
        *)
            echo "ERROR: Unknown flag"
            exit 1
            ;;
    esac
done

GC_CONTENT_FILENAME="${REFS_DIR}/${GENOME}_${BIN_SIZE}_gc.bed"  


# binning
time "${SCRIPTS_DIR}/binning/bin.sh" -f $BAM -n $BIN_OUTPUT_PREFIX -s $BIN_OUTPUT_SUFFIX -w $BIN_SIZE -g $GENOME

echo "Binning COMPLETE"
BED_FILENAME=${BIN_OUTPUT_PREFIX}_${BIN_OUTPUT_SUFFIX}.bed
echo $BED_FILENAME

# rescaling and signal to noise ratio
RESCALED_OUTPUT=${BIN_OUTPUT_PREFIX}_${RESCALED_OUTPUT_SUFFIX}.bed
SNR_OUTPUT=${BIN_OUTPUT_PREFIX}_${SNR_OUTPUT_SUFFIX}
time "${SCRIPTS_DIR}/rescaling/SubmitRescaleBinnedFiles.R" --bam_filename $BAM --binned_bed_filename $BED_FILENAME --genome $GENOME \
#--output_filename $RESCALED_OUTPUT --snr_output_filename $SNR_OUTPUT
echo "Mappability rescaling COMPLETE"
echo $RESCALED_OUTPUT

# cnv_rescaling
CNV_RESCALED_OUTPUT=${BIN_OUTPUT_PREFIX}_${CNV_OUTPUT_SUFFIX}.bed
CNV_FLAG_OUTPUT_FILENAME=${BIN_OUTPUT_PREFIX}_${CNV_FLAG_FILENAME_SUFFIX}
time "${SCRIPTS_DIR}/cnvRescaling/SubmitCNVRescale.R" --binned_bed_filename $RESCALED_OUTPUT --is_input_control $IS_INPUT_CONTROL --cnv_ratios_filename $CNV_RATIOS_FILENAME --assembly $GENOME --cnv_rescale_output $CNV_RESCALED_OUTPUT --saved_gc_filename $GC_CONTENT_FILENAME --cnv_flag_output_filename $CNV_FLAG_OUTPUT_FILENAME --cnv_rescale_success_output $CNV_RESCALE_SUCCESS_OUTPUT --bypass_cnv_rescaling_step $BYPASS_CNV_RESCALING_STEP
echo "CNV rescaling COMPLETE"
echo $CNV_RESCALED_OUTPUT

# fitting
time "${SCRIPTS_DIR}/fitting/SubmitFitDistributionWithCVM.R" --binned_bed_filename $CNV_RESCALED_OUTPUT --params_output $PARAMS_OUTPUT
echo "Fitting background distribution COMPLETE"

# calculate PBS
time "${SCRIPTS_DIR}/pbs/SubmitProbabilityBeingSignal.R" --binned_bed_filename $CNV_RESCALED_OUTPUT --params_df_filename $PARAMS_OUTPUT --pbs_filename $PBS_OUTPUT
echo "Calculating PBS COMPLETE"
