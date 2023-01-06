version 1.0

workflow CNVAnalysis {
  input {
    # Alignment Post Processing BAM
    File bam

    # CNV Ratios BED (absent for input controls)
    File? cnvRatiosBed

    String genomeName = 'hg19'

    Int binSize = 5000

    Boolean bypassCNVRescalingStep

    String dockerImage
  }

  String outPrefix = basename(bam, '.bam')

  call binning {
    input:
      bam = bam,
      outPrefix = outPrefix,
      genomeName = genomeName,
      binSize = binSize,
      dockerImage = dockerImage,
  }

  call rescaling {
    input:
      bam = bam,
      genomeName = genomeName,
      binnedBed = binning.outBed,
      outPrefix = outPrefix,
      dockerImage = dockerImage,
  }

  call cnv_rescaling {
    input:
      rescaledBed = rescaling.outBed,
      cnvRatiosBed = cnvRatiosBed,
      outPrefix = outPrefix,
      genomeName = genomeName,
      bypassCNVRescalingStep = bypassCNVRescalingStep,
      dockerImage = dockerImage,
  }

  call fitting {
    input:
      binnedBed = cnv_rescaling.outBinnedBed,
      outPrefix = outPrefix,
      dockerImage = dockerImage,
  }

  call pbs {
    input:
      binnedBed = cnv_rescaling.outBinnedBed,
      outPrefix = outPrefix,
      fitParams = fitting.outParams,
      dockerImage = dockerImage,
  }

  output {
    File binnedBed = cnv_rescaling.outBinnedBed
    Float signalToNoiseRatio = rescaling.signalToNoiseRatio
    File? cnvRatiosBedOut = cnv_rescaling.outCnvRatiosBed
    Boolean cnvsDetected = cnv_rescaling.cnvsDetected
    File fittingParams = fitting.outParams
    File pbsBed = pbs.outBed
  }
}

task binning {
  input {
    File bam
    String outPrefix
    String genomeName
    Int binSize

    String dockerImage
  }

  String outSuffix = 'binned'

  command <<<
    /scripts/binning/bin.sh \
      -f '~{bam}' \
      -g ~{genomeName} \
      -n ~{outPrefix} \
      -s ~{outSuffix} \
      -w ~{binSize}
  >>>

  runtime {
    docker: dockerImage
    disks: 'local-disk ' + ceil(1.1 * size(bam, 'G') + 1) + ' HDD'
    memory: '3G'
    cpu: 1
  }

  output {
    File outBed = '~{outPrefix}_~{outSuffix}.bed'
  }
}

task rescaling {
  input {
    File bam
    File binnedBed
    String outPrefix
    String genomeName

    String dockerImage
  }

  String outBedFile = outPrefix + '_map_scaled.bed'
  String outSnrFile = outPrefix + '_snr.txt'

  command <<<
    /scripts/rescaling/SubmitRescaleBinnedFiles.R \
      --bam_filename '~{bam}' \
      --binned_bed_filename '~{binnedBed}' \
      --genome ~{genomeName} \
      --output_filename '~{outBedFile}' \
      --snr_output_filename '~{outSnrFile}'
  >>>

  runtime {
    docker: dockerImage
    disks: 'local-disk ' + ceil(size(bam, 'G') + 12) + ' HDD'
    memory: '3G'
    cpu: 2
  }

  output {
    File outBed = outBedFile
    Float signalToNoiseRatio = read_float(outSnrFile)
  }
}

task cnv_rescaling {
  input {
    File rescaledBed
    File? cnvRatiosBed
    String outPrefix
    String genomeName
    Boolean bypassCNVRescalingStep

    String dockerImage
  }

  Boolean hasCnvRatios = defined(cnvRatiosBed)
  String outBinnedBedFile = outPrefix + '_binned_final.bed'
  String outCnvRatiosBedFile = outPrefix + '_cnv_ratios.bed'
  String outCnvFlagFile = outPrefix + '_cnv_flag.txt'

  command <<<
    /scripts/cnvRescaling/SubmitCNVRescale.R \
      --binned_bed_filename '~{rescaledBed}' \
      --cnv_rescale_output '~{outBinnedBedFile}' \
      --cnv_ratios_filename '~{if hasCnvRatios then cnvRatiosBed else outCnvRatiosBedFile}' \
      --cnv_flag_output_filename '~{outCnvFlagFile}' \
      --is_input_control ~{if hasCnvRatios then 'F' else 'T'} \
      --assembly ~{genomeName} \
      --bypass_cnv_rescaling_step ~{if bypassCNVRescalingStep then 'T' else 'F'}
  >>>

  runtime {
    docker: dockerImage
    disks: 'local-disk 1 HDD'
    memory: '2.25G'
    cpu: 1
  }

  output {
    File outBinnedBed = outBinnedBedFile
    File? outCnvRatiosBed = if hasCnvRatios then cnvRatiosBed else outCnvRatiosBedFile
    Boolean cnvsDetected = if bypassCNVRescalingStep then false else read_boolean(outCnvFlagFile)
    Boolean isInputControl = !hasCnvRatios
  }
}

task fitting {
  input {
    File? binnedBed
    String outPrefix

    String dockerImage
  }

  String outParamsFile = outPrefix + '_fit_params.txt'

  command <<<
    /scripts/fitting/SubmitFitDistributionWithCVM.R \
      --binned_bed_filename '~{binnedBed}' \
      --params_output '~{outParamsFile}'
  >>>

  runtime {
    docker: dockerImage
    disks: 'local-disk 1 HDD'
    memory: '1.25G'
    cpu: 1
  }

  output {
    File outParams = outParamsFile
  }
}

task pbs {
  input {
    File? binnedBed
    File fitParams
    String outPrefix

    String dockerImage
  }

  String outBedFile = outPrefix + '_pbs.bed'

  command <<<
    /scripts/pbs/SubmitProbabilityBeingSignal.R \
      --binned_bed_filename '~{binnedBed}' \
      --params_df_filename '~{fitParams}' \
      --pbs_filename '~{outBedFile}'
  >>>

  runtime {
    docker: dockerImage
    disks: 'local-disk 1 HDD'
    memory: '1G'
    cpu: 1
  }

  output {
    File outBed = outBedFile
  }
}