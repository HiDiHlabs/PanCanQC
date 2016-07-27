#!/bin/bash

export pid=$1
export sample=$2
export path=$3
#export FILENAME_COV_WINDOWS_1KB=$4
export CONFIG_FILE=$4

export aceseqOutputDirectory=$path/$pid/ACEseqQC_${sample}/

#set -xuv
[[ ! -d $aceseqOutputDirectory} ]] && mkdir -p ${aceseqOutputDirectory}

#export FILENAME_COV_WINDOWS_1KB=/icgc/dkfzlsdf/project/mmml/xp_mmml/sequencing/whole_genome_sequencing/view-by-pid/$pid/${sample}/paired/merged-alignment/.merging_0/qualitycontrol/merged/coverage/${sample}_${pid}_readCoverage_1kb_windows.txt
export FILENAME_COV_WINDOWS_1KB=/icgc/dkfzlsdf/analysis/mmml/WPN/results_per_pid/$pid/coverage/${sample}_${pid}_readCoverage_1kb_windows.txt

sh /home/kleinhei/ACEseqQCOTP/gcCorrect/analysisTools/vcfAnno.sh
sh /home/kleinhei/ACEseqQCOTP/gcCorrect/analysisTools/cnvMergeFilter.sh
sh /home/kleinhei/ACEseqQCOTP/gcCorrect/analysisTools/correct_gc_bias.sh

