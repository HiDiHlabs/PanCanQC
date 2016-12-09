#!/bin/bash

export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:${PATH}

source /etc/profile || exit 1

source /home/pcawg/data/scripts/PCAWG-QC_config.sh || exit 1

set -o pipefail
set -x

control_bam=$1
tumor_bam=$2

[[ ! -f $control_bam ]] || [[ ! -f $tumor_bam ]] && echo "Control of tumor bam file not found." && exit 2

#Tumor
NP_READBINS_IN_TUMOR=${localScratchDirectory}/np_readbins_in_merged_tumor
NP_COVERAGEQC_IN_TUMOR=${localScratchDirectory}/np_coverageqc_in_merged_tumor
NP_METRICS_IN_TUMOR=${localScratchDirectory}/np_metrics_in_merged_tumor
NP_COMBINEDANALYSIS_IN_TUMOR=${localScratchDirectory}/np_combinedanalysis_in_merged_tumor
NP_CALLABLE_BASES_TUMOR=${localScratchDirectory}/np_callable_bases_tumor
NP_BAMSTATS_TUMOR=${localScratchDirectory}/np_bamstats_tumor

#Control
NP_READBINS_IN_CONTROL=${localScratchDirectory}/np_readbins_in_merged_control
NP_COVERAGEQC_IN_CONTROL=${localScratchDirectory}/np_coverageqc_in_merged_control
NP_METRICS_IN_CONTROL=${localScratchDirectory}/np_metrics_in_merged_control
NP_COMBINEDANALYSIS_IN_CONTROL=${localScratchDirectory}/np_combinedanalysis_in_merged_control
NP_CALLABLE_BASES_CONTROL=${localScratchDirectory}/np_callable_bases_control
NP_BAMSTATS_CONTROL=${localScratchDirectory}/np_bamstats_control

bamname_tumor=`basename ${tumor_bam}`
mkfifo ${NP_READBINS_IN_TUMOR} ${NP_COVERAGEQC_IN_TUMOR} ${NP_COMBINEDANALYSIS_IN_TUMOR} ${NP_CALLABLE_BASES_TUMOR} ${NP_BAMSTATS_TUMOR}

bamname_control=`basename ${control_bam}`
mkfifo ${NP_READBINS_IN_CONTROL} ${NP_COVERAGEQC_IN_CONTROL} ${NP_COMBINEDANALYSIS_IN_CONTROL} ${NP_CALLABLE_BASES_CONTROL} ${NP_BAMSTATS_CONTROL}

# Tumor
# Create tree of input files/named pipes for the merged bam file
(cat ${tumor_bam} | tee ${NP_COVERAGEQC_IN_TUMOR} ${NP_READBINS_IN_TUMOR} ${NP_CALLABLE_BASES_TUMOR} ${NP_BAMSTATS_TUMOR} | ${SAMTOOLS_BINARY} view - > ${NP_COMBINEDANALYSIS_IN_TUMOR}) & procIDSAMpipe_TUMOR=$!
sleep 1

# SAM output is piped to perl script that calculates various QC measures
(${PERL_BINARY} ${TOOL_COMBINED_BAM_ANALYSIS} -i ${NP_COMBINEDANALYSIS_IN_TUMOR} -c ${CHROM_SIZES_FILE} -o ${FILENAME_DIFFCHROM_STATISTICS}_TUMOR.tmp ) & procIDCBA_TUMOR=$!

# genome coverage (depth of coverage and other QC measures in one file)
(${TOOL_COVERAGE_QC_D_IMPL} --alignmentFile=${NP_COVERAGEQC_IN_TUMOR} --outputFile=${FILENAME_GENOME_COVERAGE}_TUMOR.tmp --processors=1 --basequalCutoff=${BASE_QUALITY_CUTOFF-0} --ungappedSizes=${CHROM_SIZES_FILE}) & procIDGenomeCoverage_TUMOR=$!

# this part often fails with broken pipe, ?? where this comes from. The mbuffer did not help, maybe --processors=4 does?
(${TOOL_GENOME_COVERAGE_D_IMPL} --alignmentFile=${NP_READBINS_IN_TUMOR} --outputFile=/dev/stdout --processors=4 --mode=countReads --windowSize=${WINDOW_SIZE-1} | ${PERL_BINARY} ${TOOL_FILTER_READ_BINS} - ${CHROM_SIZES_FILE} > ${FILENAME_READBINS_COVERAGE}_TUMOR.tmp) & procIDReadbinsCoverage_TUMOR=$!



# Control
# Create tree of input files/named pipes for the merged bam file
(cat ${control_bam} | tee ${NP_COVERAGEQC_IN_CONTROL} ${NP_READBINS_IN_CONTROL} ${NP_CALLABLE_BASES_CONTROL} ${NP_BAMSTATS_CONTROL} | ${SAMTOOLS_BINARY} view - > ${NP_COMBINEDANALYSIS_IN_CONTROL}) & procIDSAMpipe_CONTROL=$!
sleep 1
# SAM output is piped to perl script that calculates various QC measures
(${PERL_BINARY} ${TOOL_COMBINED_BAM_ANALYSIS} -i ${NP_COMBINEDANALYSIS_IN_CONTROL} -c ${CHROM_SIZES_FILE} -o ${FILENAME_DIFFCHROM_STATISTICS}_CONTROL.tmp ) & procIDCBA_CONTROL=$!

# genome coverage (depth of coverage and other QC measures in one file)
(${TOOL_COVERAGE_QC_D_IMPL} --alignmentFile=${NP_COVERAGEQC_IN_CONTROL} --outputFile=${FILENAME_GENOME_COVERAGE}_CONTROL.tmp --processors=1 --basequalCutoff=${BASE_QUALITY_CUTOFF} --ungappedSizes=${CHROM_SIZES_FILE}) & procIDGenomeCoverage_CONTROL=$!

# this part often fails with broken pipe, ?? where this comes from. The mbuffer did not help, maybe --processors=4 does?
(${TOOL_GENOME_COVERAGE_D_IMPL} --alignmentFile=${NP_READBINS_IN_CONTROL} --outputFile=/dev/stdout --processors=4 --mode=countReads --windowSize=${WINDOW_SIZE} | ${PERL_BINARY} ${TOOL_FILTER_READ_BINS} - ${CHROM_SIZES_FILE} > ${FILENAME_READBINS_COVERAGE}_CONTROL.tmp) & procIDReadbinsCoverage_CONTROL=$!


# Determine callabel bases
(${SAMTOOLS_BINARY} depth ${NP_CALLABLE_BASES_CONTROL} ${NP_CALLABLE_BASES_TUMOR} | awk '($3>14 && $4>8)' | wc -l > ${localScratchDirectory}/callable_bases.txt.tmp) & procIDCALBASE=$!

sleep 2

# this is a new pipe for the read edits
(bam_stats -i ${NP_BAMSTATS_TUMOR} -o ${localScratchDirectory}/read_edits_tumor.tmp) & procIDReadEdits_TUMOR=$!

# this is a new pipe for the read edits
(bam_stats -i ${NP_BAMSTATS_CONTROL} -o ${localScratchDirectory}/read_edits_control.tmp) & procIDReadEdits_CONTROL=$!

# Tumor
# Some waits for parallel processes. This also depends on the used merge binary.
wait $procIDSAMpipe_TUMOR; [[ $? -gt 0 ]] && echo "Error from samtools SAM pipe" && exit 9
wait $procIDReadbinsCoverage_TUMOR; [[ $? -gt 0 ]] && echo "Error from genomeCoverage read bins" && exit 10
wait $procIDGenomeCoverage_TUMOR; [[ $? -gt 0 ]] && echo "Error from coverageQCD" && exit 11
wait $procIDCBA_TUMOR; [[ $? -gt 0 ]] && echo "Error from combined QC perl script" && exit 13
wait $procIDReadEdits_TUMOR; [[ $? -gt 0 ]] && echo "Error from tumor read edits" && exit 14

# Add here a final summary with all QC metrices rather than moving the files. The final results should be a table/json file with all metrices

mv ${FILENAME_DIFFCHROM_STATISTICS}_TUMOR.tmp ${FILENAME_DIFFCHROM_STATISTICS}_TUMOR
mv ${FILENAME_READBINS_COVERAGE}_TUMOR.tmp ${FILENAME_READBINS_COVERAGE}_TUMOR
mv ${FILENAME_GENOME_COVERAGE}_TUMOR.tmp ${FILENAME_GENOME_COVERAGE}_TUMOR
mv ${localScratchDirectory}/read_edits_tumor.tmp ${localScratchDirectory}/read_edits_tumor.txt

# Control
# Some waits for parallel processes. This also depends on the used merge binary.
wait $procIDSAMpipe_CONTROL; [[ $? -gt 0 ]] && echo "Error from samtools SAM pipe" && exit 9
wait $procIDReadbinsCoverage_CONTROL; [[ $? -gt 0 ]] && echo "Error from genomeCoverage read bins" && exit 10
wait $procIDGenomeCoverage_CONTROL; [[ $? -gt 0 ]] && echo "Error from coverageQCD" && exit 11
wait $procIDCBA_CONTROL; [[ $? -gt 0 ]] && echo "Error from combined QC perl script" && exit 13
wait $procIDReadEdits_CONTROL; [[ $? -gt 0 ]] && echo "Error from control read edits" && exit 14

wait $procIDCALBASE; [[ $? -gt 0 ]] && echo "Error in determination of callable bases" && exit 15

# Add here a final summary with all QC metrices rather than moving the files. The final results should be a table/json file with all metrices

mv ${FILENAME_DIFFCHROM_STATISTICS}_CONTROL.tmp ${FILENAME_DIFFCHROM_STATISTICS}_CONTROL
mv ${FILENAME_READBINS_COVERAGE}_CONTROL.tmp ${FILENAME_READBINS_COVERAGE}_CONTROL
mv ${FILENAME_GENOME_COVERAGE}_CONTROL.tmp ${FILENAME_GENOME_COVERAGE}_CONTROL
mv ${localScratchDirectory}/read_edits_control.tmp ${localScratchDirectory}/read_edits_control.txt
mv ${localScratchDirectory}/callable_bases.txt.tmp ${localScratchDirectory}/callable_bases.txt

# From here ACEseq QC from Kortine Kleinheinz
PID=pcawg_qc
# Tumor
SAMPLE=tumor
#reformat original coverage file so it can be annotated with annotate_vcf.pl
O_FILE_TUMOR="${FILENAME_COV_WINDOWS_1KB_ANNO}_TUMOR"
tmp_out_TUMOR="${O_FILE_TUMOR}_tmp"
cat ${FILENAME_READBINS_COVERAGE}_TUMOR | awk '{print $1,$2,$2+999,$3}' | sed 's/ /\t/g' | sed 's/^/chr/' |  sed '1i\#chr\tpos\tend\tcoverage' | \
${PERL_BINARY} "${TOOL_ANNOTATE_CNV_VCF}" \
               -a - \
               --aFileType=custom \
               --aChromColumn chr \
               --aPosColumn pos \
               --aEndColumn end \
               -b "${MAPPABILITY_FILE}" \
               --bFileType=bed \
               --reportBFeatCoord \
               --columnName map | \
${PYTHON_BINARY} ${TOOL_ADD_MAPPABILITY} \
                 -o "${tmp_out_TUMOR}"

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while adding the mappability values."
	exit 2
fi

mv $tmp_out_TUMOR ${O_FILE_TUMOR}

${PYTHON_BINARY} "${TOOL_MERGE_FILTER_CNV}" \
                 --inputfile    "${O_FILE_TUMOR}" \
                 --output       "${FILENAME_COV_WINDOWS_WG}_TUMOR" \
	               --coverage     ${cnv_min_coverage} \
                 --mappability  ${mapping_quality} \
                 --NoOfWindows  ${min_windows} 

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code in the cnv merge and filter process;" 
	exit 2
fi

tmp_corrected_windowfile_TUMOR=${FILENAME_GC_CORRECTED_WINDOWS}_TUMOR.tmp
tmp_corrected_table_TUMOR=${FILENAME_GC_CORRECTED_QUALITY}_TUMOR.tmp
corrected_table_slim_TUMOR=${FILENAME_GC_CORRECTED_QUALITY}_TUMOR.slim.txt

${RSCRIPT_BINARY} ${TOOL_CORRECT_GC_BIAS_R} \
	                --windowFile	${FILENAME_COV_WINDOWS_WG}_TUMOR \
	                --timefile	${REPLICATION_TIME_FILE} \
	                --chrLengthFile	${CHROMOSOME_LENGTH_FILE} \
	                --pid		${PID} \
	                --sample	${SAMPLE} \
	                --outfile	${tmp_corrected_windowfile_TUMOR} \
	                --corPlot	${FILENAME_GC_CORRECT_PLOT}_TUMOR \
	                --corTab	${tmp_corrected_table_TUMOR} \
	                --qcTab		${corrected_table_slim_TUMOR} \
	                --gcFile	${GC_CONTENT_FILE} \
	                --outDir	${localScratchDirectory} \
	                --lowess_f	${LOWESS_F} \
	                --scaleFactor	${SCALE_FACTOR} \
	                --coverageYlims ${COVERAGEPLOT_YLIMS}


if [[ $? != 0 ]]
then
	echo "Something went wrong during GC correction. Program had non-zero exit status, exiting pipeline...\n\n"
	exit 2
fi	

mv ${tmp_corrected_windowfile_TUMOR} ${FILENAME_GC_CORRECTED_WINDOWS}_TUMOR

# Control
SAMPLE=control
#reformat original coverage file so it can be annotated with annotate_vcf.pl
O_FILE_CONTROL="${FILENAME_COV_WINDOWS_1KB_ANNO}_CONTROL"
tmp_out_CONTROL="${O_FILE_CONTROL}_tmp"
cat ${FILENAME_READBINS_COVERAGE}_CONTROL | awk '{print $1,$2,$2+999,$3}' | sed 's/ /\t/g' | sed 's/^/chr/' |  sed '1i\#chr\tpos\tend\tcoverage' | \
${PERL_BINARY} "${TOOL_ANNOTATE_CNV_VCF}" \
               -a - \
               --aFileType=custom \
               --aChromColumn chr \
               --aPosColumn pos \
               --aEndColumn end \
               -b "${MAPPABILITY_FILE}" \
               --bFileType=bed \
               --reportBFeatCoord \
               --columnName map | \
${PYTHON_BINARY} ${TOOL_ADD_MAPPABILITY} \
                 -o "${tmp_out_CONTROL}"

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while adding the mappability values."
	exit 2
fi

mv $tmp_out_CONTROL ${O_FILE_CONTROL}

${PYTHON_BINARY} "${TOOL_MERGE_FILTER_CNV}" \
                 --inputfile    "${O_FILE_CONTROL}" \
                 --output       "${FILENAME_COV_WINDOWS_WG}_CONTROL" \
	               --coverage     ${cnv_min_coverage} \
                 --mappability  ${mapping_quality} \
                 --NoOfWindows  ${min_windows} 

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code in the cnv merge and filter process;" 
	exit 2
fi

tmp_corrected_windowfile_CONTROL=${FILENAME_GC_CORRECTED_WINDOWS}_CONTROL.tmp
tmp_corrected_table_CONTROL=${FILENAME_GC_CORRECTED_QUALITY}_CONTROL.tmp
corrected_table_slim_CONTROL=${FILENAME_GC_CORRECTED_QUALITY}_CONTROL.slim.txt

${RSCRIPT_BINARY} ${TOOL_CORRECT_GC_BIAS_R} \
	                --windowFile	${FILENAME_COV_WINDOWS_WG}_CONTROL \
	                --timefile	${REPLICATION_TIME_FILE} \
	                --chrLengthFile	${CHROMOSOME_LENGTH_FILE} \
	                --pid		${PID} \
	                --sample	${SAMPLE} \
	                --outfile	${tmp_corrected_windowfile_CONTROL} \
	                --corPlot	${FILENAME_GC_CORRECT_PLOT}_CONTROL \
	                --corTab	${tmp_corrected_table_CONTROL} \
	                --qcTab		${corrected_table_slim_CONTROL} \
	                --gcFile	${GC_CONTENT_FILE} \
	                --outDir	${localScratchDirectory} \
	                --lowess_f	${LOWESS_F} \
	                --scaleFactor	${SCALE_FACTOR} \
	                --coverageYlims ${COVERAGEPLOT_YLIMS}


if [[ $? != 0 ]]
then
	echo "Something went wrong during GC correction. Program had non-zero exit status, exiting pipeline...\n\n"
	exit 2
fi	

mv ${tmp_corrected_windowfile_CONTROL} ${FILENAME_GC_CORRECTED_WINDOWS}_CONTROL

# Gather QC values and print results:
CALLABLE_BASES=`cat ${localScratchDirectory}/callable_bases.txt`

# Control
DIFFCHROM_VALUE_CONTROL=`cat ${FILENAME_DIFFCHROM_STATISTICS}_CONTROL`
COVERAGE_CONTROL=`grep "^all" ${FILENAME_GENOME_COVERAGE}_CONTROL | cut -f 16 | sed 's|x$||'`
FWHM_CONTROL=`tail -n +2 ${corrected_table_slim_CONTROL} | cut -f 4`
READ_EDITS_CONTROL=`perl -e 'use strict; use warnings; open(IN, "<$ARGV[0]"); my $r1 = 0, my $r2 = 0; while(<IN>){chomp; next if($_ =~ /^bam_filename/); my @l = split("\t", $_); $r1 += $l[12]; $r2 += $l[13];}my $res = 0; if($r1 > $r2){$res = $r1/$r2;}else{$res = $r2/$r1;}print $res, "\n";' ${localScratchDirectory}/read_edits_control.txt`
MEDIAN_MEAN_CONTROL=`cat ${FILENAME_READBINS_COVERAGE}_CONTROL | perl -F"\t" -ae 'print join("\t", $F[0],$F[1],$F[1]+999,$F[2])' | intersectBed -a stdin -b ${MEDIAN_MEAN_TARGET_FILE} | Rscript ${TOOL_MEDIAN_MEAN}`

# Tumor
DIFFCHROM_VALUE_TUMOR=`cat ${FILENAME_DIFFCHROM_STATISTICS}_TUMOR`
COVERAGE_TUMOR=`grep "^all" ${FILENAME_GENOME_COVERAGE}_TUMOR | cut -f 16 | sed 's|x$||'`
FWHM_TUMOR=`tail -n +2 ${corrected_table_slim_TUMOR} | cut -f 4`
READ_EDITS_TUMOR=`perl -e 'use strict; use warnings; open(IN, "<$ARGV[0]"); my $r1 = 0, my $r2 = 0; while(<IN>){chomp; next if($_ =~ /^bam_filename/); my @l = split("\t", $_); $r1 += $l[12]; $r2 += $l[13];}my $res = 0; if($r1 > $r2){$res = $r1/$r2;}else{$res = $r2/$r1;}print $res, "\n";' ${localScratchDirectory}/read_edits_tumor.txt`
MEDIAN_MEAN_TUMOR=`cat ${FILENAME_READBINS_COVERAGE}_TUMOR | perl -F"\t" -ae 'print join("\t", $F[0],$F[1],$F[1]+999,$F[2])' | intersectBed -a stdin -b ${MEDIAN_MEAN_TARGET_FILE} | Rscript ${TOOL_MEDIAN_MEAN}`

# Combine
mkdir -p $HOME/results
echo -e "Callable bases\t${CALLABLE_BASES}\nReads mapping on different chromosomes control\t${DIFFCHROM_VALUE_CONTROL}\nReads mapping on different chromosomes tumor\t${DIFFCHROM_VALUE_TUMOR}\nCoverage control\t${COVERAGE_CONTROL}\nCoverage tumor\t${COVERAGE_TUMOR}\nFWHM control\t${FWHM_CONTROL}\nFWHM tumor\t${FWHM_TUMOR}\nMedian over Mean control\t${MEDIAN_MEAN_CONTROL}\nMedian over Mean tumor\t${MEDIAN_MEAN_TUMOR}\nRead edits control\t${READ_EDITS_CONTROL}\nRead edits tumor\t${READ_EDITS_TUMOR}" > $HOME/results/QC_results.tsv
