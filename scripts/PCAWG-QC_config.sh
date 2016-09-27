#!/bin/bash

scriptPath=/home/pcawg/data/scripts/
export localScratchDirectory=/home/pcawg/data/localScratchDirectory

# BINARIES
PERL_BINARY=perl
PYTHON_BINARY=python
RSCRIPT_BINARY=Rscript
SAMTOOLS_BINARY=samtools
MBUFFER_BINARY=mbuffer

# CONFIG VALUES
cnv_min_coverage=0
mapping_quality=1000
min_windows=5

LOWESS_F=0.1
SCALE_FACTOR=0.9
COVERAGEPLOT_YLIMS=4
GC_bias_json_key=gc-bias
BASE_QUALITY_CUTOFF=0
WINDOW_SIZE=1

# TOOLS
TOOL_ANNOTATE_CNV_VCF=${scriptPath}/annotate_vcf.pl
TOOL_ADD_MAPPABILITY=${scriptPath}/addMappability.py
TOOL_MERGE_FILTER_CNV=${scriptPath}/merge_and_filter_cnv.py
TOOL_CORRECT_GC_BIAS_R=${scriptPath}/correctGCBias.R
TOOL_CONVERT_TO_JSON=${scriptPath}/convertTabToJson.py
TOOL_ESTIMATE_SEX=${scriptPath}/getSex.R
TOOL_COMBINED_BAM_ANALYSIS=${scriptPath}/flags_isizes_PEaberrations.pl
TOOL_COVERAGE_QC_D_IMPL=/home/pcawg/binaries/coverageQc
TOOL_GENOME_COVERAGE_D_IMPL=/home/pcawg/binaries/genomeCoverage
TOOL_FILTER_READ_BINS=${scriptPath}/filter_readbins.pl
TOOL_MEDIAN_MEAN=${scriptPath}/median_mean.R

# FILENAMES
FILENAME_DIFFCHROM_STATISTICS=${localScratchDirectory}/diffchrom_file
FILENAME_GENOME_COVERAGE=${localScratchDirectory}/genome_coverage
FILENAME_READBINS_COVERAGE=${localScratchDirectory}/readbin_coverage
FILENAME_COV_WINDOWS_1KB_ANNO=${localScratchDirectory}/readCoverage_1kb_windows.anno.txt.gz
FILENAME_COV_WINDOWS_WG=${localScratchDirectory}/readCoverage_10kb_windows.filtered.txt.gz
FILENAME_GC_CORRECT_PLOT=${localScratchDirectory}/gc_corrected.png
FILENAME_GC_CORRECTED_WINDOWS=${localScratchDirectory}/readCoverage_10kb_windows.filtered.corrected.txt
FILENAME_QC_GC_CORRECTION_JSON=${localScratchDirectory}/qc_gc_corrected.json
FILENAME_GC_CORRECTED_QUALITY=${localScratchDirectory}/qc_gc_corrected.tsv

# FILES_FROM_DATABASES
CHROM_SIZES_FILE=/home/pcawg/data/database_files/hs37d5.fa.chrLenOnlyACGT_realChromosomes.tab
MAPPABILITY_FILE=/home/pcawg/data/database_files/wgEncodeCrgMapabilityAlign100mer_chr.bedGraph.gz
CHROMOSOME_LENGTH_FILE=/home/pcawg/data/database_files/chrlengths.txt
REPLICATION_TIME_FILE=/home/pcawg/data/database_files/ReplicationTime_10cellines_mean_10KB.Rda
GC_CONTENT_FILE=/home/pcawg/data/database_files/hg19_GRch37_100genomes_gc_content_10kb.txt
MEDIAN_MEAN_TARGET_FILE=/home/pcawg/data/database_files/HUMAN_hsapiens.hs37d5_reducedGenome.n300l5M.sorted.bed
