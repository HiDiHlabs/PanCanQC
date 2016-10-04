#!/usr/bin/Rscript

input <- file('stdin', 'r')
cov <- read.table( input, nrows=1.5e6, sep="\t")
colnames(cov) <- c("chr", "start", "end", "cov")
ratio <- median(cov$cov)/mean(cov$cov)

cat( ratio, "\n" )

