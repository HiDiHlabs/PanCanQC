#!/usr/bin/Rscript

#Kortine Kleinheinz 14.08.14
#This script can be used to correct for gc bias and replication timing starting from tumor and control read counts in 10kb windows
#

require("getopt")
# source("/home/pcawg/binaries/R-3.2.5/getopt/R/getopt.R")
# source("/home/pcawg/binaries/R-3.2.5/getopt/R/utils.R")

script_dir = "/home/pcawg/data/scripts/"
source(paste0(script_dir,"/qq.R"))
source(paste0(script_dir, "/getopt.R"))

#default Parameter
functionPath <- paste0(script_dir, "/correctGCBias_functions.R")
scale_factor <- 0.9
lowess_f <- 0.1
magnification_factor <- 1e9
coverageYlims <- 4

getopt2(matrix(c('timefile',				 't', 1, "character", 'file for replication timing',
		         'windowFile',			 'f', 1, "character", 'file with 10kb window coordinates and coverage',
		         'chrLengthFile',		 'h', 1, "character", 'file with chromosomelength',
		         'sample',      		 'a', 1, "character", 'normal|tumor',
	                 'pid',				 'p', 1, "character", 'patient ID',
                	 'outfile',			 'o', 1, "character", 'new file with corrected and raw coverage',
	                 'corPlot',			 'c', 1, "character", 'Name for plot with corrected GC bias',
        	         'corTab',			 'r', 1, "character", 'Name for table with GC bias parameter',
        	         'qcTab',			 'q', 1, "character", 'Name for table with GC bias parameters relevant for qc',
	                 'gcFile',			 'g', 1, "character", 'table with gc content per 10kb window',
        	         'outDir',			 'x', 1, "character", 'directory for outputfiles',
	        	 'scaleFactor',			 's', 2, "double"   , 'scaling factor to determine range of points to take as main cloud',
		         'lowess_f',			 'l', 2, "double"   , 'smoothing parameter of lowess function',
		         'magnification_factor',	 'm', 2, "double"   , 'factor to prevent underflow in slope and curvature',
		         'functionPath', 		 'u', 2, "character", 'path and name of script with functions for script',
		         'coverageYlims', 		 'y', 2, "numeric", 'Y range for corrected coverage plot'
                ), ncol = 5, byrow = TRUE))

cat(qq("windowFile: @{windowFile}\n\n"))
cat(qq("timefile: @{timefile}\n\n"))
cat(qq("chrLengthFile: @{chrLengthFile}"))
cat(qq("pid: @{pid}\n\n"))
cat(qq("outfile: @{outfile}\n\n"))
cat(qq("corPlot: @{corPlot}\n\n"))
cat(qq("corTab: @{corTab}\n\n"))
cat(qq("qcTab: @{qcTab}\n\n"))
cat(qq("gcFile: @{gcFile}\n\n"))
cat(qq("outDir: @{outDir}\n\n"))
cat(qq("coverageYlims: @{coverageYlims}\n\n"))
cat("\n")

coverageLims <- c( -abs(coverageYlims), abs(coverageYlims) )
source(functionPath)

plotDir <- dirname(corPlot)
outputfile_gc  <- paste0(corPlot) #includes plotDir already
outputfile_rep <- paste0(plotDir,"/", pid, "_qc_rep_corrected.png")

# define files for output of quantification of correction factors
outrepQuant_file <- paste0(plotDir,"/",pid,"_qc_repQuant.tsv")
outGCcorrectQuant_file <- corTab
plot_flag <- 1
restrict_flag <- 1
writeGCcorrect_flag <- 1


cat("reading input files\n")
file_gc  <- read.table(file=gcFile, header=TRUE ,sep="\t", check.names=FALSE)
file_gc$chromosome <- gsub('X', 23, file_gc$chromosome )
file_gc$chromosome <- gsub('Y', 24, file_gc$chromosome )

file_cov <- read.table(file=windowFile, head=FALSE ,sep="\t", check.names=FALSE)
colnames(file_cov) <- c('chromosome', 'start', 'coverage') 

#read replication time data
load(file=timefile)
colnames(time10) <- c('chromosome', 'tenkb', 'time')
time10$chromosome <- gsub('X',23, time10$chromosome)

#read and sort chromosome length file
chrLengthTab = read.table(chrLengthFile, header = FALSE, as.is = TRUE)
chrLengthTab = data.frame(chrLengthTab)
colnames(chrLengthTab)  <- c("chromosome", "length", "info")[1:dim(chrLengthTab)[2]]
chrLengthTab$chromosome <- gsub('chr','',chrLengthTab$chromosome)
chrLengthTab$chromosome <- gsub('X', 23, chrLengthTab$chromosome)
chrLengthTab$chromosome <- gsub('Y', 24, chrLengthTab$chromosome)
chrLengthTab$chromosome <- as.numeric(chrLengthTab$chromosome)
chrLengthTab		<- chrLengthTab[order(chrLengthTab$chromosome),]

#get normalized coverage and adjust coordinates
file_cov$cov <- file_cov$coverage / sum(as.numeric(file_cov$coverage))

file_cov$start <- floor(file_cov$start/10000)*10000
file_cov$chromosome <- gsub('chr', '', file_cov$chromosome)

#merge dataframes
file_comb <- merge( x = file_cov,
		y = file_gc,
		by=c('chromosome','start'),
		sort =F)
rm(file_gc)
rm(file_cov)
gc()

cat("correction for gc bias...\n")
#reorder original sample according to gc_content
order_file <- file_comb[order(file_comb$gc_content),]

#compute correction for gc bias by fits
cov_fit  <- lowess(order_file$gc_content,order_file$cov, lowess_f)
cov_mean <- mean(order_file$cov)
cov_mean_vector <- matrix(cov_mean,length(order_file$gc_content))
#compute correction for gc bias by fit from all points
order_file$corrected_cov <- order_file$cov/cov_fit$y

###isolate main cluster in normal sample

# modify call to the function defineMainCluster to return also the width of the main cluster
temp_main_cluster  <- defineMainCluster(order_file$cov, cov_fit, sample,  plotFlag=plotDir)
main_cluster_ind   <- temp_main_cluster$ind
main_cluster_width <- temp_main_cluster$width

main_cluster <- data.frame(cov = order_file$cov[main_cluster_ind], gc_content = order_file$gc_content[main_cluster_ind])

order_file$main_cluster <- 0
order_file$main_cluster[main_cluster_ind] <- 1

#compute a new correction for only main cluster of the control
small_fit <- lowess(main_cluster$gc_content,main_cluster$cov, lowess_f)
fit4 <- approx(small_fit$x, small_fit$y, order_file$gc_content, rule=2)
# make an equidistant vector for later extraction of quantitative values (slope, curvature etc.)
x_vector <- seq( 0, max(fit4$x), max(fit4$x)/length(fit4$x))
fit5     <- approx( fit4$x, fit4$y, x_vector, rule=2)
main_cluster_mean <- mean(main_cluster$cov)
main_cluster_mean_vector <- matrix( main_cluster_mean, length(order_file$gc_content) )

###isolate main cluster in tumor sample

#compute correction for gc bias by fits
order_file$corrected_cov4 <- (order_file$cov/fit4$y)*(main_cluster_mean/cov_mean)

#extract minimal y-values of lowess fit as quality measure
min_gc <- 0.33
max_gc <- 0.6

restrict_x <- which( fit4$x >= min_gc & fit4$x <= max_gc  )
minimal_coverage_gcfit  <- min( fit4$y[restrict_x] )

#now compute FWHM of the cluster by calling method
main_cluster_density            <- density(order_file$corrected_cov4[main_cluster_ind]) # , adjust=0.1,from=-0.5,to=2.5)
main_cluster_FWHM_data          <- extractFWHM(main_cluster_density)
main_cluster_FWHM               <- main_cluster_FWHM_data$FWHM
main_cluster_half_max_pos_right <- main_cluster_FWHM_data$half_max_pos_right
main_cluster_half_max_pos_left  <- main_cluster_FWHM_data$half_max_pos_left

if ( writeGCcorrect_flag == 1 ){
	sub_order_file = order_file[ ,c("chromosome", "start", "corrected_cov4") ]
	write.table( sub_order_file, paste0(plotDir, "/", pid, "_", sample, "_all_seg.gc_corrected.txt") ,sep="\t", col.names=T, row.names=F, quote=F )
}


cat("correction for replication timing...\n")
order_file$tenkb <- order_file$start %/% 10000

rdWithTime <- merge(order_file, time10, by.x = c("chromosome", "tenkb"), by.y = c("chromosome", "tenkb"))

#compute correction for replication timing  by fits
order_file_rt <- rdWithTime[order(rdWithTime$time),]

#compute a new correction for only main cluster
cat('fitting end\n')
small_fit_rt <- lowess( order_file_rt$time[order_file_rt$main_cluster==1], order_file_rt$corrected_cov4[order_file_rt$main_cluster==1], lowess_f )
cov_fit4_rt <- approx(small_fit_rt$x,small_fit_rt$y,order_file_rt$time, rule=2 )
main_cluster_mean_rt <- mean(main_cluster$cov)
main_cluster_mean_rt_vector <- matrix( main_cluster_mean_rt, length(order_file_rt$time) )

#compute correction for rc bias by fits
order_file_rt$corrected_cov4_rt <- (order_file_rt$corrected_cov4/cov_fit4_rt$y)*(main_cluster_mean_rt/cov_mean)

#do quantitative characterization of GC bias correction
cat('do quantitative evaluation of GC bias correction:\n')
GC_correct_lin <- lm(y~x,data=small_fit)
cov_table <- summary(GC_correct_lin)$coefficients

#compute slope (1st derivative) to estimate bias alternatively to linear regression (note that the x-values are equally spaced and allways the same, therefore one can compute the slope as Delta_y instead of Delta_y/Delta_x)
slope <- diff(fit5$y)*magnification_factor
mean_abs_slope <- mean(abs(slope))
mean_slope <- mean(slope)
#compute curvature (2nd derivative) to estimate convexity of GC bias
curvature <- diff(slope)
mean_abs_curvature <- mean(abs(curvature))
mean_curvature <- mean(curvature)
#prepare for output
GCcorrectQuant_string <- paste(t(cov_table), sep="\t")


#write table again, so it can be converted to a json file
write.table( data.frame(pid, main_cluster_FWHM_data, mean_slope, mean_abs_slope,
			mean_curvature, mean_abs_curvature, 
			main_cluster_width, main_cluster_FWHM,
			minimal_coverage_gcfit ), 
			sep="\t", file=qcTab, col.names=TRUE, row.names=FALSE, quote=F )

#do linear interpolation of replication time correction
cat('do linear regressions on rep time correction:\n')
#restrict fit for rep timing to 15 < reptime < 70 for fitting the slope 
restriction_ind <- which(small_fit_rt$x > 15 & small_fit_rt$x < 70)
model_fit_rt <- data.frame( x=small_fit_rt$x[ restriction_ind], y=small_fit_rt$y[ restriction_ind ] )
if(restrict_flag){
  rep_time_correct_lin <- lm( y~x, data=model_fit_rt)
} else {
  rep_time_correct_lin <- lm(y~x,data=small_fit_rt)
}

cov_table <- summary(rep_time_correct_lin)$coefficients
#prepare for output
repQuant_string <- paste(t(cov_table),sep="\t")
#write to pid specific file
write(c(pid,repQuant_string),sep="\t",file=outrepQuant_file,ncolumns=17)

cat(qq("writing results into @{outfile}\n\n"))
out_table		<- order_file_rt[,c('chromosome', 'start', 'cov', 'corrected_cov4_rt')]
colnames(out_table)	<- c('chromosome', 'start', 'cov') 
out_table		<- out_table[order(out_table$chromosome, out_table$start),]

sub_order_file <- order_file_rt[,c('chromosome', 'start', 'corrected_cov4_rt' )]
colnames(sub_order_file)  <- c('chromosome', 'start', 'covnorm') 
sub_order_file  	  <- sub_order_file[order(sub_order_file$chromosome, sub_order_file$start),]

#include Y chromosome without RT correction
sel_y_windows <- which(order_file$chromosome==24)
if(length(sel_y_windows)>0){
	y_windows 		<- order_file[sel_y_windows,]
	y_windows 		<- y_windows[,c("chromosome", "start", "cov", "corrected_cov4")]
	colnames(y_windows)	<- c('chromosome', 'start',  'cov') 
	y_windows 		<- y_windows[order(y_windows$chromosome, y_windows$start), ]
	out_table		<- rbind(out_table, y_windows)
  
	y_windows2  <- order_file[sel_y_windows,c("chromosome", "start", "corrected_cov4" )]
	colnames(y_windows2)  <- c('chromosome', 'start', 'covnorm' ) 
	#y_windows2  <- order_file[,]
	sub_order_file <- rbind(sub_order_file,y_windows2)
}

#TODO: insert pdf here
#      filenames?
pdf( paste0(plotDir, "/", pid, "_qc_coverageDensityByChroomosome.pdf") )
  corCovInd <- which( colnames(sub_order_file) =="covnorm" )
  diffPeaks <- checkControl( sub_order_file, corCovInd )
dev.off()

if ( ! length(diffPeaks) == sum( is.na(diffPeaks) ) ) {
  sel <- paste( which( ! is.na( diffPeaks ) ), collapse=", " )
  bodyText <- paste0("Warning ", pid, ": Errors found for chromosome ", sel )
}

#create coverage and gc/replication-timing plots
if (plot_flag) {
  
# if(sample == "control" | sample=="normal" | sample== "Control" | sample == "Normal"){
  for( chr in unique(sub_order_file$chromosome) ) {
    selSub <- which(sub_order_file$chromosome==chr)
    sub <- sub_order_file[selSub,]
    
    chr <- gsub( 24, "Y", chr)
    chr <- gsub( 23, "X", chr)
    
    png(paste0(plotDir,"/",sample, "_", pid, "_chr", chr,"_coverage.png"), width=2000, height=1000, type="cairo")
      plotCoverageSingle( sub, chr=chr, ylims=coverageLims )
    dev.off()
  }
  
  #plot whole genome coverage
  coordinates <- adjustCoordinates( chrLengthTab, sub_order_file )
  newCoverageTab    <- coordinates$coverageTab
  chromosomeBorders <- coordinates$chromosomeBorders
  png(paste0(plotDir,"/", sample, "_", pid, "_wholeGenome_coverage.png"), width=2000, height=1000, type='cairo')
      plotCoverageSingle( newCoverageTab, chromosomeBorders=chromosomeBorders, ylims=coverageLims )
  dev.off()
#  }

write.table(out_table, outfile, row.names=FALSE, col.names=TRUE,sep='\t', quote=F)
write.table(sub_order_file, qq("@{plotDir}/all_corrected.txt"), row.names=FALSE, col.names=TRUE,sep='\t', quote=F)

cat(qq("creating plots...\n\n"))
png(file=outputfile_gc, width=1000, height=2000, type='cairo')
#	par(mfrow=c(4,3))
	par(mfrow=c(2,3))

	##plot raw data
	#control raw
	ylims=c(median(order_file$cov)-10e-6, median(order_file$cov)+10e-6)
	matplot(order_file$gc_content,order_file$cov,type="p",pch=20,col="#0000FF22",xlim=c(0.28,0.71), xlab="%GC",ylab="normalized reads",main=paste(pid, sample, sep=" "),cex=c(0.2,0.2), ylim=ylims)
	points(main_cluster$gc_content,main_cluster$cov, type='p', pch=20,col="#FF000022",ylim=c(0.1,2.2),xlim=c(0.28,0.71),xlab="%GC",ylab="normalized reads",main=paste(pid, "tumor", sep=" "),cex=c(0.2,0.2))
	lines(cov_fit$x,cov_fit$y)
	lines(fit4$x,fit4$y, col='blue')
	lines(cov_fit$x,0.5*cov_fit$y)
	lines(cov_fit$x,cov_mean_vector,col="red")

	#plot corrected graphs (curves fitted to main cloud)
	matplot(order_file$gc_content,order_file$corrected_cov4,type="p",pch=20,col="#0000FF22",ylim=c(0.1,2.2),xlim=c(0.28,0.71),xlab="%GC",ylab="normalized reads",main=paste(pid, sample, sep=" "),cex=c(0.2,0.2))
  points(order_file$gc_content[main_cluster_ind],order_file$corrected_cov4[main_cluster_ind], type='p', pch=20,col="#FF000022",ylim=c(0.1,2.2),xlim=c(0.28,0.71),xlab="%GC",ylab="normalized reads",main=paste(pid, sample, sep=" "),cex=c(0.2,0.2) )

	#plot corrected plots with replication timing correction
	matplot(order_file_rt$gc_content, order_file_rt$corrected_cov4_rt, type="p",pch=20,col="#0000FF22",ylim=c(0.1,2.2),xlim=c(0.28,0.71),xlab="%GC",ylab="normalized reads",main=paste(pid, sample, sep=" "),cex=c(0.2,0.2))
	abline(h=c(0.5,1,1.5))

	#plot uncorrected curves replication timing
	matplot(order_file_rt[,"time"],order_file_rt$cov,type="p",pch=20,col="#0000FF22", xlim=c(0,80),xlab="time",ylab="normalized reads",main=paste(pid, sample, sep=" "),cex=c(0.2,0.2), ylim=ylims)

	#replication timing and gc corrected values
	#control
	matplot(order_file_rt[,"time"],order_file_rt$corrected_cov4,type="p",pch=20,col="#0000FF22",xlim=c(0,80), ylim=c(0.1,2.2), xlab="time",ylab="normalized reads",main=paste(pid, sample, sep=" "),cex=c(0.2,0.2))
	points(order_file_rt[order_file_rt$main_cluster==1,"time"], order_file_rt$corrected_cov4[order_file_rt$main_cluster==1],pch=20,col="#FF000022",cex=c(0.2,0.2))
	lines(cov_fit4_rt$x,cov_fit4_rt$y, col='blue')
	lines(cov_fit4_rt$x,0.5*cov_fit4_rt$y)
	lines(cov_fit4_rt$x,1.5*cov_fit4_rt$y)
	lines(cov_fit4_rt$x,main_cluster_mean_rt_vector, col="red")      

	#plot corrected plots      
	matplot(order_file_rt[,"time"],order_file_rt$corrected_cov4_rt,type="p",pch=20,col="#0000FF22", ylim=c(0.1,2.2), xlim=c(0,80), xlab="time",ylab="normalized reads",main=paste(pid, sample, sep=" "),cex=c(0.2,0.2))
	abline(h=c(0.5,1,1.5))

dev.off()

}

