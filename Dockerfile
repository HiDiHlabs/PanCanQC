FROM r-base:3.2.5
MAINTAINER Ivo Buchhalter <i.buchhalter@dkfz-heidelberg.de>

# Install build-dependencies and needed tools
RUN \
    apt-get -y update && \
    apt-get -y install \
      wget git unzip \
      build-essential autoconf gfortran fort77 \
      zlibc zlib1g-dev libncurses5-dev xorg-dev libx11-dev libcurl4-openssl-dev nettle-dev \
      libtry-tiny-perl \
      mbuffer \
      samtools tabix bedtools

# Install R package getopt in version 1.20.0
RUN \
    cd /tmp && \
    wget https://cran.r-project.org/src/contrib/getopt_1.20.0.tar.gz && \
    R CMD INSTALL /tmp/getopt_1.20.0.tar.gz && \
    rm -f /tmp/getopt_1.20.0.tar.gz

# Install database files
RUN \
    mkdir -p /home/pcawg/data/localScratchDirectory /home/pcawg/data/bams /home/pcawg/data/results && \
    cd /home/pcawg/data/ && \
    wget https://powerfolder.dkfz.de/dl/fiJGFePRxoTVNk47XWbHYnXK/database_files.tar.gz && \
    tar -xf database_files.tar.gz && \
    rm database_files.tar.gz && \
	chmod -R 777 /home/pcawg/data

# Add external binaries and scripts
ADD binaries/ /home/pcawg/binaries
RUN \
    ln -s /home/pcawg/binaries/* /usr/local/bin/
ADD scripts /home/pcawg/data/scripts

VOLUME /home/pcawg/data/localScratchDirectory

CMD /home/pcawg/data/scripts/PCAWG_QC.sh /home/pcawg/data/bams/control.bam /home/pcawg/data/bams/tumor.bam
