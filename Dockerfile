FROM ubuntu:xenial
MAINTAINER Ivo Buchhalter @ DKFZ

RUN \
    umask 000 && \
    apt-get -y update && \
    apt-get -y upgrade && \
    apt-get -y install apt-utils && \
    apt-get -y install apt-utils apt-transport-https && \
    apt-get -y install autoconf make && \
    apt-get -y install build-essential && \
    apt-get -y install zlibc zlib1g zlib1g-dev && \
    apt-get -y install libncurses5-dev && \
    apt-get -y install sudo && \
    apt-get -y install mbuffer && \
    apt-get -y install gfortran fort77 xorg-dev libx11-dev && \
    apt-get -y install wget && \
    apt-get -y install libcurl4-openssl-dev && \
    apt-get -y install nettle-dev && \
    apt-get -y install git && \
    apt-get -y install unzip && \
    apt-get -y install libtry-tiny-perl

RUN \
    echo "deb https://cran.r-project.org/bin/linux/ubuntu xenial/" >> /etc/apt/sources.list && \
    apt-get -y update && \
    apt-get -y build-dep r-base && \
    mkdir /program_files && \
    cd /program_files && \
    wget https://cran.r-project.org/src/base/R-3/R-3.2.5.tar.gz && \
    tar -xzf R-3.2.5.tar.gz && \
    rm R-3.2.5.tar.gz && \
    cd R-3.2.5/ && \
    ./configure && \
    make && \
    make check && \
    make install && \
    cd /program_files/ && \
    wget https://cran.r-project.org/src/contrib/getopt_1.20.0.tar.gz && \
    R CMD INSTALL /program_files/getopt_1.20.0.tar.gz && \
    chmod -R 777 /program_files && \
    mkdir -p /home/pcawg

RUN \
    mkdir -p /home/pcawg/data/localScratchDirectory && \
    mkdir /home/pcawg/data/bams && \
    cd /home/pcawg/data/ && \
    wget https://powerfolder.dkfz.de/dl/fiJGFePRxoTVNk47XWbHYnXK/database_files.tar.gz && \
    tar -xf database_files.tar.gz && \
    rm database_files.tar.gz && \
    mkdir /home/pcawg/data/results && \
    cp /etc/profile /home/pcawg/.bashrc && \
    chmod -R 777 /home/pcawg/data

RUN \
    mkdir /home/pcawg/temp/ && \
    cd /home/pcawg/temp/ && \
    git clone https://github.com/samtools/htslib.git && \
    cd htslib/ && \
    autoheader && \
    autoconf && \
    ./configure && \
    make && \
    make install && \
    cd ../ && \
    git clone https://github.com/samtools/samtools.git && \
    cd samtools && \
    make && \
    make install && \
    cd ../ && \
    ln -s samtools/samtools /bin/ && \
    ln -s htslib/tabix /bin/ && \
    ln -s htslib/bgzip /bin/

ADD binaries/ /home/pcawg/binaries

RUN \
    ln -s /home/pcawg/binaries/* /bin/
    
ADD scripts /home/pcawg/data/scripts   

ENTRYPOINT bash /home/pcawg/data/scripts/PCAWG_QC.sh

# docker run -v /ibios/co01/buchhalt/gits/mixed_projects/PanCanQC/database_files/:/home/pcawg/data/database_files/ -v /ibios/co01/buchhalt/temp/tumor_MD-207_merged.mdup.bam:/home/pcawg/data/bams/tumor.bam -v /ibios/co01/buchhalt/temp/tumor_MD-207_merged.mdup.bam:/home/pcawg/data/bams/control.bam -v /ibios/co01/buchhalt/gits/mixed_projects/PanCanQC/results/:/home/pcawg/data/results/ pcawg-qc
