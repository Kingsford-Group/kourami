# Base image
FROM ubuntu:16.04
WORKDIR /opt

# Set description and maintainer.
LABEL description="Dockerfile for installing Kourami" \
      maintainer="Vladimir Kovacevic, Seven Bridges Genomics, <vladimir.kovacevic@sbgenomics.com>"


RUN apt-get update && apt-get install -y \
    autoconf \
    automake \
    make \
    gcc \
    zlib1g-dev \
    libgsl0-dev \
    wget \
    tar \
    unzip

RUN apt-get update && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get clean;

# Fix certificate issues
RUN apt-get update && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;
#Setup JAVA_HOME -- useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME

# Maven
WORKDIR /opt
RUN wget http://www-eu.apache.org/dist/maven/maven-3/3.5.3/binaries/apache-maven-3.5.3-bin.tar.gz
RUN tar xzvf apache-maven-3.5.3-bin.tar.gz
RUN cp -rs /opt/apache-maven-3.5.3/bin/* /bin

# BWA
RUN wget --no-check-certificate https://github.com/lh3/bwa/releases/download/v0.7.15/bwakit-0.7.15_x64-linux.tar.bz2
RUN apt-get update && apt-get install -y bzip2
RUN tar -jxvf bwakit-0.7.15_x64-linux.tar.bz2
RUN cp -rs /opt/bwa.kit/bwa /bin
# Kourami
RUN wget https://github.com/Kingsford-Group/kourami/releases/download/v0.9.6/kourami-0.9.6_bin.zip
RUN unzip kourami-0.9.6_bin.zip
RUN /bin/bash -c kourami-0.9.6/scripts/download_panel.sh


WORKDIR /opt/kourami-0.9.6
RUN /opt/apache-maven-3.5.3/bin/mvn install

# Download and index reference: Primary assembly+decoy+alt+hla
RUN /bin/bash -c "scripts/download_grch38.sh hs38NoAltDH"
RUN bwa index resources/hs38NoAltDH.fa

# bamUtil
WORKDIR /opt
RUN wget https://github.com/statgen/bamUtil/archive/v1.0.14.zip
RUN unzip v1.0.14.zip
WORKDIR /opt/bamUtil-1.0.14
RUN apt-get update && apt-get install -y git g++
RUN make cloneLib
RUN make
RUN make install INSTALLDIR=/opt/bamUtil
RUN cp -rs /opt/bamUtil/bam /bin

# Samtools
WORKDIR /opt
RUN wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
RUN tar -jxvf samtools-1.8.tar.bz2
WORKDIR /opt/samtools-1.8
RUN /bin/bash -c "/opt/samtools-1.8/configure --prefix=/opt/samtools --without-curses --disable-bz2 --disable-lzma"
RUN make
RUN make install
RUN cp -rs /opt/samtools/bin/samtools /bin


# Copy Dockerfile
COPY Dockerfile /opt/
