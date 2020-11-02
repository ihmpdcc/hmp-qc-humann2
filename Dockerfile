###############################################################
# Dockerfile to build container mWGS processing pipeline image
# This pipeline runs KneadData (QC) and HUMAnN2
###############################################################

FROM ubuntu:18.04

MAINTAINER Kemi Ifeonu <kifeonu@som.umaryland.edu>

#############
## General ##
#############

ARG SRATOOLKIT_DOWNLOAD_URL=http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
ARG METAPHLAN_DOWNLOAD_URL=https://github.com/biobakery/MetaPhlAn/archive/2.7.8.tar.gz
ARG METAPHLAN_DB_URL=https://www.dropbox.com/sh/7qze7m7g9fe2xjg/AAASBOj-2gAbA53cV1bXBULYa/mpa_v20_m200.tar

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update && apt-get install -y build-essential autoconf libtool pkg-config

RUN apt-get -qq update && apt-get -qq install -y --no-install-recommends \
    bowtie2 \
    curl \
    default-jre \
    jq \
    python-pip \
    python-dev \
    python-setuptools \
    python3-dev \
    python3-pip \
    python3-setuptools \
    python3.6 \
    python2.7 \
    python3-setuptools \
    unzip \
    wget

RUN pip3 install --upgrade pip
RUN pip3 install --upgrade awscli
RUN pip2 install numpy cython
RUN pip install numpy cython
RUN pip install biopython cython matplotlib scipy==1.3.1 biom-format h5py
RUN pip2 install biopython biom-format==2.1.7
RUN pip3 install boto3

# For reference databases needed by tools
RUN mkdir /dbs



##################
## fasterq-dump ##
##################

RUN wget -O /opt/sratoolkit.tar.gz $SRATOOLKIT_DOWNLOAD_URL \
    && tar -xzf /opt/sratoolkit.tar.gz --directory /opt \
    && rm /opt/sratoolkit.tar.gz \
    && mv /opt/sratoolkit* /opt/sratoolkit

ENV PATH $PATH:/opt/sratoolkit/bin


###############
## KneadData ##
###############

RUN pip3 install kneaddata

RUN mkdir /dbs/kneaddata /kneaddata_output

#Download database
RUN kneaddata_database --download human_genome bowtie2 /dbs/kneaddata



#############
## HUMAnN2 ##
#############

RUN pip2 install humann2

#MetaPhlAn2 is required for humann2
RUN wget -O /opt/metaphlan2.tar.gz $METAPHLAN_DOWNLOAD_URL \
    && tar -xzf /opt/metaphlan2.tar.gz --directory /opt \
    && rm /opt/metaphlan2.tar.gz \
    && ln -s /opt/MetaPhlAn-2.7.8/* /usr/local/bin/

#Download databases
RUN mkdir /dbs/humann2
RUN humann2_databases --download chocophlan full /dbs/humann2/
RUN humann2_databases --download uniref uniref90_diamond /dbs/humann2/
RUN humann2_databases --download utility_mapping full /dbs/humann2/

RUN mkdir /dbs/humann2/metaphlan
RUN mkdir /dbs/humann2/metaphlan/mpa_v20_m200
RUN wget $METAPHLAN_DB_URL \
    && tar -xvf mpa_v20_m200.tar --directory /dbs/humann2/metaphlan \
    && bunzip2 /dbs/humann2/metaphlan/mpa_v20_m200.fna.bz2 \ 
    && rm mpa_v20_m200.tar \
    && bowtie2-build /dbs/humann2/metaphlan/mpa_v20_m200.fna /dbs/humann2/metaphlan/mpa_v20_m200 
RUN rm /dbs/humann2/metaphlan/mpa_v20_m200.fna
RUN mv /dbs/humann2/metaphlan/mpa_v20_m200.* /dbs/humann2/metaphlan/mpa_v20_m200/ 



##############
## PIPELINE ##
##############
COPY execute_pipeline.py /opt/scripts/execute_pipeline.py
RUN chmod 755 /opt/scripts/execute_pipeline.py

# Change to root directory
WORKDIR /
RUN mkdir /root/.ncbi
COPY user-settings.mkfg /root/.ncbi
RUN chmod 755 /root/.ncbi/user-settings.mkfg


ENTRYPOINT ["python3","-u","/opt/scripts/execute_pipeline.py"]
#ENTRYPOINT ["bash","/opt/scripts/execute_pipeline.sh"]

