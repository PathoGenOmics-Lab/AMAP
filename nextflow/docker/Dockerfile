FROM continuumio/miniconda3

RUN apt-get update && \
    apt-get install -y git build-essential zlib1g-dev && \
    conda install mamba -c conda-forge && \
    mamba install python=3.10 -c conda-forge && \
    mamba install -c bioconda -c conda-forge bakta=1.9.3 && \
    git clone https://github.com/lh3/miniprot && \
    cd miniprot && \
    make && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* 


RUN mkdir -p /data/bakta_db && \
    bakta_db download --output /data/bakta_db --type light

# Crear directorio para los archivos
RUN mkdir -p /reference

# Descargar archivos
RUN wget -P /reference https://zenodo.org/record/8024447/files/nt_correctedH37Rv.fa && \
    wget -P /reference https://zenodo.org/record/8024447/files/new_refH37Rv.fasta && \
    wget -P /reference https://zenodo.org/record/8024447/files/bad.bakta.fa && \
    wget -P /reference https://zenodo.org/record/8024447/files/sequence.gb && \
    wget -P /reference https://zenodo.org/record/8232692/files/snps_lineages.tsv && \
    wget -P /reference https://zenodo.org/record/3497110/files/MTB_ancestor_reference.fasta

ENV PATH="${PATH}:/miniprot:/data"


ENV PATH="${PATH}:/miniprot:/data/bakta_db:/reference"
