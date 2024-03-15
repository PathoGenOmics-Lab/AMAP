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

ENV PATH="${PATH}:/miniprot"