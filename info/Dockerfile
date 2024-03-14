FROM continuumio/miniconda3

RUN apt-get update && \
    apt-get install -y git build-essential zlib1g-dev && \
    conda install -c conda-forge -c bioconda bakta && \
    git clone https://github.com/lh3/miniprot && \
    cd miniprot && \
    make && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV PATH="${PATH}:/miniprot"