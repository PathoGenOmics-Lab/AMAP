FROM continuumio/miniconda3

RUN apt-get update && \
    apt-get install -y git build-essential zlib1g-dev r-base && \
    conda install -c conda-forge python=3.9 pandas plotly biopython python-louvain pyvis && \
    conda install -c bioconda -c conda-forge blast=2.15.0 && \ 
    pip install intervaltree  && \
    conda install -c bioconda -c conda-forge bakta=1.9.3 && \
    git clone https://github.com/lh3/miniprot && \
    cd miniprot && \
    make && \
    R -e "install.packages(c('ggplot2', 'dplyr'), repos='http://cran.rstudio.com/')" && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /data/bakta_db && \
    bakta_db download --output /data/bakta_db --type light

RUN mkdir -p /reference && \
    wget -P /reference https://zenodo.org/record/3497110/files/MTB_ancestor_reference.fasta && \
    wget -O /reference/aa_H37Rv.translated.fas 'https://docs.google.com/uc?export=download&id=1ibHC4VgZnIbu7EsVZ-VFmcUMCt1ACKuc'  && \
    wget -O /reference/bad.bakta.fa 'https://docs.google.com/uc?export=download&id=1sNUuEkgq9Et4sRP3nr-xsMktx1G5rOjn' && \
    wget -O /reference/bad.bakta.translated.fas 'https://docs.google.com/uc?export=download&id=1lDrQPe7S4FJmgejsKyov5jcVbckLqxj9' && \
    wget -O /reference/new_refH37Rv.fasta 'https://docs.google.com/uc?export=download&id=1ZvlTVwp3RVPqZPgIaTONg1IQEm6c98P8' && \
    wget -O /reference/nt_correctedH37Rv.fa 'https://docs.google.com/uc?export=download&id=1_wWCduMeJW4PHV084K8X3BzivaSAmbo8' && \
    wget -O /reference/sequence.gb 'https://docs.google.com/uc?export=download&id=1l0MnOZUWM6t7gYjWgfbt-1SuXOEMi7MQ' && \
    wget -O /reference/snps_lineages.tsv 'https://docs.google.com/uc?export=download&id=1C7BcWzpLnbFwMA-uTyiD_FSi-LkVSHAj' && \
    wget -O /reference/Mycobacterium_tuberculosis_H37Rv_txt_v4.txt 'https://docs.google.com/uc?export=download&id=1_yWmTKd3FRykoT7RYI74ksyDm9gyHmpM'

ENV PATH="${PATH}:/miniprot:/data/bakta_db:/reference"