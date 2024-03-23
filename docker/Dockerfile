FROM continuumio/miniconda3

RUN apt-get update && \
    apt-get install -y git build-essential zlib1g-dev r-base && \
    conda install -c conda-forge python=3.9 pandas plotly && \
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

RUN mkdir -p /reference

RUN wget -P /reference https://zenodo.org/record/3497110/files/MTB_ancestor_reference.fasta

RUN apt-get update && \
    apt-get install -y python3-pip && \
    pip install gdown && \
    gdown --id 1sNUuEkgq9Et4sRP3nr-xsMktx1G5rOjn -O /reference/bad.bakta.fa && \
    gdown --id 1ibHC4VgZnIbu7EsVZ-VFmcUMCt1ACKuc -O /reference/aa_H37Rv.translated.fas && \
    gdown --id 1ZvlTVwp3RVPqZPgIaTONg1IQEm6c98P8 -O /reference/new_refH37Rv.fasta && \
    gdown --id 1_wWCduMeJW4PHV084K8X3BzivaSAmbo8 -O /reference/nt_correctedH37Rv.fa && \
    gdown --id 1l0MnOZUWM6t7gYjWgfbt-1SuXOEMi7MQ -O /reference/sequence.gb && \
    gdown --id 1C7BcWzpLnbFwMA-uTyiD_FSi-LkVSHAj -O /reference/snps_lineages.tsv && \
    gdown --id 1lDrQPe7S4FJmgejsKyov5jcVbckLqxj9 -O /reference/bad.bakta.translated.fas

ENV PATH="${PATH}:/miniprot:/data/bakta_db:/reference"