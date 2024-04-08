FROM continuumio/miniconda3

# Establecer variables de entorno para Quarto
ENV XDG_CACHE_HOME=/tmp/quarto_cache_home
ENV XDG_DATA_HOME=/tmp/quarto_data_home

# Configurar la variable de entorno para Matplotlib
ENV MPLCONFIGDIR=/app/matplotlib_config
ENV PATH="${PATH}:/miniprot:/data/bakta_db:/reference"

# Crear directorio de configuración para Matplotlib y otros directorios necesarios
RUN mkdir -p /app/matplotlib_config /data/bakta_db /reference && \
    chmod -R 777 /app/matplotlib_config

# Instalación de paquetes y limpieza
RUN apt-get update && \
    apt-get install -y git build-essential zlib1g-dev r-base pandoc libxml2-dev libcurl4-openssl-dev libssl-dev && \
    conda install -c conda-forge python=3.9 pandas plotly biopython  matplotlib=3.8.3 && \
    conda install -c bioconda -c conda-forge eggNOG-mapper=2.1.12 bakta=1.9.3 && \
    pip install intervaltree python-louvain pyvis==0.3.2 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Instalación de paquetes de R
RUN R -e "install.packages(c('ggplot2', 'dplyr', 'rmarkdown', 'prettydoc', 'remotes'), repos='http://cran.rstudio.com/')" && \
    R -e "remotes::install_github('ropensci/plotly')"
    
RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.553/quarto-1.4.553-linux-amd64.deb && \
    dpkg -i quarto-1.4.553-linux-amd64.deb && \
    rm quarto-1.4.553-linux-amd64.deb
    
# Clonar y construir miniprot
RUN git clone https://github.com/lh3/miniprot && \
    cd miniprot && \
    make
RUN pip install --upgrade setuptools networkx
# Descargar datos de referencia
RUN wget -P /reference https://zenodo.org/record/3497110/files/MTB_ancestor_reference.fasta && \
    wget -O /reference/aa_H37Rv.translated.fas 'https://docs.google.com/uc?export=download&id=1ibHC4VgZnIbu7EsVZ-VFmcUMCt1ACKuc' && \
    wget -O /reference/bad.bakta.fa 'https://docs.google.com/uc?export=download&id=1sNUuEkgq9Et4sRP3nr-xsMktx1G5rOjn' && \
    wget -O /reference/bad.bakta.translated.fas 'https://docs.google.com/uc?export=download&id=1lDrQPe7S4FJmgejsKyov5jcVbckLqxj9' && \
    wget -O /reference/new_refH37Rv.fasta 'https://docs.google.com/uc?export=download&id=1ZvlTVwp3RVPqZPgIaTONg1IQEm6c98P8' && \
    wget -O /reference/nt_correctedH37Rv.fa 'https://docs.google.com/uc?export=download&id=1_wWCduMeJW4PHV084K8X3BzivaSAmbo8' && \
    wget -O /reference/sequence.gb 'https://docs.google.com/uc?export=download&id=1l0MnOZUWM6t7gYjWgfbt-1SuXOEMi7MQ' && \
    wget -O /reference/snps_lineages.tsv 'https://docs.google.com/uc?export=download&id=1C7BcWzpLnbFwMA-uTyiD_FSi-LkVSHAj' && \
    wget -O /reference/Mycobacterium_tuberculosis_H37Rv_txt_v4.txt 'https://docs.google.com/uc?export=download&id=1_yWmTKd3FRykoT7RYI74ksyDm9gyHmpM'

# Descargar base de datos bakta
RUN bakta_db download --output /data/bakta_db --type light
