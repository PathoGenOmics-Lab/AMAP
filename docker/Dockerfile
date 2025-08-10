# STAGE 1: Build tools that need compilation
FROM continuumio/miniconda3 AS builder

# Install build dependencies, including zlib1g-dev for miniprot compilation.
# We install only what's needed for this stage to keep it light.
RUN apt-get update && apt-get install -y --no-install-recommends \
        git \
        build-essential \
        wget \
        tar \
        bzip2 \
        zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# Set a working directory for the build artifacts
WORKDIR /build

# Compile miniprot
RUN git clone https://github.com/lh3/miniprot && \
    cd miniprot && \
    make

# Download and extract pangene
RUN wget https://github.com/lh3/pangene/releases/download/v1.1/pangene-1.1-bin.tar.bz2 -O- | tar -jxvf -


# STAGE 2: Download large datasets
# This stage is dedicated to downloading data, so it can be cached effectively.
FROM continuumio/miniconda3 AS data_downloader

# Install bakta and wget to download data, specifying a compatible Python version.
RUN conda install -c conda-forge -c bioconda -y python=3.9 bakta=1.11.3 && \
    apt-get update && apt-get install -y --no-install-recommends wget && \
    conda clean -afy && \
    rm -rf /var/lib/apt/lists/*

# Download bakta database
RUN bakta_db download --output /data/bakta_db --type light

# Download reference data
RUN mkdir /reference && \
    wget -P /reference https://zenodo.org/record/3497110/files/MTB_ancestor_reference.fasta && \
    wget -O /reference/aa_H37Rv.translated.fas 'https://docs.google.com/uc?export=download&id=1ibHC4VgZnIbu7EsVZ-VFmcUMCt1ACKuc' && \
    wget -O /reference/bad.bakta.fa 'https://docs.google.com/uc?export=download&id=1sNUuEkgq9Et4sRP3nr-xsMktx1G5rOjn' && \
    wget -O /reference/bad.bakta.translated.fas 'https://docs.google.com/uc?export=download&id=1lDrQPe7S4FJmgejsKyov5jcVbckLqxj9' && \
    wget -O /reference/new_refH37Rv.fasta 'https://docs.google.com/uc?export=download&id=1ZvlTVwp3RVPqZPgIaTONg1IQEm6c98P8' && \
    wget -O /reference/nt_correctedH37Rv.fa 'https://docs.google.com/uc?export=download&id=1_wWCduMeJW4PHV084K8X3BzivaSAmbo8' && \
    wget -O /reference/sequence.gb 'https://docs.google.com/uc?export=download&id=1l0MnOZUWM6t7gYjWgfbt-1SuXOEMi7MQ' && \
    wget -O /reference/snps_lineages.tsv 'https://docs.google.com/uc?export=download&id=1C7BcWzpLnbFwMA-uTyiD_FSi-LkVSHAj' && \
    wget -O /reference/Mycobacterium_tuberculosis_H37Rv_txt_v4.txt 'https://docs.google.com/uc?export=download&id=1_yWmTKd3FRykoT7RYI74ksyDm9gyHmpM'


# -----------------------------------------------------------------------------
# STAGE 3: Final production image (REVISED VERSION)
# -----------------------------------------------------------------------------
FROM continuumio/miniconda3

ENV MPLCONFIGDIR=/app/matplotlib_config
ENV PATH="/opt/conda/bin:/miniprot:/data/bakta_db:/reference:${PATH}"

RUN mkdir -p /miniprot /data/bakta_db /reference /app/matplotlib_config && \
    chmod -R 777 /app/matplotlib_config

# Copy artifacts from previous stages
COPY --from=builder /build/miniprot/miniprot /miniprot/
COPY --from=builder /build/pangene-1.1-bin/bin_x64-linux/* /usr/local/bin/
COPY --from=builder /build/pangene-1.1-bin/scripts/pangene.js /usr/local/bin/
RUN chmod +x /miniprot/miniprot /usr/local/bin/*
COPY --from=data_downloader /data/bakta_db /data/bakta_db
COPY --from=data_downloader /reference /reference

# --- UNIFIED PACKAGE INSTALLATION WITH CONDA ---
# We use Conda to install EVERYTHING: Python, R, and their respective packages.
# This ensures compatibility and proper dependency resolution.
RUN \
    # 1. Configure Conda channels
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict && \
    \
    # 2. Install Python and R packages with Conda in a single step.
    #    Conda will handle system dependencies (like cairo, libxml2, etc.) for us.
    conda install -y \
      # Base environment
      python=3.9 \
      r-base=4.2 \
      r-essentials \
      # Python packages
      pandas \
      biopython \
      matplotlib=3.8.3 \
      eggnog-mapper=2.1.13 \
      bakta=1.11.3 \
      # R packages (using the r- prefix)
      r-ggplot2 \
      r-dplyr \
      r-rmarkdown \
      r-prettydoc \
      r-remotes \
      r-plotly \
      r-ggpubr \
      r-dt \
      r-stringr \
      r-knitr \
      r-rstatix \
      r-ggrepel \
      r-scales \
      pandoc && \
    \
    # 3. Install Python packages with Pip that are not on Conda
    pip install --no-cache-dir \
      intervaltree \
      python-louvain \
      pyvis==0.3.2 \
      setuptools \
      networkx && \
    \
    # 4. Clean up caches to reduce the final image size
    conda clean -afy && \
    rm -rf /root/.cache/pip

# --- R PACKAGE INSTALLATION FROM GITHUB (IF NEEDED) ---
# This runs after the Conda R environment is fully set up.
# Uncomment the line below only if you truly need the development version of plotly.
 RUN R -e "remotes::install_github('maximilian-heeg/ggsimplestats', force = TRUE)"
# RUN R -e "remotes::install_github('ropensci/plotly', force = TRUE)"

# Set the final working directory
WORKDIR /app
