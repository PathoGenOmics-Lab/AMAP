# =============================================================================
# AMAP v1.0.0 — Multi-stage Docker build
# =============================================================================
# Stage 1: Compile native tools (miniprot)
# Stage 2: Download databases (Bakta DB, reference files)
# Stage 3: Final image with all runtime dependencies
#
# NOTE on reference files: Several files are hosted on Google Drive. If the
# build fails on those wget calls, the links may have expired. Replace with
# Zenodo/Figshare URLs when available. The MTB_ancestor_reference.fasta is
# already on Zenodo and should always work.
# =============================================================================

# ---------------------------------------------------------------------------
# STAGE 1: Build tools that need compilation
# ---------------------------------------------------------------------------
FROM --platform=linux/amd64 continuumio/miniconda3:24.1.2-0 AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        git build-essential wget tar bzip2 zlib1g-dev \
        curl pkg-config libssl-dev && \
    rm -rf /var/lib/apt/lists/*

# Install Rust toolchain
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

WORKDIR /build

# Compile pathotypr
RUN git clone --depth 1 https://github.com/PathoGenOmics-Lab/pathotypr && \
    cd pathotypr && cargo build --release && \
    cp target/release/pathotypr /usr/local/bin/

# Compile miniprot v0.18
RUN git clone --branch v0.18 --depth 1 https://github.com/lh3/miniprot && \
    cd miniprot && make

# Download pangene v1.1 pre-built binaries
RUN wget -q https://github.com/lh3/pangene/releases/download/v1.1/pangene-1.1-bin.tar.bz2 -O- | tar -jxf -

# ---------------------------------------------------------------------------
# STAGE 2: Download large datasets (cached independently)
# ---------------------------------------------------------------------------
FROM --platform=linux/amd64 continuumio/miniconda3:24.1.2-0 AS data_downloader

RUN conda install -c conda-forge -c bioconda -y python=3.9 bakta=1.12.0 && \
    apt-get update && apt-get install -y --no-install-recommends wget && \
    conda clean -afy && \
    rm -rf /var/lib/apt/lists/*

# Bakta light database
RUN bakta_db download --output /data/bakta_db --type light

# Reference files for MTBC mode
RUN mkdir /reference && \
    wget -q -P /reference https://zenodo.org/record/3497110/files/MTB_ancestor_reference.fasta && \
    wget -q -O /reference/aa_H37Rv.translated.fas 'https://docs.google.com/uc?export=download&id=1ibHC4VgZnIbu7EsVZ-VFmcUMCt1ACKuc' && \
    wget -q -O /reference/bad.bakta.fa 'https://docs.google.com/uc?export=download&id=1sNUuEkgq9Et4sRP3nr-xsMktx1G5rOjn' && \
    wget -q -O /reference/bad.bakta.translated.fas 'https://docs.google.com/uc?export=download&id=1lDrQPe7S4FJmgejsKyov5jcVbckLqxj9' && \
    wget -q -O /reference/new_refH37Rv.fasta 'https://docs.google.com/uc?export=download&id=1ZvlTVwp3RVPqZPgIaTONg1IQEm6c98P8' && \
    wget -q -O /reference/nt_correctedH37Rv.fa 'https://docs.google.com/uc?export=download&id=1_wWCduMeJW4PHV084K8X3BzivaSAmbo8' && \
    wget -q -O /reference/sequence.gb 'https://docs.google.com/uc?export=download&id=1l0MnOZUWM6t7gYjWgfbt-1SuXOEMi7MQ' && \
    wget -q -O /reference/snps_lineages.tsv 'https://docs.google.com/uc?export=download&id=1C7BcWzpLnbFwMA-uTyiD_FSi-LkVSHAj' && \
    wget -q -O /reference/Mycobacterium_tuberculosis_H37Rv_txt_v4.txt 'https://docs.google.com/uc?export=download&id=1_yWmTKd3FRykoT7RYI74ksyDm9gyHmpM'

# ---------------------------------------------------------------------------
# STAGE 3: Final production image
# ---------------------------------------------------------------------------
FROM continuumio/miniconda3:24.1.2-0

LABEL maintainer="Paula Ruiz-Rodriguez <paururo@github>" \
      version="1.0.0" \
      description="AMAP: Assisted Microbial Annotator Pipeline"

ENV MPLCONFIGDIR=/app/matplotlib_config
ENV PATH="/opt/conda/bin:/miniprot:/data/bakta_db:/reference:${PATH}"

RUN mkdir -p /miniprot /data/bakta_db /reference /app/matplotlib_config && \
    chmod -R 777 /app/matplotlib_config

# Copy compiled tools and databases from previous stages
COPY --from=builder /build/miniprot/miniprot /miniprot/
COPY --from=builder /usr/local/bin/pathotypr /usr/local/bin/pathotypr
COPY --from=builder /build/pangene-1.1-bin/bin_x64-linux/* /usr/local/bin/
COPY --from=builder /build/pangene-1.1-bin/scripts/pangene.js /usr/local/bin/
RUN chmod +x /miniprot/miniprot /usr/local/bin/*
COPY --from=data_downloader /data/bakta_db /data/bakta_db
COPY --from=data_downloader /reference /reference

# --- CONDA: Python + R + bioinformatics tools ---
# Single solve step for reproducibility. No r-essentials (pulls ~150 unused pkgs).
RUN \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority flexible && \
    \
    conda install -y \
      # --- Python runtime ---
      python=3.9 \
      numpy \
      pandas \
      biopython \
      matplotlib \
      plotly \
      # --- Bioinformatics tools ---
      bakta=1.12.0 \
      eggnog-mapper=2.1.13 \
      diamond=2.1.24 \
      gfatools \
      # --- R runtime (only packages used by report_AMAP.Rmd) ---
      r-base=4.4 \
      r-ggplot2 \
      r-dplyr \
      r-tidyr \
      r-stringr \
      r-scales \
      r-rmarkdown \
      r-knitr \
      r-prettydoc \
      r-plotly \
      r-ggpubr \
      r-ggrepel \
      r-rstatix \
      r-dt \
      r-pheatmap \
      r-rcolorbrewer \
      r-remotes \
      pandoc && \
    \
    # --- Pip: packages not on Conda ---
    pip install --no-cache-dir \
      intervaltree \
      python-louvain \
      pyvis==0.3.2 \
      networkx && \
    \
    # --- Cleanup ---
    conda clean -afy && \
    rm -rf /root/.cache/pip /tmp/*

# --- R: GitHub-only packages ---
RUN R -e "remotes::install_github('maximilian-heeg/ggsimplestats', force = TRUE)"

WORKDIR /app
