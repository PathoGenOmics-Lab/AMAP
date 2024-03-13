FROM ubuntu:22.04

RUN apt-get update && apt-get install -y git python3 python3-pip build-essential && \
    pip3 install numpy pandas scikit-learn && \
    git clone https://github.com/lh3/miniprot && \
    cd miniprot && \
    make && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV PATH="${PATH}:/miniprot"