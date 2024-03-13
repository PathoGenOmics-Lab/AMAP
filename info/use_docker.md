
1. Write dockerimage as Dockerfile
```
FROM ubuntu:22.04

RUN apt-get update && apt-get install -y git python3 python3-pip build-essential && \
    pip3 install numpy pandas scikit-learn && \
    git clone https://github.com/lh3/miniprot && \
    cd miniprot && \
    make && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV PATH="${PATH}:/miniprot"
```
2. Start docker
```
sudo service docker start
```
3. Build docker image
```
docker run -it amap /bin/bash
```
4. Use docker image
```
docker run -it amap /bin/bash
```

