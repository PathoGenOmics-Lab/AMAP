
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
docker build -t amap .
```
4. Use docker image
```
docker run -it amap /bin/bash
```

To upload your Docker image to Docker Hub, follow these steps:

Log in to Docker Hub:
Open the terminal and run docker login.
Enter your Docker Hub username and password.

Tag Your Image:
Before uploading, tag your local image with your Docker Hub username and the name of the image. For example:

```
docker tag amap paururo/amap:latest
```

Upload the Image:

Use the following command to upload your image to Docker Hub:
```
docker push paururo/amap:latest
```
Again, replace with your username on Docker Hub.

Once the image has been successfully uploaded, you will be able to view and manage it in your Docker Hub repository.