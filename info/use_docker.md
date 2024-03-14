
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

Para subir tu imagen de Docker a Docker Hub, sigue estos pasos:

Inicia Sesión en Docker Hub:
Abre la terminal y ejecuta docker login.
Ingresa tu nombre de usuario y contraseña de Docker Hub.

Etiqueta tu Imagen:
Antes de subir, etiqueta tu imagen localmente con tu nombre de usuario de Docker Hub y el nombre de la imagen. Por ejemplo:

```
docker tag amap tu_usuario/amap:latest
```

Sube la Imagen:

Usa el siguiente comando para subir tu imagen a Docker Hub:
```
docker push tu_usuario/amap:latest
```
De nuevo, reemplaza tu_usuario con tu nombre de usuario en Docker Hub.

Una vez que la imagen se haya subido correctamente, podrás verla y gestionarla en tu repositorio de Docker Hub.