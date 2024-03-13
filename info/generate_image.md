For generating a singularity image

1. Generate the def file
```
Bootstrap: docker
From: ubuntu:20.04

%post
    apt-get update
    apt-get install -y python3 python3-pip
    pip3 install numpy pandas scikit-learn
    apt-get clean
    rm -rf /var/lib/apt/lists/*
```

2. Generate .sif image
```
sudo singularity build AMAP.sif AMAP.def
```

3. Try the image
```
singularity exec AMAP.sif python3 -c "import numpy; print(numpy.__version__)"
```