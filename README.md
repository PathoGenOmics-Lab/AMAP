# AMAP: Assisted MTBC annotator pipeline

Execute pipeline:  
``
snakemake --use-conda -c8
``

## DAG (Directed Acyclic Graph)

<p align="center">
  <img src="https://github.com/PathoGenOmics-Lab/AMAP/blob/main/images/dag.png" alt="DAGPLOT" width="400">
</p>


### Generate Rule graph

``
snakemake --rulegraph | dot -Tpng > rulegraph.png
``

### Generate DAG (Directed Acyclic Graph)

``
snakemake --dag | dot -Tpng > dag.png
``

### Generate Filegraph

``
snakemake --filegraph | dot -Tpng > filegraph.png
``
