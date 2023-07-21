# AMAP: Assisted MTBC annotator pipeline


## DAG (Directed Acyclic Graph)


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
