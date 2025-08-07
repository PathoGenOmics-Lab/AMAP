[![PGO](https://img.shields.io/badge/PathoGenOmics-Lab-grey?style=plastic&labelColor=red)](https://github.com/PathoGenOmics-Lab)
[![AMAP](https://img.shields.io/badge/AMAP-v0.1.0-green?style=plastic)](https://github.com/PathoGenOmics-Lab/AMAP)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A5%20v23.10.1-%2330ba91?style=plastic)](https://github.com/nextflow-io/nextflow/releases/tag/v23.10.1)
[![License: GPL v3](https://img.shields.io/badge/License%3A-GPL%20v3-orange?style=plastic)](https://github.com/PathoGenOmics-Lab/AMAP/blob/main/LICENSE)
[![Docker Image](https://img.shields.io/badge/Docker%20Image-AMAP-blue?style=plastic&logo=docker)](https://hub.docker.com/repository/docker/paururo/amap/general)
[![Zenodo](https://img.shields.io/badge/doi%3A-10.5281%2Fzenodo.8024447-%23c9ac02?style=plastic)](https://doi.org/10.5281/zenodo.8024447)

# AMAP: Assisted MTBC annotator pipeline  
__Paula Ruiz-Rodriguez<sup>1</sup>__ 
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. I<sup>2</sup>SysBio, University of Valencia-CSIC, FISABIO Joint Research Unit Infection and Public Health, Valencia, Spain </sub>  

<p align="center">
  <img src=".github/AMAP.png" title="AMAP logo" style="width:650px; height: auto;">
</p>

AMAP is a bioinformatics pipeline written in Nextflow for processing, annotating and correct annotations of _Mycobacterium tuberculosis_ complex (MTBC) genomes.
## First Steps for Using AMAP

**Why Choose Nextflow with Conda and Singularity?** 
Nextflow excels in managing complex workflows, offering reproducibility, scalability, and portability. It's particularly effective for bioinformatics and data analysis tasks. Nextflow's support for containerization with Docker and Singularity ensures consistent and isolated execution environments. Conda complements this by simplifying the package management process.

### Installation Guide

To set up your environment with Nextflow and Singularity, follow these straightforward steps:

1. **Install Nextflow Using Conda**:
   Begin by downloading Miniconda, a lightweight version of Conda, which simplifies package management. You can find it at [Miniconda's official site](https://docs.conda.io/en/latest/miniconda.html).
2. **Create and Activate a Conda Environment**:
   Organize your packages and dependencies in an isolated environment to avoid conflicts.
   ```bash
   conda create -n nextflow-env
   conda activate nextflow-env
3. **Install Nextflow in Conda Environment**:
   With the Conda environment ready, install Nextflow seamlessly.
   ```bash
   conda install -c bioconda nextflow
   ```
4. **Install Singularity**:
   To handle containerization, install Singularity following the detailed instructions in the [Singularity Installation Guide](https://docs.sylabs.io/guides/4.1/user-guide/quick_start.html).
5. **Verify Installation**:
   Confirm that Nextflow and Singularity are correctly installed and ready to use.
   ```bash
   nextflow -version
   singularity --version
   ```
These steps ensure that you have a robust setup for running reproducible data analysis workflows using Nextflow, Conda, and Singularity.

### Running AMAP
Nextflow facilitates the direct execution of pipelines from GitHub repositories. This approach simplifies the process and does not require cloning the repository locally.

To run the AMAP pipeline hosted on GitHub by PathoGenOmics-Lab, follow these simple steps:

1. **Prepare Your Environment**:
   Ensure that Nextflow and Singularity are installed and properly configured in your Conda environment as previously detailed.

2. **Run the Pipeline**:
   Execute the AMAP pipeline directly from GitHub with Nextflow. Replace `[parameters]` with the specific parameters needed for your analysis.
   ```bash
   nextflow run PathoGenOmics-Lab/AMAP --tsv test_file.tsv
   ```

## Parameters

## Inputs

## Outputs
### MTB_TYPER
MTB typer uses a k-mear aproach to use robust phylogenetic markers:
- /fastas/*.fasta
- /typing/*.
### BAKTA
All inputs that Bakta software offers to annotate a genome, generates a directory by genome:
- /genome_name/*.faa
- /genome_name/*.gbff
- /genome_name/*.gff3
- /genome_name/*.hypotheticals.faa
- /genome_name/*.hypotheticals.tsv
- /genome_name/*.png
- /genome_name/*.tsv
- /genome_name/*.txt
### CDS_CASES
- /genome_name/*-case1.tsv:
- /genome_name/*-case2.tsv:
- /genome_name/*-case3.tsv:
- /genome_name/*-case4.tsv:
- /genome_name/*-rvmissing.tsv:
- /genome_name/*-single_overlap.html:
- /genome_name/*-multiple_overlap.html:
- /summary/*-summary.tsv
### AMAP
- /genome_name/*.AMAP.tsv
- /genome_name/*.AMAP.COG.tsv
- /genome_name/*.AMAP.faa
- /genome_name/*.AMAP.fna
- /genome_name/*.AMAP.MycoBrowser.tsv
- /genome_name/*.AMAP.gbk
- /genome_name/*.AMAP.gff3
- /genome_name/*.AMAPmastertable.tsv
### DUPLICATIONS
- /duplications_*
- /duplications_*_groups.tsv
- /duplications_*_network.png
- /duplications_*_interactive.html
- /duplications_*.AMAP.faa.renamed.dup
### FINAL REPORT
- homologous_communities.html
- homologous_communities.tsv
### pipeline_info
- execution_report_*.html
- execition_timeline_*.html
- pipeline_dag_*.html
- execution_trace_*.tsv

#### Network measures of Duplications
1. **Degree Centrality** ``(AverageDegree1-2)``: Measures the number of direct connections a node has. In networks based on pairwise relationships, this measure is crucial for identifying highly connected nodes. Simply put,it is the average number of connections each node in the community has.
2. **Closeness Centrality** ``(AverageCloseness1-2)``: Calculates the average distance of a node to all others in the network, identifying nodes capable of efficiently propagating information. This metric measures how'close' nodes are to each other in terms of steps or hops. A high value indicates that nodes are, on average, closer to each other.
3. **Betweenness Centrality** ``(AverageBetweenness1-2)``: Indicates how often a node acts as a bridge on the shortest path between two other nodes, thus controlling the flow of information. A high value may indicate nodesthat act as 'bridges' between different parts of the network.
4. **Eigen Centrality** ``(Eigen1-2)``: Values the number of connections of a node and the importance of the connected nodes, useful for identifying influential nodes in pairwise networks.
5. **Katz Centrality** ``(Katz1-2)``: Similar to Eigen Centrality, but it considers nodes at greater distances with diminishing weight, useful for identifying influential nodes considering indirect connections.
6. **Network Density** ``(Density1-2)``: The ratio of possible connections that are actual, indicating the degree of interconnection of the network. A density close to 1 indicates a highly connected network, while a lowdensity indicates a more dispersed network.
7. **Clustering Coefficient** ``(ClusteringCoefficient1-2)``: Measures the degree to which nodes cluster together, indicating a well-defined community structure. This metric measures the extent to which nodes in a networktend to cluster or form closed groups. A high value suggests that if two nodes are connected to a third one, they are likely also connected to each other.
8. **Average Path Length** ``(AveragePathLength1-2)``: The average distance between all pairs of nodes, with short lengths indicating a "small world" distance between any pair of nodes. It is the average number of stepsalong the shortest paths for all possible pairs of nodes. A low average path length means that most nodes can be reached from other nodes in a few steps.
9. **Transitivity** ``(transitivity1-2)``: Measures the likelihood that the neighbors of a node are also connected to each other. It is similar to the clustering coefficient but applies to the entire network. It indicatesthe general tendency of the network to form triangles or "three connected friends."
10. **Modularity** ``(modularity1-2)``: Evaluates the community structure within a network. High modularity indicates that the network has densely interconnected communities, but with few connections between these distinct communities. It helps to identify how well defined the communities are within the network.

## AMAP Workflow
```mermaid
%%{init: { 'theme': 'dark' } }%%
%%{init: { 'theme': 'dark' } }%%
flowchart TB
    subgraph "INPUT"
    v0["GENOMES"]
    end
    v3([PATHOTYPR])
    v4([BAKTA ANNOTATOR])
    v6([MINIPROT])
    v8([TYPE CASES])
    v14([FIX ANNOTATION])
    v20([GET 
    DUPLICATIONS])
    v17([EGGNOG DB])
    v18([EGGNOG 
    ANNOTATION])
    v19([FUNCTIONAL
    ANNOTATION])
    v16["FINAL RUN 
    REPORT"]
    v22([PANGENE])
    subgraph "RESULTS"
    v12["OG Bakta annotation 
    /BAKTA"]
    v10["Report cases
    /CDS_CASES"]
    v13["Corrected annotation
    /AMAP"]
    v23["Duplicated genes
    /DUPLICATIONS"]
    v21["Pangenome"
    /PANGENE]
    v15["Lineage Typing
    /PATHOTYPR"]
    
    end
    v3 --> v4
    v0 --> v3
    v4 --> v12
    v4 --> v6
    v3 --> v6
    v4 --> v8
    v6 --> v8
    v8 --> v10
    v8 --> v14
    v14 --> v13
    v20 --> v13
    v19 --> v13
    v3 --> v15
    v10 --> v16
    v13 --> v16
    v15 --> v16
    v17 --> v18
    v14 --> v19
    v18 --> v19
    v14 --> v20
    v14 --> v22
    v22 --> v21
    v20 --> v23


```

## AMAP Algorithm

## Citation
>Ruiz-Rodriguez, P., & Coscolla, M. (2025). AMAP (Version 0.1.0) [Computer software]. https://doi.org/10.5281/zenodo.1234
```bibtex
@software{Ruiz-Rodriguez_AMAP_2025,
  author = {Ruiz-Rodriguez, Paula and Coscolla, Mireia},
  doi = {10.5281/zenodo.1234},
  license = {GPL-3.0},
  month = jul,
  title = {{AMAP}},
  url = {https://github.com/PathoGenOmics-Lab/AMAP},
  version = {0.1.0},
  year = {2025}}
```
---
<h2 id="contributors" align="center">

✨ [Contributors]((https://github.com/PathoGenOmics-Lab/AMAP/graphs/contributors))
</h2>

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<div align="center">
AMAP is developed with ❤️ by:
<table>
  <tr>
    <td align="center">
      <a href="https://github.com/paururo">
        <img src="https://avatars.githubusercontent.com/u/50167687?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Paula Ruiz-Rodriguez</b></sub>
      </a>
      <br />
      <a href="" title="Code">💻</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Data">🔣</a>
      <a href="" title="Desing">🎨</a>
      <a href="" title="Tool">🔧</a>
    </td> 
    <td align="center">
      <a href="https://github.com/mireiacoscolla">
        <img src="https://avatars.githubusercontent.com/u/29301737?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Mireia Coscolla</b></sub>
      </a>
      <br />
      <a href="https://www.uv.es/instituto-biologia-integrativa-sistemas-i2sysbio/es/investigacion/proyectos/proyectos-actuales/mol-tb-host-1286169137294/ProjecteInves.html?id=1286289780236" title="Funding/Grant Finders">🔍</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Mentoring">🧑‍🏫</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="User Testing">📓</a>
    </td> 
  </tr>
</table>

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification ([emoji key](https://allcontributors.org/docs/en/emoji-key)).

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
---  
<h2 id="contributors" align="center">
💻 Software Used
</h2>

| Software                                                  | Version | Type     | Usage                           |
|-----------------------------------------------------------|---------|----------|---------------------------------|
| [Nextflow](https://github.com/nextflow-io/nextflow)       | 23.10.1 | Nextflow | Workflow Management             |
| [Bakta](https://github.com/oschwengers/bakta)             | 1.9.3   | Other    | Annotation Tools                |
| [miniprot](https://github.com/lh3/miniprot)               | 0.13    | Other    | Mapping Tools                   |
| [Python](https://www.python.org/)                         | 3.9     | Python   | Programming Language            |
| [R](https://www.R-project.org/)                           | 4.3.3   | R        | Statistical Computing           |
| [pandas](https://github.com/pandas-dev/pandas)            | 2.2.1   | Python   | Data Analysis                   |
| [intervaltree](https://github.com/chaimleib/intervaltree) | 3.1.0   | Python   | Data Structures                 |
| [plotly.py](https://github.com/plotly/plotly.py)          | 5.20.0  | Python   | Data Visualization              |
| [ggplot2](https://ggplot2.tidyverse.org/index.html)       | 3.5.0   | R        | Data Visualization              |
| [dplyr](https://dplyr.tidyverse.org/)                     | 1.1.4   | R        | Data Manipulation               |
| [Singularity](https://sylabs.io/singularity/)             | 3.11.1  | Other    | Containerization                |
| [Docker](https://www.docker.com/)                         | 25.0.4  | Other    | Containerization                |
| [Conda](https://docs.conda.io/en/latest/)                 | 24.1.2  | Other    | Package/Environment Management  |
| [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper)| 2.1.12  | Python   | Annotation Tools                |

<div align="justify">
