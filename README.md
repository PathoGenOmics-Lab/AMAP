<h0 id="badges" align="center">

![AMAP](https://img.shields.io/badge/AMAP-v0.1.0-green?style=plastic)
![Nextflow](https://img.shields.io/badge/Nextflow-v23.10.1-%2335ccae?style=plastic)
[![License: GPL v3](https://img.shields.io/badge/License%3A-GPL%20v3-orange?style=plastic)](https://github.com/PathoGenOmics-Lab/AMAP/blob/main/LICENSE)
![Docker Image](https://img.shields.io/badge/Docker%20Image-AMAP-blue?style=plastic&logo=docker)
![Zenodo](https://img.shields.io/badge/doi%3A-10.5281%2Fzenodo.8024447-%23c9ac02?style=plastic)
![PGO](https://img.shields.io/badge/PathoGenOmics-Lab-yellow?style=plastic&labelColor=red)
[![Twitter Follow](https://img.shields.io/twitter/follow/paururo)](https://twitter.com/intent/follow?screen_name=paururo)
</h0>

<h0 id="contributors" align="justify">

# AMAP: Assisted MTBC annotator pipeline  

__Paula Ruiz-Rodriguez<sup>1</sup>__ 
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. I<sup>2</sup>SysBio, University of Valencia-CSIC, FISABIO Joint Research Unit Infection and Public Health, Valencia, Spain </sub>  

AMAP is a bioinformatics pipeline in nextflow for processing, annotating and correct annotations of _Mycobacterium tuberculosis_ complex (MTBC) genomes.
## Installation  

## Run AMAP in Nextflow  

## AMAP flowchart
```mermaid
%%{init: { 'theme': 'forest' } }%%
flowchart TB
    subgraph "INPUT"
    v0["Genomes"]
    end
    v3([MTB TYPER])
    v4([BAKTA ANNOTATOR])
    v6([MINIPROT])
    v8([TYPE CASES])
    v14([FIX ANNOTATION])
    subgraph "RESULTS"
    v12["OG Bakta 
    annotation"]
    v10["Report cases"]
    v13["Corrected 
    annotation"]
    v15["Lineage Typing"]
    v16["Final Stats"]
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
    v3 --> v15
    v12 --> v16
    v10 --> v16
    v13 --> v16
    v15 --> v16
    
```
---
<h2 id="contributors" align="center">
✨ Contributors
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

| Software          | Version | Type   | Usage                 |
|-------------------|---------|--------|-----------------------|
| [Nextflow](https://github.com/nextflow-io/nextflow)     | 23.10.1 | Other  | Workflow Management   |
| [Bakta](https://github.com/oschwengers/bakta)           | 1.9.3   | Other  | Annotation Tools  |
| [miniprot](https://github.com/lh3/miniprot)             | 0.13    | Other  | Mapping Tools  |
| [Python](https://www.python.org/)                        | 3.9     | Python | Programming Language  |
| [R](https://www.R-project.org/)                          | 4.3.3   | R      | Statistical Computing |
| [pandas](https://github.com/pandas-dev/pandas)           | 2.2.1   | Python | Data Analysis         |
| [intervaltree](https://github.com/chaimleib/intervaltree) | 3.1.0 | Other  | Data Structures       |
| [plotly.py](https://github.com/plotly/plotly.py)         | 5.20.0  | Python | Data Visualization    |
| [ggplot2](https://ggplot2.tidyverse.org/index.html)      | 3.5.0   | R      | Data Visualization    |
| [dplyr](https://dplyr.tidyverse.org/)                    | 1.1.4   | R      | Data Manipulation     |
| [Singularity](https://sylabs.io/singularity/)           | 3.11.1   | Other  | Containerization      |
| [Docker](https://www.docker.com/)                        | 25.0.4   | Other  | Containerization      |


<div align="justify">