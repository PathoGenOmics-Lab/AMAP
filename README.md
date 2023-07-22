# AMAP: Assisted MTBC annotator pipeline

AMAP is a bioinformatics pipeline for processing and annotating _Mycobacterium tuberculosis_ complex (MTBC) data.

## Installation

You can easily install the AMAP pipeline using [snakemake](https://snakemake.readthedocs.io/en/stable/), a workflow management system that aims to reduce the complexity of creating workflows by providing a fast and comfortable execution environment. 

To install and execute the pipeline, use the following command:
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

---
<h3 id="contributors" align="center">
✨ Contributors
</h3>

AMAP is developed with ❤️ by:
<p align="center">
  <p align="center">
    <a href="https://github.com/paururo">
      <img src="https://avatars.githubusercontent.com/u/50167687?v=4&s=100" width="100px;" alt=""/>
    </a> 
  </p>
  <p align="center">
    <p align="center">
      <a href="https://github.com/paururo">
      Paula Ruiz-Rodriguez
      </a>
    </p>
    <p align="center">
      <small>Institute for Integrative Systems Biology</small></br>
    </p>         
  </p>
</p>

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification ([emoji key](https://allcontributors.org/docs/en/emoji-key)). Contributions of any kind welcome!

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<div align="center">
<table>
  <tr>
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
    </td> 
  </tr>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
---
