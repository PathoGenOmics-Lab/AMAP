# AMAP: Assisted MTBC annotator pipeline  
AMAP is a bioinformatics pipeline for processing and annotating _Mycobacterium tuberculosis_ complex (MTBC) data.

## Installation  
You can easily install the AMAP pipeline using [snakemake](https://snakemake.readthedocs.io/en/stable/), a workflow management system that aims to reduce the complexity of creating workflows by providing a fast and comfortable execution environment.

### Snakemake installation  
For further directions go to [snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html):

1. First you will need to have mamba installed:  
``
conda install -n base -c conda-forge mamba
``
2. Activate conda:  
``
conda activate base
``
3. Install snakemake:  
``
mamba create -c conda-forge -c bioconda -n snakemake snakemake
`` 
### Execute AMAP in Snakemake  
1. Activate snakemake:  
   ```
   conda activate snakemake
   ```
2. Clone AMAP into your system:  
   ```
   git clone https://github.com/PathoGenOmics-Lab/AMAP.git
   ```
3. Go to workflow directory:  
   ```
   cd AMAP/workflow
   ```  
4. Modify config file with your paths:
   ```
   fasta_dir: myawesome/path/with/fastas
   output_dir: whereI/store/the/results
   ```  
5. To install and execute the pipeline, use the following command:  
   ```
   snakemake --use-conda -c8
   ```

## DAG (Directed Acyclic Graph)

<p align="center">
  <img src="https://github.com/PathoGenOmics-Lab/AMAP/blob/main/images/dag.png" alt="DAGPLOT" width="900">
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
<h2 id="contributors" align="center">
✨ Contributors
</h2>

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

<h2 id="contributors" align="center">
💻 Software Used
</h2>
<div align="justify">

1. [Snakemake](https://github.com/snakemake/snakemake):  
    >Mölder F, Jablonski KP, Letcher B et al. Sustainable data analysis with Snakemake [version 1; peer review: 1 approved, 1 approved with reservations]. F1000Research 2021, 10:33 (https://doi.org/10.12688/f1000research.29032.1)
      
2. [Bakta](https://github.com/oschwengers/bakta):  
    >Schwengers O., Jelonek L., Dieckmann M. A., Beyvers S., Blom J., Goesmann A. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. Microbial Genomics, 7(11). https://doi.org/10.1099/mgen.0.000685
      
