# Citation

If AMAP is part of your work, cite it in whichever form your reference manager
takes. All three say the same thing; the DOI is the concept DOI, which always
resolves to the newest version.

**APA**

```text
Ruiz-Rodriguez, P., & Coscolla, M. (2026). AMAP (Version 1.0.1) [Computer software]. https://doi.org/10.5281/zenodo.8024446
```

**BibTeX**

```bibtex
@software{Ruiz-Rodriguez_AMAP_2026,
  author = {Ruiz-Rodriguez, Paula and Coscolla, Mireia},
  doi = {10.5281/zenodo.8024446},
  license = {GPL-3.0},
  month = jul,
  title = {{AMAP}},
  url = {https://github.com/PathoGenOmics-Lab/AMAP},
  version = {1.0.1},
  year = {2026}}
```

**CITATION.cff**

```yaml
cff-version: 1.2.0
title: AMAP
message: 'If you use this software, please cite it as below.'
type: software
authors:
  - family-names: Ruiz-Rodriguez
    given-names: Paula
    orcid: 'https://orcid.org/0000-0003-0727-5974'
    affiliation: >-
      Institute for Integrative Systems Biology (I2SysBio,
      University of Valencia - CSIC), FISABIO Joint Research
      Unit “Infection and Public Health”, Paterna, Spain
    email: paula.ruiz.rodriguez@csic.es
  - family-names: Coscolla
    given-names: Mireia
    orcid: 'https://orcid.org/0000-0003-0752-0538'
    affiliation: >-
      Institute for Integrative Systems Biology (I2SysBio,
      University of Valencia - CSIC), FISABIO Joint Research
      Unit “Infection and Public Health”, Paterna, Spain
    email: mireia.coscolla@csic.es
url: 'https://github.com/PathoGenOmics-Lab/AMAP'
license: GPL-3.0
version: 1.0.1
date-released: '2026-07-25'
doi: 10.5281/zenodo.8024446
```

The CITATION.cff tab is the repository's own file, which GitHub reads to offer
a formatted citation from the sidebar.

## Authors

**Paula Ruiz-Rodriguez** and **Mireia Coscolla**

Institute for Integrative Systems Biology (I<sup>2</sup>SysBio, University of
Valencia-CSIC), FISABIO Joint Research Unit "Infection and Public Health",
Paterna, Spain.

## Cite the tools too

AMAP orchestrates other people's work. If a result depends on one of these,
cite it as well.

| Tool | What it does here | Paper |
| :--- | :--- | :--- |
| [Nextflow](https://www.nextflow.io/) | Runs the pipeline | Di Tommaso et al., *Nature Biotechnology* 35, 316-319 (2017). [10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820) |
| [Bakta](https://github.com/oschwengers/bakta) | Baseline annotation | Schwengers et al., *Microbial Genomics* 7, 000685 (2021). [10.1099/mgen.0.000685](https://doi.org/10.1099/mgen.0.000685) |
| [miniprot](https://github.com/lh3/miniprot) | Reference proteome alignment | Li, *Bioinformatics* 39, btad014 (2023). [10.1093/bioinformatics/btad014](https://doi.org/10.1093/bioinformatics/btad014) |
| [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) | Functional annotation | Cantalapiedra et al., *Molecular Biology and Evolution* 38, 5825-5829 (2021). [10.1093/molbev/msab293](https://doi.org/10.1093/molbev/msab293) |
| [Diamond](https://github.com/bbuchfink/diamond) | Duplication searches, and eggNOG's | Buchfink et al., *Nature Methods* 18, 366-368 (2021). [10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x) |
| [pangene](https://github.com/lh3/pangene) | Pangenome graph | Li et al., *Bioinformatics* 40, btae456 (2024). [10.1093/bioinformatics/btae456](https://doi.org/10.1093/bioinformatics/btae456) |
| [pathotypr](https://github.com/PathoGenOmics-Lab/pathotypr) | MTBC lineage typing | Its repository, and its lineage markers below |

All of them in one block, for a reference manager:

```bibtex
@article{DiTommaso_Nextflow_2017,
  author = {Di Tommaso, Paolo and Chatzou, Maria and Floden, Evan W. and Barja, Pablo Prieto and Palumbo, Emilio and Notredame, Cedric},
  title = {Nextflow enables reproducible computational workflows},
  journal = {Nature Biotechnology},
  volume = {35},
  number = {4},
  pages = {316--319},
  year = {2017},
  doi = {10.1038/nbt.3820}}

@article{Schwengers_Bakta_2021,
  author = {Schwengers, Oliver and Jelonek, Lukas and Dieckmann, Marius Alfred and Beyvers, Sebastian and Blom, Jochen and Goesmann, Alexander},
  title = {Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification},
  journal = {Microbial Genomics},
  volume = {7},
  number = {11},
  pages = {000685},
  year = {2021},
  doi = {10.1099/mgen.0.000685}}

@article{Li_miniprot_2023,
  author = {Li, Heng},
  title = {Protein-to-genome alignment with miniprot},
  journal = {Bioinformatics},
  volume = {39},
  number = {1},
  pages = {btad014},
  year = {2023},
  doi = {10.1093/bioinformatics/btad014}}

@article{Cantalapiedra_eggNOGmapper_2021,
  author = {Cantalapiedra, Carlos P. and Hern{\'a}ndez-Plaza, Ana and Letunic, Ivica and Bork, Peer and Huerta-Cepas, Jaime},
  title = {{eggNOG-mapper} v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale},
  journal = {Molecular Biology and Evolution},
  volume = {38},
  number = {12},
  pages = {5825--5829},
  year = {2021},
  doi = {10.1093/molbev/msab293}}

@article{Buchfink_DIAMOND_2021,
  author = {Buchfink, Benjamin and Reuter, Klaus and Drost, Hajk-Georg},
  title = {Sensitive protein alignments at tree-of-life scale using {DIAMOND}},
  journal = {Nature Methods},
  volume = {18},
  number = {4},
  pages = {366--368},
  year = {2021},
  doi = {10.1038/s41592-021-01101-x}}

@article{Li_pangene_2024,
  author = {Li, Heng and Marin, Maximillian and Farhat, Maha R.},
  title = {Exploring gene content with pangene graphs},
  journal = {Bioinformatics},
  volume = {40},
  number = {7},
  pages = {btae456},
  year = {2024},
  doi = {10.1093/bioinformatics/btae456}}
```

The exact version of each one used by a given run is recorded in
`results/pipeline_info/software_versions.yml`.

## Reference data

The MTBC reference files bundled in the container come from public
repositories. The pathotypr lineage markers are archived on Zenodo under
[10.5281/zenodo.19210044](https://doi.org/10.5281/zenodo.19210044), and
Mycobrowser functional categories come from
[mycobrowser.epfl.ch](https://mycobrowser.epfl.ch/).

## Licence

AMAP is released under the GPL-3.0 licence. See
[LICENSE](https://github.com/PathoGenOmics-Lab/AMAP/blob/main/LICENSE).
