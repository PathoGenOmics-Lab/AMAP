---
title: Citation
---

# Citing AMAP

> Ruiz-Rodriguez, P., & Coscolla, M. (2026). AMAP (Version 1.0.0)
> [Computer software]. <https://doi.org/10.5281/zenodo.8024447>

```bibtex
@software{Ruiz-Rodriguez_AMAP_2026,
  author = {Ruiz-Rodriguez, Paula and Coscolla, Mireia},
  doi = {10.5281/zenodo.8024447},
  license = {GPL-3.0},
  month = apr,
  title = {{AMAP}},
  url = {https://github.com/PathoGenOmics-Lab/AMAP},
  version = {1.0.0},
  year = {2026}}
```

The repository also carries a `CITATION.cff`, which GitHub reads to offer a
formatted citation from the sidebar.

## Authors

**Paula Ruiz-Rodriguez** and **Mireia Coscolla**

I<sup>2</sup>SysBio, University of Valencia-CSIC, FISABIO Joint Research Unit
Infection and Public Health, Valencia, Spain.

## Cite the tools too

AMAP orchestrates other people's work. If a result depends on one of these,
cite it as well.

| Tool | What it does here |
| :--- | :--- |
| [Nextflow](https://www.nextflow.io/) | Runs the pipeline |
| [Bakta](https://github.com/oschwengers/bakta) | Baseline annotation |
| [Miniprot](https://github.com/lh3/miniprot) | Reference proteome alignment |
| [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) | Functional annotation |
| [Diamond](https://github.com/bbuchfink/diamond) | Duplication searches |
| [pangene](https://github.com/lh3/pangene) | Pangenome graph |
| [pathotypr](https://github.com/PathoGenOmics-Lab/pathotypr) | MTBC lineage typing |

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
