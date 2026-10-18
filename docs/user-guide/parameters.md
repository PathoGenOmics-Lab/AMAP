---
title: Parameters
---

# Parameters

Every AMAP parameter, what it defaults to, and when you would change it. Run
`nextflow run PathoGenOmics-Lab/AMAP --help` for the same list at the terminal.

Pipeline parameters take two dashes (`--tsv`). Nextflow's own options take one
(`-profile`, `-resume`).

## Input and output

| Parameter | Description | Type | Default | Required |
| :--- | :--- | :--- | :--- | :--- |
| `--tsv` | Sample sheet with `sampleId` and `fasta` columns, tab separated. Not needed with `-profile test`. | `file` | `samples.tsv` | **Yes** |
| `--outdir` | Directory for results. | `string` | `results` | No |

## Skip options

| Parameter | Description | Type | Default |
| :--- | :--- | :--- | :--- |
| `--skip_typing` | Skip pathotypr lineage typing. Required for non-MTBC organisms. | `boolean` | `false` |
| `--skip_functional_annotation` | Stop after annotation correction. Skips eggNOG, duplications, pangenome and the report. | `boolean` | `false` |
| `--skip_mycobrowser` | Skip Mycobrowser functional categories. Auto-enabled when `--genus` is not Mycobacterium and a custom reference is set. | `boolean` | `false` |

Their effect on the output tree is covered under
[outputs](outputs.md#effect-of-the-skip-options).

## Annotation

| Parameter | Description | Type | Default |
| :--- | :--- | :--- | :--- |
| `--genus` | Genus passed to Bakta. | `string` | `"Mycobacterium"` |
| `--species` | Species passed to Bakta. | `string` | `"Mycobacterium tuberculosis"` |

Bakta uses these to pick naming conventions and reference sets, so setting them
correctly improves the baseline annotation that AMAP then corrects.

## Resources

| Parameter | Description | Type | Default |
| :--- | :--- | :--- | :--- |
| `--max_cpus` | Ceiling on CPUs for any single process. | `integer` | `24` |
| `--max_memory` | Ceiling on memory for any single process. | `string` | `50.GB` |
| `--max_time` | Ceiling on wall-clock time for any single process. | `string` | `72.h` |
| `--scratch` | Use node-local scratch for I/O-heavy processes. On by default in the HPC profiles. | `boolean` | `false` |

These are ceilings, not requests. Each process asks for what it needs and is
then clamped to these values, so lowering them makes AMAP fit a smaller machine
rather than making it run faster.

Memory scales on retry: a process killed by the out-of-memory killer is
resubmitted with more, up to the ceiling.

## Validation and information

| Parameter | Description | Type | Default |
| :--- | :--- | :--- | :--- |
| `--validate_only` | Run input validation and stop. | `boolean` | `false` |
| `--strict_validation` | Stop on validation warnings too, not only on errors. | `boolean` | `false` |
| `--help` | Print the help text and exit. | `boolean` | `false` |
| `--version` | Print pipeline and Nextflow versions and exit. | `boolean` | `false` |

## Custom reference

Supply both to annotate any prokaryote against your own trusted annotation.

| Parameter | Description | Type | Default |
| :--- | :--- | :--- | :--- |
| `--ref_annotation` | Trusted reference annotation, GFF3 or GenBank. | `file` | `null` |
| `--ref_genome` | Trusted reference genome, FASTA. | `file` | `null` |

They must be given together; setting only one is an error. See
[custom reference](custom-reference.md) for what AMAP does with them.

## Reference and database arguments

!!! note "You almost certainly do not need these"

    Every file below ships inside the Docker and Singularity images, under
    `/reference/`. Running with `-profile docker` or `-profile singularity`
    picks them up automatically.

    They exist for running without a container, for swapping in updated
    reference files, or for adapting the pipeline by hand.

| Parameter | Description | Default inside the container |
| :--- | :--- | :--- |
| `--bakta_db` | Bakta database directory, inside the container or on the host. | `/data/bakta_db/db-light` |
| `--eggnog_db` | Pre-formatted eggNOG database directory. | `eggnog_db` |
| `--pathotyper_snps` | Lineage marker SNPs for pathotypr. | `/reference/snps_lineages.tsv` |
| `--pathotyper_ref` | Reference FASTA for pathotypr. | `/reference/MTB_ancestor_reference.fasta` |
| `--miniprot_ref_aa` | Reference proteome aligned onto each genome. | `/reference/aa_H37Rv.translated.fas` |
| `--miniprot_ref_bad_bakta` | Known problematic annotations, used to filter novel-gene candidates. | `/reference/bad.bakta.translated.fas` |
| `--miniprot_ref_h37rv` | Reference genome FASTA. | `/reference/new_refH37Rv.fasta` |
| `--bakta_proteins` | Reference proteins GenBank file for Bakta. | `/reference/sequence.gb` |
| `--mycobrowser_ref` | Mycobrowser functional categories table. | `/reference/Mycobacterium_tuberculosis_H37Rv_txt_v4.txt` |

A database or a file of your own on the host, the full Bakta database or a
newer Mycobrowser table for instance, is mounted into the tasks that use it, so
it works with the container profiles too. A path that is not on the host is
read inside the container, where the defaults are, and validation stops a run
whose file is in neither place.

The eggNOG database is the one exception: it is not bundled, because of its
size. AMAP downloads and builds it on the first run and caches it at
`--eggnog_db`, so point that at somewhere persistent and shared if you run on a
cluster.

## Nextflow options worth knowing

These belong to Nextflow rather than AMAP, so they take a single dash.

| Option | What it does |
| :--- | :--- |
| `-profile` | Selects execution profiles, comma separated: `docker`, `singularity`, `slurm`, `sge`, `pbs`, `test`. |
| `-resume` | Reuse cached results from the previous run. |
| `-work-dir` | Where intermediate files go. Defaults to `./work`, and it gets large. |
| `-process.queue` | Override the scheduler queue or partition. |
| `-process.clusterOptions` | Extra flags passed to the scheduler, such as `'--account=mygroup'`. |

See [running AMAP](running.md) for how the profiles combine.
