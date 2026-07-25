# Changelog

All notable changes to AMAP are documented in this file.

## [1.0.1] - 2026-07-25

Bug fixes, most of them found by running the pipeline on real genomes rather
than by reading it. Several changed what the pipeline outputs, so results from
1.0.0 are not directly comparable.

### 🔬 Correctness
- **Genes no longer vanish during overlap resolution**: `_insert_with_dedup` queued entries for eviction while scanning for conflicts, and applied them even when it then abandoned the insert, so genes were deleted with nothing put in their place
- **Reference pseudogenes are recognised for any organism**: the check was gated on the identifier starting with `Rv`, so on every other organism every reference pseudogene was extended to the next stop and lost its flag. The same gene came out differently depending on what it was called
- **The extracted reference proteome carries its stop codon**: `extract_ref_proteins.py` wrote coordinates spanning the stop but a protein without it, so miniprot's alignment ended three bases short of every Bakta CDS and no gene could ever reach case 2. On *E. coli* O157:H7 this moved 3325 genes out of conflict resolution
- **Mycobrowser join key**: built with `rvs.split('_')[0]`, which returned the whole string for a gene with several references (`Rv0094c, Rv3467`) and truncated any name containing an underscore (`PE_PGRS1` became `PE`). 44 more genes get a functional category on a real MTBC sample, none lose one
- **Annotation table and FASTA files agree**: `check_cds.py` wrote the FASTA files from dictionaries filled before its final deduplication, and `deduplicate_rvs.py` matched FASTA headers on the gene id alone, so removing one copy of a multi-copy gene deleted all of them
- **Untranslatable genes are dropped from the table**, not left in it with no sequence anywhere

### 🐛 Crashes
- **`--skip_mycobrowser` and organism-agnostic mode**: the placeholder frame was built with a duplicate `Locus` column, so `pd.merge` refused to run and `GET_OUTPUTS` failed on every run without Mycobrowser
- **Miniprot PAF parsing**: pandas inferred the column count from the first row, and miniprot varies its trailing tags per alignment, so any file whose first row was narrower aborted the parse
- **Biopython 1.85 and newer**: `SeqFeature.strand` was removed; `to_panaroo.py` and `extract_ref_proteins.py` now read it off the location
- **`-profile test`**: pointed at a sample sheet naming paths that are not in this repository. It now annotates the reference genome bundled in the container, so it needs no data and nothing downloaded

### 📉 Silent data loss
- **First Diamond hit**: `read_and_create` left header inference on while stripping `#` comments, so the first alignment became the column names and disappeared from the duplication analysis
- **Multi-part gene names in `rename_fastas.py`**: headers were cut at the third underscore, collapsing every PE/PPE gene of a sample onto one identifier
- **eggNOG annotation rate**: counted with `notna()`, and emapper writes `-` for fields it could not fill, so the QC metric read close to 100% whatever the run had done

### ⚡ Performance
- **`reconcile_coordinates` is linear rather than quadratic**: it scanned every protein and every reference for each gene, calling a regex on each candidate. Indexed by gene name up front, a 4000-gene genome goes from 18.3 s to 0.011 s with identical output over 200 randomised genomes
- **`deduplicate`**: `list.count()` per element and a full scan per duplicated gene, replaced with a Counter and a single grouping pass

### 📊 Report
- **Layout**: the stylesheet was losing to the one rmarkdown inlines afterwards, so the article sat underneath the table of contents and the left logo on top of the title
- **Self-contained without a network at render time**: jQuery came from a CDN and the second logo from GitHub. pandoc inlines those when the render machine has a network and leaves them dangling when it does not, so a report rendered on an isolated node arrived with a broken image
- **Half the size**: 8.2 MB to 4.4 MB, by shipping only the plotly traces actually used and quantising the figures

### 📚 Documentation and infrastructure
- **Documentation site** built with MkDocs Material, published to GitHub Pages, including a tutorial with the numbers from real runs and a rendered example report
- **Docker image workflow**: builds, smoke tests and publishes on a `v*` tag or on demand
- **Tests workflow enabled**: it had been sitting disabled and would have failed, since its dependency list was missing plotly and its syntax check had drifted five files behind
- **`requirements-test.txt`**: one dependency list for CI, the docs and a local checkout. `environment.yml` was missing the Python `plotly`, so a conda environment built from it could not start `classify_cases.py`
- **92 tests**, up from 58

## [1.0.0] — 2026-04-01

### 🔬 Correctness
- **Strand-aware gene dictionary**: `combined_dict` keys changed from `(Start, End)` to `(Start, End, Strand)` — genes on opposite strands at the same coordinates no longer collide silently
- **Consistent 1-based lengths**: all length calculations now use `End - Start + 1` (15 sites fixed)
- **Robust ID cleaning**: regex `_\d+-\d+$` replaces brittle 4-part underscore rule — PE_PGRS genes now handled correctly
- **Exact gene name matching**: prevents false matches (Rv3508 ≠ Rv3508c, Rv1 ≠ Rv10)
- **IntervalTree half-open intervals**: all 9 operations corrected from `[Start:End]` to `[Start:End+1]`
- **Cross-strand gene preservation**: merge operations now include Strand to prevent loss
- **Paralog preservation**: multi-copy PE/PPE genes no longer collapsed during groupby
- **Pseudogene propagation**: Is_Pseudogene column now reaches master table and report

### 🛡️ Robustness
- **Pre-flight validation**: new `VALIDATE_INPUTS` process checks FASTA format, N-content, GC%, genome size, duplicate headers, database paths — runs before any computation
- **`--validate_only` mode**: check all inputs without executing the pipeline
- **eggNOG DB integrity**: `bacteria.dmnd` <100 MB detected as corrupted → automatic rebuild
- **Bakta DB pre-check**: actionable error message with download command
- **OOM retry**: memory scales with `task.attempt` (e.g., 50 GB → 100 GB → 150 GB)
- **0-protein handling**: eggNOG batch handles empty inputs gracefully
- **Completion handler**: prints duration, succeeded/failed/ignored counts

### 📋 Reproducibility
- **`nextflow.manifest`**: version 1.0.0, Nextflow ≥23.04.0, DOI, homepage
- **`SOFTWARE_VERSIONS`**: YAML + TSV with exact versions of 12 tools
- **`DUMP_PARAMS`**: complete parameter dump per run
- **Diamond sorted output**: deterministic results regardless of thread count
- **`environment.yml`**: Conda environment with pinned versions
- **Louvain seed**: fixed to 42 for reproducible community detection

### ⚡ Performance & HPC
- **Dynamic resource scaling**: memory grows with each retry attempt
- **Node-local scratch**: `--scratch true` (default on SLURM/SGE/PBS/LSF)
- **Diamond block size tuning**: auto-calculated from available RAM
- **Bakta local temp**: `--tmp-dir` on node-local disk (avoids NFS metadata storms)
- **eggNOG batch mode**: single run for all samples (10× faster, 50× less filesystem I/O)
- **Thread utilization**: Bakta 12, Diamond 6-24, Miniprot 2×3, Pangene 2×N
- **SLURM/SGE/PBS/LSF profiles**: with time limits, queue sizes, rate limiting
- **`storeDir` for eggNOG DB**: persists across runs, survives `nextflow clean`

### 📐 Scalability
- **`find -print0 | xargs -0`**: replaces `echo|tr|grep` in REPORT_INPUTS (handles 500+ files)
- **`publishDir mode: 'link'`**: hardlinks instead of copies for per-sample outputs
- **`cache 'lenient'`**: tolerates NFS timestamp drift between cluster nodes
- **awk `close()`**: prevents file handle exhaustion on 500+ sample eggNOG splits
- **JVM tuning docs**: NXF_OPTS guidance for large runs

### 🧬 Organism-Agnostic
- **Custom reference**: `--ref_annotation` + `--ref_genome` replaces built-in H37Rv
- **Auto-detect mode**: `--skip_mycobrowser` auto-enabled when genus ≠ Mycobacterium
- **Dynamic reference comparison**: report uses median distribution (not hardcoded H37Rv)
- **Auto-extending color palette**: non-MTBC lineages get colors automatically
- **Generic labels**: "Ratio to Median" instead of "Ratio to H37Rv"
- **Manifest updated**: "Assisted Microbial Annotator Pipeline"

### 📊 Output Quality
- **Clean proteome**: `{sample}-AMAP.clean.faa` without pseudogenes or internal stops
- **Annotation QC**: per-sample metrics (genome size, GC%, coding density, hypothetical %, annotation rate)
- **GFF3 spec compliance**: `##sequence-region` pragmas, attribute escaping, sorted features
- **QC report section**: interactive table + coding density / hypothetical % / genome size plots

### 🧪 Code Quality
- **58 unit tests**: covering IDs, overlaps, reconciliation, validation, clean proteome, QC, merge
- **GitHub Actions CI**: Python 3.9 + 3.11 matrix (currently disabled, ready to activate)
- **Type hints**: core modules with `Coord`, `Coord3`, `CoordDict` aliases
- **Module docstrings**: coordinate conventions, naming schemes documented

### 🖥️ Usability
- **Structured startup banner**: shows all configuration in a formatted box
- **`--version` flag**: prints pipeline + Nextflow version
- **Progress logging**: `▸ [BAKTA] Annotating G26689...` for every major step
- **Actionable error messages**: boxed errors with cause + fix + re-run command
- **Sample count at startup**: confirms TSV was parsed correctly

### 📝 Documentation
- **README**: updated for v1.0.0 with all new parameters, outputs, and organism-agnostic usage
- **CHANGELOG**: this file
- **Output tree**: documents `.clean.faa`, `.qc.tsv`, `pipeline_info/` new files

### 🐛 Bug Fixes (85 total)
- 20 critical (🔴) + 65 improvements (🟡)
- See commit history on `fix` branch for full details

## [0.1.0] — 2025

- Initial release
- MTBC-focused annotation pipeline
- Bakta + Miniprot + eggNOG integration
- HTML report with COG, Mycobrowser, pangenome analysis
