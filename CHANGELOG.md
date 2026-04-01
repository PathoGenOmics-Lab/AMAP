# Changelog

All notable changes to AMAP are documented in this file.

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
