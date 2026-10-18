# Changelog

All notable changes to AMAP are documented in this file.

## [Unreleased]

Annotating H37Rv against its own reference now gives back every CDS of the
AL123456.3 GenBank record at exactly its coordinates, strand and name: 4031 of
4031, against 4020 before. The changes apply to every genome, so results are not
directly comparable with 1.0.1.

### 🔬 Correctness
- **Paralog copies no longer take another gene's place**: a protein with a paralog was placed at both loci, and which copy survived came down to insertion order. On H37Rv itself that lost Rv1150, Rv1945 and Rv2828c; on M. bovis it reported genes the strain lacks (the mce3 operon of RD7, the phiRv2 prophage of RD11, plcBC of RD5) as copies sitting on their paralogs. Each placement is now checked against where its protein belongs, predicted from the sample's own offset to the reference so that indels do not throw it off
- **Strand of a resolved conflict**: taken from the protein's first placement anywhere in the genome instead of the one being resolved, so esxM (Rv1792), Rv3466 and two IS6110 copies came out at their own coordinates on the wrong strand; a Bakta call chosen over an alignment also took the alignment's strand
- **Reference pseudogenes keep their extent**: PE21 (Rv2099c) and PPE48 (Rv3022c), which the reference annotates without a stop, were extended to the next one, and oxyR' (Rv2427A) and Rv3216 lost the partial final codon the reference gives them
- **The 1.0.1 fix for reference pseudogenes never fired**: it looked for reference coordinates in a column that has not carried them since March 2026, and its tests wrote them in by hand. `check_cds.py` now reads the reference proteome instead
- **Alignments that stop short are completed**: miniprot can miss a protein's N-terminus (Rv3611 on CP110619) or leave a repetitive C-terminus unaligned further from its stop than the 300 bp searched; the reference protein now says how much is missing
- **Ties between copies of a gene** go to the one Bakta also called, then to the one where the protein belongs by synteny, not the one that comes first in the genome (E. coli K-12 b2635)
- **Draft assemblies**: genes were compared on coordinates alone, and every contig starts at 1, so cut into 40 contigs CDC1551 had 87% of its proteins translated from the wrong contig. Every gene is now matched, translated and written on its own contig, and Bakta runs with `--keep-contig-headers` so that its contig names match the FASTA's
- **Case 4 is grouped strand by strand**: about 8% of its groups (101 of 1255 on H37Rv) put an alignment on one strand together with a Bakta call on the other
- **A complete alignment is not displaced by a longer call in its frame**: such a call ends at the same stop and only starts further upstream (plcD on CDC1551). Nor does a call that spans several reference genes win against them (Rv3022c and Rv3021c, one reading frame in CDC1551)
- **Reference genes are protected again**: the check looked for coordinates in identifiers that no longer carry them, so any call of similar length could displace a reference gene (Rv1614)
- **Missing genes are recovered by what covers them**: genes stored under Bakta locus tags looked missing by name, so a placement a codon off its gene went back in as a second copy. With the case 4 changes, ORFs annotated twice fall from 10 to 3 on CDC1551 and from 16 to 6 on M. bovis, where the reference has none
- **Real copies of a gene are kept**: `deduplicate_rvs.py` deleted every copy of a gene but one, 34 genes on E. coli O157:H7 Sakai, which were then missing from the annotation and from the duplication analysis. Further copies get identifiers of their own
- **Genes that lost their stop no longer swallow the next gene**: `check_cds.py` extended a gene without a stop to the next one in its frame, which after an indel is the next gene's, so that gene was annotated a second time in the same frame (M. bovis Rv0590 over Rv0590A, Rv3373 over Rv3374, Rv2293c over Rv2292c). Such a gene now keeps its extent as a pseudogene, unless the extension only restores residues its reference protein has. ORFs annotated twice fall from 3 to 2 on CDC1551 and from 6 to 3 on M. bovis, and none of those left comes from an extension
- **Small proteins**: miniprot placed none of the 33 proteins of 8 to 29 residues in E. coli K-12, and put the ldr toxin b4665 only on its near-identical neighbour. Every protein without an exact placement is now looked for in the six-frame translation of the genome, which takes K-12 against itself from 4252 to all 4286 of its genes, its alternative-start isoforms counted with their gene, and Sakai from 198 to 220 of 220
- **Frameshifts and lost stops**: a gene was extended to the next stop in its frame however far that was, and came out as a complete gene with a C-terminus its reference lacks; 18 in CDC1551 and 20 in M. bovis. It now keeps its extent and is flagged as a pseudogene when that stop lies more than ten residues past the end of the reference protein, as PGAP and Bakta treat a frameshift, which also covers a frame that runs into the next gene
- **Genes a contig end cuts**: in a draft assembly they came out as pseudogenes, or as complete genes, without a word. A gene whose frame runs off the contig, or that opens at its edge without a start codon, is now marked partial at that end. In CDC1551 cut into 40 contigs, the larger piece of each of the 25 genes a cut runs through is marked so
- **The other pieces of a cut gene**: miniprot reports only the placements within 5% of a protein's best, so only the larger piece of a gene a contig end cuts was annotated. What that placement lacks is now aligned again as a protein of its own, and its placement at the end of another contig is annotated under the gene's name and marked partial: 47 of the 50 pieces in the CDC1551 draft instead of 25, the other three being five or six codons long. A gene that opens at a contig's edge with a codon that can also start one, as one in six codons inside H37Rv genes can, is marked partial too when its protein begins inside its reference protein
- **Disrupted copies of a gene**: deduplication kept further copies only when intact, so a whole copy that is a pseudogene was deleted with the pieces of split genes; it now stays when nearly as long as the reference protein
- **GFF3 references**: each CDS line was translated on its own, so a CDS of several lines came out in fragments, and proteins under 11 residues were dropped to hide them. The lines are joined per CDS and translated with the bacterial code, reading `transl_except`, and E. coli K-12's GFF3 gives the same 4318 proteins as its GenBank record, 24 of which were missing

### 🐛 Pipeline
- **Pseudogenes in the GFF3 and GenBank**: the GFF3 never said which genes were pseudogenes, and the GenBank only for genes with an internal stop, though the master table knew. Both now carry AMAP's flag with its reason (`pseudo=true` and `Note`, `/pseudo` and `/note`), and mark the ends a contig cut off (`partial=true` with `start_range` or `end_range`, `<` or `>`); the master table and `*-AMAP.tsv` have them as their last columns
- **The default H37Rv reference under Docker**: its files live inside the container and were handed to processes as paths, so Nextflow mounted `/reference` from the host, which Docker Desktop refuses and Docker on Linux mounts empty over the container's own. `BUNDLED_REFERENCE` copies them out first
- **`-profile test`**: `VALIDATE_INPUTS` still required the sample sheet the test profile does without. 1.0.1 had only been checked with `-preview`, which runs no task
- **`-resume` kept results made by older code**: Nextflow hashes a task's inputs and script, not the files it runs by path, so after an edit to the classifier or any script under `bin/`, `-resume` reused the old results. Every process that runs AMAP code now lists those files in `ext code`, which puts them in its cache key; edit `bin/classify/cases.py` and TYPE_CASES and FIX_CDS run again while PATHOTYPR and Bakta stay cached
- **PATHOTYPR no longer depends on the classifier**: its script said `pathotypr classify`, and Nextflow adds any `bin/` entry a script names as a bare word to the cache key, so a change in `bin/classify/` re-ran PATHOTYPR, Bakta and everything after them. The subcommand is quoted
- **Nextflow 25.10 and later**: the pipeline stopped before running anything, because the strict parser those versions use by default rejects the helper functions in `nextflow.config`, the statements outside the workflow in `main.nf`, the bare `env(SAMPLE_NAME)` output and `publishDir` paths that name a task input. Nextflow 26.10 also refuses the `return` that ended a `--validate_only` run. The resource caps are now the `resourceLimits` directive, which still takes `--max_cpus`, `--max_memory` and `--max_time` from the command line, a profile or a `-c` config. Nextflow 24.04, where that directive first appeared, is the new minimum
- **`--skip_mycobrowser` in organism-agnostic mode**: a custom reference from another genus turned it on by reassigning the parameter, which current Nextflow ignores with a warning, so the run went ahead with Mycobrowser enabled. It is now worked out in the parameter's default
- **Tasks leave the pipeline untouched**: the first time a copy of AMAP ran, Python cached the bytecode of AMAP's own modules inside it, in `bin/classify/__pycache__`, and with `pathotypr classify` naming that directory in PATHOTYPR's script, the first `-resume` re-ran PATHOTYPR, Bakta and everything after them. Tasks now run with `PYTHONDONTWRITEBYTECODE=1`
- **Samples with no Bakta call outside cases 2 and 3**: when every Bakta call had an alignment at exactly its coordinates, which a small input can do, TYPE_CASES stopped with a KeyError, because the empty column of overlaps was read as a list of column names. Case 1 and case 4 now come out empty and the rest is classified as usual. An empty miniprot table stopped it the same way
- **Samples Bakta calls no CDS on**: with no contig names from Bakta to compare, the contig check took them for different names and stopped the run, and a Bakta table with nothing but its header did not read. The genes miniprot places now come out as missing from Bakta, like any gene it did not call
- **`--version`**: `--help` and the parameter docs have listed it since 1.0.0, but nothing read it, so `nextflow run PathoGenOmics-Lab/AMAP --version` started a run. It now prints the pipeline and Nextflow versions and exits before any task
- **Drafts named after one of their contigs**: pathotypr writes a summary row for each contig that carries lineage markers, and PATHOTYPR took the sample's name from the first, so CDC1551 cut into 40 contigs, sample `CDC1551d`, came out as `L4_contig_14`, and so did every file after it. Two drafts whose first row named the same contig became one sample, the results of one overwriting the other's. A genome is now named `<lineage>_<sampleId>` from the sample sheet, the lineage of a draft being the one the markers of all its contigs add up to, and the contigs of a draft keep their names instead of the first one taking the sample's. A closed genome whose header is its sample ID comes out as before
- **Headers with a description**: the name was the whole header line, so the tutorial's CDC1551 (`NC_002755.2 Mycobacterium tuberculosis CDC1551, complete sequence`) was named with spaces in it, Bakta failed on it three times and the run finished without it, and without failing. Such a genome now takes its sample ID like any other, and so does the test profile's: `L4_H37Rv_test`, where it was `L4_CP110619_H37Rv`
- **An input could be emptied**: PATHOTYPR wrote the renamed genome next to the link Nextflow stages its input as, so when the new name was the file's own, sample `G26689` in `L4_G26689.fasta` for instance, the input was overwritten with nothing and the run went on without it. The input is now staged in a directory of its own
- **Sample IDs with underscores**: the steps after the correction told a sample's lineage from its sample ID by splitting its name at underscores, though typing only puts the lineage in front and sample IDs may hold underscores of their own. `L4_H37Rv_test` was counted as genome `H37Rv`, the tutorial's `Beijing_1` and a `Beijing_2` as one genome `Beijing` in the category counts, the report and the duplication network, an untyped `H37Rv_test` as lineage `H37Rv`, and each gene of an untyped sample with no underscore as a genome of its own in the network. PATHOTYPR now says which lineage it named a genome after and every step is handed the sample ID and the lineage. The case summary gains a `Sample` column and names an untyped sample as it is, not `Untyped_<sampleId>`, and the report no longer reads the first gene not in the reference as a header, nor stops when every gene has one
- **Genomes kept outside the run's folder**: `VALIDATE_INPUTS` staged only the sample sheet and opened the paths written in it from inside its container, which sees only the folder that the work directory and the sheet share. Genomes anywhere else, `/data/genomes` while running from `~/project` for instance, were "FASTA file not found" and the run stopped before anything ran, and a relative path, which the sheet may use, was read from the task directory with or without a container. The genomes are now inputs of the process, each in a directory of its own so that two with the same file name do not collide, and each is read where Nextflow staged it
- **Runs with functional annotation stopped at validation**: `VALIDATE_INPUTS` looked for `--eggnog_db` in its own task directory, where the default `eggnog_db`, relative to the launch directory, never is, and took a database not built yet for an error, though `EGGNOG_DB` builds it. It now looks where `EGGNOG_DB` builds it, staged when it exists, and stops the run only when `eggnog_db.ok` is there without a complete database, since Nextflow then never runs `EGGNOG_DB` again
- **A Bakta database on the host**: `--bakta_db` reached the tasks as a path the container does not mount. With `-profile docker`, validation called the database missing and stopped the run, and where nothing validated it, as under `-profile test`, every sample failed three times and was ignored, and the run ended with nothing annotated and exit 0. A database on the host is now staged into `BAKTA_ANNOTATOR`, `PREPARE_REFERENCE` and `VALIDATE_INPUTS` like any other input; the default inside the container is used as before
- **A warning on every run**: validation expected a `db.sqlite3` that no Bakta database has. It now checks for `version.json` and `bakta.db`, which Bakta cannot run without, so their absence stops the run, and names the database's type and version
- **`--strict_validation`**: the documentation said validation warnings could be made to stop the run, but nothing passed the script's `--strict` on. The new parameter does
- **eggNOG's database under Docker**: `EGGNOG_DB` built it in a folder it worked out for itself from `--eggnog_db`, and having no inputs, its container saw that folder only when Nextflow happened to mount it with the work directory and the pipeline's `bin/`, through a parent folder they share. Otherwise, as when running from `/data/project` with the pipeline in `~/.nextflow`, the database was downloaded and built inside the container and lost with it, and only `eggnog_db.ok` reached the host: later runs skipped `EGGNOG_DB`, and `ANNOT_EGGNOGMAP_BATCH` found no `bacteria.dmnd`. A database already there without `eggnog_db.ok` was built again. `main.nf` now makes the folder on the host and hands it to `EGGNOG_DB` as an input. A folder left with only `eggnog_db.ok` needs that file deleted, as validation says
- **eggNOG's database in a folder of any name**: `ANNOT_EGGNOGMAP_BATCH` asked eggNOG-mapper for `eggnog_db/bacteria.dmnd` and `--data_dir eggnog_db`, but Nextflow stages a folder under its own name, so with a `--eggnog_db` called anything else, as the HPC guide's `/shared/databases/eggnog` or the tutorial's `~/databases/eggnog`, `EGGNOG_DB` built the database and the annotation stopped at "DIAMOND database eggnog_db/bacteria.dmnd not present", with or without a container. The folder is now staged under that one name, whatever it is called. The task's tag and first log line also gave a single sample's FAA size in bytes as the number of samples (`batch_1519647_samples` under `-profile test`)
- **Reference files on the host**: `BAKTA_ANNOTATOR` gave Bakta `--ref_annotation`, or `--bakta_proteins`, as a path written into its script, which a container does not see, nor a task directory when it is relative, so the organism-agnostic mode failed under `-profile docker`: with K-12 as the reference for Sakai, Bakta called the "user proteins file" not valid three times, the sample was ignored and the run ended with exit 0. Likewise `--pathotyper_snps` and `--pathotyper_ref` on the host left every genome untyped, as `NA_<sampleId>`, `--mycobrowser_ref` left the Mycobrowser categories out, and one `--miniprot_ref_*` on the host while the others kept their defaults made `BUNDLED_REFERENCE` fail, so nothing was corrected. A file on the host is now staged into the tasks that read it, each in a directory of its own, and a default is read inside the container as before
- **Reference files that are nowhere**: a typo in one of those paths, or a default on a run without the container, went unnoticed until the samples were ignored or left untyped. Validation now checks each reference file the run reads and that is not on the host, inside the container where the tasks read it, and stops the run when it is not there. A missing Mycobrowser table only warns, since `GET_OUTPUTS` goes on without it
- **The pangenome without lineage typing**: PANGENE aligned each sample's proteins onto the genomes PATHOTYPR renamed, so with `--skip_typing`, which every custom reference from another genus needs, it had none and wrote empty files for every sample. It now aligns them onto the sample sheet's genomes, each copied to `<sampleId>.fasta` first, so that two files of the same name (`assembly.fasta` from two assemblies) stay apart and a `.fna` is not passed over
- **Genes not in the reference, sample by sample**: the report counted only the samples that had some, so a sample whose genes all matched the reference dropped out of its lineage's box instead of counting 0, and when no sample had any (H37Rv against itself, as in `-profile test`) the section was an empty plot. Every sample is now counted, the section says so when none has any, and each point is drawn at its count, on an axis from 0, where it was jittered up and down by as much as 0.4
- **Genes across the origin of a circular sequence**: Bakta reads a sequence as a complete replicon when its FASTA description says `complete` (NCBI's "complete genome" and "complete sequence") or `circular=true` (Unicycler), joins its ends and calls a gene across its origin from its start to its stop, start after end. TYPE_CASES stopped on such a gene: E. coli O157:H7 Sakai with NCBI's headers, whose stcE runs across the origin of plasmid pO157 (92563..2502), failed three times and was left out, and the run finished without failing. Past it, `check_cds.py` dropped such a gene as invalid coordinates, deduplication and the QC table gave it a negative length, the GFF3 wrote it backwards, the GenBank writer stopped on it and the Panaroo GFF3 spanned it over the whole sequence. It now goes through whole and is written as Bakta writes it: start after end in the tables, `join()` across the origin in a circular GenBank record, and in the GFF3 its end plus the sequence's length on a region marked `Is_circular=true`. A gene near the origin of such a sequence is read across it too, so Rop on plasmid pOSAK1, of which Bakta called the first 141 bp, comes out whole. Annotated against K-12, Sakai with NCBI's headers now has 5319 genes in common with the same genome under bare names, which Bakta reads as linear, and differs in stcE and Rop, whole where the linear reading has 804 of stcE's 886 residues and no Rop, and in two phage genes Bakta starts elsewhere on a complete genome

### 🧹 Code
- **`reconcile_coordinates` is gone**: it read reference coordinates from protein identifiers that have not carried them since March 2026, so it changed nothing. The paralog check does its job with the sample's own offsets

### 🧪 Testing
- **`tests/selftest/h37rv_selftest.sh`** annotates H37Rv with itself inside the container and fails unless every GenBank CDS comes back exactly
- **`tests/test_cache_keys.py`** fails when a process runs AMAP code that is not in its cache key, or when its script names a `bin/` entry as a bare word without declaring it
- **`tests/benchmark/benchmark.sh`** scores the correction algorithm on seven genomes whose right answer is known (H37Rv twice, CDC1551 complete and in 40 contigs, M. bovis, E. coli K-12 and Sakai), downloading them and running Bakta and miniprot in the container the first time, and lists every gene that changed since an earlier run
- **`tests/test_name_by_lineage.py`** runs PATHOTYPR's own script, rendered as Nextflow renders it, with a stand-in for pathotypr that writes the summaries pathotypr wrote for CDC1551 complete and in 40 contigs
- **`tests/test_sample_ids.py`** runs TYPE_CASES, GET_OUTPUTS, GET_COUNT, GET_DUPLICATIONS, REPORT_INPUTS and GET_ALLDUP as Nextflow renders them, chained as the pipeline chains them, for a typed and an untyped run whose sample IDs hold underscores, and then the report's own R code where R is installed, and FINAL_REPORT's own script on a run's tables and on one that leaves every table it can empty
- **`tests/test_validate_inputs.py`** runs `VALIDATE_INPUTS`'s own script, rendered as Nextflow renders it, in a task directory laid out as Nextflow lays one out: genomes in a folder the task cannot see, relative paths, two genomes with the same file name and a sheet of 500 samples; and checks both databases, eggNOG's in each state `EGGNOG_DB` can leave it in
- **`tests/test_eggnog_db.py`** runs `EGGNOG_DB`'s own script, rendered as Nextflow renders it, with its folder staged as Nextflow stages an input and stand-ins for eggNOG-mapper's download and build: the database is built in that folder whatever it is called, one already there is used as it is, and a short `bacteria.dmnd` is built again
- **`tests/test_eggnog_annotation.py`** runs `ANNOT_EGGNOGMAP_BATCH`'s own script the same way, with a stand-in for `emapper.py` that looks for the database as eggNOG-mapper does: the database is found whatever its folder is called, a folder called `split`, where the script writes each sample's annotations, is left as it was, and the tag and log count a single sample as one
- **`tests/test_reference_files.py`** runs the scripts of `BAKTA_ANNOTATOR`, `PATHOTYPR`, `GET_OUTPUTS` and `BUNDLED_REFERENCE`, rendered as Nextflow renders them, their own Groovy included, with stand-ins for Bakta and pathotypr and each reference file on the host, under a relative path or inside the container; and fails when a process reads a reference file's param other than as what it falls back on when nothing is staged
- **`tests/test_pangene_genomes.py`** runs PANGENE_GENOME and PANGENE's scripts as Nextflow renders them, with stand-ins for miniprot and pangene, and fails unless each sample is aligned onto every other sample's genome and never onto its own, typed or not
- **465 tests**, up from 92

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
- **The version in the report is the version that produced it**: it was written by hand in three places in the Rmd, so any release that forgot to update them published a report naming the wrong version. It now comes from the manifest in `nextflow.config`

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
