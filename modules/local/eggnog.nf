process EGGNOG_DB {
    tag 'eggnog_db_setup'

    // storeDir caches the DB outside the work directory.
    // Once built, it's reused across ALL pipeline runs — never rebuilt.
    storeDir "${params.eggnog_db}"

    output:
    path "eggnog_db.ok", emit: db_path

    script:
    """
    set -euo pipefail

    # Resolve the user-provided eggnog_db path
    # Use launchDir (where user ran nextflow) for relative paths, NOT projectDir
    # (which points to the pipeline cache and is wrong for user-supplied paths)
    EGGNOG_DB_DIR="${params.eggnog_db}"
    if [[ ! "\$EGGNOG_DB_DIR" = /* ]]; then
        EGGNOG_DB_DIR="${launchDir}/${params.eggnog_db}"
    fi

    mkdir -p "\$EGGNOG_DB_DIR"

    # ------------------------------------------------------------
    # 0) Quick exit if the bacteria DB already exists
    # ------------------------------------------------------------
    if [ -f "\$EGGNOG_DB_DIR/bacteria.dmnd" ]; then
        # Integrity check: Diamond DB should be >100 MB
        DB_SIZE=\$(stat -f%z "\$EGGNOG_DB_DIR/bacteria.dmnd" 2>/dev/null || stat -c%s "\$EGGNOG_DB_DIR/bacteria.dmnd" 2>/dev/null || echo 0)
        if [ "\$DB_SIZE" -lt 100000000 ]; then
            echo "WARNING: bacteria.dmnd is only \${DB_SIZE} bytes (expected >100 MB) — may be corrupted, rebuilding..." >&2
            rm -f "\$EGGNOG_DB_DIR/bacteria.dmnd"
        else
            echo "==> bacteria.dmnd found (\${DB_SIZE} bytes) – skipping download and build."
            touch eggnog_db.ok
            exit 0
        fi
    fi

    # ------------------------------------------------------------
    # 1) Download metadata only (skip 60-GB Diamond core)
    # ------------------------------------------------------------
    echo "==> Downloading metadata only (skipping Diamond core)…"
    download_eggnog_data.py -D -y --data_dir "\$EGGNOG_DB_DIR"

    # ------------------------------------------------------------
    # 2) Build the bacteria-only Diamond index (~3 GB)
    # ------------------------------------------------------------
    echo "==> Building bacteria.dmnd (taxid 2)…"
    create_dbs.py -m diamond \
                  --dbname bacteria \
                  --taxa Bacteria \
                  -y \
                  --data_dir "\$EGGNOG_DB_DIR"

    # ------------------------------------------------------------
    # 3) Clean up build-only files
    # ------------------------------------------------------------
    echo "==> Removing build-only FASTA and taxid table…"
    rm -f "\$EGGNOG_DB_DIR"/e5.proteomes.faa  || true
    rm -f "\$EGGNOG_DB_DIR"/e5.taxid_info.tsv || true

    touch eggnog_db.ok
    """
}

/*
 * BATCH eggNOG annotation — runs eggNOG-mapper ONCE on all samples combined.
 *
 * WHY: On HPC clusters with shared filesystems (NFS/Lustre), running N
 * parallel eggNOG-mapper jobs that each load the ~3 GB Diamond database
 * causes I/O contention. With 200+ concurrent jobs, the filesystem becomes
 * the bottleneck and annotation is slower than serial.
 *
 * HOW: Concatenate all per-sample FAA files with sample-prefixed headers,
 * run eggNOG-mapper once with all CPUs, then split the output back by sample.
 * The DB is loaded into memory once (--dbmem), and Diamond uses all cores.
 */
process ANNOT_EGGNOGMAP_BATCH {
  tag "batch_${faa_files.size()}_samples"

  input:
  path faa_files    // collected FAA files from all samples
  val  sample_ids   // list of sampleIds in matching order
  path eggnog_db

  output:
  path "split/*.emapper.annotations", emit: split_annotations

  script:
  """
  set -euo pipefail
  mkdir -p tmp scr_tmp split

  echo "▸ [eggNOG] Preparing batch annotation for ${faa_files.size()} samples..." >&2

  # ================================================================
  # 1. Concatenate all FAAs with sample-prefixed headers
  # ================================================================
  #    Delimiter: @@SMPL@@ (safe — no protein ID contains this)
  #    Original:  >L8_G26689_Rv0001;+;0-1524
  #    Tagged:    >SAMPLE@@SMPL@@L8_G26689@@SMPL@@L8_G26689_Rv0001;+;0-1524
  for faa in ${faa_files}; do
      sid=\$(basename "\$faa" .AMAP.faa)
      awk -v sid="\$sid" '/^>/{sub(/^>/, ""); print ">SAMPLE@@SMPL@@" sid "@@SMPL@@" \$0; next} {print}' "\$faa"
  done > all_proteins.faa

  N_PROTS=\$(grep -c '^>' all_proteins.faa || echo 0)
  if [ "\$N_PROTS" -eq 0 ]; then
      echo "WARNING: No proteins to annotate (all samples have 0 CDS)" >&2
      echo "# emapper: no input proteins" > batch.emapper.annotations
  else
      # ================================================================
      # 2. Deduplicate by sequence — run eggNOG only on unique proteins
      # ================================================================
      #    With N closely-related genomes (e.g. 480 MTB strains), most
      #    proteins are identical across samples. Deduplicating reduces
      #    the Diamond search from N_total to N_unique (typically 5-20×
      #    fewer queries).
      #
      #    Strategy:
      #    - Hash each protein sequence → pick one representative header
      #    - Build a mapping: representative → all duplicate headers
      #    - Run eggNOG on representatives only
      #    - Propagate annotations to all duplicates

      echo "==> Deduplicating \$N_PROTS proteins by sequence..." >&2

      # dedup_proteins.py: reads FASTA, groups by sequence, writes
      #   unique.faa (representatives) and header_map.tsv (rep → dups)
      python3 - <<'PYEOF'
import sys
from collections import OrderedDict

seq_to_headers = OrderedDict()   # sequence → [header1, header2, ...]
header_to_seq = {}               # header → sequence (for ordering)

current_header = None
chunks = []

with open("all_proteins.faa") as fh:
    for line in fh:
        if line.startswith(">"):
            if current_header is not None:
                seq = "".join(chunks)
                header_to_seq[current_header] = seq
                seq_to_headers.setdefault(seq, []).append(current_header)
            current_header = line[1:].rstrip()
            chunks = []
        else:
            chunks.append(line.strip())
    # last record
    if current_header is not None:
        seq = "".join(chunks)
        header_to_seq[current_header] = seq
        seq_to_headers.setdefault(seq, []).append(current_header)

# Write unique representatives
with open("unique.faa", "w") as out_faa:
    for seq, headers in seq_to_headers.items():
        rep = headers[0]
        out_faa.write(f">{rep}\n{seq}\n")

# Write mapping: rep<TAB>dup1<TAB>dup2<TAB>...
with open("header_map.tsv", "w") as out_map:
    for seq, headers in seq_to_headers.items():
        rep = headers[0]
        out_map.write(rep + "\t" + "\t".join(headers[1:]) + "\n")

n_unique = len(seq_to_headers)
n_total = sum(len(h) for h in seq_to_headers.values())
ratio = n_total / n_unique if n_unique > 0 else 0
print(f"    Unique: {n_unique} / {n_total} ({ratio:.1f}× redundancy)", file=sys.stderr)
PYEOF

      N_UNIQUE=\$(grep -c '^>' unique.faa || echo 0)
      echo "==> Annotating \$N_UNIQUE unique proteins (from \$N_PROTS total)..." >&2

      # ================================================================
      # 3. Run eggNOG-mapper on unique sequences only
      # ================================================================
      MEM_GB=\$(echo "${task.memory}" | grep -oP '\\d+' | head -1)
      BLOCK_SIZE=\$(( MEM_GB / 8 ))
      [ \$BLOCK_SIZE -lt 2 ] && BLOCK_SIZE=2
      [ \$BLOCK_SIZE -gt 12 ] && BLOCK_SIZE=12

      emapper.py \\
        -i unique.faa \\
        -m diamond \\
        --dmnd_db eggnog_db/bacteria.dmnd \\
        --cpu ${task.cpus} \\
        --usemem \\
        --dbmem \\
        --sensmode fast \\
        --block_size \$BLOCK_SIZE \\
        --index_chunks 1 \\
        --tax_scope 2 \\
        --target_orthologs one2one \\
        -o batch \\
        --data_dir eggnog_db \\
        --temp_dir tmp \\
        --scratch_dir scr_tmp

      # ================================================================
      # 4. Expand annotations back to all duplicate headers
      # ================================================================
      #    For each annotation line, check if the query is a representative
      #    that has duplicates → emit one line per duplicate with its
      #    original header restored.
      echo "==> Expanding annotations to all \$N_PROTS proteins..." >&2

      python3 - <<'PYEOF'
import sys

# Load representative → duplicates mapping
rep_to_dups = {}  # rep_header → [dup1, dup2, ...]
with open("header_map.tsv") as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        rep = parts[0]
        dups = [d for d in parts[1:] if d]  # skip empty fields
        if dups:
            rep_to_dups[rep] = dups

# Read batch annotations, expand duplicates
with open("batch.emapper.annotations") as fin, \
     open("batch_expanded.emapper.annotations", "w") as fout:
    for line in fin:
        # Pass through header/comment lines
        if line.startswith("#"):
            fout.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        query = fields[0]
        # Write the original line (representative)
        fout.write(line)
        # Write duplicates with their own headers
        if query in rep_to_dups:
            for dup in rep_to_dups[query]:
                fields[0] = dup
                fout.write("\t".join(fields) + "\n")
PYEOF

      # Replace batch output with expanded version
      mv batch_expanded.emapper.annotations batch.emapper.annotations
  fi

  # ================================================================
  # 5. Split results back into per-sample annotation files
  # ================================================================
  #    All in one awk pass — O(N) instead of O(N×M) subprocesses
  grep '^#' batch.emapper.annotations > header.txt || true

  # Create per-sample files with headers
  for faa in ${faa_files}; do
      sid=\$(basename "\$faa" .AMAP.faa)
      cp header.txt "split/\${sid}.emapper.annotations"
  done

  # Split annotation lines by sample — single awk pass over the file
  # NOTE: close() after each print >> to avoid hitting ulimit -n (open files)
  grep -v '^#' batch.emapper.annotations 2>/dev/null | awk -F'\\t' '
  {
      query = \$1
      split(query, parts, "@@SMPL@@")
      sid = parts[2]
      orig_query = parts[3]
      \$1 = orig_query
      outfile = "split/" sid ".emapper.annotations"
      print >> outfile
      close(outfile)
  }
  ' OFS='\\t' || true
  """
}



