#!/bin/bash
# Annotate H37Rv with itself and check the correction algorithm gets every gene
# right: each CDS of AL123456.3 at exactly its coordinates, strand and name.
#
# Runs inside the AMAP container, which holds the reference, the GenBank record
# and the Bakta database, so nothing large has to live in the repository:
#
#   docker run --rm --platform linux/amd64 -v "$PWD:/amap" -w /amap \
#       paururo/amap:latest bash tests/selftest/h37rv_selftest.sh [workdir]
#
# Bakta and miniprot run once, with the arguments BAKTA_ANNOTATOR and MINIPROT
# use, and are reused afterwards. After that a rerun only exercises the
# algorithm (the last step of MINIPROT, TYPE_CASES and FIX_CDS) and takes
# seconds, which is what makes it usable while working on bin/classify or
# bin/check_cds.py.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
WORK=${1:-$REPO/selftest_h37rv}
SAMPLE=L4_H37Rv
REF=/reference
mkdir -p "$WORK/inputs"
cd "$WORK/inputs"

if [ ! -s "$REF/sequence.gb" ]; then
    echo "ERROR: $REF/sequence.gb not found. Run this inside the AMAP container." >&2
    exit 2
fi

# AL123456.3 is the assembly the reference proteome's coordinates describe.
# The H37Rv the test profile annotates, CP110619, is another one, 6.4 kb longer.
if [ ! -s "$SAMPLE.fasta" ]; then
    python3 - "$REF/sequence.gb" "$SAMPLE" > "$SAMPLE.fasta" <<'EOF'
import sys
from Bio import SeqIO
record = next(SeqIO.parse(sys.argv[1], "genbank"))
seq = str(record.seq)
print(f">{sys.argv[2]}")
print("\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)))
EOF
fi

if [ ! -s "$SAMPLE/$SAMPLE.tsv" ]; then
    echo "==> Bakta (once)" >&2
    export TMPDIR=$(mktemp -d) MPLCONFIGDIR=$(mktemp -d)
    bakta --db /data/bakta_db/db-light --keep-contig-headers --skip-trna --threads "${THREADS:-4}" \
        --genus "Mycobacterium" --species "Mycobacterium tuberculosis" \
        --force --output "$SAMPLE" --prefix "$SAMPLE" --locus-tag "$SAMPLE" \
        --proteins "$REF/sequence.gb" --tmp-dir "$TMPDIR" "$SAMPLE.fasta" > bakta.log 2>&1
fi

if [ ! -s "$SAMPLE-ref.paf" ]; then
    echo "==> miniprot (once)" >&2
    miniprot --outs=0.95 -Iut4 -S "$SAMPLE.fasta" "$REF/aa_H37Rv.translated.fas" > "$SAMPLE.paf" 2>/dev/null
    miniprot --outs=0.95 -Iut4 -S "$SAMPLE.fasta" "$REF/bad.bakta.translated.fas" > "$SAMPLE-badbakta.paf" 2>/dev/null
    miniprot --outs=0.95 -Iut4 -S "$REF/new_refH37Rv.fasta" "$SAMPLE/$SAMPLE.faa" > "$SAMPLE-ref.paf" 2>/dev/null
fi

echo "==> correction algorithm" >&2
IN=$WORK/inputs
rm -rf "$WORK/run"
mkdir -p "$WORK/run"
cd "$WORK/run"
# MINIPROT's pass over what contig ends cut off runs AMAP code, so it is never reused.
python3 "$REPO/bin/contig_end_pieces.py" "$IN/$SAMPLE.paf" "$REF/aa_H37Rv.translated.fas" \
    > "$SAMPLE-pieces.faa" 2> pieces.log
miniprot --outs=0.95 -L 10 -n 2 -Iut4 -S "$IN/$SAMPLE.fasta" "$SAMPLE-pieces.faa" \
    > "$SAMPLE-pieces.paf" 2>> pieces.log
python3 "$REPO/bin/classify_cases.py" -b "$IN/$SAMPLE/$SAMPLE.tsv" -p "$IN/$SAMPLE.paf" \
    -bb "$IN/$SAMPLE-badbakta.paf" -c "$IN/$SAMPLE-ref.paf" \
    --genome "$IN/$SAMPLE.fasta" --ref-proteome "$REF/aa_H37Rv.translated.fas" \
    --pieces "$SAMPLE-pieces.paf" > classify.log 2>&1
python3 "$REPO/bin/check_cds.py" "$IN/$SAMPLE.fasta" "${SAMPLE}_corrected.tsv" \
    --ref-proteome "$REF/aa_H37Rv.translated.fas" > check.log 2>&1
python3 "$REPO/bin/deduplicate_rvs.py" "$SAMPLE.AMAP.tsv" --faa "$SAMPLE.AMAP.faa" \
    --fna "$SAMPLE.AMAP.fna" --bakta "$IN/$SAMPLE/$SAMPLE.tsv" \
    --ref-proteome "$REF/aa_H37Rv.translated.fas" --genome "$IN/$SAMPLE.fasta" > dedup.log 2>&1

echo "==> score against the GenBank record" >&2
python3 "$REPO/tests/selftest/score_selftest.py" "$REF/sequence.gb" \
    "$REF/aa_H37Rv.translated.fas" "$SAMPLE.AMAP.tsv"
