#!/bin/bash
# Measure how well the correction algorithm annotates seven genomes whose right
# answer is known, and compare a run with an earlier one.
#
# Runs inside the AMAP container, which holds Bakta, miniprot and the H37Rv
# reference; the other genomes are downloaded from NCBI the first time:
#
#   docker run --rm --platform linux/amd64 -v "$PWD:/amap" -w /amap \
#       paururo/amap:latest bash tests/benchmark/benchmark.sh [workdir] [run] [compare-with]
#
# workdir       where inputs, truths and runs are kept (default: benchmark_work)
# run           name of this run, kept under workdir/runs (default: latest)
# compare-with  name of an earlier run, to list the genes that changed
#
# The first run downloads four genomes, prepares the E. coli reference and
# runs Bakta eight times, a quarter of an hour with 12 threads.
# Everything up to the algorithm is cached in the workdir, so later runs only
# redo classify_cases.py, check_cds.py and deduplicate_rvs.py and score them,
# a few minutes for all seven genomes. See tests/benchmark/README.md.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
WORK=$(mkdir -p "${1:-$REPO/benchmark_work}" && cd "${1:-$REPO/benchmark_work}" && pwd)
RUN=${2:-latest}
COMPARE=${3:-}
THREADS=${THREADS:-8}
REF=/reference
BAKTA_DB=${BAKTA_DB:-/data/bakta_db/db-light}
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF
MINIMAP2_VERSION=2.28
BENCH="python3 $REPO/tests/benchmark/benchmark.py"

if [ ! -s "$REF/sequence.gb" ]; then
    echo "ERROR: $REF/sequence.gb not found. Run this inside the AMAP container." >&2
    exit 2
fi
mkdir -p "$WORK/inputs" "$WORK/truth" "$WORK/tools"
log() { echo "[$(date +%H:%M:%S)] $*" >&2; }

# A FASTA whose records are called NAME, NAME_2, NAME_3 and so on, the way a
# sample's contigs reach the pipeline.
rename_fasta() {
    awk -v n="$2" '/^>/ {c++; print ">" n (c > 1 ? "_" c : ""); next} {print}' "$1" > "$3"
}

fetch() {
    [ -s "$2" ] && return 0
    log "download $(basename "$1")"
    curl -sSfL --retry 3 -o "$2.gz" "$1"
    gunzip -f "$2.gz"
}

# --- inputs ------------------------------------------------------------------
I=$WORK/inputs
mkdir -p "$I/mtb_al" "$I/mtb_cp" "$I/mtb_cdc1551" "$I/mtb_mbovis" "$I/mtb_draft" "$I/ecoli"

if [ ! -s "$I/mtb_al/L4_H37Rv.fasta" ]; then
    python3 - "$REF/sequence.gb" > "$I/mtb_al/L4_H37Rv.fasta" <<'EOF'
import sys
from Bio import SeqIO
seq = str(next(SeqIO.parse(sys.argv[1], "genbank")).seq)
print(">L4_H37Rv")
print("\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)))
EOF
fi
[ -s "$I/mtb_cp/L4_H37Rv.fasta" ] || rename_fasta "$REF/new_refH37Rv.fasta" L4_H37Rv "$I/mtb_cp/L4_H37Rv.fasta"

fetch "$NCBI/000/008/585/GCF_000008585.1_ASM858v1/GCF_000008585.1_ASM858v1_genomic.fna.gz" "$I/mtb_cdc1551/raw.fna"
[ -s "$I/mtb_cdc1551/L4_CDC1551.fasta" ] || rename_fasta "$I/mtb_cdc1551/raw.fna" L4_CDC1551 "$I/mtb_cdc1551/L4_CDC1551.fasta"
fetch "$NCBI/000/195/835/GCF_000195835.3_ASM19583v2/GCF_000195835.3_ASM19583v2_genomic.fna.gz" "$I/mtb_mbovis/raw.fna"
[ -s "$I/mtb_mbovis/La1_Mbovis.fasta" ] || rename_fasta "$I/mtb_mbovis/raw.fna" La1_Mbovis "$I/mtb_mbovis/La1_Mbovis.fasta"
[ -s "$I/mtb_draft/L4_CDC1551d.fasta" ] || $BENCH draft "$I/mtb_cdc1551/L4_CDC1551.fasta" \
    "$I/mtb_draft/L4_CDC1551d.fasta" "$I/mtb_draft/contig_map.tsv" L4_CDC1551d 40

fetch "$NCBI/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz" "$I/ecoli/k12.gbff"
fetch "$NCBI/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz" "$I/ecoli/k12.fna"
fetch "$NCBI/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz" "$I/ecoli/sakai.fna"
[ -s "$I/ecoli/K12self.fasta" ] || rename_fasta "$I/ecoli/k12.fna" K12self "$I/ecoli/K12self.fasta"
[ -s "$I/ecoli/Sakai.fasta" ] || rename_fasta "$I/ecoli/sakai.fna" Sakai "$I/ecoli/Sakai.fasta"

# minimap2 for the truth liftover; the container has none, and its release
# binary is static. Python unpacks it, since the container has no tar.
MINIMAP2=$WORK/tools/minimap2-${MINIMAP2_VERSION}_x64-linux/minimap2
if [ ! -x "$MINIMAP2" ]; then
    log "download minimap2 $MINIMAP2_VERSION"
    curl -sSfL --retry 3 -o "$WORK/tools/mm2.tar.bz2" \
        "https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2"
    python3 -c "import tarfile, sys; tarfile.open(sys.argv[1]).extractall(sys.argv[2])" \
        "$WORK/tools/mm2.tar.bz2" "$WORK/tools"
    rm -f "$WORK/tools/mm2.tar.bz2"
fi

# --- what the pipeline does before the algorithm --------------------------------
export MPLCONFIGDIR=$(mktemp -d)

# The custom-reference path, PREPARE_REFERENCE, with E. coli K-12 as the reference.
E=$I/ecoli
if [ ! -s "$E/custom_ref_bad_bakta.faa" ]; then
    log "prepare the E. coli K-12 reference (Bakta on the reference)"
    (
        cd "$E"
        python3 "$REPO/bin/extract_ref_proteins.py" --ref-annotation k12.gbff \
            --ref-genome k12.fna --prefix custom_ref
        T=$(mktemp -d)
        bakta --db "$BAKTA_DB" --skip-trna --threads "$THREADS" --genus Escherichia \
            --species "Escherichia coli" --force --output bakta_ref --prefix bakta_ref \
            --locus-tag bakta_ref --tmp-dir "$T" k12.fna
        python3 "$REPO/bin/find_bakta_false_positives.py" --bakta-tsv bakta_ref/bakta_ref.tsv \
            --curated-faa custom_ref_ref_proteins.faa --bakta-faa bakta_ref/bakta_ref.faa \
            --prefix custom_ref
    ) > "$E/prepare_reference.log" 2>&1
fi

# BAKTA_ANNOTATOR and MINIPROT, with the arguments the pipeline gives them.
# Arguments: dir sample genus species bakta-proteins proteome bad-bakta reference-genome
annotate() {
    local dir=$1 s=$2
    [ -s "$dir/$s-ref.paf" ] && return 0
    log "Bakta and miniprot: $s"
    (
        cd "$dir"
        T=$(mktemp -d)
        bakta --db "$BAKTA_DB" --keep-contig-headers --skip-trna --threads "$THREADS" \
            --genus "$3" --species "$4" --force --output "$s" --prefix "$s" --locus-tag "$s" \
            --proteins "$5" --tmp-dir "$T" "$s.fasta"
        miniprot --outs=0.95 -Iut"$THREADS" -S "$s.fasta" "$6" > "$s.paf"
        miniprot --outs=0.95 -Iut"$THREADS" -S "$s.fasta" "$7" > "$s-badbakta.paf"
        miniprot --outs=0.95 -Iut"$THREADS" -S "$8" "$s/$s.faa" > "$s-ref.paf.tmp"
        mv "$s-ref.paf.tmp" "$s-ref.paf"
    ) > "$dir/annotate_$s.log" 2>&1
}

MTB_AA=$REF/aa_H37Rv.translated.fas
MTB=("Mycobacterium" "Mycobacterium tuberculosis" "$REF/sequence.gb" "$MTB_AA"
     "$REF/bad.bakta.translated.fas" "$REF/new_refH37Rv.fasta")
ECO=("Escherichia" "Escherichia coli" "$E/k12.gbff" "$E/custom_ref_ref_proteins.faa"
     "$E/custom_ref_bad_bakta.faa" "$E/k12.fna")
annotate "$I/mtb_al" L4_H37Rv "${MTB[@]}"
annotate "$I/mtb_cp" L4_H37Rv "${MTB[@]}"
annotate "$I/mtb_cdc1551" L4_CDC1551 "${MTB[@]}"
annotate "$I/mtb_mbovis" La1_Mbovis "${MTB[@]}"
annotate "$I/mtb_draft" L4_CDC1551d "${MTB[@]}"
annotate "$E" K12self "${ECO[@]}"
annotate "$E" Sakai "${ECO[@]}"

# --- the truth ------------------------------------------------------------------
T=$WORK/truth
# Tables from before the 'aligned' column (October 2026) are made again.
if [ ! -s "$T/sakai.tsv" ] || ! head -1 "$T/sakai.tsv" | grep -qw aligned; then
    log "truth tables"
    $BENCH truth-genbank "$REF/sequence.gb" "$MTB_AA" "$T/al.tsv"
    $BENCH truth-liftover "$REF/sequence.gb" "$I/mtb_cp/L4_H37Rv.fasta" "$T/al.tsv" "$T/cp_lifted.tsv" \
        --minimap2 "$MINIMAP2" --threads "$THREADS"
    $BENCH truth-curated "$T/cp_lifted.tsv" "$I/mtb_cp/L4_H37Rv.fasta" "$MTB_AA" "$T/cp.tsv"
    $BENCH truth-liftover "$REF/sequence.gb" "$I/mtb_cdc1551/L4_CDC1551.fasta" "$T/al.tsv" "$T/cdc1551.tsv" \
        --minimap2 "$MINIMAP2" --threads "$THREADS"
    $BENCH truth-liftover "$REF/sequence.gb" "$I/mtb_mbovis/La1_Mbovis.fasta" "$T/al.tsv" "$T/mbovis.tsv" \
        --minimap2 "$MINIMAP2" --threads "$THREADS"
    $BENCH truth-genbank "$E/k12.gbff" "$E/custom_ref_ref_proteins.faa" "$T/k12.tsv"
    $BENCH truth-liftover "$E/k12.gbff" "$E/Sakai.fasta" "$T/k12.tsv" "$T/sakai.tsv" \
        --minimap2 "$MINIMAP2" --threads "$THREADS"
fi

# --- the algorithm ----------------------------------------------------------------
# TYPE_CASES and FIX_CDS as modules/local/cases.nf runs them.
R=$WORK/runs/$RUN
rm -rf "${R:?}"
mkdir -p "$R"
correct() {
    local label=$1 dir=$2 s=$3 proteome=$4 out=$R/$1
    mkdir -p "$out"
    (
        cd "$out"
        python3 "$REPO/bin/classify_cases.py" -b "$dir/$s/$s.tsv" -p "$dir/$s.paf" \
            -bb "$dir/$s-badbakta.paf" -c "$dir/$s-ref.paf" \
            --genome "$dir/$s.fasta" --ref-proteome "$proteome" > classify.log 2>&1
        python3 "$REPO/bin/check_cds.py" "$dir/$s.fasta" "${s}_corrected.tsv" \
            --ref-proteome "$proteome" > check.log 2>&1
        python3 "$REPO/bin/deduplicate_rvs.py" "$s.AMAP.tsv" --faa "$s.AMAP.faa" --fna "$s.AMAP.fna" \
            --bakta "$dir/$s/$s.tsv" --ref-proteome "$proteome" > dedup.log 2>&1
    ) || { echo "ERROR: the algorithm failed on $label, see $out" >&2; exit 1; }
}

log "correction algorithm, run '$RUN'"
correct H37Rv-AL "$I/mtb_al" L4_H37Rv "$MTB_AA"
correct H37Rv-CP "$I/mtb_cp" L4_H37Rv "$MTB_AA"
correct CDC1551 "$I/mtb_cdc1551" L4_CDC1551 "$MTB_AA"
correct M.bovis "$I/mtb_mbovis" La1_Mbovis "$MTB_AA"
correct CDC1551-draft "$I/mtb_draft" L4_CDC1551d "$MTB_AA"
correct K-12 "$E" K12self "$E/custom_ref_ref_proteins.faa"
correct Sakai "$E" Sakai "$E/custom_ref_ref_proteins.faa"

# --- scores ------------------------------------------------------------------------
echo
echo "Run '$RUN' of $(cd "$REPO" && git rev-parse --short HEAD 2>/dev/null || echo 'this tree'):"
S=$R/scores.tsv
$BENCH score H37Rv-AL "$T/al.tsv" "$R/H37Rv-AL/L4_H37Rv.AMAP.tsv" --tsv "$S"
$BENCH score H37Rv-CP "$T/cp.tsv" "$R/H37Rv-CP/L4_H37Rv.AMAP.tsv" --tsv "$S"
$BENCH score CDC1551 "$T/cdc1551.tsv" "$R/CDC1551/L4_CDC1551.AMAP.tsv" --tsv "$S"
$BENCH score M.bovis "$T/mbovis.tsv" "$R/M.bovis/La1_Mbovis.AMAP.tsv" --tsv "$S"
$BENCH score CDC1551-draft "$T/cdc1551.tsv" "$R/CDC1551-draft/L4_CDC1551d.AMAP.tsv" --tsv "$S" \
    --contig-map "$I/mtb_draft/contig_map.tsv" --fasta "$I/mtb_draft/L4_CDC1551d.fasta" \
    --faa "$R/CDC1551-draft/L4_CDC1551d.AMAP.faa"
$BENCH score K-12 "$T/k12.tsv" "$R/K-12/K12self.AMAP.tsv" --tsv "$S"
$BENCH score Sakai "$T/sakai.tsv" "$R/Sakai/Sakai.AMAP.tsv" --tsv "$S"
echo "Scores in $S"

if [ -n "$COMPARE" ]; then
    echo
    echo "Changes from run '$COMPARE':"
    for g in H37Rv-AL:L4_H37Rv H37Rv-CP:L4_H37Rv CDC1551:L4_CDC1551 M.bovis:La1_Mbovis \
             CDC1551-draft:L4_CDC1551d K-12:K12self Sakai:Sakai; do
        label=${g%%:*} s=${g##*:}
        old=$WORK/runs/$COMPARE/$label/$s.AMAP.tsv
        [ -s "$old" ] || { echo "$label: no run '$COMPARE' to compare with"; continue; }
        echo "$label:"
        $BENCH compare "$old" "$R/$label/$s.AMAP.tsv"
    done
fi
