#!/bin/bash
# GET_DUPLICATIONS on each corrected proteome, then GET_ALLDUP on what they
# found, with the commands of modules/local/get.nf (tools/check_site.py checks
# that the Diamond options are still the same). rename_fastas.py is left out:
# it only appends each protein's eggNOG COG letter to its name, and needs
# eggNOG's output. Run in the AMAP container, for tools/dup_data.py:
#
#   bash tools/dup_run.sh WORK L4_CDC1551.AMAP.faa La1_Mbovis.AMAP.faa
set -euo pipefail
work=$1
shift
mkdir -p "$work/all"
bin=$(cd "$(dirname "$0")/../bin" && pwd)
for faa in "$@"; do
    sampleId=$(basename "$faa" .AMAP.faa)
    cp "$faa" "$work/$sampleId.AMAP.faa.renamed"
    (
        cd "$work"
        diamond makedb --in "$sampleId.AMAP.faa.renamed" --db renamed_db --threads 6
        diamond blastp \
            --db renamed_db \
            -q "$sampleId.AMAP.faa.renamed" \
            --threads 6 \
            --tmpdir "${TMPDIR:-/tmp}" \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen nident \
            --very-sensitive \
            --no-self-hits \
            --header \
            > "duplications_$sampleId.raw"
        grep '^#' "duplications_$sampleId.raw" > "duplications_$sampleId" || true
        grep -v '^#' "duplications_$sampleId.raw" | sort -k1,1 -k2,2 -k12,12gr >> "duplications_$sampleId"
        rm "duplications_$sampleId.raw"
        python3 "$bin/get_duplications.py" "duplications_$sampleId"
        python3 "$bin/dup_ex_multifasta.py" "duplications_${sampleId}_groups.tsv" "$sampleId.AMAP.faa.renamed"
        mv "$sampleId.AMAP.faa.renamed.dup" "$sampleId.dup"
    )
done
cd "$work/all"
cat ../*.dup > all_fastas.fa
# Which genome each sequence comes from: the sample its file is named after.
for f in ../*.dup; do
    genome=$(basename "$f" .dup)
    awk -v genome="$genome" '/^>/ { print substr($1, 2) "\t" genome }' "$f"
done > genomes.tsv
diamond makedb --in all_fastas.fa --db all_fastas_db --threads 6
diamond blastp \
    --db all_fastas_db \
    -q all_fastas.fa \
    --threads 6 \
    --tmpdir "${TMPDIR:-/tmp}" \
    --block-size 2 \
    --index-chunks 1 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen nident \
    --very-sensitive \
    --no-self-hits \
    --header \
    > all_vs_all.raw
grep '^#' all_vs_all.raw > all_vs_all.diamond || true
grep -v '^#' all_vs_all.raw | sort -k1,1 -k2,2 -k12,12gr >> all_vs_all.diamond
rm all_vs_all.raw
python3 "$bin/network_genomes.py" all_vs_all.diamond --genomes genomes.tsv
