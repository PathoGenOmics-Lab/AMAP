#!/usr/bin/env bash
# Check that an AMAP image is actually usable before it gets published.
#
# Verifies the command line tools, the Python and R libraries the pipeline
# imports, the bundled reference data, and then runs the repository test suite
# with the image interpreter. That last part is the point of the whole script:
# importing a module only proves its imports resolve, so it would happily pass
# on an image whose Biopython no longer has SeqFeature.strand. Running the
# tests exercises the code paths and catches that kind of drift.
#
# Usage:
#   docker/smoke_test.sh [image] [repo_dir]
#
#   image     image to test        (default: amap:ci)
#   repo_dir  AMAP checkout to test (default: the repo this script lives in)
#
# Run it against a local build with:
#   docker build --platform linux/amd64 -f docker/Dockerfile -t amap:ci .
#   docker/smoke_test.sh amap:ci

set -euo pipefail

IMAGE="${1:-amap:ci}"
REPO_DIR="${2:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

# The checkout goes in read-only and is copied to a scratch dir inside the
# container, because py_compile and pytest both want to write next to it.
run() { docker run --rm -v "${REPO_DIR}:/src:ro" "$IMAGE" "$@"; }

echo "Testing image: $IMAGE"
echo "Against repo:  $REPO_DIR"
echo

echo "== command line tools =="
run bash -lc '
  set -e
  for t in miniprot pathotypr diamond bakta emapper.py pangene gfatools pandoc; do
    command -v "$t" >/dev/null || { echo "  MISSING: $t"; exit 1; }
    echo "  ok  $t  ($(command -v "$t"))"
  done
  echo "  ok  $(R --version | head -1)"
'

echo
echo "== python packages used by bin/ =="
run python3 -c '
import importlib, sys
mods = ["Bio", "community", "intervaltree", "matplotlib", "networkx",
        "numpy", "pandas", "plotly", "pyvis"]
bad = []
for m in mods:
    try:
        importlib.import_module(m)
        print("  ok ", m)
    except Exception as exc:
        bad.append("%s: %s" % (m, exc))
if bad:
    print("  MISSING:")
    for b in bad:
        print("   ", b)
    sys.exit(1)
'

echo
echo "== R packages used by the report =="
run Rscript -e '
  pkgs <- c("DT","RColorBrewer","dplyr","ggplot2","ggpubr","ggrepel","knitr",
            "pheatmap","plotly","rstatix","scales","stringr","tidyr","rmarkdown")
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing)) { cat("  MISSING:", missing, "\n"); quit(status = 1) }
  cat("  ok ", length(pkgs), "packages\n")
'

echo
echo "== bundled reference data =="
run bash -lc '
  set -e
  for f in /reference/MTB_ancestor_reference.fasta \
           /reference/aa_H37Rv.translated.fas \
           /reference/bad.bakta.translated.fas \
           /reference/new_refH37Rv.fasta \
           /reference/sequence.gb \
           /reference/Mycobacterium_tuberculosis_H37Rv_txt_v4.txt \
           /reference/snps_lineages.tsv; do
    [ -s "$f" ] || { echo "  MISSING or empty: $f"; exit 1; }
    echo "  ok  $f ($(stat -c%s "$f") bytes)"
  done
  [ -d /data/bakta_db ] && [ -n "$(ls -A /data/bakta_db)" ] \
    || { echo "  MISSING: Bakta database at /data/bakta_db"; exit 1; }
  echo "  ok  /data/bakta_db"
'

echo
echo "== repository code against the image interpreter =="
# The scripts live in the repo, not the image, and Nextflow runs them from
# projectDir. This is where the image and the code get checked together.
run bash -lc '
  set -e
  cp -r /src /work
  cd /work

  for s in bin/*.py; do
    python3 -m py_compile "$s" || { echo "  FAILED to compile: $s"; exit 1; }
  done
  echo "  ok  all bin/*.py compile"

  PYTHONPATH=/work/bin python3 -c "
import importlib, pathlib, sys
names = [p.stem for p in sorted(pathlib.Path(\"/work/bin\").glob(\"*.py\"))]
names += [\"classify.\" + p.stem
          for p in sorted(pathlib.Path(\"/work/bin/classify\").glob(\"*.py\"))
          if p.stem != \"__init__\"]
bad = []
for n in names:
    try:
        importlib.import_module(n)
    except Exception as exc:
        bad.append(\"%s: %s\" % (n, exc))
if bad:
    print(\"  FAILED to import:\")
    for b in bad:
        print(\"   \", b)
    sys.exit(1)
print(\"  ok  %d pipeline modules import\" % len(names))
"

  # pytest is not part of the runtime image, so pull it in just for this.
  pip install --quiet --no-cache-dir pytest >/dev/null 2>&1
  python3 -m pytest tests/ -q -p no:cacheprovider
'

echo
echo "All smoke tests passed."
