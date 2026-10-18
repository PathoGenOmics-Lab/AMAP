# Contributing

Bug reports, fixes and suggestions are all welcome. The full guidelines live in
[CONTRIBUTING.md](https://github.com/PathoGenOmics-Lab/AMAP/blob/main/.github/CONTRIBUTING.md),
and the project follows a
[code of conduct](https://github.com/PathoGenOmics-Lab/AMAP/blob/main/.github/CODE_OF_CONDUCT.md).

## Reporting a bug

Open an issue at
[github.com/PathoGenOmics-Lab/AMAP/issues](https://github.com/PathoGenOmics-Lab/AMAP/issues)
with the command you ran, the failing task's `.command.err`, and both
`software_versions.yml` and `params.yml` from `pipeline_info/`. Those last two
pin the exact versions and parameters, which is usually what the answer turns
on. [Troubleshooting](../deployment/troubleshooting.md) covers the common cases
first.

## Working on the code

```bash
git clone https://github.com/PathoGenOmics-Lab/AMAP
cd AMAP
pip install -r requirements-test.txt
python -m pytest tests/ -v
```

`requirements-test.txt` covers everything `bin/` imports, not only what the
tests import directly, and CI installs from the same file. A shorter list
appears to work until a test touches `classify_cases.py`, which pulls in
plotly through `classify/plotting.py`.

The suite covers ID handling, overlap arithmetic, case classification, the
rules that resolve and recover genes, and the parsers, and it runs in about a second. Anything
touching `bin/classify/` should come with a test, since the correction logic is
where a silent change does the most damage.

### Annotating H37Rv with itself

A change to the correction algorithm should also pass the self-test, which has
exactly one right answer: annotating H37Rv against its own reference must give
back every CDS of the AL123456.3 GenBank record, at exactly its coordinates and
strand, under its own name.

```bash
docker run --rm --platform linux/amd64 -v "$PWD:/amap" -w /amap \
    paururo/amap:latest bash tests/selftest/h37rv_selftest.sh
```

```text
GenBank CDS: 4031   AMAP genes: 4031
  exact and named right: 4031
  not at exactly their coordinates: 0
  named wrong: 0
  AMAP genes that are no GenBank CDS: 0
```

The first run spends a couple of minutes on Bakta and miniprot and keeps their
output under `selftest_h37rv/`. After that it only runs `classify_cases.py`,
`check_cds.py` and `deduplicate_rvs.py`, in seconds, and exits non-zero with
the offending genes listed if anything is missing, misplaced or misnamed.

H37Rv is one genome, and the one the defaults were tuned on. The accuracy
benchmark in `tests/benchmark/` scores six more whose right answer is known:
another H37Rv assembly, CDC1551 complete and cut into 40 contigs, M. bovis, and
E. coli K-12 and Sakai against K-12 as a custom reference. Give a run a name,
and the next one lists every gene that changed:

```bash
docker run --rm --platform linux/amd64 -v "$PWD:/amap" -w /amap \
    paururo/amap:latest bash tests/benchmark/benchmark.sh benchmark_work before
# change the algorithm, then
docker run --rm --platform linux/amd64 -v "$PWD:/amap" -w /amap \
    paururo/amap:latest bash tests/benchmark/benchmark.sh benchmark_work after before
```

`tests/benchmark/README.md` explains how the truth is built and what to expect.

The quickest way to run the whole pipeline is the test profile:

```bash
nextflow run . -profile docker,test
```

### Pipeline code and `-resume`

Nextflow hashes a task's inputs and script, not the files the script calls by
path, so `python3 ${projectDir}/bin/check_cds.py` alone would let `-resume` hand
back what the old `check_cds.py` produced. Every process that runs AMAP code
therefore lists those files in `ext code`, and its script names `task.ext.code`
in a comment, which puts them in the cache key (Nextflow 24.04 and later):

```groovy
ext code: { files("${projectDir}/bin/{check_cds.py,deduplicate_rvs.py,classify/*.py}") }
```

A script that imports `classify` needs `classify/*.py` in the list. The reverse
also bites: Nextflow adds to a task's cache key any `bin/` entry its script names
as a bare word, which is why PATHOTYPR quotes its `'classify'` subcommand.
`tests/test_cache_keys.py` fails when a process runs a file it does not list,
lists one that does not exist, or names a `bin/` entry it does not declare.

## Working on the container

```bash
docker build --platform linux/amd64 -f docker/Dockerfile -t amap:dev .
docker/smoke_test.sh amap:dev
```

The smoke test checks the tools, the Python and R libraries, the bundled
reference data, and runs the test suite with the image interpreter. Worth doing
before proposing a change to the Dockerfile, since a library update can break
the pipeline scripts without the image looking any different.

## Working on these docs

The website is plain HTML in `docs/`, with nothing to build on the server:
GitHub Pages serves the folder as it is, and it opens from a clone too. The long
pages are the markdown files, which GitHub also renders, and
`tools/build_docs.py` writes each one's `.html` next to it, with the site's
header, page list, contents and search.

```bash
python3 tools/build_docs.py          # after editing any .md under docs/
python3 tools/check_site.py          # what CI runs
python3 -m http.server -d docs 8000  # then open http://localhost:8000
```

The front page (`docs/index.html`) and the illustrated walk through the
correction (`docs/how-it-works.html`) are written by hand. Their drawings come
from `docs/assets/figures.js` and `how.js`, and every gene in them from the
windows of real runs in `docs/assets/figure-data.js`. A line
`<!-- figure: NAME -->` in a markdown page puts one of those drawings there, and
GitHub shows nothing in its place.

`check_site.py` fails when a page no longer matches its markdown, a link or an
anchor is broken, a page uses an option the pipeline does not have, a threshold
on a page differs from the code, the metro map's stations are out of the order
`main.nf` runs them in, or a drawing shows something `bin/` would not do: it
runs every window of `figure-data.js` through the classification,
`check_cds.py` and `deduplicate_rvs.py`. When a change to the algorithm changes
a drawn window, rerun the benchmark, write the windows again, and read the
words of the drawing against what changed:

```bash
python3 tools/figure_data.py benchmark_work final --reference /reference
```

The workflow in `.github/workflows/docs.yml` runs both checks on every pull
request that touches the docs, `bin/`, `modules/` or `main.nf`.

## Style

Keep to what is already there: short paragraphs, GitHub's alerts (`> [!WARNING]`)
for the things that actually bite, and concrete numbers rather than "large" or
"slow". Where a threshold appears in the docs it should match the constant in
the code, so if you change one, change the other; `check_site.py` holds the ones
it knows to the code.
