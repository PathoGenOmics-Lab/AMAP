---
title: Contributing
---

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

```bash
pip install -r requirements-docs.txt
mkdocs serve
```

That serves the site at `http://127.0.0.1:8000` and reloads as you edit. Pages
live in `docs/` and the navigation is the `nav` block in `mkdocs.yml`. Every
page has an edit link in its top right corner, which goes straight to the file
on GitHub.

The site is built with `--strict`, so a broken internal link fails the build
rather than shipping. Run `mkdocs build --strict` before opening a pull request.

Publishing is automatic: a push to `main` that touches `docs/`, `mkdocs.yml`,
`CHANGELOG.md` or the workflow itself rebuilds and deploys to GitHub Pages.

## Style

Keep to what is already there: short paragraphs, admonitions for the things that
actually bite, and concrete numbers rather than "large" or "slow". Where a
threshold appears in the docs it should match the constant in the code, so if
you change one, change the other.
