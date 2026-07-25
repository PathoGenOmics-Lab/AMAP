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
reconciliation rules, and the parsers, and it runs in about a second. Anything
touching `bin/classify/` should come with a test, since the correction logic is
where a silent change does the most damage.

The quickest way to run the whole pipeline is the test profile:

```bash
nextflow run . -profile docker,test
```

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
