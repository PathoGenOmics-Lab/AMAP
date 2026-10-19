"""Everything that cites AMAP cites one DOI: the concept DOI of the manifest.

Zenodo gives AMAP a concept DOI, which always resolves to the latest version,
and a DOI for each version. main.nf, the manifest and the new report cited the
concept DOI, while the README, CITATION.cff and every page of the website cited
the DOI of a single version.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEXT = {'.md', '.cff', '.nf', '.config', '.py', '.js', '.html', '.css', '.Rmd', '.txt', '.yml'}
DOI = re.compile(r'10\.5281(?:/|%2F)zenodo\.(\d+)')
# Other works on Zenodo the docs cite, which are not AMAP
OTHER_WORKS = {
    '19210044': 'the pathotypr lineage markers',
}


def _concept():
    m = re.search(r"doi\s*=\s*'https://doi\.org/10\.5281/zenodo\.(\d+)'", (ROOT / 'nextflow.config').read_text())
    assert m, 'nextflow.config has no doi in its manifest'
    return m.group(1)


def _sources():
    yield from (ROOT / name for name in ('README.md', 'CITATION.cff', 'main.nf', 'nextflow.config'))
    for folder in ('bin', 'docs', 'tools', 'modules'):
        yield from (p for p in sorted((ROOT / folder).rglob('*')) if p.is_file() and p.suffix in TEXT)


def test_one_doi_everywhere():
    cited = {}
    for path in _sources():
        for record in DOI.findall(path.read_text(encoding='utf-8', errors='replace')):
            cited.setdefault(record, set()).add(str(path.relative_to(ROOT)))
    concept = _concept()
    assert concept in cited
    others = set(cited) - {concept} - set(OTHER_WORKS)
    assert not others, {r: sorted(cited[r])[:5] for r in others}
