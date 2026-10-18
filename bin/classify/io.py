"""I/O: readers for Bakta TSV, Miniprot PAF, and summary writers.

Every frame carries the sequence each feature sits on as 'Contig'. A draft
assembly numbers each contig from 1, so coordinates alone say nothing about
which contig a gene is on, and anything that compares them has to compare the
contig too. Bakta must run with --keep-contig-headers for its names to match
the FASTA and miniprot's: by default it renames every sequence to contig_N.
"""

import pandas as pd
from .ids import clean_id

# The columns -rvmissing.tsv has always had; the R report reads it by position.
PAF_COLS = ["Prot Id", "Prot len", "Strand", "Start", "End", "Match Nt", "Num nt"]


def read_bakta(bakta_path):
    """Read and process a Bakta TSV file, returning CDS-only entries."""
    bakta_df = pd.read_csv(bakta_path, comment='#', sep='\t', header=None)
    bakta_df.columns = [
        "Sequence Id", "Type", "Start", "End", "Strand",
        "Locus Tag", "Gene", "Product", "DbXrefs"
    ]
    bakta_df['Locus Tag'] = bakta_df['Locus Tag'].apply(clean_id)
    bakta_cds = bakta_df[bakta_df['Type'] == 'cds'].copy()
    bakta_cds['Length_nt'] = bakta_cds['End'] - bakta_cds['Start'] + 1
    bakta_cds['Contig'] = bakta_cds['Sequence Id'].astype(str)
    return bakta_cds[["Locus Tag", "Start", "End", "Strand", 'Length_nt', 'Contig']]


def _custom_filter(group):
    """Keep only perfect-match rows if any exist; otherwise keep all."""
    perfect = group[group['Num nt'] - group['Match Nt'] == 0]
    return perfect if len(perfect) > 0 else group


def _read_paf(path):
    """Read a PAF file robustly, handling variable numbers of tag columns.
    
    PAF has 12 fixed columns; columns 13+ are optional SAM-style tags
    whose count varies by aligner version and parameters. We read all
    columns and rename only the first 12 that we actually use.

    The tag count also varies BETWEEN ROWS of a single file: miniprot emits
    fs:i only for frameshifted alignments, and cg:Z/cs:Z depend on the run
    flags. pandas infers the column count from the first row, so a later row
    with more tags aborts the parse with a ParserError. We therefore scan the
    file for the widest row first and read with that many explicit columns;
    shorter rows are padded with NaN.
    """
    paf_cols = [
        "Prot Id", "Prot len", "Prot Start", "ProtStop", "Strand",
        "fasta Tag", "fasta len", "Start", "End", "Match Nt", "Num nt",
        "mapping qual"
    ]
    try:
        with open(path) as fh:
            n_cols = max(
                (line.count('\t') + 1 for line in fh
                 if line.strip() and not line.startswith('#')),
                default=0
            )
    except OSError:
        n_cols = 0

    if n_cols == 0:
        # Empty PAF (e.g., empty bad_bakta.faa → no alignments)
        return pd.DataFrame(columns=paf_cols)

    try:
        # index_col=False keeps pandas from promoting leading columns to the
        # index when `names` is wider than a given row.
        df = pd.read_csv(path, sep='\t', header=None, comment='#',
                         names=range(n_cols), index_col=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=paf_cols)
    # Rename first 12 columns; keep remaining as-is (tag columns)
    rename_map = {i: name for i, name in enumerate(paf_cols) if i < len(df.columns)}
    df = df.rename(columns=rename_map)
    return df


def _tag(df, name):
    """A numeric SAM-style tag ('st:i:1') of every row of a PAF; 0 where absent.

    Tags are not at fixed columns: miniprot writes fs:i only for alignments
    with a frameshift, and none at all for a protein it could not place.
    """
    value = pd.Series(0, index=df.index, dtype=float)
    for col in [c for c in df.columns if isinstance(c, int)]:
        text = df[col].astype(str)
        hit = text.str.startswith(name + ':')
        value[hit] = pd.to_numeric(text[hit].str.split(':').str[2], errors='coerce')
    return value.fillna(0)


def read_miniprot(miniprot_path, ref_coords=False):
    """Read Miniprot PAF output, apply best-match filtering per protein.

    For multi-copy genes (PE/PPE, IS elements), the same protein can map
    to multiple genomic positions. We deduplicate per (Prot Id, Contig, Start,
    End), keeping the row with the highest Match Nt, and prefer perfect matches.

    clean_id strips the reference coordinates from each protein ID. With
    ref_coords=True they are kept as 'Ref Start' and 'Ref End' for
    synteny.drop_paralog_placements, which drops them again.
    """
    miniprot_df = _read_paf(miniprot_path)
    extra = []
    if ref_coords:
        coords = miniprot_df['Prot Id'].astype(str).str.extract(r'_(\d+)-(\d+)$')
        miniprot_df['Ref Start'] = pd.to_numeric(coords[0])
        miniprot_df['Ref End'] = pd.to_numeric(coords[1])
        extra = ['Ref Start', 'Ref End']
    miniprot_df['Contig'] = miniprot_df['fasta Tag'].astype(str)
    # Whether the alignment is the whole gene: every residue of the reference
    # protein, its final '*' included, so the sample's gene ends in a stop
    # codon right where the reference one does. A protein the reference
    # annotates without a stop (Rv3022c runs on into Rv3021c by a frameshift)
    # aligns end to end with no stop codon (st:i:0), and the sample's reading
    # frame may go on well past it: in CDC1551 an 81-residue copy of Rv3022c
    # opens a 933 bp PPE gene.
    miniprot_df['Complete'] = ((pd.to_numeric(miniprot_df['Prot Start'], errors='coerce') == 0)
                               & (pd.to_numeric(miniprot_df['ProtStop'], errors='coerce')
                                  >= pd.to_numeric(miniprot_df['Prot len'], errors='coerce'))
                               & (_tag(miniprot_df, 'st') > 0))
    miniprot_df['Prot Id'] = miniprot_df['Prot Id'].apply(clean_id)
    miniprot_df['Start'] = miniprot_df['Start'] + 1  # 0-based → 1-based

    last = miniprot_df[PAF_COLS + ['Contig', 'Complete'] + extra].copy()

    # Mark perfect matches (aligned nt == total nt)
    last['_is_perfect'] = (last['Num nt'] - last['Match Nt'] == 0).astype(int)

    # Sort so perfect matches and highest Match Nt come first
    last = last.sort_values(['Prot Id', 'Contig', 'Start', 'End', '_is_perfect', 'Match Nt'],
                            ascending=[True, True, True, True, False, False])

    # Keep the best hit per placement. That preserves multi-copy placements, including
    # copies at the same coordinates of different contigs.
    last = last.drop_duplicates(subset=['Prot Id', 'Contig', 'Start', 'End'], keep='first')
    last = last.drop(columns=['_is_perfect']).reset_index(drop=True)

    return last


def read_miniprot_unfiltered(miniprot_path):
    """Read Miniprot PAF without best-match filtering (for 'bad bakta' alignments)."""
    miniprot_df = _read_paf(miniprot_path)
    miniprot_df['Contig'] = miniprot_df['fasta Tag'].astype(str)
    miniprot_df['Prot Id'] = miniprot_df['Prot Id'].apply(clean_id)
    miniprot_df['Start'] = miniprot_df['Start'] + 1
    return miniprot_df[PAF_COLS + ['Contig']]


def write_summary(base_name, data_frames, case_types, sample=None, lineage=None):
    """Write a summary TSV with counts per case type.

    Sequence is the name the sample's files carry, Sample its ID in the sample
    sheet and Lineage what lineage typing called it. The pipeline passes the
    last two: typing names a genome <lineage>_<sampleId>, and a sample ID can
    hold underscores of its own, so neither can be read back out of the name
    (L4_H37Rv_test, or H37Rv_test under --skip_typing). Without them the sample
    is its name and nothing typed it.
    """
    summary_data = {
        "Type": case_types,
        "Number": [len(df) for df in data_frames],
        "Sequence": base_name,
        "Lineage": lineage or "Untyped",
        "Sample": sample or base_name,
    }
    pd.DataFrame(summary_data).to_csv(f"{base_name}-summary.tsv", sep='\t', index=False)
