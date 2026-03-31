"""I/O: readers for Bakta TSV, Miniprot PAF, and summary writers."""

import pandas as pd
from .ids import clean_id


def read_bakta(bakta_path):
    """Read and process a Bakta TSV file, returning CDS-only entries."""
    bakta_df = pd.read_csv(bakta_path, comment='#', sep='\t', header=None)
    bakta_df.columns = [
        "Sequence Id", "Type", "Start", "End", "Strand",
        "Locus Tag", "Gene", "Product", "DbXrefs"
    ]
    bakta_df['Locus Tag'] = bakta_df['Locus Tag'].apply(clean_id)
    bakta_cds = bakta_df[bakta_df['Type'] == 'cds'].copy()
    bakta_cds['Length_nt'] = bakta_cds['End'] - bakta_cds['Start']
    return bakta_cds[["Locus Tag", "Start", "End", "Strand", 'Length_nt']]


def _custom_filter(group):
    """Keep only perfect-match rows if any exist; otherwise keep all."""
    perfect = group[group['Num nt'] - group['Match Nt'] == 0]
    return perfect if len(perfect) > 0 else group


def read_miniprot(miniprot_path):
    """Read Miniprot PAF output, apply best-match filtering per protein."""
    miniprot_df = pd.read_csv(miniprot_path, sep='\t', header=None)
    miniprot_df.columns = [
        "Prot Id", "Prot len", "Prot Start", "ProtStop", "Strand",
        "fasta Tag", "fasta len", "Start", "End", "Match Nt", "Num nt",
        "mapping qual", "AS", "MS", "NP", "FS", "ST", "DA", "DO", "CG", "CS"
    ]
    miniprot_df['Prot Id'] = miniprot_df['Prot Id'].apply(clean_id)
    miniprot_df['Start'] = miniprot_df['Start'] + 1  # 0-based → 1-based

    last = miniprot_df[["Prot Id", "Prot len", "Strand", "Start", "End", "Match Nt", "Num nt"]]

    def _apply_filter(group):
        result = _custom_filter(group)
        return result.drop(columns=['Prot Id'], errors='ignore')

    try:
        return last.groupby('Prot Id').apply(_apply_filter, include_groups=False).reset_index()
    except TypeError:
        # pandas < 2.2 compatibility
        return last.groupby('Prot Id').apply(_apply_filter).reset_index()


def read_miniprot_unfiltered(miniprot_path):
    """Read Miniprot PAF without best-match filtering (for 'bad bakta' alignments)."""
    miniprot_df = pd.read_csv(miniprot_path, sep='\t', header=None)
    miniprot_df.columns = [
        "Prot Id", "Prot len", "Prot Start", "ProtStop", "Strand",
        "fasta Tag", "fasta len", "Start", "End", "Match Nt", "Num nt",
        "mapping qual", "AS", "MS", "NP", "FS", "ST", "DA", "DO", "CG", "CS"
    ]
    miniprot_df['Prot Id'] = miniprot_df['Prot Id'].apply(clean_id)
    miniprot_df['Start'] = miniprot_df['Start'] + 1
    return miniprot_df[["Prot Id", "Prot len", "Strand", "Start", "End", "Match Nt", "Num nt"]]


def write_summary(base_name, data_frames, case_types):
    """Write a summary TSV with counts per case type."""
    summary_data = {
        "Type": case_types,
        "Number": [len(df) for df in data_frames],
        "Sequence": base_name,
        "Lineage": base_name.split("_")[0]
    }
    pd.DataFrame(summary_data).to_csv(f"{base_name}-summary.tsv", sep='\t', index=False)
