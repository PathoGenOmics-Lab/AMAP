#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import logging
import re
import warnings

__author__ = 'Paula Ruiz-Rodriguez'
__credits__ = ['Paula Ruiz-Rodriguez', 'PathoGenOmics']
__license__ = ['GPL 3']
__version__ = '1.0.0'
__maintainer__ = '@paururo'
__email__ = 'paula.ruiz.rodriguez@csic.es'
__status__ = 'finished'
__lab__ = 'PathoGenOmics (I2SysBio)'

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Redirect warnings to logging
def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return f'{filename}:{lineno}: {category.__name__}: {message}\n'

warnings.formatwarning = warning_on_one_line
logging.captureWarnings(True)

def find_next_stop_codon(fasta_seq, start_index, strand, max_distance=300):
    """Find the next in-frame stop codon in the DNA sequence.
    
    For '+' strand: scans forward from start_index in steps of 3.
    For '-' strand: scans backward (toward lower coordinates) from start_index
                    in steps of 3, reading codons on the reverse complement.
    
    Args:
        fasta_seq: Biopython SeqRecord of the genome.
        start_index: 0-based position to start scanning from.
        strand: '+' or '-'.
        max_distance: Maximum distance (bp) to search before giving up.
    
    Returns the 0-based position of the first base of the stop codon, or None.
    """
    sequence = fasta_seq.seq
    stop_codons = {"TAA", "TAG", "TGA"}

    logging.info(f"Finding next stop codon from index {start_index} on strand {strand}")

    if strand == "+":
        end_limit = min(start_index + max_distance, len(sequence) - 2)
        for i in range(start_index, end_limit, 3):
            codon = str(sequence[i:i + 3])
            if codon in stop_codons:
                logging.info(f"Found stop codon '{codon}' at position {i}")
                return i
    else:
        # On minus strand, codons are read as reverse complement.
        # We scan toward lower coordinates (3' direction on minus strand).
        end_limit = max(start_index - max_distance, 0)
        for i in range(start_index, end_limit, -3):
            if i < 0 or i + 3 > len(sequence):
                continue
            codon = str(sequence[i:i + 3].reverse_complement())
            if codon in stop_codons:
                logging.info(f"Found stop codon '{codon}' (rc) at position {i}")
                return i

    logging.warning(f"No stop codon found within {max_distance}bp")
    return None

def get_goodseq(sequence, start, end, strand):
    """Extend a CDS to include its stop codon, adjusting frame if needed.
    
    Only extends if the current sequence does NOT already end with a stop codon.
    This prevents extending genes that already have the stop included in their
    coordinates (as is the convention in Sanger/RefSeq annotations).
    
    Args:
        sequence: Biopython SeqRecord of the genome.
        start: 0-based start coordinate (always < end).
        end: 0-based end coordinate (exclusive).
        strand: '+' or '-'.
    
    Returns:
        Tuple of (start, end, strand) with adjusted coordinates.
    """
    logging.info(f"Getting good sequence from {start} to {end} on strand {strand}")
    seq = sequence.seq

    if strand == "+":
        # Ensure length is divisible by 3
        resto = (end - start) % 3
        if resto != 0:
            end = end - resto
        # Check if gene already ends with stop codon
        if end >= start + 3:
            last_codon = str(seq[end - 3:end])
            if last_codon in {"TAA", "TAG", "TGA"}:
                logging.info(f"Gene already ends with stop codon '{last_codon}' — no extension needed")
                return start, end, strand
        # Search for stop codon starting at current end
        new_end = find_next_stop_codon(sequence, end, strand)
        if new_end is not None:
            end = new_end + 3  # Include the stop codon
    else:
        # For minus strand, the 3' end is at the lower coordinate (start).
        length = end - start
        resto = length % 3
        if resto != 0:
            start = start + resto  # Trim from the 3' end (lower coord)
        # Check if gene already has stop codon at 3' end (lower coordinate on minus strand)
        if start + 3 <= end:
            first_codon_rc = str(seq[start:start + 3].reverse_complement())
            if first_codon_rc in {"TAA", "TAG", "TGA"}:
                logging.info(f"Minus-strand gene already ends with stop codon '{first_codon_rc}' (rc) — no extension needed")
                return start, end, strand
        # Search for stop codon scanning toward lower coordinates
        new_start = find_next_stop_codon(sequence, max(start - 3, 0), strand)
        if new_start is not None:
            start = new_start  # Include the stop codon

    logging.info(f"Adjusted sequence from {start} to {end} with length {end - start}")
    return start, end, strand

def read_to_pandas(file):
    logging.info(f"Reading annotation file {file}")
    df = pd.read_csv(file, sep='\t', header=None)
    df.columns = ["id", "strand", "start", "end", "rvs"]
    return df

def extract_and_translate_cds(fasta_file, annotation_df):
    logging.info(f"Extracting and translating CDS from FASTA file {fasta_file}")
    # Load all contigs into a dict for multi-contig support
    fasta_records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    # Sort contigs by length (longest first) for coordinate-based matching
    sorted_contigs = sorted(fasta_records.values(), key=lambda r: len(r.seq), reverse=True)
    is_multi_contig = len(fasta_records) > 1
    if is_multi_contig:
        logging.warning(f"Multi-contig assembly: {len(fasta_records)} contigs — matching genes to contigs by coordinate")
    default_contig = sorted_contigs[0]
    cds_amino_acids = {}
    cds_nts = {}

    indices_to_drop = []
    for index, row in annotation_df.iterrows():
        start = row['start'] - 1  
        end = row['end']
        strand = str(row['strand']).strip()
        ide = row['id']
        rvs = row['rvs']

        # Skip genes with invalid strand
        if strand not in ('+', '-'):
            logging.warning(f"Skipping {ide}: invalid strand '{strand}'")
            indices_to_drop.append(index)
            continue

        # Select the correct contig for this gene
        if is_multi_contig:
            fasta_seq = default_contig
            for contig in sorted_contigs:
                if end <= len(contig.seq):
                    fasta_seq = contig
        else:
            fasta_seq = default_contig
        
        # Validate coordinates
        if start < 0 or end <= start or end > len(fasta_seq.seq):
            logging.warning(f"Skipping {ide}: invalid coordinates {start+1}-{end}")
            indices_to_drop.append(index)
            continue

        sequence = fasta_seq.seq[start:end]
        
        if strand == '-':
            sequence = sequence.reverse_complement()

        logging.info(f"Translating sequence for ide {ide}")

        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                translated_seq = sequence.translate()
                for warning in w:
                    logging.warning(f"Warning for ide {ide}: {warning.message}")
                    pass

        except Exception as e:

            logging.error(f"Error translating sequence for ide {ide}: {e}")
            continue

        # Track original coordinates for fallback
        orig_start = start
        orig_end = end

        if len(translated_seq) > 0 and translated_seq[-1] == "*":
            # Gene already ends with a stop codon — good as is
            cds_id = f"{ide};{strand};{start}-{end}"
            cds_amino_acids[cds_id] = translated_seq
            cds_nts[cds_id] = sequence
        elif len(translated_seq) > 0:
            # No terminal stop codon — try to extend to find one
            # Check if this is a reference gene whose coordinates match the ID-embedded
            # coordinates. If so, it's a known pseudogene — don't extend.
            # If coords DON'T match the ID, it's a partial alignment that needs extension.
            is_ref_pseudogene = False
            ref_id = str(rvs) if pd.notna(rvs) and str(rvs).startswith('Rv') else str(ide)
            m = re.search(r'_(\d+)-(\d+)$', ref_id)
            if m:
                id_start = int(m.group(1))
                id_end = int(m.group(2))
                # Coords match the ID → this is a reference gene with its true coordinates
                if abs((start + 1) - id_start) <= 3 and abs(end - id_end) <= 3:
                    is_ref_pseudogene = True
            
            if is_ref_pseudogene:
                # Reference gene without stop = known pseudogene. Keep original coords.
                logging.info(f"Reference gene {ide} has no stop codon (likely pseudogene), keeping original coordinates")
                cds_id = f"{ide};{strand};{start}-{end}"
                cds_amino_acids[cds_id] = translated_seq
                cds_nts[cds_id] = sequence
            else:
                new_start, new_end, strand = get_goodseq(fasta_seq, start, end, strand)
                
                # Validate that we got valid coordinates
                if new_start is not None and new_end is not None and new_end > new_start:
                    start, end = new_start, new_end
                    if strand == "-":
                        sequence_nt = fasta_seq.seq[start:end].reverse_complement()
                    else:
                        sequence_nt = fasta_seq.seq[start:end]

                    new_translated = sequence_nt.translate()
                    
                    # Check that extension actually added a stop and didn't introduce internal stops
                    if len(new_translated) > 0 and new_translated[-1] == "*" and "*" not in str(new_translated)[:-1]:
                        cds_id = f"{ide};{strand};{start}-{end}"
                        cds_amino_acids[cds_id] = new_translated
                        cds_nts[cds_id] = sequence_nt
                    else:
                        # Extension failed — keep original coordinates
                        logging.warning(f"Extension for {ide} produced errors, keeping original")
                        start, end = orig_start, orig_end
                        if strand == "-":
                            sequence_nt = fasta_seq.seq[start:end].reverse_complement()
                        else:
                            sequence_nt = fasta_seq.seq[start:end]
                        cds_id = f"{ide};{strand};{start}-{end}"
                        cds_amino_acids[cds_id] = sequence_nt.translate()
                        cds_nts[cds_id] = sequence_nt
                else:
                    logging.warning(f"Could not find valid stop codon for {ide}, keeping original")
                    start, end = orig_start, orig_end
                    cds_id = f"{ide};{strand};{start}-{end}"
                    cds_amino_acids[cds_id] = translated_seq
                    cds_nts[cds_id] = sequence
        else:
            logging.warning(f"Empty translation for {ide}, skipping")
            continue

        # Update annotation coordinates (they may have changed from stop extension)
        annotation_df.at[index, 'start'] = start + 1  
        annotation_df.at[index, 'end'] = end

    # Drop invalid rows collected during iteration (avoids mutation during iterrows)
    annotation_df.drop(indices_to_drop, inplace=True)
    logging.info(f"CDS extraction and translation complete")
    return cds_amino_acids, cds_nts, annotation_df

def main():
    parser = argparse.ArgumentParser(description="Extract and translate CDS from a FASTA file to amino acids")
    parser.add_argument("fasta_file", help="FASTA file containing the sequence")
    parser.add_argument("annotation_file", help="Annotation file in tabulated format")
    args = parser.parse_args()
    base_name = args.annotation_file.replace('_corrected.tsv', '')
    annotation_df = read_to_pandas(args.annotation_file)

    annotation_df = annotation_df[annotation_df['end'] != 0]
    
    # Remove duplicate (start, end) entries BEFORE translation to avoid wasted work
    n_before = len(annotation_df)
    annotation_df = annotation_df.drop_duplicates(subset=['start', 'end'], keep='first')
    n_after = len(annotation_df)
    if n_before != n_after:
        logging.info(f"Removed {n_before - n_after} duplicate coordinate entries before translation")
    
    logging.info("Annotation DataFrame loaded and filtered")
    cds_amino_acids, cds_nts, annotation_df = extract_and_translate_cds(args.fasta_file, annotation_df)
    
    annotation_df['rvs'] = annotation_df['rvs'].fillna('')

    # Find duplicates based on 'start' and 'end' columns
    duplicados = annotation_df[annotation_df[['start', 'end']].duplicated(keep=False)]
    logging.info("Found duplicates")
    logging.info(duplicados)

    # Keep only one entry per duplicated (start, end) pair — remove all but the first
    if not duplicados.empty:
        indices_para_mantener = duplicados.groupby(['start', 'end']).apply(lambda x: x.index[0])
        indices_para_eliminar = duplicados.index.difference(indices_para_mantener)
        logging.info(f"Keeping {len(indices_para_mantener)} unique entries, removing {len(indices_para_eliminar)} duplicates")

        annotation_df = annotation_df.drop(indices_para_eliminar)
    
    logging.info(f"Saving filtered annotation DataFrame to {base_name + '.AMAP.tsv'}")
    annotation_df.to_csv(base_name + '.AMAP.tsv', sep='\t', index=False)

    logging.info(f"Saving translated sequences to {base_name + '.AMAP.faa'} and {base_name + '.AMAP.fna'}")
    with open(base_name + '.AMAP.faa', 'w') as file:
        for id, amino_acid_seq in cds_amino_acids.items():
            # Strip trailing stop codon '*' for tool compatibility (HMMER, InterProScan)
            seq_str = str(amino_acid_seq).rstrip('*')
            file.write(f">{id}\n{seq_str}\n")
    with open(base_name + '.AMAP.fna', 'w') as file:
        for id, amino_acid_seq in cds_nts.items():
            file.write(f">{id}\n{amino_acid_seq}\n")

if __name__ == "__main__":
    main()
