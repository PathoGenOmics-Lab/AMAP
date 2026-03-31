#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import logging
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

def find_next_stop_codon(fasta_seq, start_index, strand):
    """Find the next in-frame stop codon in the DNA sequence.
    
    For '+' strand: scans forward from start_index in steps of 3.
    For '-' strand: scans backward (toward lower coordinates) from start_index
                    in steps of 3, reading codons on the reverse complement.
    
    Returns the 0-based position of the first base of the stop codon, or None.
    """
    sequence = fasta_seq.seq
    stop_codons = {"TAA", "TAG", "TGA"}

    logging.info(f"Finding next stop codon from index {start_index} on strand {strand}")

    if strand == "+":
        for i in range(start_index, len(sequence) - 2, 3):
            codon = str(sequence[i:i + 3])
            if codon in stop_codons:
                logging.info(f"Found stop codon '{codon}' at position {i}")
                return i
    else:
        # On minus strand, codons are read as reverse complement.
        # We scan toward lower coordinates.
        for i in range(start_index, 1, -3):
            if i + 3 > len(sequence):
                continue
            codon = str(sequence[i:i + 3].reverse_complement())
            if codon in stop_codons:
                logging.info(f"Found stop codon '{codon}' (rc) at position {i}")
                return i

    logging.warning("No stop codon found")
    return None

def get_goodseq(sequence, start, end, strand):
    """Extend a CDS to include its stop codon, adjusting frame if needed.
    
    Args:
        sequence: Biopython SeqRecord of the genome.
        start: 0-based start coordinate (always < end).
        end: 0-based end coordinate (exclusive).
        strand: '+' or '-'.
    
    Returns:
        Tuple of (start, end, strand) with adjusted coordinates.
    """
    logging.info(f"Getting good sequence from {start} to {end} on strand {strand}")

    if strand == "+":
        # Ensure length is divisible by 3
        resto = (end - start) % 3
        if resto != 0:
            end = end - resto
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
        # Search for stop codon scanning toward lower coordinates
        new_start = find_next_stop_codon(sequence, start - 3, strand)
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
    if len(fasta_records) == 1:
        fasta_seq = next(iter(fasta_records.values()))
    else:
        # For single-contig genomes (most MTBC), use the first
        fasta_seq = next(iter(fasta_records.values()))
        logging.info(f"Multi-contig assembly: {len(fasta_records)} contigs, using first for CDS extraction")
    cds_amino_acids = {}
    cds_nts = {}

    for index, row in annotation_df.iterrows():
        start = row['start'] - 1  
        end = row['end']
        strand = row['strand']
        ide = row['id']
        rvs = row['rvs']

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

        if len(translated_seq) > 0 and translated_seq[-1] == "*":
            cds_id = f"{ide};{strand};{start}-{end}"
            cds_amino_acids[cds_id] = translated_seq
            cds_nts[cds_id] = sequence
        else:
            if len(translated_seq) > 0:
                # Use get_goodseq to find the stop codon — it handles both strands
                new_start, new_end, strand = get_goodseq(fasta_seq, start, end, strand)
                
                # Validate that we got valid coordinates
                if new_start is not None and new_end is not None and new_end > new_start:
                    start, end = new_start, new_end
                    if strand == "-":
                        sequence_nt = fasta_seq.seq[start:end].reverse_complement()
                    else:
                        sequence_nt = fasta_seq.seq[start:end]

                    cds_id = f"{ide};{strand};{start}-{end}"
                    cds_amino_acids[cds_id] = sequence_nt.translate()
                    cds_nts[cds_id] = sequence_nt
                else:
                    logging.warning(f"Could not find valid stop codon for {ide}, skipping")
                    continue

            annotation_df.at[index, 'start'] = start + 1  
            annotation_df.at[index, 'end'] = end   

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
            file.write(f">{id}\n{amino_acid_seq}\n")
    with open(base_name + '.AMAP.fna', 'w') as file:
        for id, amino_acid_seq in cds_nts.items():
            file.write(f">{id}\n{amino_acid_seq}\n")

if __name__ == "__main__":
    main()
