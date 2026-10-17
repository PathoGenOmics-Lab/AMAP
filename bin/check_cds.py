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
__version__ = '1.0.1'
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

def get_goodseq(sequence, start, end, strand, max_distance=300):
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
        new_end = find_next_stop_codon(sequence, end, strand, max_distance)
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
        new_start = find_next_stop_codon(sequence, max(start - 3, 0), strand, max_distance)
        if new_start is not None:
            start = new_start  # Include the stop codon

    logging.info(f"Adjusted sequence from {start} to {end} with length {end - start}")
    return start, end, strand

def read_to_pandas(file):
    logging.info(f"Reading annotation file {file}")
    df = pd.read_csv(file, sep='\t', header=None)
    df.columns = ["id", "strand", "start", "end", "rvs"]
    return df

# Start codons of translation table 11, the bacterial code.
START_CODONS = {"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"}


def read_reference_proteins(path):
    """Reference proteins by name, with the length of the span each was cut from.

    Headers are NAME_START-END (Rv0001_1-1524), as aa_H37Rv.translated.fas and
    extract_ref_proteins.py write them, and NAME is what classify_cases.py puts
    in the rvs column. The span is None for a header without coordinates.
    """
    refs = {}
    for rec in SeqIO.parse(path, "fasta"):
        m = re.match(r'^(.+)_(\d+)-(\d+)$', rec.id)
        if m:
            name, span = m.group(1), int(m.group(3)) - int(m.group(2)) + 1
        else:
            name, span = rec.id, None
        refs.setdefault(name, (str(rec.seq).upper(), span))
    return refs


def reference_for(rvs, refs):
    """The (protein, span) a gene was matched to, or None.

    A gene matched to several references at once carries all of them in rvs;
    they are copies of one protein, so the first one known is as good as any.
    """
    if not refs or pd.isna(rvs):
        return None
    for name in str(rvs).split(','):
        if name.strip() in refs:
            return refs[name.strip()]
    return None


def lacks_stop(ref):
    """Whether the reference annotates this gene without a stop codon.

    Its protein then spans the whole annotated stretch codon for codon, with
    nothing left over for a stop (Rv2099c: 58 residues over 174 bp). A protein
    without a final '*' is not enough on its own: for a frameshifted (join) CDS
    extract_ref_proteins.py leaves the '*' off because the outer span is not
    3 x (residues + 1), stop or no stop (E. coli b0240: 133 residues, 1178 bp).
    """
    protein, span = ref
    return (not protein.endswith('*') and span is not None and span % 3 == 0
            and len(protein) == span // 3)


def ends_like(translated, ref_protein, n=10):
    """Whether a gene's last residues are the reference protein's last ones."""
    got, ref = str(translated).rstrip('*'), ref_protein.rstrip('*')
    n = min(n, len(got), len(ref))
    return n > 0 and sum(a == b for a, b in zip(got[-n:], ref[-n:])) >= 0.8 * n


def complete_n_terminus(fasta_seq, start, end, strand, ref_protein, translated):
    """Extend a gene whose protein lacks the reference's first residues.

    Miniprot can stop an alignment short of the N-terminus: on the CP110619
    assembly of H37Rv it aligned Rv3611 from residue 30, over sequence identical
    to AL123456.3 where it aligned the whole protein, and the gene came out 90 bp
    short. When a gene ends in a stop codon and is k residues shorter than its
    reference, the k codons upstream are tried: they must open with a start
    codon, contain no stop, and translate to the reference's first residues.

    A gene k residues short is not necessarily missing its first k: M. bovis
    Rv3159c (PPE53) is one residue short through an in-frame deletion and
    starts at the same ATG as H37Rv. So the gene must line up with the
    reference from residue k, and not already from residue 0, before anything
    is added; growing its 5' end would only put a CTG in front of a correct start.

    Coordinates are 0-based, end exclusive. Returns them changed or not.
    """
    ref = ref_protein.rstrip('*')
    got = str(translated).rstrip('*')
    k = len(ref) - len(got)
    if k <= 0 or not str(translated).endswith('*'):
        return start, end

    def identity(a, b):
        n = min(len(a), len(b))
        return sum(x == y for x, y in zip(a, b)) / n if n else 0.0

    # Residue 0 depends on which start codon a naive translation met, so the
    # line-up is judged from residue 1.
    w = min(20, len(got))
    if identity(got[1:w], ref[1:w]) >= 0.9:
        return start, end
    if identity(got[1:w], ref[k + 1:k + w]) < 0.9 or identity(got, ref[k:]) < 0.9:
        return start, end

    if strand == '+':
        new_start, new_end = start - 3 * k, end
        if new_start < 0:
            return start, end
        ext = fasta_seq.seq[new_start:start]
    else:
        new_start, new_end = start, end + 3 * k
        if new_end > len(fasta_seq.seq):
            return start, end
        ext = fasta_seq.seq[end:new_end].reverse_complement()

    ext_aa = str(ext.translate())
    if str(ext[:3]).upper() not in START_CODONS or '*' in ext_aa:
        return start, end
    # The reference proteome is a naive translation too, so the same start
    # codon gives the same first residue and every one of the k can be held to it.
    if identity(ext_aa, ref[:k]) < 0.9:
        return start, end
    logging.info(f"Restored {k} N-terminal residue(s) missing from the alignment")
    return new_start, new_end

def extract_and_translate_cds(fasta_file, annotation_df, ref_proteins=None):
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
    pseudogene_ids = set()  # Track genes without terminal stop / with internal stops

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

        # Select the correct contig for this gene.
        # For multi-contig: find the smallest contig where the gene fits.
        # Sorted desc by length, so we iterate and keep the LAST match
        # (= smallest fitting contig, tightest fit).
        # WARNING: This is ambiguous when contigs overlap in coordinate ranges.
        if is_multi_contig:
            fasta_seq = default_contig
            best_len = float('inf')
            for contig in sorted_contigs:
                clen = len(contig.seq)
                if start >= 0 and end <= clen and clen < best_len:
                    fasta_seq = contig
                    best_len = clen
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
            # Drop the row too: keeping it would emit a gene in the .AMAP.tsv
            # (and hence the master table and GFF) with no matching record in
            # the .AMAP.faa / .AMAP.fna.
            indices_to_drop.append(index)
            continue

        # Track original coordinates for fallback
        orig_start = start
        orig_end = end
        ref = reference_for(rvs, ref_proteins)

        if len(translated_seq) > 0 and translated_seq[-1] == "*" and ref is not None:
            new_start, new_end = complete_n_terminus(fasta_seq, start, end, strand,
                                                     ref[0], translated_seq)
            if (new_start, new_end) != (start, end):
                start, end = new_start, new_end
                sequence = fasta_seq.seq[start:end]
                if strand == '-':
                    sequence = sequence.reverse_complement()
                translated_seq = sequence.translate()

        if len(translated_seq) > 0 and translated_seq[-1] == "*":
            # Gene ends with a stop codon — check for internal stops
            cds_id = f"{ide};{strand};{start}-{end}"
            cds_amino_acids[cds_id] = translated_seq
            cds_nts[cds_id] = sequence
            if "*" in str(translated_seq)[:-1]:
                n_internal = str(translated_seq)[:-1].count("*")
                logging.warning(f"Gene {ide} has {n_internal} internal stop codon(s) — likely pseudogene or frameshift")
                pseudogene_ids.add(ide)
        elif len(translated_seq) > 0:
            # No terminal stop codon. Normally the gene is cut short and gets
            # extended to the next one, but not when the reference has no stop
            # either. This used to be decided from coordinates embedded in the
            # reference IDs, which clean_id has stripped from rvs since March
            # 2026, so it never fired. The reference proteome answers it
            # directly, for any organism.
            partial = (ref is not None and ref[1] and ref[1] % 3
                       and (end - start) == ref[1] - ref[1] % 3
                       # the bases go on the 3' end, which is the low end on '-'
                       and (end + ref[1] % 3 <= len(fasta_seq.seq) if strand == '+'
                            else start - ref[1] % 3 >= 0))

            if ref is not None and lacks_stop(ref) and ends_like(translated_seq, ref[0]):
                # A pseudogene the reference annotates without a stop codon
                # (PE21, Rv2099c; PPE48, Rv3022c). Running on to the next stop
                # would invent a C-terminus the reference does not have. Only
                # when the gene ends where the reference protein does, though:
                # a frameshifted reference pseudogene translated straight
                # through also lacks a final '*', because its real stop sits in
                # another frame, and the C-terminal fragment a sample carries
                # in that frame does have to run on to it (E. coli b0240).
                logging.info(f"Reference protein of {ide} has no stop codon either (reference pseudogene), keeping original coordinates")
                cds_id = f"{ide};{strand};{start}-{end}"
                cds_amino_acids[cds_id] = translated_seq
                cds_nts[cds_id] = sequence
                pseudogene_ids.add(ide)
            elif partial:
                # A frameshifted reference pseudogene annotated with a partial
                # final codon (oxyR', Rv2427A; Rv3216). An alignment only covers
                # whole codons, so put back the trailing bases the reference has.
                if strand == '+':
                    end += ref[1] % 3
                    sequence = fasta_seq.seq[start:end]
                else:
                    start -= ref[1] % 3
                    sequence = fasta_seq.seq[start:end].reverse_complement()
                logging.info(f"Reference pseudogene {ide} ends in a partial codon, restoring it")
                cds_id = f"{ide};{strand};{start}-{end}"
                cds_amino_acids[cds_id] = sequence[:len(sequence) - len(sequence) % 3].translate()
                cds_nts[cds_id] = sequence
                pseudogene_ids.add(ide)
            else:
                # 300 bp is enough for an alignment that stopped a few codons
                # short. When the reference says how many residues are missing,
                # look that far: miniprot left the last 128 residues of a
                # repetitive PE_PGRS (Rv1003720new) unaligned, 387 bp from its stop.
                reach = 300
                if ref is not None and ref[0].endswith('*'):
                    missing = len(ref[0].rstrip('*')) - len(str(translated_seq).rstrip('*'))
                    reach = max(reach, 3 * missing + 30)
                new_start, new_end, strand = get_goodseq(fasta_seq, start, end, strand, reach)
                
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
                        pseudogene_ids.add(ide)
                else:
                    logging.warning(f"Could not find valid stop codon for {ide}, keeping original")
                    start, end = orig_start, orig_end
                    cds_id = f"{ide};{strand};{start}-{end}"
                    cds_amino_acids[cds_id] = translated_seq
                    cds_nts[cds_id] = sequence
                    pseudogene_ids.add(ide)
        else:
            logging.warning(f"Empty translation for {ide}, skipping")
            # Same as above: no sequence was stored, so the row must not
            # survive into the annotation table.
            indices_to_drop.append(index)
            continue

        # Update annotation coordinates (they may have changed from stop extension)
        annotation_df.at[index, 'start'] = start + 1  
        annotation_df.at[index, 'end'] = end

    # Drop invalid rows collected during iteration (avoids mutation during iterrows)
    annotation_df.drop(indices_to_drop, inplace=True)
    # Mark pseudogenes in annotation
    annotation_df['pseudogene'] = annotation_df['id'].apply(
        lambda x: 'yes' if str(x).split(';')[0] in pseudogene_ids else 'no'
    )

    logging.info(f"CDS extraction and translation complete ({len(pseudogene_ids)} pseudogenes detected)")
    return cds_amino_acids, cds_nts, annotation_df

def main():
    parser = argparse.ArgumentParser(description="Extract and translate CDS from a FASTA file to amino acids")
    parser.add_argument("fasta_file", help="FASTA file containing the sequence")
    parser.add_argument("annotation_file", help="Annotation file in tabulated format")
    parser.add_argument("--ref-proteome", default=None,
                        help="Reference proteome (NAME_START-END headers). Lets reference "
                             "pseudogenes keep their extent and restores N-termini an "
                             "alignment clipped. Without it, every gene missing a stop "
                             "codon is extended to the next one.")
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
    ref_proteins = read_reference_proteins(args.ref_proteome) if args.ref_proteome else None
    cds_amino_acids, cds_nts, annotation_df = extract_and_translate_cds(args.fasta_file, annotation_df,
                                                                        ref_proteins)
    
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

    # The sequence dictionaries were filled while iterating, before the rows
    # above were dropped, and the FASTA files below are written straight from
    # them. Without this they keep proteins whose row no longer exists, so the
    # .faa ends up larger than the .tsv and the two disagree about which genes
    # the sample has. Rebuild the keys the same way extract_and_translate_cds
    # does: it stores the 0-based start, the table carries it 1-based.
    surviving = {
        f"{row['id']};{row['strand']};{int(row['start']) - 1}-{int(row['end'])}"
        for _, row in annotation_df.iterrows()
    }
    n_orphans = len(cds_amino_acids) - len(surviving & set(cds_amino_acids))
    if n_orphans:
        logging.info(f"Dropping {n_orphans} sequence(s) whose annotation row was removed")
    cds_amino_acids = {k: v for k, v in cds_amino_acids.items() if k in surviving}
    cds_nts = {k: v for k, v in cds_nts.items() if k in surviving}

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
