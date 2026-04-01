#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script converts a GenBank file (GBK) into a GFF3 file formatted
specifically for compatibility with Panaroo, mimicking the Prokka GFF3 style.
It processes 'gene' and 'CDS' features, creating appropriate attributes
and parent-child relationships, and appends the FASTA sequence to the end.
"""

import argparse
from pathlib import Path
from typing import Dict, List, Any
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

__author__ = 'Paula Ruiz-Rodriguez'
__credits__ = ['Paula Ruiz-Rodriguez', 'PathoGenOmics']
__license__ = ['GPL 3']
__version__ = '0.1.0'
__maintainer__ = '@paururo'
__email__ = 'paula.ruiz.rodriguez@csic.es'
__status__ = 'finished'
__lab__ = 'PathoGenOmics (I2SysBio)'

def format_gff_attributes(feature: SeqFeature) -> str:
    """
    Formats the GFF3 attributes column for a given feature.

    Args:
        feature: A Biopython SeqFeature object.

    Returns:
        A semicolon-separated string of GFF3 attributes.
    """
    attributes = []
    locus_tag = feature.qualifiers.get("locus_tag", [None])[0]

    if not locus_tag:
        return ""

    if feature.type == 'gene':
        attributes.append(f"ID={locus_tag}_gene")
        attributes.append(f"Name={locus_tag}")
        attributes.append(f"gene={locus_tag}")
    
    elif feature.type == 'CDS':
        attributes.append(f"ID={locus_tag}")
        attributes.append(f"Parent={locus_tag}_gene")
        if "product" in feature.qualifiers:
            attributes.append(f"product={feature.qualifiers['product'][0]}")
        if "db_xref" in feature.qualifiers:
            # Join multiple db_xrefs if they exist
            db_xrefs = ",".join(feature.qualifiers['db_xref'])
            attributes.append(f"db_xref={db_xrefs}")
        if "inference" in feature.qualifiers:
            attributes.append(f"inference={feature.qualifiers['inference'][0]}")

    return ";".join(attributes)

def convert_gbk_to_panaroo_gff3(gbk_path: Path, gff3_path: Path):
    """
    Reads a GenBank file and writes a Panaroo-compatible GFF3 file.

    Args:
        gbk_path: The path to the input GenBank file.
        gff3_path: The path to the output GFF3 file.
    """
    try:
        gbk_records = list(SeqIO.parse(gbk_path, "genbank"))
    except FileNotFoundError:
        print(f"Error: Input file not found at {gbk_path}")
        return
    except Exception as e:
        print(f"An error occurred while reading the GenBank file: {e}")
        return

    if not gbk_records:
        print(f"Error: No records found in {gbk_path}")
        return

    with open(gff3_path, "w") as gff3_file:
        # Write GFF3 header
        gff3_file.write("##gff-version 3\n")
        for gbk_record in gbk_records:
            gff3_file.write(f"##sequence-region {gbk_record.id} 1 {len(gbk_record.seq)}\n")

        # Process and write features for all records/contigs
        # Panaroo requires both gene and CDS features with proper Parent references.
        # tsv_to_gbk only creates CDS features, so we emit a gene feature for
        # every CDS to satisfy the Panaroo GFF3 requirement.
        for gbk_record in gbk_records:
            for feature in gbk_record.features:
                if feature.type not in ["gene", "CDS"]:
                    continue

                start = feature.location.start + 1
                end = feature.location.end
                strand = '+' if feature.strand == 1 else '-'
                score = '.'
                source = "AMAP"

                locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                if not locus_tag:
                    continue

                if feature.type == 'CDS':
                    # Emit a parent gene feature first (required by Panaroo)
                    gene_attrs = f"ID={locus_tag}_gene;Name={locus_tag};gene={locus_tag}"
                    gff3_file.write(
                        f"{gbk_record.id}\t{source}\tgene\t"
                        f"{start}\t{end}\t{score}\t{strand}\t.\t{gene_attrs}\n"
                    )
                    # Then emit the CDS with Parent reference
                    cds_attrs = format_gff_attributes(feature)
                    gff3_file.write(
                        f"{gbk_record.id}\t{source}\tCDS\t"
                        f"{start}\t{end}\t{score}\t{strand}\t0\t{cds_attrs}\n"
                    )
                elif feature.type == 'gene':
                    # Gene feature already exists in GBK — emit as-is
                    attribute_str = format_gff_attributes(feature)
                    if attribute_str:
                        gff3_file.write(
                            f"{gbk_record.id}\t{source}\tgene\t"
                            f"{start}\t{end}\t{score}\t{strand}\t.\t{attribute_str}\n"
                        )

        # Append FASTA sequences for all records at the end
        gff3_file.write("##FASTA\n")
        for gbk_record in gbk_records:
            gff3_file.write(f">{gbk_record.id}\n")
            for i in range(0, len(gbk_record.seq), 60):
                gff3_file.write(str(gbk_record.seq[i:i+60]) + "\n")

    print(f"Successfully converted {gbk_path.name} to {gff3_path.name}")

def main():
    """
    Main function to parse command-line arguments and initiate the conversion.
    """
    parser = argparse.ArgumentParser(
        description='Convert a GenBank (.gbk) file to a Panaroo-compatible GFF3 file.'
    )
    parser.add_argument(
        '--id',
        type=str,
        required=True,
        help='The base ID for the input and output filenames (e.g., "A1_G48905").'
    )
    args = parser.parse_args()
    
    # Construct file paths using pathlib for better cross-platform compatibility
    gbk_filename = Path(f"{args.id}-AMAP.gbk")
    gff3_filename = Path(f"{args.id}.panaroo.gff")

    convert_gbk_to_panaroo_gff3(gbk_filename, gff3_filename)

if __name__ == "__main__":
    main()
