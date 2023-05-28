#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from concurrent.futures import ThreadPoolExecutor
import subprocess
from tempfile import NamedTemporaryFile
from io import StringIO
import os
from collections import defaultdict

def get_seq(chck_rds: list, fasta_file: str, fasta_name: str) -> List[Tuple[str, int]]:
    chck_rds_set = set(chck_rds)  # Convert the list to a set for faster searching
    selected_seqs = [seq_record for seq_record in SeqIO.parse(fasta_file, "fasta") if seq_record.id in chck_rds_set] # Use list comprehension for faster, more readable code
    SeqIO.write(selected_seqs, fasta_name, "fasta") # Write the selected sequences to a fasta file
    seq_lengths = [(seq_record.id, len(seq_record.seq)) for seq_record in selected_seqs] # Use list comprehension to create the list of sequence lengths
    return seq_lengths

def do_ranges_overlap(a, b):
    a_start, a_end = sorted(a)
    b_start, b_end = sorted(b)

    return (a_start < b_end) and (a_end > b_start)

def create_database(fasta: str, basename: str) -> None:
    subprocess.run(['makeblastdb', '-in', fasta, '-dbtype', 'nucl', '-out', basename])

def blastn(in_query: str, in_subject: str) -> pd.DataFrame:
    blastn_cline_tab = NcbiblastnCommandline(query=in_query, db=in_subject, outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen', dust="no")
    blastn_out_tab, _ = blastn_cline_tab()
    blastn_cline_xml = NcbiblastnCommandline(query=in_query, db=in_subject, outfmt='5', dust="no")
    blastn_out_xml, _ = blastn_cline_xml()

    tab_data = [
        dict(zip(['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs', 'qcovhsp', 'qlen', 'slen'], line.split("\t")))
        for line in blastn_out_tab.strip().split("\n")
    ] if len(blastn_out_tab.strip()) > 0 else []

    blastn_out_xml = StringIO(blastn_out_xml)
    if blastn_out_xml.getvalue().strip():
        blast_records = list(NCBIXML.parse(blastn_out_xml))
    else:
        blast_records = []

    combined_results = []

    for tab in tab_data:
        query_id = tab['qseqid']
        subject_id = tab['sseqid']
        xml_record = next((record for record in blast_records if record.query == query_id and any(subject_id in alignment.hit_def for alignment in record.alignments)), None)

        if xml_record is not None:
            combined_result = {
                "tab_data": tab,
                "xml_data": xml_record,
                "subject_name": xml_record.descriptions[0].title.split()[1] if len(xml_record.descriptions) > 0 else None
            }
            combined_results.append(combined_result)

    df = pd.DataFrame.from_records(
        [
            {
                **{"" + k: v for k, v in result["tab_data"].items()},
                "xml_data": result["xml_data"],
                "subject_name": result["subject_name"]
            }
            for result in combined_results
        ]
    )

    return df
