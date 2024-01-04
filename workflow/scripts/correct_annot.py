#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbitblastnCommandline
import subprocess
from tempfile import NamedTemporaryFile
from io import StringIO
import re
import argparse
import glob
from functools import partial
import os
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from typing import Tuple, Dict
from concurrent.futures import ThreadPoolExecutor, as_completed

def run_blast(in_query: str, in_subject: str, outfmt: str) -> str:
    '''Runs blastn with the given parameters and returns the output.'''
    cline = NcbitblastnCommandline(query=in_query, db=in_subject, outfmt=outfmt, seg="no")
    out, _ = cline()
    return out

def run_blast_tab(in_query: str, in_subject: str): 
    '''Runs blastn with the given parameters and returns the output in tabular format.'''
    return run_blast(in_query, in_subject, '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen')

def run_blast_xml(in_query: str, in_subject: str): 
    '''Runs blastn with the given parameters and returns the output in xml format.'''
    return run_blast(in_query, in_subject, '5')

def parse_tab_data(blastn_out_tab: str) -> list:
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs', 'qcovhsp', 'qlen', 'slen']
    return [dict(zip(columns, line.split("\t"))) for line in blastn_out_tab.strip().split("\n") if line]

def parse_xml_data(blastn_out_xml: str) -> dict:
    return {record.query: record for record in NCBIXML.parse(StringIO(blastn_out_xml)) if blastn_out_xml.strip()}

def blastn(in_query: str, in_subject: str) -> pd.DataFrame:
    '''Runs blastn with the given parameters and returns the output in a pandas DataFrame.'''
    
    with ThreadPoolExecutor(max_workers=2) as executor:
        future_tab = executor.submit(run_blast_tab, in_query, in_subject)
        future_xml = executor.submit(run_blast_xml, in_query, in_subject)
        
        blastn_out_tab = future_tab.result()
        blastn_out_xml = future_xml.result()

    tab_data = parse_tab_data(blastn_out_tab)
    blast_records = parse_xml_data(blastn_out_xml)

    combined_data = []
    for tab in tab_data:
        xml_record = blast_records.get(tab['qseqid'])
        if xml_record and any(tab['sseqid'] in alignment.hit_def for alignment in xml_record.alignments):
            combined_data.append({**tab,
                                  "xml_data": xml_record,
                                  "subject_name": xml_record.descriptions[0].title.split()[1] if xml_record.descriptions else None})
    
    return pd.DataFrame(combined_data)