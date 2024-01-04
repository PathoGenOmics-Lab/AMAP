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
