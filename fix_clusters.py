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
