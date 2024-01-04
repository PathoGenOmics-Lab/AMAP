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