from lib2to3.pytree import BasePattern
import anndata
import numpy as np
from scipy.sparse import csr_matrix
import glob
import ntpath
import scanpy as sc
import pandas as pd
import os

'''
Find the files, rename them and 
1. With BC file remove the extra columns and the top row as there should not be any col names
2. With Genes file replace , with tabs, remove the top colum, take the first 2 columns
'''

import sys
files=sys.argv[0]
#files = "/media/data/AtteR/scifi-analysis/parsebio/test_data/donors"

samples_list = [dir for dir in sorted(glob.glob(files +"/*"))] 
samples_list
import subprocess

def bash_command(cmd):
    subprocess.Popen(cmd, shell=True, executable='/bin/bash')

for f in samples_list:
    #print(f)
    if "DGE.mtx.gz" in f:
        orig = os.path.join(files + "/"+ ntpath.basename(f))
        print(orig)
        new = os.path.join(files + "/" + "matrix.mtx.gz")
        os.rename(orig, new)
    if "cell_metadata.csv" in f:
        os.chdir(os.path.dirname(f))
        print(os.path.dirname(f))
        bash_command("ls")
        bash_command("cat cell_metadata.csv | tail -n +2 | cut -d ',' -f1 > barcodes.tsv.gz")
    if "genes.csv" in f:
        bash_command("cat genes.csv | tail -n +2 | cut -d ',' -f1,2 | sed 's:,:\t:g' > features.tsv.gz")
