from email.mime import base
from turtle import color
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
import anndata as ad
import os
import plotly 
import plotly.graph_objs as go   #numba needs a numpy of 1.21 or less 
import pickle
import glob
import ntpath
import seaborn as sns
import scipy
import scipy.io as scio
import re
import multiprocessing as mp
import sys
import time
from multiprocessing import Process, Pool, Queue, Manager
import gzip
import random
from functools import reduce
import csv 
import tempfile
np1 = np.arange(1,10)
np1
np2 = np.array([[1,4,3], [4,10,11], [1,3,2,15,6]], dtype=object)
np2


outer_product = np.outer(np1, np2)
outer_product
'''
Vectorisation approach - take the seq ids as a np 1 and the clusters with seq ids as np 2. 
np1 = 1d np array
np2 = 2d np array

#each element of this np array corresponds to the seq id when counting from 0 to N
[[1,4,3], [4,10,11] [1,3,2,15,6]]

We go over each of them and multiply by the element of np1 and then 

seq id: 4

So products when each element of np2 are calculated are:
4
16
12


Pool and map() arent suited to situations that need to maintain state over time or, especially, situations where there are 
two or more different operations that need to run and interact with each other in some way. For these kinds of problems, one
needs to make use of somewhat more complex features of the multiprocessing module, such as Process, Queue and Event. Using 
these features introduces a number of complicated issues that now need to be managed, especially with regard to cleanly 
starting and stopping subprocesses, as well as coordinating between them.

The Queue class includes the put method for depositing data and the get method for retrieving data

“A process pool object which controls a pool of worker processes to which jobs can be submitted. It supports asynchronous 
results with timeouts and callbacks and has a parallel map implementation”.
Using Pool we can assign as many parallel processes as we like, but only the `processes` number of threads will be 
active at any given moment.
'''

basePath = 'SoloNova10x_correctAAV/'
basePath= "/media/data/AtteR/scifi-analysis/scifi6/2nd_try/output_enriched/"

r21_files = sorted(glob.glob(basePath+'*R21.fastq.gz'))
r21_files= [item for item in r21_files if 'oDT' not in item]  
r21_files= [item for item in r21_files if 'SciFi_Pool1-2' in item]  

r21_files
r21_sc_files = {}

for f in r21_files:
    #print(f)
    r21_sc_files[f]=f.replace("R21.fastq.gz", "starcode.tsv")
for a,b in r21_sc_files.items():
    print(a)
#tmp_loc.cleanup()

for r21, sc_f in r21_sc_files.items():
    #print(BBC_path)
    #print(starcode_outp)
    if os.path.exists(r21) == False:
        continue
    if os.path.exists(sc_f) == False:
        continue
    #(r21 + "---" + sc_f)
    tmp_loc = tempfile.TemporaryDirectory(dir = "/media/data/AtteR/scifi-analysis/scifi6/2nd_try")
    tmp_dir = str(tmp_loc).split("'", 1)[1].split("'")[0]

    #Number the seqs from fastq file, make them as dict keys with empty values. Then go over the starcode file, saving cluster as key and 
    #the locations of seqs belonging to the cluster into the value. go over each each line in fastq
    #file and compare the results to the values in the dict. If a match is found, add a count to
    #the value part of the dict
    sc_f
    BC_clusters = pd.read_csv(sc_f, sep='\t', header=None, index_col=0) #dtype={0: str, 1: 'Int64'})
    clusters = list(BC_clusters.index)
    print(BC_clusters)
    print("Cluster found")
    #Get the BC cluster name and seqs 
    BC_clusters_dict = {}
    for x in range(len(BC_clusters.index)):
        BC_clusters.iloc[0,:1]
        BC_clusters_dict[clusters[x]]=BC_clusters.iloc[x,1]
    ###########

    if len(BC_clusters_dict)<=2:
        print("Cluster length too short, skipping!")
        continue
    #get the 21 seq file, count the number of seqs and save as dict 
    with open(r21, 'rb') as fastq:
        gzip_fastq = gzip.GzipFile(fileobj=fastq)
        fastq = gzip_fastq.read().splitlines() 
        print("number of sequences: " + str(len(fastq)/4))
        uniq_bcs_cell = {}
        for n in range(int(len(fastq)/4)):
            uniq_bcs_cell[n] = 0  #number of seqs as key, as value we will save the count of unique BCs
    seqs = list(uniq_bcs_cell.keys())


    nprocs = mp.cpu_count()
    p_size = nprocs - 5
    #############

    def progressbar(it, prefix="", size=100, file=sys.stdout):
        count = len(it)
        try:
            def show(j):
                x = int(size*j/count)
                file.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, count))
                file.flush()        
            show(0)
            for i, item in enumerate(it):
                yield item
                show(i+1)
            file.write("\n")
            file.flush()
        except ZeroDivisionError:
            print("Division by zero. The range length was " + str(count))


    #since we are not importing the entire seq range as a whole into the function below (to risk of too high I/O stream), 
    #we will name the files based on the range instead and save into tmp dir and afterwards merge them all into same 
    def unique_bcs_per_seq(BC_clusters, seq_no, filepath, q): #take in the contents of starcode file as a dict (BC_clusters) and the seq number
        uniqe_bcs_cell = {}
        fname = os.path.basename(filepath)
        f_id = fname.split(".")[0]
        save_path = os.path.dirname(filepath)
        #print("splitting the " + str(seq_no))
        m = str(seq_no).split('(', 1)[1].split(')')[0]
        #print(type(m))
        a = m.split(",")[0]
        b = m.split(" ")[1]
        file_sub_id = "_r_" + a + "-" + b
        for seq in progressbar(seq_no, "Computing: ", 40):  #when we split the numbers, we put them into ranges of certain size so seq_no is now in range so need to iterate over it
        #for seq in seq_no:  #when we split the numbers, we put them into ranges of certain size so seq_no is now in range so need to iterate over it
            counts=0
            for key, value in BC_clusters.items():
                #print("entered cluster")
                cluster_seq_ids = value.split(",")
                if str(seq) in cluster_seq_ids:
                    counts+=1
                    uniqe_bcs_cell[seq]=counts
            #print(str(seq) + " count: " + str(counts))
                #q.put(uniqe_bcs_cell)   #q object to save the result
            #print("Processed Seq " + str(seq) + ". Unique BCs found: " + str(counts))
        final_df = pd.DataFrame(list(uniqe_bcs_cell.items()),columns = ['Cell','N(Unique BCs/cell)']) 
        #final_df=final_df.T
        #print(final_df.head())
        with open(tmp_dir + "/" +  fname + file_sub_id +".csv", "w") as out: #write the dataframe into an output file
            #for key in uniqe_bcs_cell.keys():
            #    out.write("%s,%s\n"%(key,uniqe_bcs_cell[key]))
            final_df.to_csv(out)
            # # df.to_string(out, index=None)
        print("Saved as " + tmp_dir + "/" + fname + file_sub_id + ".csv")
        return(final_df)
    #############

def reader(i,q):
    message = q.get() #retrieves the queue whereas the writer function puts its output into the queue
    print(message)
def slice_data(data, nprocs):
    aver, res = divmod(len(data), nprocs)
    nums = []
    for proc in range(nprocs):
        if proc < res:
            nums.append(aver + 1)
        else:
            nums.append(aver)
    count = 0
    slices = []
    for proc in range(nprocs):
        slices.append(data[count: count+nums[proc]])
        count += nums[proc]
    return slices

#I can either run this script on the file prior to demultiplexing. This means that when I demultiplex the files, I must index
#the lines... So instead run this script on the demultiplexed samples (so first starcode, then this), then map and then import into anndata

#I must iterate over the files and get the R21 read and then the corresponding starcode tsv

#multiprocess the files
p = Pool(p_size)
m = Manager()

q = m.Queue()
seq_lists = slice_data(range(len(seqs)), nprocs)
print(seq_lists)
#multi_result = [pool.apply_async(power_n_list, (inp, 2)) for inp in inp_lists] uniq_bcs_cell
multi_result = [p.apply_async(unique_bcs_per_seq, (BC_clusters_dict, seq,r21,q,)) for seq in seq_lists]
final_result = [result.get() for result in multi_result]
'''
#final_result
vect_func = np.vectorize(unique_bcs_per_seq)
for seq in seq_lists
    vect_func(BC_clusters_dict, seq)
'''
#get the files from the tmp dir
count_files = [i for i in glob.glob(tmp_dir +"/*.csv")]
def common_name(sa, sb):
    """ returns the longest common substring from the beginning of sa and sb """
    def _iter():
        for a, b in zip(sa, sb):
            if a == b:
                yield a
            else:
                return

    return ''.join(_iter())
shared_path = common_name(count_files[0], count_files[1])

shared_name = '_'.join(shared_path.split("/")[-1].split("_")[:-3])
#combine all files in the list
combined_csv = pd.concat([pd.read_csv(f) for f in count_files ])
#export to csv
combined_csv.to_csv(basePath + "/" + shared_name + "_Uniq_BCs_per_cell.csv", index=False)
df_counts_full = pd.read_csv(basePath + "/" + shared_name + "_Uniq_BCs_per_cell.csv")
df_counts_full =df_counts_full.sort_values(by = "Cell", ascending=True)
del df_counts_full["Unnamed: 0"]
#df_counts_full
df_counts_full.to_csv(basePath + "/" + shared_name + "_Uniq_BCs_per_cell.csv", index=False)
print("Done")
tmp_loc.cleanup()
