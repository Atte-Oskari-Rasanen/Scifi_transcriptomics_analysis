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

from multiprocessing import Process, Pool, Queue, Manager
import gzip
import random
from functools import reduce
import csv 
import tempfile


'''
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

r21_sc_files = {}

for f in r21_files:
    print(f)
    r21_sc_files[f]=f.replace("R21.fastq.gz", "starcode.tsv")
r21_sc_files

#tmp_loc.cleanup()
r21_sc_files
for r21, sc_f in r21_sc_files.items():
    BBC_path = r21
    starcode_outp=sc_f
    #print(BBC_path)
    #print(starcode_outp)
    if os.path.exists(BBC_path) == False:
        continue
    if os.path.exists(starcode_outp) == False:
        continue
    print(BBC_path)
    tmp_loc = tempfile.TemporaryDirectory(dir = "/media/data/AtteR/scifi-analysis/scifi6/2nd_try")
    #tmp_loc
    tmp_dir = str(tmp_loc).split("'", 1)[1].split("'")[0]
    tmp_dir


    #Number the seqs from fastq file, make them as dict keys with empty values. Then go over the starcode file, saving cluster as key and 
    #the locations of seqs belonging to the cluster into the value. go over each each line in fastq
    #file and compare the results to the values in the dict. If a match is found, add a count to
    #the value part of the dict

    BC_clusters = pd.read_csv(starcode_outp, sep='\t', header=None, index_col=0, dtype={0: str, 1: 'Int64'})
    #len(BC_clusters.index)
    clusters = list(BC_clusters.index)
    #clusters


    #Get the BC cluster name and seqs 
    BC_clusters_dict = {}
    for x in range(len(BC_clusters.index)):
        BC_clusters.iloc[0,:1]
        BC_clusters_dict[clusters[x]]=BC_clusters.iloc[x,1]
    ###########

    #get the 21 seq file, count the number of seqs and save as dict 
    with open(BBC_path, 'rb') as fastq:
        gzip_fastq = gzip.GzipFile(fileobj=fastq)
        fastq = gzip_fastq.read().splitlines() 
        type(fastq)
        len(fastq)
        print("number of sequences: " + str(len(fastq)/4))
        uniq_bcs_cell = {}
        for n in range(int(len(fastq)/4)):
            uniq_bcs_cell[n] = 0  #number of seqs as key, as value we will save the count of unique BCs
    seqs = list(uniq_bcs_cell.keys())


    nprocs = mp.cpu_count()
    p_size = nprocs - 25
    #############


    #since we are not importing the entire seq range as a whole into the function below (to risk of too high I/O stream), 
    #we will name the files based on the range instead and save into tmp dir and afterwards merge them all into same 
    def unique_bcs_per_seq(BC_clusters, seq_no, filepath, q): #take in the contents of starcode file as a dict (BC_clusters) and the seq number
        uniqe_bcs_cell = {}
        fname = os.path.basename(filepath)
        f_id = fname.split(".")[0]
        save_path = os.path.dirname(filepath)
        m = str(seq_no).split('(', 1)[1].split(')')[0]
        print(m)
        a = m.split(",")[0]
        b = m.split(" ")[1]
        file_sub_id = "_r_" + a + "-" + b

        for seq in seq_no:  #when we split the numbers, we put them into ranges of certain size so seq_no is now in range so need to iterate over it
            counts=0

            for key, value in BC_clusters.items():
                cluster_seq_ids = value.split(",")
                if str(seq) in cluster_seq_ids:
                    counts+=1
                    uniqe_bcs_cell[seq]=counts
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

    #multi_result = [pool.apply_async(power_n_list, (inp, 2)) for inp in inp_lists] uniq_bcs_cell
    multi_result = [p.apply_async(unique_bcs_per_seq, (BC_clusters_dict, seq,BBC_path,q,)) for seq in seq_lists]
    final_result = [result.get() for result in multi_result]
    #final_result



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
    df_counts_full = df_counts_full.drop('Unnamed: 0', inplace=True, axis=1)
    df_counts_full =df_counts_full.sort_values(by = "Cell", ascending=True)
    df_counts_full.to_csv(basePath + "/" + shared_name + "_Uniq_BCs_per_cell.csv", index=False)

    tmp_loc.cleanup()








#child_r, parent_w = os.pipe()
#OSError: [Errno 24] Too many open files
#######################################
#######################################
#######################################

#for n in test_seqs:
#    uniq_bcs_cell[n] = 0  #number of seqs as key, as value we will save the count of unique BCs
#uniq_bcs_cell
for seq_no in uniq_bcs_cell.keys(): #go over each seq in read 21 file
    uniq_bc_count=sum(any(str(seq_no) in s for s in subList) for subList in BC_clusters_dict.values())
    print(uniq_bc_count)
    uniq_bcs_cell[seq_no] = uniq_bc_count
uniq_bcs_cell
for seq_no in uniq_bcs_cell.keys(): #go over each seq in read 21 file
    unique_bs = 0
    round_count = 0
    print(seq_no)
    for cluster in BC_clusters_dict.keys(): #go over each cluster. if only 1 match found this mean that this seq no found only in one cluster BC
        cluster_seq_ids = BC_clusters_dict[cluster].split(",")  #split each cluster'd id seqs into separate elements
        #cluster_seq_ids 
        round_count +=1
        if str(seq_no) in cluster_seq_ids:  #you are counting this as a single occurence so the match is always 1!!!!
            #if seq_no is found in the cluster_seq_ids, then add a count and update the location in the dict of uniq_bcs_cell
            unique_bs+=1
            uniq_bcs_cell[seq_no] = unique_bs
            print("match found:" + str(unique_bs) +", which is:" + str(seq_no))
    # print("Round count: " + str(round_count))
print("Done")
uniq_bcs_cell
'''

############################################################
import random
from functools import reduce
import multiprocessing as mp
print("Number of processors: ", mp.cpu_count())

# Step 1: Init multiprocessing.Pool()
#pool = mp.Pool(mp.cpu_count())

def calculate_uniq_bcs_cell(BBC_path):
    with open(BBC_path, 'rb') as fastq:
        gzip_fastq = gzip.GzipFile(fileobj=fastq)
        fastq = gzip_fastq.read().splitlines() 
        print("number of sequences: " + str(len(fastq)/4))
        uniq_bcs_cell = {}
        for n in range(int(len(fastq)/4)):
            uniq_bcs_cell[n] = 0  #number of seqs as key, as value we will save the count of unique BCs

        for seq_no in uniq_bcs_cell.keys():
            print(seq_no)
            uniq_bcs_cell[seq_no] = unique_bcs_per_seq(BC_clusters_dict, str(seq_no))
            print(uniq_bcs_cell[seq_no])
    return(uniq_bcs_cell)
# Step 2: `pool.apply` the `howmany_within_range()`
#results = [pool.apply(uniq_bcs_cell, args=(seq_no, 4, 8)) for seq_no in uniq_bcs_cell]
# Step 3: Don't forget to close
#pool.close()    

###########
def unique_bcs_per_seq1(BC_clusters, seqs, q): #take in the contents of starcode file as a dict (BC_clusters) and the seq number
    for seq_no in seqs: #go over each seq in read 21 file
        uniqe_bcs_cell = {}
        counts=0
        for key, value in BC_clusters.items():
            cluster_seq_ids = value.split(",")
            if seq_no in cluster_seq_ids:
                counts+=1
                uniqe_bcs_cell[seq_no]=counts
            q.put(uniqe_bcs_cell)   #q object to save the result

        print("Processed Seq " + str(seq_no) + ". Unique BCs found: " + str(counts))
    return(uniqe_bcs_cell)
#############
def unique_bcs_per_seq(BC_clusters, seq_no,q): #take in the contents of starcode file as a dict (BC_clusters) and the seq number
    uniqe_bcs_cell = {}
    counts=0
    for key, value in BC_clusters.items():
        cluster_seq_ids = value.split(",")
        if seq_no in cluster_seq_ids:
            counts+=1
            uniqe_bcs_cell[seq_no]=counts
        q.put(uniqe_bcs_cell)   #q object to save the result
    print("Processed Seq " + str(seq_no) + ". Unique BCs found: " + str(counts))
    return(uniqe_bcs_cell)
#############

with open(BBC_path, 'rb') as fastq:
    gzip_fastq = gzip.GzipFile(fileobj=fastq)
    fastq = gzip_fastq.read().splitlines() 
    type(fastq)
    len(fastq)
    print("number of sequences: " + str(len(fastq)/4))
    uniq_bcs_cell = {}
    for n in range(int(len(fastq)/4)):
        uniq_bcs_cell[n] = 0  #number of seqs as key, as value we will save the count of unique BCs
seqs = list(uniq_bcs_cell.keys())
seqs
from multiprocessing import Pool, Process
nprocs = mp.cpu_count()
print(f"Number of CPU cores: {nprocs}")
nprocs=48-5

###################
# '''
# Pool and map() arent suited to situations that need to maintain state over time or, especially, situations where there are 
# two or more different operations that need to run and interact with each other in some way. For these kinds of problems, one
# needs to make use of somewhat more complex features of the multiprocessing module, such as Process, Queue and Event. Using 
# these features introduces a number of complicated issues that now need to be managed, especially with regard to cleanly 
# starting and stopping subprocesses, as well as coordinating between them.


# The Queue class includes the put method for depositing data and the get method for retrieving data

###################
#QUEUING method
# qout = mp.Queue()
# # Setup a list of processes that we want to run
# #processes = [mp.Process(target=unique_bcs_per_seq1, args=(BC_clusters_dict, seq,qout)) for seq in range(len(seqs))]
# processes = [mp.Process(target=unique_bcs_per_seq1, args=(BC_clusters_dict, seqs,qout))]

# processes.start()
# processes.join()

# result = qout.get()
# print(result)

######################################
from multiprocessing import Process, Queue
    
if __name__ ==  '__main__':
    q = Queue()
    for seq_no in range(len(seqs)):
        Process(target=unique_bcs_per_seq, args=(BC_clusters_dict, seq_no,q,)).start()
    for i in range(10):
        message = q.get()
        print(message)

######################################
# Run processes
processes = [mp.Process(target=unique_bcs_per_seq1, args=(BC_clusters_dict, seq,qout)) for seq in range(len(seqs))]
for p in processes:
    p.start()

# Exit the completed processes
for p in processes:
    p.join()

# Get process results from the output queue
results = [output.get() for p in processes]

print(results)

pool = mp.Pool(processes=40)
results = [pool.apply(unique_bcs_per_seq, args=(BC_clusters_dict, seq)) for seq in range(len(seqs))]
print(results)

###################
import multiprocessing

ret = {'foo': False}

def worker(queue):
    ret = queue.get()
    ret['foo'] = True
    queue.put(ret)

if __name__ == '__main__':
    queue = multiprocessing.Queue()
    queue.put(BC_clusters_dict)
    p = multiprocessing.Process(target=worker, args=(queue,))
    p.start()
    p.join()
    print(queue.get())  # Prints {"foo": True}
###################
import multiprocessing as mp

num_workers = mp.cpu_count()  

pool = mp.Pool(num_workers)
for task in tasks:
    pool.apply_async(func, args = (task,))

pool.close()
pool.join()



pool = mp.Pool(processes=nprocs)
#result = pool.starmap(unique_bcs_per_seq, [(BC_clusters_dict, seqs)])

#input list range in our case is the seqs and the function is BC_clusters
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
'''
Queue is a blocking, thread-safe queue that you can use to store the return values 
from the child processes. So you have to pass the queue to each process. Something less 
obvious here is that you have to get() from the queue before you join the Processes or else the
queue fills up and blocks everything.
'''
inp_lists = slice_data(seqs, nprocs)
multi_result = [pool.apply_async(unique_bcs_per_seq, (inp, 2)) for inp in inp_lists]

result = [x for p in multi_result for x in p.get()]
print(result)
unique_bcs_cell = pd.DataFrame([result.T])

p = Process(target=unique_bcs_per_seq, args=(BC_clusters_dict, seqs,))
p.start()
p.join()


p =Pool(5)

print(p.map(unique_bcs_per_seq, uniq_bcs_cell, seqs))


'''
pool = mp.Pool(mp.cpu_count())

pool = mp.Pool(processes=4)
processes = [pool.apply(uniq_bcs_cell, args=(BC_clusters, seq_no)) for seq_no in uniq_bcs_cell.keys()]
print(processes)
'''

# Define an output queue
output = mp.Queue()
#processes = [pool.apply(uniq_bcs_cell, args=(seq_no, 4, 8)) for seq_no in uniq_bcs_cell.keys()]


processes = [mp.Process(target=uniq_bcs_cell, args=(BC_clusters, seq_no)) for seq_no in seqs]
# Run processes
for p in processes:
    p.start()

# Exit the completed processes
for p in processes:
    p.join()

# Get process results from the output queue
results = [output.get() for p in processes]

print(results)
pool.close()    

print(results[:10])

    #Maybe save the number of seqs as list or use it more like as a range?
    #We 
    # Number of sublists containing 'Mary'
    for seq_no in uniq_bcs_cell.keys(): #go over each seq in read 21 file
        uniq_bc_count=sum(any(str(seq_no) in s for s in subList) for subList in BC_clusters_dict.values())
        print(uniq_bc_count)
        uniq_bcs_cell[seq_no] = uniq_bc_count
    uniq_bcs_cell[2]
    for seq_no in uniq_bcs_cell.keys():
        print(seq_no)
        uniq_bcs_cell[seq_no] = unique_bcs_per_seq(BC_clusters_dict, str(seq_no))
        print(uniq_bcs_cell[seq_no])
    # Number of strings containing 'Mary'
   # print sum(sum('Mary' in s for s in subList) for subList in myDict.values())
    #uniq_bcs_cell
    len(uniq_bcs_cell.keys())
    for seq_no in uniq_bcs_cell.keys(): #go over each seq in read 21 file
        unique_bs = 0
        for cluster in BC_clusters_dict.keys(): 
            print(cluster)
            cluster_seq_ids = BC_clusters_dict[cluster].split(",")  #split each cluster'd id seqs into separate elements
            if str(seq_no) in cluster_seq_ids:
                unique_bs+=1
                uniq_bcs_cell[seq_no] = unique_bs
                print("Match found:" + str(unique_bs)+ " - Seq nro: "+ str(seq_no) )
    print("Done")
uniq_bcs_cell


##################################
##################################
##################################

a = '347380'
for cluster in BC_clusters_dict.keys(): 
        cluster_seq_ids = BC_clusters_dict[cluster].split(",")  #split each cluster'd id seqs into separate elements
        #print(cluster_seq_ids)

        if a in cluster_seq_ids :
            print(True)
            '''
            with each read, iterate over all starcode cluster bcs, check if a matching seq id 
            found inside them, if it is, then include this into the count in uniq_bcs_cell values
            '''    

gzip_fastq = gzip.GzipFile(fileobj=fastq)
BBC = list(SeqIO.parse(gzip_fastq, "fasta"))
print(BBC[0].id)  # first record
print(BBC[-1].id)  # last record

BBC = list(SeqIO.parse(BBC_path, "fasta"))
print(BBC[0].id)  # first record
print(BBC[-1].id)  # last record
BC_clusters = pd.read_csv(basePath + "/SciFi_Pool1-2_EnrichedBCs_All_AAV-tagBFP_starcode.tsv", sep='\t')

#What you need to do is to extract the bead barcode from the R21 read 
#and pair it with the number of unique barcodes from the starcode file 
#for that bead
'''
But starcode runs it for all R3 seqs. 
You pair it based on the lines with R21 and count the number of Clusters that match to it?

So from starcode file, for example lets say the first cluster contains seqs in places 1,2,3,
this means that i take the cluster and match it with the lines 1,2,3 in the r21 file and after 
matching all the clusters, count their numbers corresponding to each seq?
I would save this data into a table and transfer it to the anndata object later on?



We have read 21 which we trim. then we extract the bead bc. presumably the trimmed 
version is the bc + umi right? We pair this with the corresponding samples unique BCs extracted
from read 3

So do I dont demultiplex the reads but rather take aav_asyn_R3's barcodes and 
match them to r21? but these are not the same. how do i match them? 

Do i demultiplex all the reads after trimming r21, THEN run starcode to read 3. starcode
output is the cluster seq, number of seqs in it and then shows the read positions IN read3 input file right?
I need to pair this bead bc with number of unique barcodes.

I take the R21 fastq file. I count the number of BCs this file has


each R21 seq in lets say A9 file are bead bcs that had this A9 primer identifier. 
'''
