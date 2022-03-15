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



basePath = 'SoloNova10x_correctAAV/'
basePath= "/media/data/AtteR/scifi-analysis/scifi6/2nd_try/output_enriched/"

#count the number of seqs in the original file, label based on this

#from starcode merge the seq to its respective seq. i.e. col1 content is added 
#to all sequences based on the positions in col2
from Bio import SeqIO
import gzip

'''
for f in os.listdir(basePath):
    if ".tsv" in f and "AAV" in f and not "Log" in f:
        print(f)
        BC_clusters = pd.read_csv(f, sep='\t')
        orig_file = f.replace("oDT", "WP")
'''
# starcode_outp = basePath + "/SciFi_Pool1-2_EnrichedBCs_All_AAV-tagBFP_starcode.tsv"
# BBC_path = basePath + "/SciFi_Pool1-2_EnrichedBCs_All_AAV-tagBFP_R21.fastq.gz"

# starcode_outp = basePath + "/SciFi_Pool1-2_EnrichedBCs_All_AAV-aSyn_starcode.tsv"
# BBC_path = basePath + "/SciFi_Pool1-2_EnrichedBCs_All_AAV-aSyn_R21.fastq.gz" 


starcode_outp = basePath + "/SciFi_Pool1-2_EnrichedBCs_All_AAV-aSyn_trimmed-WP_A9_16_starcode.tsv"
BBC_path = basePath + "/SciFi_Pool1-2_EnrichedBCs_All_AAV-aSyn_trimmed-WP_A9_16_R21.fastq.gz" 
#Number the seqs from fastq file, make them as dict keys with empty values. Then go over the starcode file, saving cluster as key and 
#the locations of seqs belonging to the cluster into the value. go over each each line in fastq
#file and compare the results to the values in the dict. If a match is found, add a count to
#the value part of the dict

BC_clusters = pd.read_csv(starcode_outp, sep='\t', header=None, index_col=0, dtype={0: str, 1: 'Int64'})
len(BC_clusters.index)
BC_clusters.iloc[3,1]
clusters = list(BC_clusters.index)
clusters

BC_clusters_dict = {}
for x in range(len(BC_clusters.index)):
    BC_clusters.iloc[0,:1]
    BC_clusters_dict[clusters[x]]=BC_clusters.iloc[x,1]
#BC_clusters_dict['CCGGGAAGACATAATTCAGGAAGC']
#BC_clusters_dict.keys()
#BC_clusters_dict
#len(BC_clusters_dict.keys())

############################################################
############################################################
############################################################
'''
BC_clusters_dict = {'BBC1':"1,2,3,4", 'BBC2': "1,3,5", 'BBC3': "2,3,4,5", 'BBC4':"1,2", 'BBC5':"1,3"}
BC_clusters_dict
test_seqs = random.sample(range(1, 6), 5)
test_seqs
#test_seqs.append(975423)
uniq_bcs_cell = {}uniq_bcs_cell = {}
for n in test_seqs:
    uniq_bcs_cell[n] = 0  #number of seqs as key, as value we will save the count of unique BCs
uniq_bcs_cell
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
            if str(seq_no) in cluster_seq_ids:
                counts+=1
                uniqe_bcs_cell[seq_no]=counts
            q.put(uniqe_bcs_cell)   #q object to save the result

        print("Processed Seq " + str(seq_no) + ". Unique BCs found: " + str(counts))
    return(uniqe_bcs_cell)
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
type(seqs[0])
from multiprocessing import Pool, Process
nprocs = mp.cpu_count()
print(f"Number of CPU cores: {nprocs}")
nprocs=nprocs-5

###################
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

p_size = nprocs - 25
from multiprocessing import Process, Queue, Manager
#############


def unique_bcs_per_seq(BC_clusters, seq_no, filepath, q): #take in the contents of starcode file as a dict (BC_clusters) and the seq number
    uniqe_bcs_cell = {}
    fname = os.path.basename(filepath)
    f_id = fname.split(".")[0]
    save_path = os.path.dirname(filepath)

    counts=0
    for seq in seq_no:  #when we split the numbers, we put them into ranges of certain size so seq_no is now in range so need to iterate over it
        for key, value in BC_clusters.items():
            cluster_seq_ids = value.split(",")
            if str(seq) in cluster_seq_ids:
                counts+=1
                uniqe_bcs_cell[seq]=counts
            #q.put(uniqe_bcs_cell)   #q object to save the result
        print("Processed Seq " + str(seq) + ". Unique BCs found: " + str(counts))
    final_df = pd.DataFrame([uniqe_bcs_cell])
    print(save_path)
    with open(save_path + "/" +  fname +".csv", "w") as out: #write the dataframe into an output file
        final_df.to_csv(out, sep='\t')
        # df.to_string(out, index=None)
        print(fname + ' output info file saved!')

    return(uniqe_bcs_cell)
#############


def reader(i,q):
    message = q.get() #retrieves the queue
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
#the lines... So instead run this script on the demultiplexed samples (so first starcode, then this), then map and then import into
#anndata
if __name__ ==  '__main__':
    p = Pool(p_size)
    m = Manager()

    q = m.Queue()
    # with mp.Manager() as manager:
    #     q = manager.Queue()
    #     d = manager.dict()
    #     with manager.Pool(p_size) as p:
    #         #pool.map(f, repeat(d, 10))
    #         seq_lists = slice_data(range(len(seqs)), nprocs)
    #         #multi_result = [pool.apply_async(power_n_list, (inp, 2)) for inp in inp_lists]
    #         multi_result = [p.apply_async(unique_bcs_per_seq, (BC_clusters_dict, seq,q,)) for seq in seq_lists]
    #         result = [x for p in multi_result for x in p.get()]
###############################################################################################################
###############################################################################################################
    seq_lists = slice_data(range(len(seqs)), nprocs)
    #multi_result = [pool.apply_async(power_n_list, (inp, 2)) for inp in inp_lists] uniq_bcs_cell
    multi_result = [p.apply_async(unique_bcs_per_seq, (BC_clusters_dict, seq,BBC_path,q,)) for seq in seq_lists]
#workers = [pool.apply_async(toParallel, args=(ht,token,)) for ht in hashtag['hashtag']]
    #multi_result = [p.apply_async(unique_bcs_per_seq, (BC_clusters_dict, seq,q,)) for seq in seq_lists]
    #result = [x for p in multi_result for x in p.get()]
    final_result = [result.get() for result in multi_result]
    final_result
    result
    type(result)
    print(result['8970'])

    df = pd.DataFrame.from_dict(data)
    # #manager enables us to manage the queue and access the different workers (unique_bcs_per_seq functions). It also closes pools when they are done.
    # m = Manager()

    # q = m.Queue()
    # # Create a group of parallel writers and start them
    # for seq_no in range(len(seqs)):
    #     Process(target=unique_bcs_per_seq, args=(BC_clusters_dict, seq_no,q,)).start()
    #     #results = pool.map(myLevenshteinFunction, stringList)
    # # Create multiprocessing pool
    # p = Pool(p_size)

    # # Create a group of parallel readers and start them
    # # Number of readers is matching the number of writers (unique_bcs_per_seq)
    # # However, the number of simultaneously running
    # # readers is constrained to the pool size
    # readers = []
    # for i in range(seqs):
    #     readers.append(p.apply_async(reader, (seqs,q,)))
    #     readers.append(p.apply_async(unique_bcs_per_seq, (inp, 2)) for inp in inp_lists])

    #     p.close()
    #     p.join()


    # [r.get() for r in readers]
    #for i in range(10):
    #    message = q.get()
    #    print(message)
    '''
    Queue is a blocking, thread-safe queue that you can use to store the return values 
    from the child processes. So you have to pass the queue to each process. Something less 
    obvious here is that you have to get() from the queue before you join the Processes or else the
    queue fills up and blocks everything.
    '''
    inp_lists = slice_data(seqs, nprocs)
    multi_result = [pool.apply_async(unique_bcs_per_seq, (inp, 2)) for inp in inp_lists]

    result = [x for p in multi_result for x in p.get()]

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
