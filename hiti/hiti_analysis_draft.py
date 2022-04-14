from curses import window
from fnmatch import translate
import ntpath
from re import L
from wsgiref import headers
from nbdev.showdoc import *
import os
import warnings
import matplotlib.pyplot as plt
from subprocess import call
import subprocess
import glob
import pandas as pd
import tempfile
from skbio.alignment import local_pairwise_align_ssw
from skbio import DNA

import os
os.getcwd()
os.chdir("/media/data/AtteR/projects/hiti")
from alignment_scripts import *
import Bio.Align.Applications; dir(Bio.Align.Applications)
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline#Read in unfiltered data
from Bio import AlignIO
from Bio import pairwise2
from Bio import SeqIO
from Bio.pairwise2 import *
from alignment_scripts import *

#into the id save the row number, perc.seqs and perc.match to the ref

#we wont be able to find 100% match as the scar region seems to vary 

#now we have found the longest seq and trimmed based on that and calculated match % to the reference
#next put this into function and compare results. If more than one match of max length, then report how many seqs have this
#we should then take each one of these seqs, make them as seqrecords and save as a fasta file which then can be visualised
#-----try this with the msa file!

#we take the consensus seqs and we trim the seqs as in the upstream. We do NT MSA but we will find out that at 5' end we will get perfect matches at times (this applies with old data?)
#but not true with new data. on 5' we will see that 90%ish is a perfect alignment etc. 
#we may get single base reading errors in some areas of the seq but they are irrelevant. we care about the scar. We line all them up and find the perfect matches or close to it


#at 3' end we link the mcherry to the arc which is why we need to make sure that this area is translated into AAs properly
#Start with 3' end ------- do MSA ------ ok lets say the first one doesnt have a perfect match but rather lets say 93%, then 30% etc and on row 12 we then have a perfect match - matches the 6 bps to the 6bps in scar
#we can make a cutoff the way that we show the rows that have lets say 90% match. or we can count the number of rows

#we must make sure that the seq is A)in frame and B)SCAR exists, then the seq is fine
#MSA will align the full consensus seq. we trim based on the perfect match, but this is relevant when visualising the seqs

#we go from 3' end (where we dont have errors) till the 5'end (where we do have erros) in terms of the seq rows - we want to represent the sides in a balanced way
#so either we represent the same number of rows or we represent into the same number of frequency as the accurate alignment here (used when visualising)
#We get the row N from the 3' end 

#look at the 3' end seq only and plot it till the perfect match ------- so this is the frequency of correct edits (seq cluster percentage-number of seqs in this cluster?) respective to the specific PAM edit
#we take all the clusters (with respective info regarding the percentage)

#we count the Nrow we needed to get to the perfect match at 5' end and this is the same number of consensus seqs from 5' end that we would also align via MSA


#when translating, you need to know the start codon


#with aa we do the msa, then we can tell which of the seqs are in reading frame and thus generate a protein ---- this will tell us how many of the red cells can we trust
#to be the fusion between mcherry and arc

#then we translate to AAs as we want to check that things are in reading frame which is done by first figuring out the start codon


#5' will be accurate and thus not very descriptive in terms of whether our tech worked or not
    #the function returns a df of certain animal, e.g. 8. and we merge the dfs into a complete one later. 
    #now we have h and s groups though which we want to keep as separate (before normalising the counts, then merge?). Normalisation should happen before inputting the seqs
    #into the starcode or after we have input them, after which we find the reads belonging to the certain clusters. 
#we need to normalise so that one group will not skew the results but rather so that both have same effect to the downstream analysis. thus, BEFORE starsolo normalise them?

#currently when we pool the 3p samples, we pool them (s and h groups) together as raw counts. they dont match by percentage. the raw counts can make things skewed. 
#first we must normalise the samples prior to pooling them. 

#lets say we have 1m reads (s group) and the other has 100 000 (h group). easiest is to divide all read counts by the total reads of the sample. 
#each will contribute equally.
#---- so we take the first sample's first group 7s, take its reads. then we take the second group 7h and its reads. we cluster the reads based on starsolo and get the seq clusters.
#or NO! Out of the 1m reads of 7s a certain number goes into lets say cluster 1 of starsolo seq etc. so cluster 1 has 400 000 reads, cluster 2 100 000 etc. we divide all of these 
#results with the total number and we get percentages instead of RAW counts. THEN we pool them and the toal will be 2? So we take the matching reads from column 7s and 7h and merge them, the rest as extra.


#we will get 7h and 7s reads so lets say 1 000 000 and 100 000. we then normalise them against themselves. then when we sum them we get 1+1=2.
#


"""
understand how the pairwise alignment was done and put into df
use different ones and compare results
visualise
maybe apply same with trimming methods
visualise
find the frame
translate to AAs
visualise
"""
transgene = 'mCherry'
assay_end = '3p'
read_fwd = True
animal_list = [7, 8, 9, 10, 11, 12] 
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' #region prior to r2 primer
lliteral = ' literal=GGCGGCATGGACGAG' #to check that its on target with mcherry
rliteral = ' literal=CATATGACCACCGG'
#base_path = '/home/lcadmin/mnm-lts/SequencingFiles/arc_hiti_asyn/HITI1-8977973/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
#export_path = '/home/lcadmin/mnm-lts/HITI-analysis/'
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path = '/media/data/AtteR/projects/hiti/output/'
target_sequence = "CGGCGGCATGGACGAGCTGTACAAGGTCGGTGCTGCGGCTCCGCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAGGTCTTCTACCGTCTGGAGAGG"

animal_nr = str(9)
"Filters and trims the reads"
search_path = base_path+animal_nr+'*'+transgene+'*'+assay_end+'*/'
search_path

hip=animal_nr + "_" + transgene + "*h_"
str=animal_nr + transgene + "*s_"
    

from functools import reduce


#you could trim and starcode all individual lane files of certain animal first. after this you sum these based on matches

def trimRead_hiti(data_dict, transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd):
    complete_df = pd.DataFrame({'sequence': ['CTGTACAAGGTCGGTGCTGCGGCTCCGCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGAC']})
    complete_df
    for animal in data_dict.keys():
        animal_group_name=animal.split("_")[0] + "_" + animal.split("_")[2]

        dfs_lane=[]
        for search_path in data_dict[animal]:
            animal_p5_cat = tempfile.NamedTemporaryFile(suffix = '.fastq.gz').name
            animal_p7_cat = tempfile.NamedTemporaryFile(suffix = '.fastq.gz').name
            test_file_p5_out = tempfile.NamedTemporaryFile(suffix = '.fastq').name
            test_file_p7_out = tempfile.NamedTemporaryFile(suffix = '.fastq').name
            test_file_p5_filter = tempfile.NamedTemporaryFile(suffix = '.fastq').name

            if read_fwd:
                animal_p5 = glob.glob(search_path+'/*R1*')
                animal_p7 = glob.glob(search_path+'/*R2*')
                #display('Forward run Animal: '+animal_nr)
            else:
                animal_p5 = glob.glob(search_path+'/*R2*')
                animal_p7 = glob.glob(search_path+'/*R1*')
                #display('Reverse run Animal: '+animal_nr)
            animal_p7

            cat_p5= "cat "+" ".join(animal_p5)+" > "+animal_p5_cat
            print(cat_p5)
            #os.system(cat_p5)
            call([cat_p5], shell=True) #call caused the terminal to freeze so switched to os
            cat_p7= "cat "+" ".join(animal_p7)+" > "+animal_p7_cat
            call([cat_p7], shell=True)
            #os.system(cat_p7)

            stats_out = export_path+group+'_'+transgene+'_'+assay_end+'_stats-filter.txt'

            kmer = '20'
            hdist = '3'
            param=" k="+kmer+" hdist="+hdist+" rcomp=f skipr2=t threads=32 overwrite=true"

            #to check if the read is an amplicon
            call_sequence = "/media/data/AtteR/Attes_bin/bbmap/bbduk.sh in="+animal_p7_cat+" in2="+animal_p5_cat+" outm1="+test_file_p7_out+" outm2="+test_file_p5_out+" literal="+filterlitteral+" stats="+stats_out + param
            call([call_sequence], shell=True)

            #actual trimming
            call_sequence = "/media/data/AtteR/Attes_bin/bbmap/bbduk.sh in="+test_file_p5_out+" out="+test_file_p5_filter+ " literal=AAAAAAAAA,CCCCCCCCC,GGGGGGGGG,TTTTTTTTT k=9 mm=f overwrite=true minlength=40"
            call([call_sequence], shell=True)

            test_file_p5_out_starcode = tempfile.NamedTemporaryFile(suffix = '.tsv').name
            starcode_call= "/media/data/AtteR/Attes_bin/starcode/starcode -i "+test_file_p5_filter+" -t 32 -o "+test_file_p5_out_starcode
            call([starcode_call], shell=True)

            df=pd.read_csv(test_file_p5_out_starcode, sep='\t', header=None)
            df = df.rename(columns={0: 'sequence', 1:'count'})
            dfs_lane.append(df)
            print(animal_group_name + " done!")
        #we iterate over all the individual dfs and merge them by taking the seq column of all dfs and placing them under the new dfs seq and do the same with counts
        df_all_lanes=reduce(lambda  left,right: pd.merge(left,right,on='sequence', how='outer'), dfs_lane)
        #reduce is useful when you need to apply a function to an iterable and reduce it to a single cumulative value.
        df_all_lanes["count"]=df_all_lanes.sum(axis=1) #make a column with total count sum of reads and remove the rest. This gives a df that has the seqs and the total counts from all lanes
        df_all_lanes.drop(df_all_lanes.iloc[:, 1:((len(df_all_lanes.columns)-1))], inplace = True, axis = 1)

        #Once you have combined all the lane dfs, then you take the percentage
        total_counts = int(df_all_lanes[['count']].sum())
        df_all_lanes['percent'] = (df_all_lanes['count'] / total_counts)
        df_all_lanes = df_all_lanes.rename(columns={'percent':animal_group_name+'_percent','count':animal_group_name+'_count',})
        complete_df = pd.merge(complete_df, df_all_lanes, on="sequence", how='outer')
        print("A full df containing the sum from all lanes of " + animal_group_name + " is done!")

    return complete_df


#Pool the reads based on striatum and hippocampus but prior to pooling them, normalise against oneself as otherwise one will contribute 
#more than the other. 

#hip_folders = [folder for folder in os.listdir(base_path) if "mCherry" in folder and "h_" in folder or "s_" in folder]
group_folders = [folder for folder in os.listdir(base_path) if "mCherry" in folder]
#str_folders = [folder for folder in os.listdir(base_path) if "mCherry" in folder and "s_" in folder]
group_folders
search_paths_groups = []
search_paths_s = []


def animal_names(group_folders):
    animals=[]
    animal_list = [*range(13)]

    for s in group_folders:
        animal_name="_".join(s.split("_")[:3])
        if int(animal_name.split("_")[0]) in animal_list:
            print(animal_name)
            animals.append("_".join(s.split("_")[:3]))
    #animals=list(set(animals))
    return(sorted(list(((set(animals))))))

animals = animal_names(group_folders)
animals

#key: animal number and whether it comes from striatum or hippocampus, value: the paths to all lane subdirs

data_dict = dict()
for animal_group in animals:
    lanes=[]
    for g_f in group_folders:
        if animal_group in g_f:
            g_p=base_path+g_f
            lanes.append(g_p)
            data_dict[animal_group]=lanes

data_dict

for animal in data_dict.keys():
    print(animal)
    animal_group_name=animal.split("_")[0] + "_" + animal.split("_")[2]
    group=search_path.split("/")[-1].split("_")[2]
    dfs_lane=[]
    for search_path in data_dict[animal]:
        print(search_path)
        group=search_path.split("/")[-1].split("_")[2]
        print(group)

#take subdirs as same group (dic key) if they come from the same animal and group but different lanes
#2_mCherry_4h_5p_L


# for h_f in hip_folders:
#     s_p=base_path+h_f
#     search_paths_h.append(s_p)

# for s_f in str_folders:
#     s_p=base_path+s_f
#     search_paths_s.append(s_p)
# search_paths_s

#this line caused the perfect match to occur since we had the target seq included!
#complete_df = pd.DataFrame({'sequence': [target_sequence]})

# go over each individual animal subfolder i.e. s and h, process them and generate a starcode file, then a df
# so unlike with last function which gave individual df8, df9 etc., now we get df8s, df8h etc. where they have been normalised against the total number 

#after which we merge the ones that have matching numbers
complete_df = pd.DataFrame({'sequence': ['CTGTACAAGGTCGGTGCTGCGGCTCCGCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGAC']})
complete_df

data_dict
#you could trim and starcode all individual lane files of certain animal first. after this you sum these based on matches
full_df = trimRead_hiti(data_dict,transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd)
#the function returns you a complete dataframe containing animals_x_brain_area (h,s). since each animal_x_brain_area contains data from 4 different subdirs, these have been
# summed into the same ones, i.e. column 12_6h contains the sequences from lanes 1-4 all summed up  
#now take the percentage values of each animal (so brain area s and h), merge into same column so that the percentages will be a sum of the two.

#each animal contains data either from striatum or hippocampus after all.

a = [f for f in full_df.columns if "h_percent" in f]
a
#animal group contains all lanes of the certain data
#sum the percentages of each seq cluster (animals 7-12)      --------------- I did not use this downstream. Percent sum the numbers and thus use these instead of animal spec ones
full_df = full_df.fillna(value=0)
perc_cols = [col for col in full_df.columns if 'percent' in col]
perc_cols
full_df['percent_sum'] = full_df[perc_cols].sum(axis=1)
full_df.sort_values(by=['percent_sum'], ascending=False, inplace=True)
full_df.head()

full_df.index[0]
full_df.head()
#trim seqs that contribute less than 0.01% percentage
rows_drop=[]
for i, perc in enumerate(full_df.iloc[:,-1]):
    if perc<0.001:
        rows_drop.append(i)
full_df_trim=full_df.drop(rows_drop, axis=0, inplace=False)

full_df
rows_drop
#generate a column with seq alignment where you take the seq cluster from sequence column and map it to the target --- part of Thomas' script
complete_df.loc[:,'sequence_align'] = complete_df.loc[:,'sequence'].apply(lambda x: align_to_ref(x, target_sequence))
export_csv = export_path+transgene+'_'+assay_end+'.csv'
complete_df.to_csv(export_csv, index=False)

complete_df.loc[:,'sequence_align']

from Bio.Seq import Seq 

complete_df.iloc[:,0]
s = Seq(complete_df.iloc[1,0])
s

#transforms the aligned seqs into seqrecord objects
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


#extract the aligned seq part, ignoring the scores etc. 
#to process and save the files downstream, must convert the result into seq record
#first take specific alignment from the alignment object via indexing, then format it into a string, split based on line changes
#and take the 3rd element, i.e the NTs matching the template one, i.e. the alignment. Transform this into seq object, save
#objects into a list and then iteratively save into fasta file. visualise them
#try different alignment methods


#we need to know what percentage of the clusters mapped with high accuracy, ideally 100%, to the reference genome. Thus
#we must save the percentages from the starcluster data
############################
############################
seq_and_perc = {}

full_df.iloc[:, [0,-1]]


#We take the seq and percent sum for each column and then align. the full_df contains ALL the clustered seqs from all the animals. we have taken the percentage amount of raw counts
#for each animal and then we have summed them together, giving as a unit providing information how much of the certain seq cluster is when summed over all the animals.
#why is this done? after we map the sequences to the reference and find the alignments, we can use this as information regarding how many of the amplicons were of certain kind and
#how well it maps to the ref. 

#take the percentage group (8_percent etc) as the key and sequence as value
for cols in full_df.columns:
    #print(cols)
    if "percent" in cols and not "sum" in cols:
        print(cols)
        seq_and_perc[cols]=complete_df.loc[:,["sequence", cols]]
seq_and_perc["8_percent"].iloc[:,0]
############################
############################
seq_and_perc
#if we align using muscle we will need to convert our seq files into fasta first
from Bio.Align.Applications import MuscleCommandline

result="/media/data/AtteR/projects/hiti/mcherry_p3_seq_cluster.fasta"
id_f = 1
with open(result, "w") as handle:
#write ref template to the top of the file
    seq_obj = SeqRecord(Seq(target_sequence), id=str(0), description="mcherry_p3_seq_clusters")
    count = SeqIO.write(seq_obj, handle, "fasta")
    for seq in complete_df.iloc[:,0]:
        seq_obj = SeqRecord(Seq(seq), id=str(id_f), description="mcherry_p3_seq_clusters")
        count = SeqIO.write(seq_obj, handle, "fasta")
        id_f+=1
############################
############################

complete_df
#make a file for especially muscle alignment which you can then visualise. This file contains the 
#percentages of seq clusters

#save seq dict as a fasta file

'''
Hi Tomas I have gone over the data and reprocessed it. a thing that i noticed was that each animal 
actually only contained data from hippocampus or striatum, not both (given that I have read it right).
You sure you hadnt already had excluded the brain area that had poor number of reads in the past?
For example animal 8:
8_mCherry_4h_3p_L001-ds.23808003b8304ec79a33a3d33e2afa46:
8_mCherry_4h_3p_L002-ds.9ab8f311ca9044bdb4ba14a98707ceff:
8_mCherry_4h_3p_L003-ds.f49838d03d6445429f6d7ef3f2b99167:
8_mCherry_4h_3p_L004-ds.1f5f77a9c4b0416b837d6256cf022335:

Anyways, I rewrote the function and reprocessed the data, summed the lanes of each animal into one df and in the end merged all the animals as one big df. summed percentages across animals for a given sequence cluster.
I wrote these sequences into fasta files, giving as IDs the summed percentage values and the line numbers. I am running the mapping using different algorithms for pairwise and MS. its taking a while though.

Also, the data contained animals from 1-12 but your script had only taken the ones from 8-12 before. I was wondering what the reasoning behind this was? were the others poor quality?
'''

full_df.iloc[0,-1]
def save_fasta(filename, full_df):
    id_f = 1
    with open(result, "w") as handle:
        seq_obj = SeqRecord(Seq(target_sequence), id=str(0), description="mcherry_p3_seq_ref")
        count = SeqIO.write(seq_obj, handle, "fasta")
        for seq_i in range(len(full_df.iloc[:,-1])):
            descr="CluSeq_%: " + str((full_df.iloc[seq_i,-1]))
            print(descr)
            seq_obj = SeqRecord(Seq(full_df.iloc[seq_i,0]), id=str(id_f), description=descr)
            count = SeqIO.write(seq_obj, handle, "fasta")
            id_f+=1
    print("Saved!")


for seq in range(len(full_df.iloc[:,-1])):
    print(seq)
#this saves them via seq record formatting. However, when importing this content back for translation,
#there are new line characters in some of the longer seqs so need to create another approach too
#---or open them via seqio and it will retain the correct formatting

def save_fasta_seqrec(filename, seq_and_perc):
    id_f=1
    with open(filename, "w") as handle:
        seq_obj = SeqRecord(Seq(target_sequence), id=str(0), description="mcherry_p3_seq_ref")
        header=">0"+" mcherry_p3_seq_ref"
        handle.write(header + "\n" + target_sequence + "\n")

        for group in seq_and_perc.keys():
            for seq_i in range(len(seq_and_perc[group])):
                header=">"+ str(id_f)+" CluSeq_%: " + str(round((seq_and_perc[group].iloc[seq_i,1]*100),4))
                handle.write(header + "\n" + seq_and_perc[group].iloc[seq_i,0] + "\n")
                id_f+=1

result="/media/data/AtteR/projects/hiti/mcherry_p3_seq_cluster_all.fasta"
result="/media/data/AtteR/projects/hiti/mcherry_p3_seq_cluster_all_redo.fasta"

save_fasta_seqrec(result, seq_and_perc)
save_fasta(result, full_df_trim)
aligned_seqs = []


###############
#ALIGNMENTS
def align_global(amplicon,target_sequence):
    alignments = pairwise2.align.globalms(target_sequence, amplicon,  2, -1, -.5, -.1)
    #alignments = pairwise2.align.globalms(target_sequence, seq_and_perc[group][:,0],  2, -1, -.5, -.1)

    alignm=format_alignment(*alignments[0])
    #make string objects, save into a list. then count which one has least ----, cut off the rest of the reference based on this?
    seq_align = alignm.split("\n")[2]
    #nt_count=count_nts(seq_align)
    #seq_and_perc[group]["match"]
    return(seq_align)
def align_global2(amplicon,target_sequence):
    alignments = pairwise2.align.globalxx(target_sequence, amplicon)

    alignm=format_alignment(*alignments[0])
    #make string objects, save into a list. then count which one has least ----, cut off the rest of the reference based on this?
    seq_align = alignm.split("\n")[2]
    #nt_count=count_nts(seq_align)
    #seq_and_perc[group]["match"]
    return(seq_align)
def align_local(amplicon,target_sequence):
    alignments = pairwise2.align.localxx(target_sequence, amplicon)
    #alignments = pairwise2.align.globalms(target_sequence, seq_and_perc[group][:,0],  2, -1, -.5, -.1)

    alignm=format_alignment(*alignments[0])
    #make string objects, save into a list. then count which one has least ----, cut off the rest of the reference based on this?
    seq_align = alignm.split("\n")[2]
    #nt_count=count_nts(seq_align)
    #seq_and_perc[group]["match"]
    return(seq_align)


def align_local2(query_sequence, target_sequence):
    alignment, score, start_end_positions = local_pairwise_align_ssw(DNA(target_sequence),DNA(query_sequence),gap_open_penalty = 3,gap_extend_penalty = 1)
    out_align = ('-'*start_end_positions[0][0])+str(alignment[1])+('-'*(len(target_sequence)-start_end_positions[0][1]-1))
    return out_align

def align_local3(amplicon,target_sequence):
    alignments = pairwise2.align.localms(target_sequence, amplicon, 2, -1, -.5, -.1)
    #alignments = pairwise2.align.globalms(target_sequence, seq_and_perc[group][:,0],  2, -1, -.5, -.1)
    alignm=format_alignment(*alignments[0])
    #make string objects, save into a list. then count which one has least ----, cut off the rest of the reference based on this?
    seq_align = alignm.split("\n")[2]
    #nt_count=count_nts(seq_align)
    #seq_and_perc[group]["match"]
    return(seq_align)


test = align_local2(seq_and_perc["8_percent"].iloc[1,0], target_sequence)
tests = pairwise2.align.globalmx(seq_and_perc["8_percent"].iloc[1,0], target_sequence,  2, -1, -.5, -.1)
tests
test_a = format_alignment(*tests[0])
test_a = pairwise2.align.globalms(target_sequence, seq_and_perc["8_percent"].iloc[1,0],  2, -1, -.5, -.1)
test_a[0]


seq_align = test_a.split("\n")[2]
seq_align
def longest_align(seqs):
    return max(seqs, key=len)


id_f=1
aligned_data=dict()
#align all the data, save into dict, then ensure that all the seqs are same length (take the longest seq). IF not, then make them equal length ny adding Nx"-"
for seq_i in range(len(full_df_trim.iloc[:,-1])):
        header=">"+ str(id_f)+" CluSeq: " + str((full_df_trim.iloc[seq_i,-1]))
        seq_obj_1= align_local3(full_df_trim.iloc[seq_i,0], target_sequence)
        aligned_data[header]=seq_obj_1
        id_f+=1

longest_seq=longest_align(aligned_data.values())
aligned_data
longest_seq
len(longest_align(aligned_data.values()))
aligned_data['>582 CluSeq: 0.0']
N_dashes=len(longest_align(aligned_data.values()))-len(aligned_data['>582 CluSeq: 0.0'])
N_dashes
aligned_data['>582 CluSeq: 0.0']=aligned_data['>582 CluSeq: 0.0']+N_dashes*"-"
for id in aligned_data.keys():
    if len(aligned_data[id])==len(longest_align(aligned_data.values())):
        continue
    else:
        N_dashes=len(longest_align(aligned_data.values()))-len(aligned_data[id])
        aligned_data[id]=aligned_data[id]+N_dashes*"-"
aligned_data['>579 CluSeq: 1.1998809718075966e-06']
test_a = pairwise2.align.globalms(target_sequence, aligned_data['>579 CluSeq: 1.1998809718075966e-06'],  2, -1, -.5, -.1)
alignment, score, start_end_positions = local_pairwise_align_ssw(DNA(target_sequence),DNA(aligned_data['>579 CluSeq: 1.1998809718075966e-06']),gap_open_penalty = 3,gap_extend_penalty = 1)
out_align = ('-'*start_end_positions[0][0])+str(alignment[1])+('-'*(len(target_sequence)-start_end_positions[0][1]-1))

test_a

import re

pattern = r'[(\d|\s]'
s="2 C-GCGGCATGGACGAGCTGTACAAGGTCGGTG----GG-"
# Remove characters 's', 'a' and 'i' from a string
mod_string = re.sub(r'[(\d|\s]', '', s)
mod_string

len(aligned_data['>579 CluSeq: 1.1998809718075966e-06'])
#downstream an issue with visualising the seqs using mview is that all the seqs are not SAME length. thus, needs to fix this
def align_and_save(filename, full_df):
    id_f=1
    aligned_data=dict()
    #align all the data, save into dict, then ensure that all the seqs are same length (take the longest seq). IF not, then make them equal length ny adding Nx"-"
    for seq_i in range(len(full_df.iloc[:,-1])):
            header=">"+ str(id_f)+" CluSeq: " + str((full_df.iloc[seq_i,-1]))
            seq_obj_1= align_local(full_df.iloc[seq_i,0], target_sequence)
            seq_obj_1 = re.sub(r'[(\d|\s]', '', seq_obj_1) #remove digits from the string caused by the alignment and empty spaces from the start
            aligned_data[header]=seq_obj_1
            id_f+=1

    #longest_seq=longest_align(aligned_data.values())
    for id in aligned_data.keys():
        if len(aligned_data[id])==len(target_sequence):
            continue
        else:
            N_dashes=len(target_sequence)-len(aligned_data[id]) 
            aligned_data[id]=aligned_data[id]+N_dashes*"-"

    with open(filename, "w") as handle:
        header=">0"+" mcherry_p3_seq_ref"
        handle.write(header + "\n" + target_sequence + "\n")
        for seq_i in aligned_data.keys():
            handle.write(seq_i + "\n" + aligned_data[seq_i] + "\n")

def align_and_save2(filename, seq_and_perc):
    id_f=1
    with open(filename, "w") as handle:
        seq_obj = SeqRecord(Seq(target_sequence), id=str(0), description="mcherry_p3_seq_ref")
        header=">0"+" mcherry_p3_seq_ref"
        handle.write(header + "\n" + target_sequence + "\n")

        for group in seq_and_perc.keys():
            for seq_i in range(len(seq_and_perc[group])):
            header=">"+ str(id_f)+" CluSeq: " + str((full_df.iloc[seq_i,-1]))
                seq_obj_1 = align_global(seq_and_perc[group].iloc[seq_i,0], target_sequence)
                handle.write(header + "\n" + seq_obj_1 + "\n")
                id_f+=1


    print("Saved into " + str(filename))


DNA(target_sequence)
target_sequence
seq_obj_1= align_local2(full_df.iloc[1,0], target_sequence)
alignment, score, start_end_positions = local_pairwise_align_ssw(target_sequence,full_df.iloc[1,0],gap_open_penalty = 3,gap_extend_penalty = 1)
seq_obj_1
full_df.iloc[1,0]

align_pairwise_loc1="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_aligned_pairwise_local_1.fasta"
align_pairwise_loc2="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_aligned_pairwise_local_2sk_redo.fasta"
align_pairwise_loc3="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_aligned_pairwise_local_3.fasta"

align_pairwise_glob1="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_aligned_pairwise_glob.fasta"
align_pairwise_glob2="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_aligned_pairwise_glob2.fasta"

musc = "/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_muscle_redo.fasta"
from Bio import SeqIO

for record in SeqIO.parse(musc, "fasta"):
    print(record.id)
    print(len(record.seq))
    print("==================")

#maybe add progress bars?
align_and_save(align_pairwise_loc1, full_df_trim)
align_and_save(align_pairwise_loc3, full_df_trim)

align_and_save(align_pairwise_glob1, full_df_trim)
align_and_save(align_pairwise_glob2, full_df_trim)

align_and_save2(align_pairwise_glob1, full_df_trim)

test
#returns a dict with percentage value of the certain cluster seq and the aligned seq as the value

outp="/media/data/AtteR/projects/hiti/mcherry_p3_seq_clusters_all2.fasta"

id_f=1
with open(outp, "w") as handle:
    seq_obj = SeqRecord(Seq(target_sequence), id=str(0), description="mcherry_p3_seq_ref")
    header=">0"+" mcherry_p3_seq_ref"
    handle.write(header + "\n" + target_sequence + "\n")

    for group in seq_and_perc.keys():
        for seq_i in range(len(seq_and_perc[group])):
            header=">"+ str(id_f)+" CluSeq_%: " + str(round((seq_and_perc[group].iloc[seq_i,1]*100),4))
            seq_obj_1 = seq_and_perc[group].iloc[seq_i,0]
            handle.write(header + "\n" + seq_obj_1 + "\n")
            id_f+=1



###############
#MUSCLE alignment
import subprocess
input = "/media/data/AtteR/projects/hiti/mcherry_p3_seq_clusters.fasta"
output = "/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_clusters_aligned.fasta"

subprocess.call('muscle -align %s -output %s'%(input,output)) #does not work, need to run from terminal 
###############



###############
#TRANSLATION
###############

from Bio.Data import CodonTable
from Bio.Seq import Seq
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
standard_table.start_codons
nt_seq=SeqIO.parse(result, "fasta")

s="CGGCGGCATGGACGAGCTGTACAAGGTCGGTGCTGCGGCTCCGAGATGGAGCTGGTCGAGATGAGCACCGCCGG"
ac=str(Seq(s).translate())
ac
def find_overlap(amplicon, mcherry_full):
    longest_seq = None
    n=0
    for i in range(1, min(len(amplicon), len(mcherry_full))):
        # Store the current iteration's intersection
        current_seq = intersects(amplicon, mcherry_full, i)
        
        # If this was one character too long, return.
        # Else, a longer intersection was found. Store it.
        if current_seq == None:
            print("No overlap found!")
            n+=1
            return 0
        else:
            longest_seq = current_seq
    # If we get here, the strings were the same.
    # For consistency, return longest_seq and its length.
    print("no overlap found between " + str(n) + " amplicons!")
    return longest_seq


#Translate and save as dictionary containing header as the key and the translated AA seq as value


#go over the two seqs in a frame of lets 5

amplicon = "TGTTCAAATAAG"
amplicon=seq_and_perc["8_percent"].iloc[1,0]
from itertools import islice

def window(seq, n):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result
window_prev=0
mcherry_chunks = list(window(mcherry_full,10))


from difflib import SequenceMatcher

def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()
#iterate over based on blocks of certain size, also considers cases where the seq length is not divisible by the step
def blocks(mcherry_full, step):
    size = len(mcherry_full)

    for pos in range(0, size - step, step): 
        #print(pos, step)
        #print(pos)
        #print(window)
        yield mcherry_full[pos:pos+step] 

    else: 
        if size % step != 0: 
            # print("not an even divide, the last block:")
            # print(pos)
            # print(pos+size%step)
            yield mcherry_full[pos:pos+size%step+step]

'''
The overlap part is useless as mcherry will map to amplicons mcherry with 100% accuracy! Instead we want to align them, then look at the reference and calculate
in codon pairs till the primer binding site and see if the rest of the seq is in frame or not. if not, then when translating the protein, we will add +1 or +2 etc 
and then translate into the AAs to SEE if the produced Arc after the scar is a functional protein. mcherry will be fine nevertheless!
CGGCGGCA
CGGCGGCA

TGGACGAGCTGTACAAGGt
TGGACGAGCTGTACAAGGTC

GAG ... GGT
'''
#a = (2891-2210) / 3
a = (2919-2210) % 3

prim_ma03_frame=(2906-2210) % 3 #ma03 found in the amplicon
prim_ma03_frame=(2909-2210) % 3 #based on the short 

prim_ma03_frame
a #out of frame by 1
def hash_sequence(string, k):

    dictionary={}
    for i in range(len(string)-(k-1)):
        sequence = string[i:i+k]
        dictionary.setdefault(sequence,[]).append(i)
    return dictionary



def intersects(string1, string2, k): #what if k=0?
    dictionary = hash_sequence(string1, k)

    for i in range(len(string2)-1): #O(n) for sybstring in string

        if string2[i:i+k] in dictionary: #O(n
            return string2[i:i+k]
    return None

def find_overlap(amplicon, mcherry_full):
    longest_seq = None
    n=0
    for i in range(1, min(len(amplicon), len(mcherry_full))):
        # Store the current iteration's intersection
        current_seq = intersects(amplicon, mcherry_full, i)
        
        # If this was one character too long, return.
        # Else, a longer intersection was found. Store it.
        if current_seq == None:
            print("No overlap found!")
            n+=1
            return 0
        else:
            longest_seq = current_seq
    # If we get here, the strings were the same.
    # For consistency, return longest_seq and its length.
    print("no overlap found between " + str(n) + " amplicons!")
    return longest_seq


def find_frame(overlap_bases): #rewrite this based on taking the aligned sequence and then calculate frame by taking the full mcherry
    start_codon_i=mcherry_full.index("ATG")
    if overlap_bases==0:
        return 0
    else:
        seq=mcherry_full[start_codon_i:mcherry_full.index(overlap_bases)]
        if len(seq)%3==0:
            print("In frame!")
            frame_N=0
        else:
            print("out of frame:" + str(len(seq)%3))
            frame_N=len(seq)%3
            return frame_N

aa_and_perc={}

#no intersection found for most of the sequences....how to deal with these cases? this would imply that the cluster seq
#exists but downstream
#issues with aligning with muscle downstream

result="/media/data/AtteR/projects/hiti/mcherry_p3_seq_cluster_all.fasta"

#translate NT sequences to the AAs based on the right codon order. Thus, we need the full mcherry sequence from which we find the overlap between the amplicon and mcherry.
#from this, we calculate the position from full mcherry's start codon till the start of the amplicon where the overlap has occured. If this is divisible by 3, the read is 
#in frame. If not, then start counting from the index where it is. Then ,translate the amplicon seq into AAs  


def translate_nt_aa(input, output_aa, frame):
    with open(input) as nt_seq:
        header=[]
        for i, line in enumerate(nt_seq):
            if line.startswith(">"):
                header.append(line.strip())
            else:
                amplicon=line.strip()
                aa_and_perc[header[0]]=seq=Seq(amplicon[frame:]).translate()
                # overlap=find_overlap(line, mcherry_full)
                # print(overlap)
                # frameN=find_frame(overlap)
                # aa_and_perc[header[0]]=str(Seq(amplicon[frameN:]).translate())
                header=[]
                #print(str(Seq(next(nt_seq)).translate()))
    with open(output_aa, "w") as aa_fasta:
        for id in aa_and_perc.keys():
            aa_fasta.write(id + "\n" + str(aa_and_perc[id]) + "\n")

    return(aa_and_perc)  #save as fasta which you will then align with muscle

    print("done")
output_aa="/media/data/AtteR/projects/hiti/mcherry_p3_seq_clusters_all_AA_inframe.fasta"

str(aa_and_perc[">348 CluSeq_%: 0.0135"])
aa_dict=translate_nt_aa(result,output_aa, 0)
aa_dict


outp="/media/data/AtteR/projects/hiti/mcherry_p3_seq_clusters_all_AA.fasta"
with open(outp, "w") as handle:
    for id in aa_and_perc.keys():
        handle.write(id + "\n" + aa_and_perc[id] + "\n")


#RRHGRAVQGRCCGSAEP--QHRRPDGAGPYDHRRPPRLPCPAGWAGRQTQCDPADW*VPS*DAGTRTEDPPASVDRSVQAGGARAERVAQVGGQAGEQLGRLRAHRRLTALEEVHQGLSLPLPGDHRQPGALGQA*DARVEGGLLPSGE
#5'-3'
#bottom strand
#mcherry_full="CTTGTACAGCTCGTCCATGCCGCCGGTGGAGTGGCGGCCCTCGGCGCGTTCGTACTGTTCCACGATGGTGTAGTCCTCGTTGTGGGAGGTGATGTCCAACTTGATGTTGACGTTGTAGGCGCCGGGCAGCTGCACGGGCTTCTTGGCCTTGTAGGTGGTCTTGACCTCAGCGTCGTAGTGGCCGCCGTCCTTCAGCTTCAGCCTCTGCTTGATCTCGCCCTTCAGGGCGCCGTCCTCGGGGTACATCCGCTCGGAGGAGGCCTCCCAGCCCATGGTCTTCTTCTGCATTACGGGGCCGTCGGAGGGGAAGTTGGTGCCGCGCAGCTTCACCTTGTAGATGAACTCGCCGTCCTGCAGGGAGGAGTCCTGGGTCACGGTCACCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAAGCCCTCGGGGAAGGACAGCTTCAAGTAGTCGGGGATGTCGGCGGGGTGCTTCACGTAGGCCTTGGAGCCGTACATGAACTGAGGGGACAGGATGTCCCAGGCGAAGGGCAGGGGGCCACCCTTGGTCACCTTCAGCTTGGCGGTCTGGGTGCCCTCGTAGGGGCGGCCCTCGCCCTCGCCCTCGATCTCGAACTCGTGGCCGTTCACGGAGCCCTCCATGTGCACCTTGAAGCGCATGAACTCCTTGATGATGGCCATGTTATCCTCCTCGCCCTTGCTCACCAT"
#top
#mcherry_full="CTTGTACAGCTCGTCCATGCCGCCGGTGGAGTGGCGGCCCTCGGCGCGTTCGTACTGTTCCACGATGGTGTAGTCCTCGTTGTGGGAGGTGATGTCCAACTTGATGTTGACGTTGTAGGCGCCGGGCAGCTGCACGGGCTTCTTGGCCTTGTAGGTGGTCTTGACCTCAGCGTCGTAGTGGCCGCCGTCCTTCAGCTTCAGCCTCTGCTTGATCTCGCCCTTCAGGGCGCCGTCCTCGGGGTACATCCGCTCGGAGGAGGCCTCCCAGCCCATGGTCTTCTTCTGCATTACGGGGCCGTCGGAGGGGAAGTTGGTGCCGCGCAGCTTCACCTTGTAGATGAACTCGCCGTCCTGCAGGGAGGAGTCCTGGGTCACGGTCACCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAAGCCCTCGGGGAAGGACAGCTTCAAGTAGTCGGGGATGTCGGCGGGGTGCTTCACGTAGGCCTTGGAGCCGTACATGAACTGAGGGGACAGGATGTCCCAGGCGAAGGGCAGGGGGCCACCCTTGGTCACCTTCAGCTTGGCGGTCTGGGTGCCCTCGTAGGGGCGGCCCTCGCCCTCGCCCTCGATCTCGAACTCGTGGCCGTTCACGGAGCCCTCCATGTGCACCTTGAAGCGCATGAACTCCTTGATGATGGCCATGTTATCCTCCTCGCCCTTGCTCACCAT"

#3'-5'
#bottom
mcherry_full="TACCACTCGTTCCCGCTCCTCCTATTGTACCGGTAGTAGTTCCTCAAGTACGCGAAGTTCCACGTGTACCTCCCGAGGCACTTGCCGGTGCTCAAGCTCTAGCTCCCGCTCCCGCTCCCGGCGGGGATGCTCCCGTGGGTCTGGCGGTTCGACTTCCACTGGTTCCCACCGGGGGACGGGAAGCGGACCCTGTAGGACAGGGGAGTCAAGTACATGCCGAGGTTCCGGATGCACTTCGTGGGGCGGCTGTAGGGGCTGATGAACTTCGACAGGAAGGGGCTCCCGAAGTTCACCCTCGCGCACTACTTGAAGCTCCTGCCGCCGCACCACTGGCACTGGGTCCTGAGGAGGGACGTCCTGCCGCTCAAGTAGATGTTCCACTTCGACGCGCCGTGGTTGAAGGGGAGGCTGCCGGGGCATTACGTCTTCTTCTGGTACCCGACCCTCCGGAGGAGGCTCGCCTACATGGGGCTCCTGCCGCGGGACTTCCCGCTCTAGTTCGTCTCCGACTTCGACTTCCTGCCGCCGGTGATGCTGCGACTCCAGTTCTGGTGGATGTTCCGGTTCTTCGGGCACGTCGACGGGCCGCGGATGTTGCAGTTGTAGTTCAACCTGTAGTGGAGGGTGTTGCTCCTGATGTGGTAGCACCTTGTCATGCTTGCGCGGCTCCCGGCGGTGAGGTGGCCGCCGTACCTGCTCGACATGTTC"
#top
#mcherry_full="ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAG"


#########################################
#here trim the extra seqs from muscle file? do later
muscle_alignment="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_clusters_all_aligned_nonrecseq.fasta"
muscle_alignment

from Bio import SeqIO
from difflib import SequenceMatcher

def similarity(a, b):
    return SequenceMatcher(None, a, b).ratio()
###########################################


#find the highest match
#parse the muscle file, put into dict. once we have specific seq, calculate the similarity score, add to the id key for later

#another option is to make a function that takes in all the seqs and finds the highest match, then trims all the seqs based on this.
#function takes in the dict of muscle alignment file

#no need to change the translation start, just report in the id whether the read is in frame or not

muscle_alignment_AA="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_clusters_all_aligned_AA.fasta"
muscle_alignment = "/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_clusters_all_aligned.fasta"
#Takes in the muscle aligned file, parses through it, calculates the match to the reference seq and saves the file
def muscle_align_postproc(muscle_f,target_sequence):
    f_name=os.path.basename(muscle_f)
    muscle_dict = {}
    for seq_record in SeqIO.parse(muscle_alignment, "fasta"):
        seq=(str(seq_record.seq))
        match_score=similarity(seq,target_sequence)
        id_head=">"+str(seq_record.description) + ", Match_%: " + str(round(match_score*100,4))
        muscle_dict[id_head]=seq
    outp=os.path.dirname(muscle_f) + "/" + f_name.split(".")[0] + "_2.fasta"
    print(outp)
    with open(outp, "w") as handle:
        for id in muscle_dict.keys():
            handle.write(id + "\n" + muscle_dict[id] + "\n")

target_sequence_AA=Seq(target_sequence).translate()
target_sequence_AA

muscle_align_postproc(muscle_alignment_AA,target_sequence_AA)
muscle_align_postproc(muscle_alignment,target_sequence)

muscle_alignment="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_clusters_all_aligned_2.fasta"
os.path.basename(muscle_alignment).split(".")[0]
m=muscle_alignment.split(".")[-2].split("/")[-1]
m

with open(muscle_alignment, "r") as align_file:
    for line in align_file:
        if not line.startswith(">"):
            print(line)
            print("#############")




AAs=nt_seq.translate(table=standard_table)



def count_nts(seq):
    nts=0
    for base in seq:
        if base.isalpha():
            nts+=1
    return nts



######
#TRIM (Do later)
######
#you need to find the highest count but you also need to keep the index so save this info into the same dic. Then we define the highest matching seq's len as 
#the cut off for the rest and thus we cut off the extra from the others as well.

#get all the results and save them into the list, find the highest one along with its index, then trim the lengths of all seqs based on this, then calculate the percentage of
#match relative to the reference template, save into list, put into the df

#if several nts with equally good match, get the one which extends longest and trim seqs based on this one
 
all_nt_counts=[]
for p in align_and_perc.keys():
    N_nts = count_nts(align_and_perc[p])
    all_nt_counts.append(N_nts)
all_nt_counts
highest = max(all_nt_counts)
highest

indices = []
i=0
for n in all_nt_counts:
    if n==highest:
        indices.append(i)
    i+=1
indices  #report this number as it tells how many seqs had this highest match number
keys_list=list(align_and_perc.keys())

#need to retain the index info so that we can get back to the original dict with the percentages and the seq clusters. 
#out of the best mathcing seqs, we take the one that goes the farthest and trim all the seqs based on this one's length
index_best = {}

for seq_i in indices:
    match_seq_perc=keys_list[seq_i] #get the key value based on the index (we retrieve the key for the value that was one of the many longest ones, e.g. when many seqs of length 75 were found)
    best_match = align_and_perc[match_seq_perc]
    index_best[match_seq_perc]=best_match
best_match
#get the longest seq
longest=max(list(index_best.values()), key = len)
longest_seq_key=list(align_and_perc.keys())[list(align_and_perc.values()).index(longest)]

#trim the seqs
align_and_perc_trim = {}
for keys, values in align_and_perc.items():
    align_and_perc_trim[keys]=values[:len(longest)]

target_sequence_trim=target_sequence[:len(longest)]

len(align_and_perc_trim.keys())
len(align_and_perc_trim.values())
df = pd.DataFrame(list(zip(align_and_perc_trim.keys(), align_and_perc_trim.values(), all_nt_counts)), columns =['Percentage_seqs', 'Aligned_seq', 'N(bases)'])

df2 = df.assign(Percent_match = lambda x: x['N(bases)']/len(target_sequence_trim))
del df2["N(bases)"]


#5'-3'
#bottom strand
mcherry_full="CTTGTACAGCTCGTCCATGCCGCCGGTGGAGTGGCGGCCCTCGGCGCGTTCGTACTGTTCCACGATGGTGTAGTCCTCGTTGTGGGAGGTGATGTCCAACTTGATGTTGACGTTGTAGGCGCCGGGCAGCTGCACGGGCTTCTTGGCCTTGTAGGTGGTCTTGACCTCAGCGTCGTAGTGGCCGCCGTCCTTCAGCTTCAGCCTCTGCTTGATCTCGCCCTTCAGGGCGCCGTCCTCGGGGTACATCCGCTCGGAGGAGGCCTCCCAGCCCATGGTCTTCTTCTGCATTACGGGGCCGTCGGAGGGGAAGTTGGTGCCGCGCAGCTTCACCTTGTAGATGAACTCGCCGTCCTGCAGGGAGGAGTCCTGGGTCACGGTCACCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAAGCCCTCGGGGAAGGACAGCTTCAAGTAGTCGGGGATGTCGGCGGGGTGCTTCACGTAGGCCTTGGAGCCGTACATGAACTGAGGGGACAGGATGTCCCAGGCGAAGGGCAGGGGGCCACCCTTGGTCACCTTCAGCTTGGCGGTCTGGGTGCCCTCGTAGGGGCGGCCCTCGCCCTCGCCCTCGATCTCGAACTCGTGGCCGTTCACGGAGCCCTCCATGTGCACCTTGAAGCGCATGAACTCCTTGATGATGGCCATGTTATCCTCCTCGCCCTTGCTCACCAT"
#top
mcherry_full="CTTGTACAGCTCGTCCATGCCGCCGGTGGAGTGGCGGCCCTCGGCGCGTTCGTACTGTTCCACGATGGTGTAGTCCTCGTTGTGGGAGGTGATGTCCAACTTGATGTTGACGTTGTAGGCGCCGGGCAGCTGCACGGGCTTCTTGGCCTTGTAGGTGGTCTTGACCTCAGCGTCGTAGTGGCCGCCGTCCTTCAGCTTCAGCCTCTGCTTGATCTCGCCCTTCAGGGCGCCGTCCTCGGGGTACATCCGCTCGGAGGAGGCCTCCCAGCCCATGGTCTTCTTCTGCATTACGGGGCCGTCGGAGGGGAAGTTGGTGCCGCGCAGCTTCACCTTGTAGATGAACTCGCCGTCCTGCAGGGAGGAGTCCTGGGTCACGGTCACCACGCCGCCGTCCTCGAAGTTCATCACGCGCTCCCACTTGAAGCCCTCGGGGAAGGACAGCTTCAAGTAGTCGGGGATGTCGGCGGGGTGCTTCACGTAGGCCTTGGAGCCGTACATGAACTGAGGGGACAGGATGTCCCAGGCGAAGGGCAGGGGGCCACCCTTGGTCACCTTCAGCTTGGCGGTCTGGGTGCCCTCGTAGGGGCGGCCCTCGCCCTCGCCCTCGATCTCGAACTCGTGGCCGTTCACGGAGCCCTCCATGTGCACCTTGAAGCGCATGAACTCCTTGATGATGGCCATGTTATCCTCCTCGCCCTTGCTCACCAT"

#3'-5'
#bottom
mcherry_full="TACCACTCGTTCCCGCTCCTCCTATTGTACCGGTAGTAGTTCCTCAAGTACGCGAAGTTCCACGTGTACCTCCCGAGGCACTTGCCGGTGCTCAAGCTCTAGCTCCCGCTCCCGCTCCCGGCGGGGATGCTCCCGTGGGTCTGGCGGTTCGACTTCCACTGGTTCCCACCGGGGGACGGGAAGCGGACCCTGTAGGACAGGGGAGTCAAGTACATGCCGAGGTTCCGGATGCACTTCGTGGGGCGGCTGTAGGGGCTGATGAACTTCGACAGGAAGGGGCTCCCGAAGTTCACCCTCGCGCACTACTTGAAGCTCCTGCCGCCGCACCACTGGCACTGGGTCCTGAGGAGGGACGTCCTGCCGCTCAAGTAGATGTTCCACTTCGACGCGCCGTGGTTGAAGGGGAGGCTGCCGGGGCATTACGTCTTCTTCTGGTACCCGACCCTCCGGAGGAGGCTCGCCTACATGGGGCTCCTGCCGCGGGACTTCCCGCTCTAGTTCGTCTCCGACTTCGACTTCCTGCCGCCGGTGATGCTGCGACTCCAGTTCTGGTGGATGTTCCGGTTCTTCGGGCACGTCGACGGGCCGCGGATGTTGCAGTTGTAGTTCAACCTGTAGTGGAGGGTGTTGCTCCTGATGTGGTAGCACCTTGTCATGCTTGCGCGGCTCCCGGCGGTGAGGTGGCCGCCGTACCTGCTCGACATGTTC"
#top
mcherry_full="ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAG"

#top

result="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_clusters_pairwise.fasta"
with open(result, "w") as handle:
#write ref template to the top of the file
    seq_obj = SeqRecord(Seq(target_sequence), id=str(0), description="mcherry_p3_seq_ref")
    count = SeqIO.write(seq_obj, handle, "fasta")
    for i, seq in enumerate(df2.iloc[:,1]):
        descr="ClusterSeq%: " + str(round((df2.iloc[i,0]*100),2)) + "% --- Match: " + str(round((df2.iloc[i,2]*100),2)) + "%"
        seq_obj = SeqRecord(Seq(seq), id=str(id_f), description=descr)
        count = SeqIO.write(seq_obj, handle, "fasta")
        id_f+=1


#are you sure we should use msa as with that we cant align strictly to the reference but it is aligned along with the others?
#







df.assign(Percentage = lambda x: float((all_nt_counts/len(target_sequence_trim), 2)))
df
#find its index. i.e. percentage value
match_seq_perc=keys_list[ind]


ind = all_nt_counts.index(highest)
ind
keys_list=list(align_and_perc.keys())
len(keys_list)
match_seq_perc=keys_list[ind]

best_match = align_and_perc[match_seq_perc]  #a perfect match found, a full length seq... should not be possible????
highest


len(best_match)
len(target_sequence)
#once we have all the aligned seqs, then we go over each seq and count the one with the highest % match based on the one which has fewest "-". Then we define this seq's len as 
#the cut off for the rest and thus we cut off the extra from the others as well.

#If more than one seq with high number of matches, then 


        #this part done later
        #need to transform into seqrecord for downstream visualisation
        seq_align = SeqRecord(Seq(alignm.split("\n")[2]), id=str(a), description="mcherry_p3_alignment")
        seq_record_list.append(seq_align)
result="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_aligned_log.fasta"

with open(result, "w") as handle:
    #place the template seq on top
    count = SeqIO.write(template_seq, handle, "fasta")

    for seq in seq_record_list:
        count = SeqIO.write(seq, handle, "fasta")

#def 
alignments = pairwise2.align.globalms(target_sequence, complete_df.iloc[3,0],  2, -1, -.5, -.1) #localxx earlier
alignments[1]  #each index gives alignment but varying start positions but the score is the same
ss = SeqRecord(Seq((format_alignment(*alignments[0]).split(",")[2])), id=str(1), description="mcherry_p3_alignment")
seq_record_list=[]
for i in range(len(complete_df.index)):
    alignments = pairwise2.align.localxx(target_sequence, complete_df.iloc[i,0], one_alignment_only=True) #generates a signle alignment!!
    for a in range(len(alignments)):  #why iterate over each letter of the single alignment? we get the length of it but it should be 1
        alignm=format_alignment(*alignments[a])
        seq_align = SeqRecord(Seq(alignm.split("\n")[2]), id=str(a), description="mcherry_p3_alignment")
        seq_record_list.append(seq_align)
    result="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_aligned_log.fasta"

    with open(result, "w") as handle:
        #place the template seq on top
        count = SeqIO.write(template_seq, handle, "fasta")

        for seq in seq_record_list:
            count = SeqIO.write(seq, handle, "fasta")

str(alignments)
#iterate over each seq_record_list sequence, transform into dic[id]=str(seq), count number of "-" and the one that has the least
#save its id for the future just in case?
c = str(seq_record_list[4].seq)
c
seq_align



############################

from Bio import Align
aligner = Align.PairwiseAligner()

seq1 = SeqRecord(
    Seq(target_sequence),
    id="1",
    name="mcherry_p3_template",
)
seq2 = SeqRecord(
    Seq(complete_df.iloc[1,0]),
    id="2",
    name="mcherry_p3_template",
)
####################
alignment, score, start_end_positions = local_pairwise_align_ssw(
    DNA(target_sequence),
    DNA(complete_df.iloc[1,0]),
    gap_open_penalty = 3,
    gap_extend_penalty = 1
)
all = local_pairwise_align_ssw(
    DNA(target_sequence),
    DNA(complete_df.iloc[1,0]),
    gap_open_penalty = 3,
    gap_extend_penalty = 1
)
out_align = ('-'*start_end_positions[0][0])+str(alignment[1])+('-'*(len(target_sequence)-start_end_positions[0][1]-1))
out_align #gives the seq areas that did align between the seqs when they have been positioned respective to the ref template?
###########

seq1.id
seq1.name
#alignments[0]
#alignments = aligner.align(target_sequence, complete_df.iloc[1,0])
SeqRecord(Seq(alignments[0]), id="1", name="a")
alignments[9]
#SeqRecord(Seq(alignments[0]), id=seq1[0].id, name=seq1[0].name)


a = format_alignment(*alignments[0])
type(a)
type(alignments[0])
alignments[0].aligned
print(str(alignments[0]))
record = SeqRecord(
    Seq(a),
    id="1",
    name="mcherry_p3",
)
record

for alignment in alignments: 
    print(format_alignment(*alignment)) 
print(format_alignment(*alignments[0]))
alignments_f = format_alignment(*alignments[0])


