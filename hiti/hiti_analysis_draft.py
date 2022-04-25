from curses import window
from email.mime import base
from fnmatch import translate
import ntpath
from re import L
from wsgiref import headers
from nbdev.showdoc import *
import os
import warnings
import matplotlib.pyplot as plt
from subprocess import call
import glob
import pandas as pd
import tempfile
import numpy as np
from skbio.alignment import *
from skbio import DNA

import os
os.getcwd()
os.chdir("/media/data/AtteR/projects/hiti")
from alignment_scripts import *
import Bio.Align.Applications; dir(Bio.Align.Applications)
from Bio.Align.Applications import MuscleCommandline#Read in unfiltered data
from Bio import AlignIO
from Bio import pairwise2
from Bio import SeqIO
from Bio.pairwise2 import *
from alignment_scripts import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
stri=animal_nr + transgene + "*s_"
    

from functools import reduce


#you could trim and starcode all individual lane files of certain animal first. after this you sum these based on matches

#Takes in all the reads based on lane, trims them, combines the lanes, sums the seq counts, calculates the percentage of 
#reads being the certain sequence (starcode cluster) out of all the reads
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

            stats_out = export_path+animal_group_name+'_'+transgene+'_'+assay_end+'_stats-filter.txt'

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


def create_datadict(base_path):

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

    return(data_dict)

# after generation of the full df containing all animals, the following function calculates the stats (sd, mean,
# perc. sum per seq cluster and raw counts sum per seq cluster)
def postprocess_data(full_df):
    import statistics as st
    #animal group contains all lanes of the certain data
    full_df = full_df.fillna(value=0)
    perc_cols = [col for col in full_df.columns if 'percent' in col]
    count_cols = [col for col in full_df.columns if 'count' in col]

    perc_cols
    #sum the percentages of each seq cluster of each animal together to see how much the certain seq is found 
    full_df['percent_sum_unit'] = full_df[perc_cols].sum(axis=1)  

    total_perc_unit=full_df.iloc[:,-1].sum()
    full_df['percent_sum'] = (full_df['percent_sum_unit'] / total_perc_unit)
    full_df.head()
    count_cols
    #full_df['sd']=full_df[count_cols].std()
    full_df['total_reads_seq'] = full_df[count_cols].sum(axis=1)  

    full_df['sd']=full_df[perc_cols].std(axis=1)
    full_df['mean']=full_df[perc_cols].mean(axis=1)

    #remove sequences that have 0-3 reads in total across groups

    #get the total number of percentages from the percent_sum column and then divide all the perc units with this
    #to get the percs
    #calculate the SD for each seq
    full_df.sort_values(by=['percent_sum_unit'], ascending=False, inplace=True)
    full_df.head()
    full_df.columns
    #discard seqs that contribute less than 0.0001% percentage
    rows_drop=[]
    full_df[count_cols]
    full_df_trim = full_df.drop(full_df[full_df["total_reads_seq"] <= 3].index)
    return(full_df_trim)


#after which we merge the ones that have matching numbers
# complete_df = pd.DataFrame({'sequence': ['CTGTACAAGGTCGGTGCTGCGGCTCCGCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGAC']})
# complete_df
data_dict=create_datadict(base_path)

#you could trim and starcode all individual lane files of certain animal first. after this you sum these based on matches
full_df = trimRead_hiti(data_dict,transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd)
#the function returns you a complete dataframe containing animals_x_brain_area (h,s). since each animal_x_brain_area contains data from 4 different subdirs, these have been
# summed into the same ones, i.e. column 12_6h contains the sequences from lanes 1-4 all summed up  
#now take the percentage values of each animal (so brain area s and h), merge into same column so that the percentages will be a sum of the two.

full_df_trim = postprocess_data(full_df)

full_df_trim.to_csv("/media/data/AtteR/projects/hiti/dfs/full_df_trim_mcherry_p3_exp1.csv")

#transforms the aligned seqs into seqrecord objects


#extract the aligned seq part, ignoring the scores etc. 
#to process and save the files downstream, must convert the result into seq record
#first take specific alignment from the alignment object via indexing, then format it into a string, split based on line changes
#and take the 3rd element, i.e the NTs matching the template one, i.e. the alignment. Transform this into seq object, save
#objects into a list and then iteratively save into fasta file. visualise them
#try different alignment methods

'''
class for df stuff

instantiate

method:gain df
method:trim df
method:save appropriate fasta file
method:import the fasta for analysis
'''


#we need to know what percentage of the clusters mapped with high accuracy, ideally 100%, to the reference genome. Thus
#we must save the percentages from the starcluster data
############################
############################
#seq_and_perc = {}

#full_df.iloc[:, [0,-1]]


#We take the seq and percent sum for each column and then align. the full_df contains ALL the clustered seqs from all the animals. we have taken the percentage amount of raw counts
#for each animal and then we have summed them together, giving as a unit providing information how much of the certain seq cluster is when summed over all the animals.
#After we map the sequences to the reference and find the alignments, we can use this as information regarding how many of the amplicons were of certain kind and
#how well it maps to the ref. 


result="/media/data/AtteR/projects/hiti/mcherry_p3_seq_cluster.fasta"

#add proper delimiters 

full_df_trim.loc[1, "percent_sum"]
id_f=1
full_df_trim.columns
full_df_trim.iloc[0,-1]
full_df_trim.iloc[0,-2]
full_df_trim.iloc[0,-3]
full_df_trim.iloc[0,-4]

full_df_trim["percent_sum"]
full_df_trim.head()
len(full_df_trim.index)



def save_fasta(filename, full_df, target_sequence):
    id_f = 1
    ref="mcherry_p3_seq_ref"
    with open(filename, "w") as handle:
        seq_obj = SeqRecord(Seq(target_sequence), id=str(0), description=ref)
        count = SeqIO.write(seq_obj, handle, "fasta")
        for seq_i in range(len(full_df.index)):
            print("seq_i:" + str(seq_i))
            descr="CluSeq:" + str(round(full_df.iloc[seq_i,-4],5)) + "_sd:" + str(round(full_df.iloc[seq_i,-2],5))
            print(descr)
            seq_obj = SeqRecord(Seq(full_df.iloc[seq_i,0]), id=str(id_f), description=descr)
            print(seq_obj)
            count = SeqIO.write(seq_obj, handle, "fasta")
            id_f+=1
    print("Saved!")

#result="/media/data/AtteR/projects/hiti/mcherry_p3_seq_cluster_all_redo.fasta"
result="/media/data/AtteR/projects/hiti/mcherry_p3_seq_cluster_all.fasta"
save_fasta(result, full_df_trim, target_sequence)
aligned_seqs = []


###############
#ALIGNMENTS
from Bio.SubsMat import MatrixInfo as matlist
Bio.Align.substitution_matrices

class align():
    def __init__(self, amplicon, target_sequence,gop=3, gep=1):
        self.amplicon=amplicon
        self.target_sequence=target_sequence
        self.gop=gop
        self.gep=gep

    def align_local(self):
        alignments = pairwise2.align.localxx(self.target_sequence, self.amplicon)
        #alignments = pairwise2.align.globalms(target_sequence, seq_and_perc[group][:,0],  2, -1, -.5, -.1)

        alignm=format_alignment(*alignments[0])
        #make string objects, save into a list. then count which one has least ----, cut off the rest of the reference based on this?
        seq_align = alignm.split("\n")[2]
        #nt_count=count_nts(seq_align)
        #seq_and_perc[group]["match"]
        return(seq_align)

    def align_local2(self):
        alignment, score, start_end_positions = local_pairwise_align_ssw(DNA(self.target_sequence),DNA(self.amplicon),gap_open_penalty =self.gop,gap_extend_penalty = self.gep)
        out_align = ('-'*start_end_positions[0][0])+str(alignment[1])+('-'*(len(target_sequence)-start_end_positions[0][1]-1))
        return out_align

    def align_local3(self):
        alignments = pairwise2.align.localms(self.target_sequence, self.amplicon, 2, -1, self.gop, self.gep)
        #alignments = pairwise2.align.localms(self.target_sequence, self.amplicon, 2, -1, -.5, -.1)
        #alignments= pairwise2.align.localds(self.target_sequence, self.amplicon, 2, -1, -.5, -.1)
        alignm=format_alignment(*alignments[0])
        seq_align = alignm.split("\n")[2]
        #nt_count=count_nts(seq_align)
        #seq_and_perc[group]["match"]
        return(seq_align)
    def align_global(self):
        alignments = pairwise2.align.globalms(self.target_sequence, self.amplicon,  2, -1, -.5, -.1)
        alignm=format_alignment(*alignments[0])
        seq_align = alignm.split("\n")[2]
        #nt_count=count_nts(seq_align)
        #seq_and_perc[group]["match"]
        return(seq_align)
    def align_global2(self):
        alignments = pairwise2.align.globalxx(self.target_sequence, self.amplicon)

        alignm=format_alignment(*alignments[0])
        #make string objects, save into a list. then count which one has least ----, cut off the rest of the reference based on this?
        seq_align = alignm.split("\n")[2]
        #nt_count=count_nts(seq_align)
        #seq_and_perc[group]["match"]
        return(seq_align)

full_df_trim.head()
full_df_trim.columns

full_df_trim.loc[1,"percent_sum"]
import re
#downstream an issue with visualising the seqs using mview is that all the seqs are not SAME length. thus, needs to fix this
def align_and_save(filename, full_df,target_sequence):
    id_f=1
    aligned_data=dict()
    #align all the data, save into dict, then ensure that all the seqs are same length (take the longest seq). IF not, then add padding!
    for seq_i in range(len(full_df.iloc[:,-1])):
            header=">"+ str(id_f)+"CluSeq:" + str((round(full_df.iloc[seq_i,-4],5))) + "_sd:" + str((round(full_df.iloc[seq_i,-2],5)))
            align_inst=align(full_df.iloc[seq_i,0], target_sequence)
            seq_obj_1=align_inst.align_local2()
            #seq_obj_1= align_local(full_df.iloc[seq_i,0], target_sequence)
            seq_obj_1 = re.sub(r'[(\d|\s]', '', seq_obj_1) #remove digits from the string caused by the alignment and empty spaces from the start
            aligned_data[header]=seq_obj_1
            id_f+=1

    #longest_seq=longest_align(aligned_data.values())
    for id in aligned_data.keys():
        if len(aligned_data[id])==len(target_sequence):
            continue
        if len(aligned_data[id])>len(target_sequence):
            N_dashes=len(target_sequence)-len(aligned_data[id])
            aligned_data[id]=aligned_data[id][:N_dashes]
            print("Seq length larger than ref by " + str(N_dashes) + " ... \n After removal length: " + str(len(aligned_data[id][:N_dashes])))
        else:
            N_dashes=len(target_sequence)-len(aligned_data[id])
            aligned_data[id]=aligned_data[id]+ N_dashes*"-"
            print("Seq length smaller than ref by " + str(N_dashes) + " ... \n After addition length: " + str(len(aligned_data[id])))
    with open(filename, "w") as handle:
        header=">0"+" mcherry_p3_seq_ref"
        handle.write(header + "\n" + target_sequence + "\n")
        for seq_i in aligned_data.keys():
            handle.write(seq_i + "\n" + aligned_data[seq_i] + "\n")
#we want to show if these are alignments that can happen via framshift but still be in frame

#we have deletions for sure, not sure about insertions

'''
RRHGRAVQGRCCGSAEPQHRRPDGAGPYDHRRPPRLPCPAGWAGRQTQCDPADW*VPS*DAGTRTEDPPASVDRSVQAGGARAERVAQVGGQAGEQLGRLRAHRRLTALEEVHQGLSLPLPGDHRQPGALGQA*DARVEGGLLPSGE


'''
#feed in the predefined aligment file 

align_pairwise_loc1="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_aligned_pairwise_local_1.fasta"
align_pairwise_loc2="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_aligned_pairwise_local_2sk.fasta"
align_pairwise_loc3="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_aligned_pairwise_local_3_go3_ge1.fasta"

align_pairwise_glob1="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_aligned_pairwise_glob.fasta"
align_pairwise_glob2="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_aligned_pairwise_glob2.fasta"

musc = "/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_muscle_redo.fasta"
from Bio import SeqIO

for record in SeqIO.parse(musc, "fasta"):
    print(record.id)
    print(record.seq)
    print("==================")

#maybe add progress bars?
align_and_save(align_pairwise_loc3, full_df_trim, target_sequence)
align_and_save(align_pairwise_loc2, full_df_trim, target_sequence)

align_and_save(align_pairwise_glob1, full_df_trim)
align_and_save(align_pairwise_glob2, full_df_trim)

#returns a dict with percentage value of the certain cluster seq and the aligned seq as the value

outp="/media/data/AtteR/projects/hiti/mcherry_p3_seq_clusters_all2.fasta"


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
result="/media/data/AtteR/projects/hiti/mcherry_p3_seq_cluster_all.fasta"

nt_seq=SeqIO.parse(result, "fasta")

#no intersection found for most of the sequences....how to deal with these cases? this would imply that the cluster seq
#exists but downstream
#issues with aligning with muscle downstream

result="/media/data/AtteR/projects/hiti/mcherry_p3_seq_cluster_all.fasta"
# result='/media/data/AtteR/projects/hiti/mcherry_p3_seq_cluster_all_redo.fasta'
#translate NT sequences to the AAs based on the right codon order. Thus, we need the full mcherry sequence from which we find the overlap between the amplicon and mcherry.
#from this, we calculate the position from full mcherry's start codon till the start of the amplicon where the overlap has occured. If this is divisible by 3, the read is 
#in frame. If not, then start counting from the index where it is. Then ,translate the amplicon seq into AAs  

'''
GGMDELYKVGAAAPRSRSTDDQMELDHMTTGGLHAYPAPRGGPAAKPNVILQIGKCRAEMLEHVRRTHRHLLTEVSKQVERELKGLHRSVGKLENNLDGYVPTGDSQRWKKSIKACLCRCQETIANLERWVKREMHVWREVFYRLER
GGMDELYKVGAAAPTRWSWTI*PP
-this one encodes just mcherry
---translate the seqs in diff frames, align to the ref, visualise (0,2). 
-the out of frame ones will have several stop codons but we should be able to see Arc when we try
different frames.(which reading frame is the arc read in!)
-the frame change only applied to the ref, then align

sds, delimiters, so that its easier to extract the percentages and delimiters

-diff reading frame 
'''

nt_seq.seq
for record in SeqIO.parse(result, "fasta"):
    print(record.id)

def translate_nt_aa(input, output_aa, frame_ampl, frame_ref):
    output_aa=output_aa.split(".")[0] + "_frameref_" + str(frame_ref) + ".fasta"
    aa_and_perc={}
    for record in SeqIO.parse(input, "fasta"):
        if record.id=="0":
            aa_and_perc[">"+str(record.description) +"_transl.frame:" + str(frame_ref)]=Seq(record.seq[frame_ref:]).translate()
        else:
            aa_and_perc[">"+str(record.description) + "_transl.frame:" + str(frame_ampl)]=Seq(record.seq[frame_ampl:]).translate()
    with open(output_aa, "w") as aa_fasta:
        for id in aa_and_perc.keys():
            aa_fasta.write(id + "\n" + str(aa_and_perc[id]) + "\n")

    return(aa_and_perc)  #save as fasta which you will then align with muscle

    print("done")
output_aa="/media/data/AtteR/projects/hiti/mcherry_p3_seq_clusters_all_AA_inframe_redo.fasta"
output_aa="/media/data/AtteR/projects/hiti/mcherry_p3_seq_clusters_all_AA"

#output_aa="/media/data/AtteR/projects/hiti/translated/mcherry_p3_seq_clusters_all_AA_inframe_aligned_pw_loc2_go5_ge3.fasta"
result
#str(aa_and_perc[">348 CluSeq_%: 0.0135"])

#ONLY translate the ref seq with different frames
aa_dict=translate_nt_aa(result,output_aa, 1, 1)
aa_dict



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

def find_overlap(amplicon, ref):
    longest_seq = None
    n=0
    for i in range(1, min(len(amplicon), len(ref))):
        # Store the current iteration's intersection
        current_seq = intersects(amplicon, ref, i)
        
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

#Get the full ref seq. then with the amplicon, go over it AA by AA and compare to the ref. If match found start another loop
#where you start counting how many similar AAs, starting from match are found. Keep count of on which AA this match was found as 
#later on you will fill the empty space on the left with it



#for each AA, save info about aa index, the name and the index of match with the ref (if any)
#after the information has been saved for each AA, then you iterate over it and find the part with stretch of matches
#then you go back to the original amplicon, remove everythin before this match and replace with "-", then save these
#as overlap
aa_info=[]
ans = 0;

ref="RRHGRAVQGRCCGSAEPQHRRPDGAGPYDHRRPPRLPCPAGWAGRQTQCDPADW*VPS*DAGTRTEDPPASVDRSVQAGGARAERVAQVGGQAGEQLGRLRAHRRLTALEEVHQGLSLPLPGDHRQPGALGQA*DARVEGGLLPSGE"
amplicon="GGMDELYKVGAAPDGAGPYDHRRP"
match = SequenceMatcher(None, ref, amplicon).find_longest_match(0, len(ref), 0, len(amplicon))
ref[:match.a+match.size]

temp = SequenceMatcher(None,ref ,amplicon)
temp
print(temp.get_matching_blocks())

from difflib import SequenceMatcher
import csv

#the first line in the input AA file must be the ref seq
def align_AAs(inputAA):
    outputAA=inputAA.split(".")[0] + "_aligned.csv"
    with open(outputAA, 'w') as f:
        with open(AAfile) as file:
            lines = [line.rstrip("\n") for line in file]
            for i, line in enumerate(lines):
                if ">0" in line:
                    ref=lines[i+1]
                if ">" in line and not ">0" in line:
                    seq_info=line
                    amplicon=lines[i+1]
                    #find all matches with size greater than 3, assign to the seq
                    
                    match = SequenceMatcher(None, ref, lines[i+1]).find_longest_match(0, len(ref), 0, len(lines[i+1]))
                    
                    matched_ampl= len(ref[:match.a])*"-" + str(lines[i+1][match.b:]) + "|"
                    # create the csv writer
                    writer = csv.writer(f)
                    writer.writerow([seq_info])
                    writer.writerow(["Ref | "+ ref[:match.a+match.size]])
                    writer.writerow(["Seq | "+ matched_ampl])
AAfile="/media/data/AtteR/projects/hiti/mcherry_p3_seq_clusters_all_AA_frameref_0.fasta"
AAfile="/media/data/AtteR/projects/hiti/mcherry_p3_seq_clusters_all_AA_frameref_2.fasta"
AAfile="/media/data/AtteR/projects/hiti/mcherry_p3_seq_clusters_all_AA_frameref_1.fasta"

align_AAs(AAfile)

#write into html instead....


'''
modify the above function to create a summary file: a script that calculates number of seqs with aligned mcherry - lets say 8AAs match from start,
then it counts the number of seqs with matching arc WHEN ref is 1) in frame (+1), 2) out of frame (0,+2)

the number of AAs
Sum up the seqs with matching mcherry, lets say if theres 8AA match, then sum up the appropriate percentages.
Do the same with Arc
Summarise 

Ilmoita kuinka suuri osa sekvenseistä sisältää arkin
--- täten voimme todeta onko yhdistelmä onnistunut vai ei, mutta jos arkki ei ole raamissaan niin sitä
ei voida tuottaa solussa. tehdäänkö arvelle mitään?

'''
