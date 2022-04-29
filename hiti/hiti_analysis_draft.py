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
lliteral = ' literal=GGCGGCATGGACGAGC' #added C at the end as its part of tje primer #to check that its on target with mcherry - the primer
rliteral = ' literal=CATATGACCACCGG'
#base_path = '/home/lcadmin/mnm-lts/SequencingFiles/arc_hiti_asyn/HITI1-8977973/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
#export_path = '/home/lcadmin/mnm-lts/HITI-analysis/'
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path = '/media/data/AtteR/projects/hiti/output/'
target_sequence = "GGCGGCATGGACGAGCTGTACAAGGTCGGTGCTGCGGCTCCGCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAGGTCTTCTACCGTCTGGAGAGG"
#removed the first C from the target seq
"Filters and trims the reads"
search_path = base_path+animal_nr+'*'+transgene+'*'+assay_end+'*/'
search_path

# hip=animal_nr + "_" + transgene + "*h_"
# stri=animal_nr + transgene + "*s_"
'''
checks with r2 that it has arc in it - fliteral
in r1 to make sure it has the amplicon - then trim it away (the primer)
the  scar is close to the pcr primer. 
we placed the primer close to the scar - to make sure that sequencing works?

compare the aliments with and without removing the primer and try different alignments

'''
    

from functools import reduce


#you could trim and starcode all individual lane files of certain animal first. after this you sum these based on matches
test_file_p5_filter2 = tempfile.NamedTemporaryFile(suffix = '.fastq').name #when cutadapt applied
test_file_p5_filter2
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
            test_file_p5_filter = tempfile.NamedTemporaryFile(suffix = '.fastq').name#when bbduk applied


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
            test_file_p5_filter2 = tempfile.NamedTemporaryFile(suffix = '.fastq').name #when cutadapt applied

            cutadapt_call="cutadapt -g "+lliteral+" -o " + test_file_p5_filter2 + " " + test_file_p5_filter
            call([cutadapt_call], shell=True)

            print("Cutadapt done! Performed on test_file_p5_filter2: "+ test_file_p5_filter2)
            test_file_p5_out_starcode = tempfile.NamedTemporaryFile(suffix = '.tsv').name
            print("test_file_p5_out_starcode: "+ test_file_p5_out_starcode)
            starcode_call= "/media/data/AtteR/Attes_bin/starcode/starcode -i "+test_file_p5_filter2+" -t 32 -o "+test_file_p5_out_starcode
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
full_df_trim.head()
#full_df_trim.to_csv("/media/data/AtteR/projects/hiti/dfs/full_df_trim_mcherry_p3_exp2_retrimC.csv")
full_df_trim.to_csv("/media/data/AtteR/projects/hiti/dfs/full_df_trim_mcherry_p3_exp2_trimPrimer.csv")

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


result="/media/data/AtteR/projects/hiti/mcherry_p3_seq_cluster_retrim.fasta"

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
            descr="CluSeq:" + str(round(full_df.iloc[seq_i,-3],5)) + "_sd:" + str(round(full_df.iloc[seq_i,-1],5))
            print(descr)
            seq_obj = SeqRecord(Seq(full_df.iloc[seq_i,0]), id=str(id_f), description=descr)
            print(seq_obj)
            count = SeqIO.write(seq_obj, handle, "fasta")
            id_f+=1
    print("Saved!")

#result="/media/data/AtteR/projects/hiti/mcherry_p3_seq_cluster_all_redo.fasta"
result="/media/data/AtteR/projects/hiti/mcherry_p3_seq_cluster_all_retrim.fasta"
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

import re
full_df_trim.columns
#downstream an issue with visualising the seqs using mview is that all the seqs are not SAME length. thus, needs to fix this
def align_and_save(filename, full_df,target_sequence):
    id_f=1
    aligned_data=dict()
    #align all the data, save into dict, then ensure that all the seqs are same length (take the longest seq). IF not, then add padding!
    for seq_i in range(len(full_df.iloc[:,-1])):
            header=">"+ str(id_f)+"CluSeq:" + str((round(full_df.iloc[seq_i,-3],5))) + "_sd:" + str((round(full_df.iloc[seq_i,-1],5)))
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
align_pairwise_loc2="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_aligned_pairwise_local_2sk_retrim.fasta"
align_pairwise_loc3="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_aligned_pairwise_local_3_go3_ge1.fasta"

align_pairwise_glob1="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_aligned_pairwise_glob.fasta"
align_pairwise_glob2="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_aligned_pairwise_glob2.fasta"

musc = "/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_muscle_redo.fasta"

#maybe add progress bars?
align_and_save(align_pairwise_loc3, full_df_trim, target_sequence)
align_and_save(align_pairwise_loc2, full_df_trim, target_sequence)

align_and_save(align_pairwise_glob1, full_df_trim)
align_and_save(align_pairwise_glob2, full_df_trim)

#returns a dict with percentage value of the certain cluster seq and the aligned seq as the value

outp="/media/data/AtteR/projects/hiti/mcherry_p3_seq_clusters_all2.fasta"

'''
Import aligned fasta file - save the header as the key and the seq as the value
'''

out_starcode = tempfile.NamedTemporaryFile(suffix = '.tsv').name
align_pairwise_loc2
starcode_call= "/media/data/AtteR/Attes_bin/starcode/starcode -i "+align_pairwise_loc2+" -t 32 -o "+out_starcode
call([starcode_call], shell=True)
df=pd.read_csv(out_starcode, sep='\t', header=None)
df = df.rename(columns={0: 'sequence', 1:'count'})



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
'''
Make a holistic summary script with colour coded AAs where:
   +1(in-frame)   0       +2
ref |AAAAAAAAAAAAAA|BBBBBBB|CCCCCC
Seq1|AAAAAAAAAAAAAA|BBBBBBB|------
Seq2|AAAAAAAAAAAAAA|-------|CCCCCC

So align the AAs with the ref in each frame, save the coordinates of the match in each case, then merge by making the changes to the ref
as well.
'''
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
from difflib import Differ, SequenceMatcher
dif = Differ()
ref="RRHGRAVQGRCCGSAEPQHRRPDGAGPYDHRRPPRLPCPAGWAGRQTQCDPADW*VPS*DAGTRTEDPPASVDRSVQAGGARAERVAQVGGQAGEQLGRLRAHRRLTALEEVHQGLSLPLPGDHRQPGALGQA*DARVEGGLLPSGE"
amplicon="GGMDELYKVGAAPDGAGPYDHRRP"
df = list(dif.compare(ref , amplicon))
df

l= [ref, amplicon]
match = get_close_matches(ref,amplicon , n=2 , cutoff=0.5)
match
match = SequenceMatcher(None, ref, amplicon).find_longest_match(0, len(ref), 0, len(amplicon))
match
ref[:match.a+match.size]

temp = SequenceMatcher(None,ref ,amplicon)
temp

"RRHGRAVQGRCCGSAEPQHRRPDGAGPYDHRRPPRLPCPAGWAGRQTQCDPADW*VPS*DAGTRTEDPPASVDRSVQAGGARAERVAQVGGQAGEQLGRLRAHRRLTALEEVHQGLSLPLPGDHRQPGALGQA*DARVEGGLLPSGE"
"---------RCCGSA"
ref[:5]
matches=[]
matched_ampls=[]

a="ABCDEFGHIJK"
b="KDCDEFIGHI"
temp = SequenceMatcher(None,a ,b)
seqs=[]
i=temp.get_matching_blocks()[1]
c = len(a[:i.a])*"-"+b[i.b:i.b+i.size] #this is how you will get the coordinates and the length of each matching block.
seqs.append(len(a[:i.a])*"-"+b[i.b:i.b+i.size])
range_line=0
len(a[:i.a])
#get the total length of the ref seq, start going along with it, for each mathcing amplicon's subseq, place them in the right
#place in the ref, the rest add just "--"
c

matched_ampl= len(ref[:match.a])*"-" + str(lines[i+1][match.b:]) + "|"

seqs=[]
range_line=0
for i in range(len(temp.get_matching_blocks())):
    match=temp.get_matching_blocks()[i]
    seqs.append(len(a[range_line:match.a])*"-"+b[match.b:match.b+match.size])
    range_line=match.a+match.size
    #based on the coordinates, extract the start and end, get the seqs between them
seqs
s=''.join(seqs)
a
ABCDEFGHIJK
--CDEFGHI--
s
matched_ampls
ref

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

#function takes in the ref as a string and the amplicon seqs as fasta file as well as the correct reading frame.
#First the ref seq is translated in diff frames and saved into dict with key telling the frame and value showing the translated seq.
#amplicons are saved as AAs into another dict as values with keys containing the header information.
#We make a dict with We find all the matches between amplicons and the ref seq and locate their coordinates and generate an aligned
#sequence this way.
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

#maybe make a DF - rows being headers, then columns for the original amplicon seqs, then columns for the alignments with refs translated in different
#frames

corr_frame=1
input=result
#def translate_nt_aa(input, output_aa, corr_frame):
out_of_frames=[0,1,2]
out_of_frames.remove(corr_frame)
refs_aa_frames={}
output_aa=output_aa.split(".")[0] + "_frameref_" + str(corr_frame) + ".fasta"
aa_and_perc={}


for record in SeqIO.parse(input, "fasta"):
    if record.id=="0":
        refs_aa_frames["Frame:" + str(corr_frame)]=str(Seq(record.seq[corr_frame:]).translate())
        for alt_frame in out_of_frames:
            refs_aa_frames["Frame:" + str(alt_frame)]=str(Seq(record.seq[alt_frame:]).translate())
    else:
        aa_and_perc[">"+str(record.description) + "_transl.frame:" + str(corr_frame)]=str(Seq(record.seq[corr_frame:]).translate())
refs_aa_frames
aa_and_perc
#you go over the ref seqs in different frames and align all the amplicons to them. save the alignment into the df's specific column. initially
#save into a list or such
ref_x_alignment={}
for frame_ref in refs_aa_frames.keys():
    alignments_per_ref=[]
    for ampl in aa_and_perc.keys():
        matches=SequenceMatcher(None,refs_aa_frames[frame_ref],aa_and_perc[ampl])
        seqs=[]
        #you use range_line so that when you fill the remnants from left side of the match, you wont keep adding from
        #beginning since in the end, we merge the seq list elements into the whole alignment of the amplicon against the ref
        range_line=0
        for i in range(len(matches.get_matching_blocks())):
            match=matches.get_matching_blocks()[i]
            seqs.append(len(refs_aa_frames[frame_ref][range_line:match.a])*"-"+str(aa_and_perc[ampl])[match.b:match.b+match.size])
            range_line=match.a+match.size
        alignments_per_ref.append(''.join(seqs))
    ref_x_alignment[frame_ref + "|Ref:" +refs_aa_frames[frame_ref]]=alignments_per_ref
len(ref_x_alignment)
ref_x_alignment.keys()
ref_x_alignment["Frame:1|Ref:AAWTSCTRSVLRLRGAAAPTTRWSWTI*PPAASTPTLPRGVGRPPNPM*SCRLVSAELRCWNTYGGPTGIC*PKCPSRWSES*KGCTGRWASWRTTWTATCPPATHSAGRSPSRPVFAAARRPSPTWSAGSSVRCTCGGRSSTVWR"]
seq_info={"Seq_info:":aa_and_perc.keys()}
keys=list(aa_and_perc.keys())

ref_x_alignment["Frame:1|Ref:AAWTSCTRSVLRLRGAAAPTTRWSWTI*PPAASTPTLPRGVGRPPNPM*SCRLVSAELRCWNTYGGPTGIC*PKCPSRWSES*KGCTGRWASWRTTWTATCPPATHSAGRSPSRPVFAAARRPSPTWSAGSSVRCTCGGRSSTVWR"][0]

seq_info=["Seq_info"]+list(aa_and_perc.keys())
ref_x_alig_list=[]
for keys, values in ref_x_alignment.items():
    print("".join(list((keys))[:]))
    #ref_x_alig_list.append([("".join(list((keys))))]+list(values))
    ref_x_alig_list.append([keys]+list(values))

#df = pd.DataFrame([seq_info, ref_x_alig_list[0],ref_x_alig_list[1],ref_x_alig_list[2]], columns =[seq_info[0],ref_x_alig_list[0][0], ref_x_alig_list[1][0], ref_x_alig_list[2][0]])
df = pd.DataFrame(data= {seq_info[0]: seq_info[1:], ref_x_alig_list[0][0]:ref_x_alig_list[0][1:], ref_x_alig_list[1][0]:ref_x_alig_list[1][1:], ref_x_alig_list[2][0]:ref_x_alig_list[2][1:]})
len(df)

#to show alignments to ref seq that has been translated in frame
df.iloc[:,0:2]
#now you have a df where first col contains the information about the specific amplicon, the cols 1-3 are
#named based on the ref codon frame + the entire translated AA seq

#Go over each alignment col by row, find the longest match and align these next to each other.
#then move onto the next col (ref seq translated using a different frame), find the longest consequtive seq
#in the amplicon and align this to the ref too. repeat with the last col as well.

#when aligning to the ref using e.g. ref1 (in frame) and amplicon 1, cut off the amplicon 1 after
#the longest conseq. seq has been found. Then take the same amplicon but when it has been aligned to the
#ref2 (out of frame), find the longest seq, remove the extra "---" from the beginning relative to the 
#amplicon 1 and then align to the ref. merge the two. 
df
df.columns[1].split(":")[-1]
df.iloc[1,3]

len(df.columns)
match_all = SequenceMatcher(None, df.columns[1].split(":")[-1], df.iloc[1,1])
match_all.get_matching_blocks()

match1 = SequenceMatcher(None, df.columns[1].split(":")[-1], df.iloc[0,1]).find_longest_match(0, len(df.columns[1].split(":")[-1]), 0, len(df.iloc[0,1]))
match1
match2 = SequenceMatcher(None, df.columns[3].split(":")[-1], df.iloc[0,3]).find_longest_match(0, len(df.columns[3].split(":")[-1]), 0, len(df.iloc[0,3]))
match2

#seq_match1 + "-"*len(match2.b-match1.b) + seq_match2
#with the out of frame cols compare which one has the longest seq
df.columns[1].split(":")[-1]
end_of_match=match1.b+match.size
"AAWTSCTRSVLRLRGAAAPTTRWSWTI*PPAASTPTLPRGVGRPPNPM*SCRLVSAELRCWNTYGGPTGIC*PKCPSRWSES*KGCTGRWASWRTTWTATCPPATHSAGRSPSRPVFAAARRPSPTWSAGSSVRCTCGGRSSTVWR"
"AAWTSCTRSVLRLRGA-------------------------------------------------GP----------R--------------------------------------------------------------------"

df.iloc[1,1][match1.b:match1.b+match1.size]
df.columns[3].split(":")[-1]
df.iloc[1,3][match2.b:match2.b+match2.size]

#make a loop where you take the matched_ampl2 for column 2 and 3 (out of frame RFs) and then one with longest
#match, select and merge with the matched_sampl1 
#######

df.columns[1].split(":")[-1]

matched_ampl1= len(df.columns[1].split(":")[-1][:match1.a])*"-" + str(df.iloc[1,1][match1.b:]) + "|"
#matched_ampl2= len(df.columns[3].split(":")[-1][:match.a])*"-" + str(df.iloc[1,3][match.b:]) + "|"
matched_ampl2= len(df.columns[3].split(":")[-1][:match2.a])*"-" + str(df.iloc[1,3][match2.b:]) + "|"

end_of_match=match1.b+match1.size
ref_a=df.columns[1].split(":")[-1][0:end_of_match] + "|"

"AAWTSCTRSVLRLRGPQHRRPDGAGPYDHRRPPRLPCPAGWAGRQTQCDPADW*VPS*DAGTRTEDPPASVDRSVQAGGARAERVAQVGGQAGEQLGRLRAHRRLTALEEVHQGLSLPLPGDHRQPGALGQA*DARVEGGLLPSGE"
"AAWTSCTRSVLRLRG----RPDGAGPYDHR--------------------------------------------------------------------------------------------------------------------|"
a=matched_ampl1[0:end_of_match] + "|"
len(matched_ampl2)
b=matched_ampl2[end_of_match:] + "|"
ref_b=df.columns[3].split(":")[-1][end_of_match:] + "|"
ref_b
a
b
c=a+b
c
#######

"-------------------RPDGAGPYDHR--------------------------------------------------------------------------------------------------------------------|"

from fpdf import FPDF
df.iloc[1,0]

#you still need to make the merged ref seq


'''
G, P, S, T		Orange
H, K, R 		Red
F, W, Y 		Blue
I, L, M, V		Green
'''

    
print any(any('Mary' in s for s in subList) for subList in myDict.values())
from colorama import Fore, Style
color_scheme={"RED": ["H","K","R"], "BLUE":["F","W","Y"], "GREEN":["I","L","M","V"], "YELLOW":["G","P","S","T"]}

def find_colour(aa,color_scheme):
    return [key for key,val in color_scheme.items() if any("P" in s for s in val)]


class write_pdf(FPDF):
    #color_scheme={"RED": ["H","K","R"], "BLUE":["F","W","Y"], "GREEN":["I","L","M","V"], "YELLOW":["G","P","S","T"]}
    #color_scheme={[255,0,0]: ["H","K","R"], [0,0,255]:["F","W","Y"], [0,128,0]:["I","L","M","V"], [255,140,0]:["G","P","S","T"], [255,255,255]: ["-"]}
    color_scheme={"RED": [[255,0,0], ["H","K","R"]], "BLUE":[[[255,0,0]],["F","W","Y"]], "GREEN":[[0,128,0], ["I","L","M","V"]], "ORANGE":[[255,140,0],["G","P","S","T"]], "BLACK": [[255,255,255],["-"]]}

    num=1
    def __init__(self, seq_infos, aligs1, aligs2):
        self.seq_infos=seq_infos
        self.aligs1 = aligs1
        self.aligs2 = aligs2

    def header(self):
        # Arial 12
        self.set_font('Arial', '', 12)
        # Background color
        self.set_fill_color(200, 220, 255)
        # Title
        self.cell(0, 6, 'Seq %d : %s' % (1, self.seq_infos), 0, 1, 'L', 1)
        # Line break
        self.ln(4)

    def find_colour(self):
        col=[key for key,val in color_scheme.items() if any(self.aa in s for s in val)]
        return(color_scheme.keys().index(col))

    def add_colour(self,color_i,aa):
        self.cell(w=0, h = 0, txt = aa, border = 0, ln = 0, align = '', fill = False, link = '')
        self.set_font('Arial', '', 10)
        self.set_text_color(color_scheme[color_i][0], color_scheme[color_i][1], color_scheme[color_i][2])

    def color_AAs(self,color_i, aa):
        self.cell(w=0, h = 0, txt = aa, border = 0, ln = 0, align = '', fill = False, link = '')
        self.set_font('Arial', '', 10)
        self.set_text_color(color_scheme[color_i][0], color_scheme[color_i][1], color_scheme[color_i][2])

    def color_AAs(self):
        pdf.add_page()
        col_seq=[]
        #go over the individual alignemnts
        for info, a1,a2 in zip(self.seq_infos, self.aligs1, self.aligs2):
            #go over the individual AAs of the alignemnts
            for h, aa1,aa2, in zip(info,a1,a2):
                color_i=find_colour(self.seq)
                print(color_i)
                self.cell(w=0, h = 0, txt = h, border = 0, ln = 0, align = '', fill = False, link = '')
                self.set_font('Arial', '', 12)
                color_AAs(color_i,aa1)
            #line break after the whole seq has been coloured
            self.ln()
AAfile="/media/data/AtteR/projects/hiti/mcherry_p3_seq_clusters_all_AA_ALIGN.csv"

match1 = SequenceMatcher(None, df.columns[1].split(":")[-1], df.iloc[0,1]).find_longest_match(0, len(df.columns[1].split(":")[-1]), 0, len(df.iloc[0,1]))
match1
match2 = SequenceMatcher(None, df.columns[3].split(":")[-1], df.iloc[0,3]).find_longest_match(0, len(df.columns[3].split(":")[-1]), 0, len(df.iloc[0,3]))
match2
matched_ampl1= len(df.columns[1].split(":")[-1][:match1.a])*"-" + str(df.iloc[1,1][match1.b:]) + "|"
matched_ampl2= len(df.columns[3].split(":")[-1][:match2.a])*"-" + str(df.iloc[1,3][match2.b:]) + "|"

seqinfos=[]
aligs_merged_seq=[]
aligs_merged_ref=[]
for ampl_row in range(round(len(df.index)*0.2)):
    match1 = SequenceMatcher(None, df.columns[1].split(":")[-1], df.iloc[ampl_row,1]).find_longest_match(0, len(df.columns[1].split(":")[-1]), 0, len(df.iloc[ampl_row,1]))
    end_of_match=match1.b+match.size
    print(match1)
    matched_ampl1= len(df.columns[1].split(":")[-1][:match1.a])*"-" + str(df.iloc[ampl_row,1][match1.b:]) + "|"
    seq_inframe=matched_ampl1[0:end_of_match] + "|"
    print(seq_inframe)
    ref_inframe1=df.columns[1].split(":")[-1][0:end_of_match] + "|"
    #get the matches with the other reading frames (i.e. the last 2 cols of the df)
    match2_a = SequenceMatcher(None, df.columns[2].split(":")[-1], df.iloc[ampl_row,2]).find_longest_match(0, len(df.columns[2].split(":")[-1]), 0, len(df.iloc[ampl_row,2]))
    match2_b = SequenceMatcher(None, df.columns[3].split(":")[-1], df.iloc[ampl_row,3]).find_longest_match(0, len(df.columns[3].split(":")[-1]), 0, len(df.iloc[ampl_row,3]))
    if match2_a.size>match2_b.size:
        matched_ampl2= len(df.columns[2].split(":")[-1][:match2_a.a])*"-" + str(df.iloc[ampl_row,2][match2_a.b:]) + "|"
        ref_outframe2=df.columns[2].split(":")[-1][0:end_of_match] + "|"
        # seq_inframe2=matched_ampl2[0:end_of_match]
        seq_outframe=matched_ampl2[end_of_match:]
        print("match2_a larger than b!")
    else:
        matched_ampl2= len(df.columns[3].split(":")[-1][:match2_b.a])*"-" + str(df.iloc[ampl_row,3][match2_b.b:]) + "|"
        ref_outframe2=df.columns[3].split(":")[-1][end_of_match:] + "|"
        seq_outframe=matched_ampl2[end_of_match:]
        # seq_inframe2=matched_ampl2[0:end_of_match]
        seq_outframe=matched_ampl2[end_of_match:]
        print("match2_b larger than a!")

    merged_align_seq=seq_inframe+seq_outframe
    # ref_inframe=matched_ampl1[0:end_of_match]
    # ref_outframe=matched_ampl2[end_of_match:]
    print("inframe1:" + seq_inframe)

    print("inframe2:" + seq_outframe)
    print("merged:" + merged_align_seq)
    merged_align_ref=ref_inframe1+ref_outframe2
    header_info=df.iloc[ampl_row,0]
    seqinfos.append(df.iloc[ampl_row,0])
    aligs_merged_seq.append(merged_align_seq)
    aligs_merged_ref.append(merged_align_ref)

# pdf=write_pdf(seqinfos,aligs_merged_seq,aligs_merged_ref)
# pdf.add_page()
# pdf.color_AAs()
# pdf.header(header_info)
# pdf.color_AAs(merged_align_ref)
# pdf.color_AAs(merged_align_seq)
# pdf.output('/media/data/AtteR/projects/hiti/tuto3.pdf', 'F')

########
########

#FIX
def find_colour(aa,color_scheme):
    col=[key for key,val in color_scheme.items() if any(aa in s for s in val)]
    return("".join(col))

color_scheme={"RED": [[255,0,0], ["H","K","R"]], "BLUE":[[255,0,0],["F","W","Y"]], "GREEN":[[0,128,0], ["I","L","M","V"]], "ORANGE":[[255,140,0],["G","P","S","T"]], "BLACK": [[255,255,255],["-", "|"]]}
ex_s="--------------GAAAPTTRWSWTI--|"

pdf=FPDF(format='letter', unit='in')
color_scheme["RED"][0][0]
color_i = find_colour("H", color_scheme)
print(color_i)

pdf.add_page()
# Set color red
pdf.set_font('Arial', '', 10)
pdf.set_text_color(color_scheme[color_i][0][0], color_scheme[color_i][0][1], color_scheme[color_i][0][2])
pdf.cell(80)

pdf.cell(0, 0, txt = "H",border = 0, ln = 0)
color_i = find_colour("F", color_scheme)
color_scheme[color_i][0][1]
pdf.set_font('Arial', '', 10)
pdf.set_text_color(color_scheme[color_i][0][0], color_scheme[color_i][0][1], color_scheme[color_i][0][2])
pdf.cell(0, 0, txt = "F",border = 3, ln = 0)

# # Colors of frame, background and text
# pdf.set_draw_color(0, 80, 180)
# pdf.set_fill_color(230, 230, 0)
out="/media/data/AtteR/projects/hiti/tuto3.pdf"
out="/media/data/AtteR/projects/hiti/tuto3.rtf"

pdf.output(out,'F')
template = '{\\rtf1\\ansi\\ansicpg1252\\deff0\\nouicompat\\deflang1033{\\fonttbl{\\f0\\fnil\\fcharset0 Times New Roman;}}\n{\\colortbl;\\red255\\green0\\blue0;\\red0\\green0\\blue255;\\red0\\green255\\blue0;\\red255\\green0\\blue255;}'
with open(out, 'w') as f:
    f.write(template)
###########
style = """<style type='text/css'>
html {
  font-family: Courier;
}
r {
  color: #ff0000;
}
g {
  color: #00ff00;
}
b {
  color: #0000ff;
}
</style>"""

RED = 'r'
GREEN = 'g'
BLUE = 'b'

def write_html(f, type, str_):
    f.write('<%(type)s>%(str)s</%(type)s>' % {
            'type': type, 'str': str_ } )

f = open(out, 'w')
f.write('<html>')
f.write(style)

write_html(f, RED, 'My name is so foo..\n')
write_html(f, BLUE, '102838183820038.028391')

f.write('</html>')
#######

'''
<html>

    <body>

        <h1>Heading</h1>

    </body>

</html>
'''
def find_colour(aa,color_scheme):
    col=[key for key,val in color_scheme.items() if any(aa in s for s in val)]
    out="".join(col)
    if out=='':
        out="BLACK"
        return
    else:
        return(out)

def aa_coloring(seq_infos, aligs1, aligs2, output):
    # color_scheme={"RED": [[255,0,0], ["H","K","R"]], "BLUE":[[[255,0,0]],["F","W","Y"]], "GREEN":[[0,128,0], ["I","L","M","V"]], "ORANGE":[[255,140,0],["G","P","S","T"]], "BLACK": [[255,255,255],["-"]]}
    #using lesk color code scheme:
    color_scheme={"RED": [["#FF0000"], ["D", "E"]], "BLUE":[["#6495ED"],["K", "R"]], "GREEN":[["#9FE2BF"], ["C", "V", "I", "L", "P", "F", "Y", "M", "W"]], "ORANGE":[["#FF7F50"],["G", "A", "S", "T"]], "MAGENTA":[["#DE3163"], ["N", "Q", "H"]], "BLACK": [["#000000"],["-","|"]]}
    f = open(output,"w")
    for info, a1,a2 in zip(seq_infos, aligs1, aligs2):
        print(info)
        a1,a2=a1[0:round(len(a1)*0.6)], a2[0:round(len(a2)*0.6)]
        f.write("<html> \n <body> <p>"+info +"</p> <p>")
    #go over the individual AAs of the alignemnts
        #write the html colouring into a string var for each aa and append into a list, then join and make a massive string,
        # then after all AAs have been iterated over for the certain seq, then write into the html 
        aa1_list=[]
        aa2_list=[]
        for aa1,aa2, in zip(a1,a2):
            aa2_list.append('<span style="color:'+ color_scheme[find_colour(aa2,color_scheme)][0][0]+ '">' + aa2 + '</span>')
            aa1_list.append('<span style="color:'+ color_scheme[find_colour(aa1,color_scheme)][0][0]+ '">' + aa1 + '</span>')
        coloured_ref="".join(aa2_list)
        coloured_seq="".join(aa1_list)
        f.write(coloured_ref +"<br>")
        f.write(coloured_seq)
        f.write("</p>")

    f.write("</body>")
    f.close()
out="/media/data/AtteR/projects/hiti/tuto3.html"
aa_coloring(seqinfos,aligs_merged_seq,aligs_merged_ref,out)

q="--------------GAAAPTTRWSWTI--|T-----------------------------------------------------------------------------RW--S----------TI-----------------------|"

q[17]
i=1
for aa in q:
    color_scheme[find_colour(aa,color_scheme)][0]
    i+=1
    i


<link rel="stylesheet" href="myCs325Style.css" type="text/css"/>
from fpdf import FPDF
 
# Create instance of FPDF class
# Letter size paper, use inches as unit of measure
pdf=FPDF(format='letter', unit='in')
# Add new page. Without this you cannot create the document.
pdf.add_page()
# Set font face to Times, size 10.0 pt
pdf.set_font('Times','',10.0)
 
# Set color red
pdf.set_text_color(255,0,0)    
pdf.cell(1.0,0.0,'Hello World!')
# Line break 0.15 inches height
pdf.ln(0.15)
 
# Set color green
pdf.set_text_color(0,255,0)    
pdf.cell(1.0,0.0,'Hello World!')
pdf.ln(0.15)
 
# Set color blue
pdf.set_text_color(0,0,255)    
pdf.cell(1.0,0.0,'Hello World!')
pdf.ln(0.15)
 
# output content into a file ('F') named 'hello3.pdf'
pdf.output(out,'F')

def aa_coloring(seq_infos, aligs1, aligs2, output):
    color_scheme={"RED": [[255,0,0], ["H","K","R"]], "BLUE":[[[255,0,0]],["F","W","Y"]], "GREEN":[[0,128,0], ["I","L","M","V"]], "ORANGE":[[255,140,0],["G","P","S","T"]], "BLACK": [[255,255,255],["-"]]}

    pdf=FPDF("P", "mm", "A4")
    # add a page 
    pdf.add_page()
    # Set font: Times, normal, size 10
    pdf.set_font('Times','', 12)
    for info, a1,a2 in zip(seq_infos, aligs1, aligs2):
    #go over the individual AAs of the alignemnts
        for h, aa1,aa2, in zip(info,a1,a2):
            #write in the header of certain seq alignment
            pdf.cell(w=0, h = 0, txt = h, border = 0, ln = 0, align = '', fill = False, link = '')
            pdf.set_font('Arial', h, 12)
            print(aa1)
            #write in the ref seq that has been merged with in-and out of frame versions 
            color_i=find_colour(aa1,color_scheme)
            print(color_i)
            pdf.cell(0, 2, txt = aa1, border = 0, ln = 0, align = '', fill = False, link = '')
            pdf.set_font('Arial', '', 10)
            pdf.set_text_color(color_scheme[color_i][0], color_scheme[color_i][1], color_scheme[color_i][2])
            #write in the aligned seq that has been merged with in-and out of frame versions 
            color_i=find_colour(aa2,color_scheme)
            print(color_i)
            pdf.cell(0, 2, txt = aa2, border = 0, ln = 0, align = '', fill = False, link = '')
            pdf.set_font('Arial', '', 10)
            pdf.set_text_color(color_scheme[color_i][0], color_scheme[color_i][1], color_scheme[color_i][2])

    #line break after the whole seq has been coloured
        pdf.ln()

    pdf.output(output,'F')
aligs_merged_seq
out="/media/data/AtteR/projects/hiti/tuto3.pdf"
aa_coloring(seqinfos,aligs_merged_seq,aligs_merged_ref, out)

class pdf():
    color_scheme={"RED": [[255,0,0], ["H","K","R"]], "BLUE":[[[255,0,0]],["F","W","Y"]], "GREEN":[[0,128,0], ["I","L","M","V"]], "ORANGE":[[255,140,0],["G","P","S","T"]], "BLACK": [[255,255,255],["-"]]}

    pdf=FPDF("P", "mm", "A4")
    # add a page 
    pdf.add_page()
    # Set font: Times, normal, size 10
    pdf.set_font('Times','', 12)
    pdf.output(out,'F')


    num=1
    def __init__(self, seq_infos, aligs1, aligs2, output):
        self.seq_infos=seq_infos
        self.aligs1 = aligs1
        self.aligs2 = aligs2
        self.output=output

    def color_AAs(self,color_i, aa):
        self.cell(0, 2, txt = aa, border = 0, ln = 0, align = '', fill = False, link = '')
        self.set_font('Arial', '', 10)
        self.set_text_color(color_scheme[color_i][0], color_scheme[color_i][1], color_scheme[color_i][2])

    def color_AAs(self):
        pdf.add_page()
        col_seq=[]
        #go over the individual alignemnts
        for info, a1,a2 in zip(self.seq_infos, self.aligs1, self.aligs2):
            #go over the individual AAs of the alignemnts
            for h, aa1,aa2, in zip(info,a1,a2):
                color_i=find_colour(self.seq)
                print(color_i)
                self.cell(w=0, h = 0, txt = h, border = 0, ln = 0, align = '', fill = False, link = '')
                self.set_font('Arial', '', 12)
                color_AAs(color_i,aa1)
            #line break after the whole seq has been coloured
            self.ln()

#make a function that calls the fpdf and the rest INSIDE it

pdf=pdf()
a=pdf._save('/media/data/AtteR/projects/hiti/tuto3.pdf')


'''
if 
    self.set_text_color(220, 50, 50)
    self.set_text_color(color_scheme[0][0], color_scheme[0][1], color_scheme[0][2])

'''
color_scheme={[220, 50, 50]: ["H","K","R"], "BLUE":["F","W","Y"], "GREEN":["I","L","M","V"], "YELLOW":["G","P","S","T"]}

def color_AAs(seq):
    col_seq=[]
    color_scheme={"RED": ["H","K","R"], "BLUE":["F","W","Y"], "GREEN":["I","L","M","V"], "YELLOW":["G","P","S","T"]}
    for aa in seq:
        col=find_colour(aa)
        if not col:  #set line colour
            coloured=f"{Style.RESET_ALL}aa"
        else:
            coloured=f"{Fore.col}" aa
        print(coloured)
colors = {'reset':'\033[0m', 'blue':'\033[34m'}
color_asci={"RED": "\033[1;31;40m"}
print ("\033[91m" + 'I know I can be red' + "\033[0m")
person = 'you'
formattext = 'How are %s%s%s' %(colors['blue'], person, colors['reset'])
my_str = f"{Fore.BLUE}Hello, {Style.RESET_ALL} guys. {Fore.RED} I should be red."
my_str
import termcolor
string = "type-name-function-location"
string = string.replace('-', termcolor.colored('-', 'red'))
string
end_of_match

matched_ampl= len(df.columns[1].split(":")[-1][:match.a])*"-" + str(df.iloc[1,1][match.b:]) + "|"

ref_col_i=1
for ref_col_i in range(len(df.columns())):
    for ampl_row in range(len(df.index)):
        match = SequenceMatcher(None, df.columns[ref_col_i].split(":")[-1], df.iloc[ampl_row,ref_col_i]).find_longest_match(0, len(df.columns[ref_col_i].split(":")[-1]), 0, len(df.iloc[ampl_row,ref_col_i]))



#ref_x_alignment key contains info about the frame used to translate the ref seq along with the translated ref seq, the values are amplicon alignments and
#these become columns in the final df
#now we have the dict containing the ref AA seqs translated in different frames as dict and the amplicon AAs with header infos as another dict
    

   # return(aa_and_perc)  #save as fasta which you will then align with muscle


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


