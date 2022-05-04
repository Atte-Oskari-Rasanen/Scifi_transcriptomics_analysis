from difflib import Differ, SequenceMatcher
from Bio.Data import CodonTable
from Bio.Seq import Seq
import re
from Bio.SubsMat import MatrixInfo as matlist
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
from functools import reduce
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
import statistics as st

def calculate_perc_sd(full_df):
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


def import_reads_process(data_dict, transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd):
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
    #full_df_trim=calculate_perc_sd(full_df)
    return(complete_df)



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


def save_fasta(filename, full_df, target_sequence):
    id_f = 1
    ref="Ref_seq"
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

import re
################
#ALIGNMENTS
from Bio.SubsMat import MatrixInfo as matlist
Bio.Align.substitution_matrices

#########
#ALIGNMENT CLASSES TO USE
#########
class align_local():
    aligned_data=dict()
    def __init__(self, amplicon, target_sequence,gop=-3, gep=-1):
        self.amplicon=amplicon
        self.target_sequence=target_sequence
        self.gop=gop
        self.gep=gep

    def align(self):
        alignments = pairwise2.align.localxx(self.target_sequence, self.amplicon)

        #alignments = pairwise2.align.globalms(target_sequence, seq_and_perc[group][:,0],  2, -1, -.5, -.1)
        alignm=format_alignment(*alignments[0])
        #make string objects, save into a list. then count which one has least ----, cut off the rest of the reference based on this?
        seq_align = alignm.split("\n")[2]
        return(seq_align)

class align_local2():
    aligned_data=dict()
    def __init__(self, amplicon,ref,gop=-3, gep=-1):
        self.amplicon=amplicon
        self.ref=ref
        self.gop=gop
        self.gep=gep

    def align(self):
        alignment, score, start_end_positions = local_pairwise_align_ssw(DNA(self.ref),DNA(self.amplicon),gap_open_penalty=self.gop,gap_extend_penalty = self.gep)
        out_align = ('-'*start_end_positions[0][0])+str(alignment[1])+('-'*(len(self.ref)-start_end_positions[0][1]-1))
        return(out_align)

class align_local3():
    aligned_data=dict()
    def __init__(self, amplicon, ref,gop=-3, gep=-1):
        self.amplicon=amplicon
        self.ref=ref
        self.gop=gop
        self.gep=gep
    def align(self):
        alignments = pairwise2.align.localms(self.ref, self.amplicon, 2, -1, self.gop, self.gep)
        #alignments = pairwise2.align.localms(self.target_sequence, self.amplicon, 2, -1, -.5, -.1)
        #alignments= pairwise2.align.localds(self.target_sequence, self.amplicon, 2, -1, -.5, -.1)
        alignm=format_alignment(*alignments[0])
        seq_align = alignm.split("\n")[2]
        return(seq_align)
class align_global():
    aligned_data=dict()

    def __init__(self, amplicon, ref,gop=-3, gep=-1):
        self.amplicon=amplicon
        self.ref=ref
        self.gop=gop
        self.gep=gep
    def align(self):
        alignments = pairwise2.align.globalms(self.ref, self.amplicon,  2, -1, -.5, -.1)
        alignm=format_alignment(*alignments[0])
        seq_align = alignm.split("\n")[2]
        return(seq_align)
    #nt_count=count_nts(seq_align)
    #seq_and_perc[group]["match"]
class align_global2():
    aligned_data=dict()

    def __init__(self, amplicon, ref,gop=-3, gep=-1):
        self.amplicon=amplicon
        self.ref=ref
        self.gop=gop
        self.gep=gep
    def align(self):
        alignments = pairwise2.align.globalxx(self.ref, self.amplicon)
        alignm=format_alignment(*alignments[0])
        #make string objects, save into a list. then count which one has least ----, cut off the rest of the reference based on this?
        seq_align = alignm.split("\n")[2]
        return(seq_align)

import inspect

def align_trimmer(aligned_data,target_sequence):
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
    return(aligned_data)
def write_align(aligned_data, filename, target_sequence):
    with open(filename, "w") as handle:
        header=">0"+" Ref_seq"
        handle.write(header + "\n" + target_sequence + "\n")
        for seq_i in aligned_data.keys():
            handle.write(seq_i + "\n" + aligned_data[seq_i] + "\n")
import subprocess
def bash_command(cmd):
    subprocess.Popen(cmd, shell=True, executable='/bin/bash')

#takes in the df and the choice of the alignment method. methods are found in class
#the class must be instantiated inside the function and the appropriate method is called
#by index passed by the user into function
def aligner(full_df, target_sequence, align_method, filename, output_path, gop=3, gep=1):
    align_class = {"align_local": align_local,
            "align_local2": align_local2,
            "align_local3":align_local3,
            "align_global":align_global,
            "align_global2":align_global2}  
    id_f=1
    aligner_init = align_class.get(str(align_method), None)  # Get the chosen class, or None if input is bad
    aligner_init
    aligned_data=dict()
    #align all the data, save into dict, then ensure that all the seqs are same length (take the longest seq). IF not, then add padding!
    for seq_i in range(len(full_df.iloc[:,-1])):
        #yield iteratively the header of the certain seq and the corresponding seq
        header=">"+ str(id_f)+"CluSeq:" + str((round(full_df.iloc[seq_i,-3],5))) + "_sd:" + str((round(full_df.iloc[seq_i,-1],5)))
        #seq_obj_align = aligner_init(full_df.iloc[seq_i,0], target_sequence, gop, gep).align()

        print("===SEQ===:" + full_df.iloc[seq_i,0])
        seq_obj_align = aligner_init(full_df.iloc[seq_i,1], target_sequence, gop, gep).align()

        seq_obj_align = re.sub(r'[(\d|\s]', '', seq_obj_align) #remove digits from the string caused by the alignment and empty spaces from the start
        aligned_data[header]=seq_obj_align
        id_f+=1
    aligned_data_trim=align_trimmer(aligned_data, target_sequence)
    write_align(aligned_data_trim, filename, target_sequence)
    #Generate a visual alignment file using mview
    mview_file=output_path + "/" + filename.split("/")[-1].split(".")[-2] + ".html"
    mview_command='/media/data/AtteR/Attes_bin/mview -in fasta -html head -css on -coloring any ' + filename + '>' + mview_file
    call([mview_command], shell=True)

    #subprocess.run(['/media/data/AtteR/Attes_bin/mview', '-in fasta', '-html head', '-css on', '-coloring any', filename, '>', mview_file])


