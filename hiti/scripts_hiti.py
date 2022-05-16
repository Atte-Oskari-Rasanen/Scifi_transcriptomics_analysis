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
import math

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
def save_fasta(filename, full_df, target_sequence):
    id_f = 1
    ref="Ref"
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

            test_file_p5_filter2 = tempfile.NamedTemporaryFile(suffix = '.fastq').name #when cutadapt applied on 5'

            #cutadapt_call="cutadapt -g "+lliteral+" -o " + test_file_p5_filter2 + " " + test_file_p5_filter
            cutadapt_call="cutadapt -g "+lliteral+" " + test_file_p5_filter + " > " + test_file_p5_filter2

            call([cutadapt_call], shell=True)
            test_file_p5_filter3 = tempfile.NamedTemporaryFile(suffix = '.fastq').name #when cutadapt applied on 3'

            #cutadapt_call="cutadapt -a "+rliteral+" -o " + test_file_p5_filter2 + " " + test_file_p5_filter2
            cutadapt_call="cutadapt -a "+rliteral+" " + test_file_p5_filter2 + " > " + test_file_p5_filter3

            call([cutadapt_call], shell=True)

            print("Cutadapt done! Performed on test_file_p5_filter2: "+ test_file_p5_filter3)
            test_file_p5_out_starcode = tempfile.NamedTemporaryFile(suffix = '.tsv').name
            print("test_file_p5_out_starcode: "+ test_file_p5_out_starcode)
            starcode_call= "/media/data/AtteR/Attes_bin/starcode/starcode -i "+test_file_p5_filter3+" -t 32 -o "+test_file_p5_out_starcode
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
    return(complete_df)
from functools import reduce

def import_reads_process_mini(base_path, ref,filterlitteral,lliteral,rliteral,read_fwd):
    complete_df = pd.DataFrame({'sequence': [ref]})
    df_animal=[]
    seq_animal=[]
    for read in os.listdir(base_path):
        animal_group_name=read.split("_")[3] + "_" + read.split("_")[4]
        print(animal_group_name)
        if "R1" in read:
            animal_p5_cat = tempfile.NamedTemporaryFile(suffix = '.fastq.gz').name
            animal_p7_cat = tempfile.NamedTemporaryFile(suffix = '.fastq.gz').name
            test_file_p5_out = tempfile.NamedTemporaryFile(suffix = '.fastq').name
            test_file_p7_out = tempfile.NamedTemporaryFile(suffix = '.fastq').name
            test_file_p5_filter = tempfile.NamedTemporaryFile(suffix = '.fastq').name#when bbduk applied


            if read_fwd:
                animal_p5 = glob.glob(base_path+'/'+read)
                animal_p7 = glob.glob(base_path+'/' + read.replace("R1","R2"))
                #display('Forward run Animal: '+animal_nr)
            else:
                print("===")
                break
                # animal_p5 = glob.glob(search_path+'/*R2*')
                # animal_p7 = glob.glob(search_path+'/*R1*')
                #display('Reverse run Animal: '+animal_nr)

            cat_p5= "cat "+" ".join(animal_p5)+" > "+animal_p5_cat
            print(cat_p5)
            #os.system(cat_p5)
            call([cat_p5], shell=True) #call caused the terminal to freeze so switched to os
            cat_p7= "cat "+" ".join(animal_p7)+" > "+animal_p7_cat
            call([cat_p7], shell=True)
            #os.system(cat_p7)

            kmer = '20'
            hdist = '3'
            param=" k="+kmer+" hdist="+hdist+" rcomp=f skipr2=t threads=32 overwrite=true"

            #to check if the read is an amplicon
            call_sequence = "/media/data/AtteR/Attes_bin/bbmap/bbduk.sh in="+animal_p7_cat+" in2="+animal_p5_cat+" outm1="+test_file_p7_out+" outm2="+test_file_p5_out+" literal="+filterlitteral+param
            call([call_sequence], shell=True)
            #actual trimming
            call_sequence = "/media/data/AtteR/Attes_bin/bbmap/bbduk.sh in="+test_file_p5_out+" out="+test_file_p5_filter+ " literal=AAAAAAAAA,CCCCCCCCC,GGGGGGGGG,TTTTTTTTT k=9 mm=f overwrite=true minlength=40"
            call([call_sequence], shell=True)
            test_file_p5_filter2 = tempfile.NamedTemporaryFile(suffix = '.fastq').name #when cutadapt applied on 5'

            #cutadapt_call="cutadapt -g "+lliteral+" -o " + test_file_p5_filter2 + " " + test_file_p5_filter
            cutadapt_call="cutadapt -g "+lliteral+" " + test_file_p5_filter + " > " + test_file_p5_filter2

            call([cutadapt_call], shell=True)
            test_file_p5_filter3 = tempfile.NamedTemporaryFile(suffix = '.fastq').name #when cutadapt applied on 3'

            #cutadapt_call="cutadapt -a "+rliteral+" -o " + test_file_p5_filter2 + " " + test_file_p5_filter2
            cutadapt_call="cutadapt -a "+rliteral+" " + test_file_p5_filter2 + " > " + test_file_p5_filter3

            call([cutadapt_call], shell=True)

            print("Cutadapt done! Performed on test_file_p5_filter2: "+ test_file_p5_filter3)
            test_file_p5_out_starcode = tempfile.NamedTemporaryFile(suffix = '.tsv').name
            print("test_file_p5_out_starcode: "+ test_file_p5_out_starcode)
            starcode_call= "/media/data/AtteR/Attes_bin/starcode/starcode -i "+test_file_p5_filter3+" -t 32 -r 5 -o "+test_file_p5_out_starcode
            call([starcode_call], shell=True)

            df=pd.read_csv(test_file_p5_out_starcode, sep='\t', header=None)
            df = df.rename(columns={0: 'sequence', 1:'count'})
            total_counts = int(df[['count']].sum())
            df['percent'] = (df['count'] / total_counts)
            df = df.rename(columns={'percent':animal_group_name+'_percent','count':animal_group_name+'_count',})
            df_animal.append(df)
    #we iterate over all the individual dfs and merge them by taking the seq column of all dfs and placing them under the new dfs seq and do the same with counts
    df_all_lanes=reduce(lambda  left,right: pd.merge(left,right,on=['sequence'], how='outer'), df_animal)
    #reduce is useful when you need to apply a function to an iterable and reduce it to a single cumulative value.
    return(df_all_lanes)

#the identical cluster seqs may be carryovers from different DFs. starcode probs find same clusters
#from the different files. Thus the script needs to merge the dfs and merge the perc and count values of
#seqs where they are identical.
'''
t_l=[]
data1 = [['tom', 10], ['juli', 14], ['atte', 2], ['santeri', 6], ['julle', 4]]
d1 = pd.DataFrame(data1, columns = ['sequence', 'Age'])
data2 = [['ee', 10], ['poo', 14], ['tom', 2], ['eiole', 6], ['slumpmässig', 4]]
d2 = pd.DataFrame(data2, columns = ['sequence', 'Age'])
t_l.append(d1)
t_l.append(d2)
t_l
dfall=reduce(lambda  left,right: pd.merge(left,right,on=['sequence'], how='outer'), t_l)
dfall=dfall.fillna(0)
dfall['Age_sum']=dfall['Age_x']+dfall['Age_y']

df2 = dfall.groupby(['sequence'],as_index=False).agg({'Age': 'sum'})
'''
def create_datadict(base_path, transgene):

    #hip_folders = [folder for folder in os.listdir(base_path) if "mCherry" in folder and "h_" in folder or "s_" in folder]
    group_folders = [folder for folder in os.listdir(base_path) if transgene in folder]
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
        alignments = pairwise2.align.localms(str(self.ref), str(self.amplicon), 2, -1, self.gop, self.gep)
        #alignments = pairwise2.align.localms(self.target_sequence, self.amplicon, 2, -1, -.5, -.1)
        #alignments= pairwise2.align.localds(self.target_sequence, self.amplicon, 2, -1, -.5, -.1)
        try:
            alignm=format_alignment(*alignments[0])
            seq_align = alignm.split("\n")[2]
            return(seq_align)

        except IndexError:
            return(None)

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


def find_remove_duplicates(aligned_data_trim):
    id_keys_dict=dict()
    algs=list(aligned_data_trim.values())

    for i, id in enumerate(aligned_data_trim.keys()):
        id_keys_dict[i]=id

    #go over the indeces where the duplicates of a seq found, two at a time and sum them all up together
    def remove_multiples(index_pos, id_keys_dict):
        perc_sum_all=0
        var_sum_all=0
        sd_sum_all=0
        for i in index_pos:
            perc_sum_all+=float(id_keys_dict[i].split("_")[1].split(":")[1])
            var_sum_all+=float(id_keys_dict[i].split("_")[2].split(":")[1])
            sd_sum_all+=float(id_keys_dict[i].split("_")[3].split(":")[1])
        #fix the sd summing!
        merged_id=">"+str(i) + "CluSeq:"+str(round(perc_sum_all,5))+ "_var:"+str(round(var_sum_all,5)) +"_sd:"+str(round(math.sqrt(var_sum_all),5))
        return(merged_id)

    #gets the indices of the matching seqs
    def get_indeces_matching_seqs(list_of_elems, element):
        index_pos_list = []
        index_pos = 0
        while True:
            try:
                # Search for item in list from indexPos to the end of list
                index_pos = list_of_elems.index(element, index_pos)
                # Add the index position in list
                index_pos_list.append(index_pos)
                index_pos += 1
            except ValueError as e:
                break
        return index_pos_list

    def remove_key_from_dict(index_pos,id_keys_dict,aligned_data_trim):
        for pos in index_pos:
            #id_keys_dict[pos] value corresponds to the key of the alignments 
            print(id_keys_dict[pos])
            print(aligned_data_trim[id_keys_dict[pos]])
            del aligned_data_trim[id_keys_dict[pos]]
        return(aligned_data_trim)
    import collections
    #merged values
    dups = [item for item, count in collections.Counter(algs).items() if count >= 2]
    merged_multiples=dict()
    for dupl_seq in dups:
        index_pos = get_indeces_matching_seqs(algs, dupl_seq)
        #take the seq from any of the matches
        #seq=algs[index_pos[0]]
        #print(dupl_seq[:40] + " found at pos " + str(index_pos))
        aligned_data_trim=remove_key_from_dict(index_pos,id_keys_dict,aligned_data_trim)
        #key includes the merged matches values and the value the seq
        merged_multiples[remove_multiples(index_pos, id_keys_dict)]=algs[index_pos[0]]
    #now we have removed the seqs that have duplicates and effectively merged them together into merged_multiples dict. This
    #is then added back to the original dict from which the duplicates had been removed.
    merged_multiples.update(aligned_data_trim)
    return(merged_multiples)
#once you get the indices of the duplicates, get the keys based on these indices, merge the key values
#Now we have a script for detecting duplicates from the data dict, removing them via merging the values
#need to remove the seqs from the original dict using the index_pos approach
def addprimers(fasta_dict, filename,lliteral,rliteral):
    with open(filename, "w") as f: 
        for id, seq_prims in fasta_dict.items():
            f.write(id + "\n")
            f.write(lliteral.split("=")[1] + "|" + seq_prims + "|" + rliteral.split("=")[1] + "\n")




#takes in the df and the choice of the alignment method. methods are found in class
#the class must be instantiated inside the function and the appropriate method is called
#by index passed by the user into function
def aligner(full_df, target_sequence, align_method, filename, output_path, lliteral, rliteral gop=3, gep=1):
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
        header=">"+ str(id_f)+"_CluSeq:" + str((round(full_df.iloc[seq_i,-4],5))) + "_var:"+str((round(full_df.iloc[seq_i,-2],5))) +"_sd:" + str((round(full_df.iloc[seq_i,-1],5)))
        #seq_obj_align = aligner_init(full_df.iloc[seq_i,0], target_sequence, gop, gep).align()

        seq_obj_align = aligner_init(full_df.iloc[seq_i,0], target_sequence, gop, gep).align()
        seq_obj_align = re.sub(r'[(\d|\s]', '', seq_obj_align) #remove digits from the string caused by the alignment and empty spaces from the start
        aligned_data[header]=seq_obj_align
        id_f+=1

    aligned_data_trim=align_trimmer(aligned_data, target_sequence)
    data_trim_nodupl=find_remove_duplicates(aligned_data_trim)
    #Add primers to both ends of the seq and save
    #write_align(data_trim_nodupl, filename, target_sequence)
    addprimers(data_trim_nodupl, filename, lliteral,rliteral)

    #Generate a visual alignment file using mview
    mview_file=output_path + "/" + filename.split("/")[-1].split(".")[-2] + ".html"
    mview_command='/media/data/AtteR/Attes_bin/mview -in fasta -html head -css on -coloring any ' + filename + '>' + mview_file
    call([mview_command], shell=True)
    print("html file created as "+ mview_file)
    #os.system('/media/data/AtteR/Attes_bin/mview -in fasta -html head -css on -coloring any {} > {}'.format(str(filename), str(mview_file))) 
    #subprocess.run(['/media/data/AtteR/Attes_bin/mview', '-in fasta', '-html head', '-css on', '-coloring any', filename, '>', mview_file])

def translate_nt_aa(result, corr_frame):
    out_of_frames=[0,1,2]
    #out_of_frames.remove(corr_frame) dont take out any of the frames as the seqs after the scar may be in frame
    refs_aa_frames={}
    aa_and_perc={}

    len(aa_and_perc)
    for record in SeqIO.parse(result, "fasta"):
        if record.id=="0":
            refs_aa_frames["Frame:" + str(corr_frame)]=str(Seq(record.seq[corr_frame:]).translate())
            for alt_frame in out_of_frames:
                refs_aa_frames["Frame:" + str(alt_frame)]=str(Seq(record.seq[alt_frame:]).translate())
        else:
            aa_and_perc[">"+str(record.description) + "_transl.frame:" + str(corr_frame)]=str(Seq(record.seq[corr_frame:]).translate())

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
    seq_info={"Seq_info:":aa_and_perc.keys()}
    keys=list(aa_and_perc.keys())


    seq_info=["Seq_info"]+list(aa_and_perc.keys())
    ref_x_alig_list=[]
    for keys, values in ref_x_alignment.items():
        #print("".join(list((keys))[:]))
        #ref_x_alig_list.append([("".join(list((keys))))]+list(values))
        ref_x_alig_list.append([keys]+list(values))

    df = pd.DataFrame(data= {seq_info[0]: seq_info[1:], ref_x_alig_list[0][0]:ref_x_alig_list[0][1:], ref_x_alig_list[1][0]:ref_x_alig_list[1][1:], ref_x_alig_list[2][0]:ref_x_alig_list[2][1:]})
    return(df)
def translate_nt_aa_hiti2(result, corr_frame, output_html):
    refs_aa_frames={}

    aa_and_perc={}
    for record in SeqIO.parse(result, "fasta"):
        if record.id=="0":
            refs_aa_frames["Frame:" + str(corr_frame)]=str(Seq(record.seq[corr_frame:]).translate())
        else:
            aa_and_perc[">"+str(record.description) + "_transl.frame:" + str(corr_frame)]=str(Seq(record.seq[corr_frame:]).translate())
    ref_x_alignment={}
    alignments_per_ref=[]
    for ampl in aa_and_perc.keys():
        matches=SequenceMatcher(None,refs_aa_frames["Frame:" + str(corr_frame)],aa_and_perc[ampl])
        seqs=[]
        #you use range_line so that when you fill the remnants from left side of the match, you wont keep adding from
        #beginning since in the end, we merge the seq list elements into the whole alignment of the amplicon against the ref
        range_line=0
        for i in range(len(matches.get_matching_blocks())):
            match=matches.get_matching_blocks()[i]
            seqs.append(len(refs_aa_frames["Frame:" + str(corr_frame)][range_line:match.a])*"-"+str(aa_and_perc[ampl])[match.b:match.b+match.size])
            range_line=match.a+match.size
        alignments_per_ref.append(''.join(seqs))
    ref_x_alignment["Frame:" + str(corr_frame) + "|Ref:" +refs_aa_frames["Frame:" + str(corr_frame)]]=alignments_per_ref
    seq_info={"Seq_info:":aa_and_perc.keys()}
    keys=list(aa_and_perc.keys())
    seq_info=["Seq_info"]+list(aa_and_perc.keys())
    ref_x_alig_list=[]
    for keys, values in ref_x_alignment.items():
        #make into a list with first being the ref, the rest being the aligned seqs. 
        ref_x_alig_list.append([keys]+list(values))
    #so we have a list containing sublist pairs of the ref seq
    #ref_x_alig_list[0][0]
    df = pd.DataFrame(data= {seq_info[0]: seq_info[1:], ref_x_alig_list[0][0]:ref_x_alig_list[0][1:]})
    
    f = open(output_html,"w")

    color_scheme={"RED": [["#FF0000"], ["D", "E"]], "BLUE":[["#6495ED"],["K", "R"]], "GREEN":[["#9FE2BF"], ["C", "V", "I", "L", "P", "F", "Y", "M", "W"]], "ORANGE":[["#FF7F50"],["G", "A", "S", "T"]], "MAGENTA":[["#DE3163"], ["N", "Q", "H"]], "BLACK": [["#000000"]]}

    for info, a1,a2 in zip(df["Seq_info"], list(df.iloc[:,1]), [df.columns[1].split(":")[2]]*len(list(df.iloc[:,1]))):
        aa1_list=[]
        aa2_list=[]

        a1,a2=a1[0:round(len(a1)*0.6)], a2[0:round(len(a2)*0.6)]
        f.write("<html> \n <body> <p>"+ "_".join(info.split("_")[:1]) +"</p> <p>")
        #test_html.append('<!DOCTYPE html><html> <head><link rel="stylesheet" href="format.css"></head> <meta charset="utf-8"> <body> <p>'+ "_".join(info.split("_")[:2]) + "<br>" + frame +'</p> <p>')
    #go over the individual AAs of the alignemnts
        #write the html colouring into a string var for each aa and append into a list, then join and make a massive string,
        # then after all AAs have been iterated over for the certain seq, then write into the html 
        for aa1,aa2, in zip(a1,a2):
            print("aa1: "+ aa1)
            print("aa2: "+ aa2)
            # try:
            #     aa2_list.append('<span style="color:'+ color_scheme[find_colour(aa2,color_scheme)][0][0]+ '">' + aa2 + '</span>')
            # except KeyError:
            #     print("Keyerror with aa2 being: " + aa2)
            # try:
            #     aa1_list.append('<span style="color:'+ color_scheme[find_colour(aa1,color_scheme)][0][0]+ '">' + aa1 + '</span>')
            # except KeyError:
            #     print("Keyerror with aa1 being: " + aa1)
            
            aa2_list.append('<span style="color:'+ color_scheme[find_colour(aa2,color_scheme)][0][0]+ '">' + aa2 + '</span>')
            aa1_list.append('<span style="color:'+ color_scheme[find_colour(aa1,color_scheme)][0][0]+ '">' + aa1 + '</span>')

            #print("============")
        coloured_ref="".join(aa2_list)

        coloured_seq="".join(aa1_list)

        f.write(coloured_ref +"<br>")
        f.write(coloured_seq)
        f.write("</p>")

    f.write("</body></html>")
    f.close()

    return(df)

#to show alignments to ref seq that has been translated in frame

def find_colour(aa,color_scheme):
    col=[key for key,val in color_scheme.items() if any(aa in s for s in val)]
    out="".join(col)
    if out==None or not col:
        out="BLACK"
        return(out)
    else:
        return(out)

def aa_coloring(seq_infos, aligs1, aligs2, frame_info, output):
    #using lesk color code scheme:
    #color_scheme={"RED": [["#FF0000"], ["D", "E"]], "BLUE":[["#6495ED"],["K", "R"]], "GREEN":[["#9FE2BF"], ["C", "V", "I", "L", "P", "F", "Y", "M", "W"]], "ORANGE":[["#FF7F50"],["G", "A", "S", "T"]], "MAGENTA":[["#DE3163"], ["N", "Q", "H"]], "BLACK": [["#000000"],["-","|", "*"]]}
    color_scheme={"RED": [["#FF0000"], ["D", "E"]], "BLUE":[["#6495ED"],["K", "R"]], "GREEN":[["#9FE2BF"], ["C", "V", "I", "L", "P", "F", "Y", "M", "W"]], "ORANGE":[["#FF7F50"],["G", "A", "S", "T"]], "MAGENTA":[["#DE3163"], ["N", "Q", "H"]], "BLACK": [["#000000"]]}
    print(output)
    f = open(output,"w")
    for info, a1,a2, frame in zip(seq_infos, aligs1, aligs2,frame_info):
        aa1_list=[]
        aa2_list=[]

        print(frame)
        a1,a2=a1[0:round(len(a1)*0.6)], a2[0:round(len(a2)*0.6)]
        #f.write("<html> \n <body> <p>"+ "_".join(info.split("_")[:1]) +"</p> <p>")
        f.write('<!DOCTYPE html><html> <head><link rel="stylesheet" href="format.css"></head> <meta charset="utf-8"> <body> <p>'+ "_".join(info.split("_")[:2]) + "<br>" + frame +'</p> <p>')

    #go over the individual AAs of the alignemnts
        #write the html colouring into a string var for each aa and append into a list, then join and make a massive string,
        # then after all AAs have been iterated over for the certain seq, then write into the html 
        for aa1,aa2, in zip(a1,a2):
            # try:
            #     aa2_list.append('<span style="color:'+ color_scheme[find_colour(aa2,color_scheme)][0][0]+ '">' + aa2 + '</span>')
            # except KeyError:
            #     print("Keyerror with aa2 being: " + aa2)
            # try:
            #     aa1_list.append('<span style="color:'+ color_scheme[find_colour(aa1,color_scheme)][0][0]+ '">' + aa1 + '</span>')
            # except KeyError:
            #     print("Keyerror with aa1 being: " + aa1)
            
            aa2_list.append('<span style="color:'+ color_scheme[find_colour(aa2,color_scheme)][0][0]+ '">' + aa2 + '</span>')
            aa1_list.append('<span style="color:'+ color_scheme[find_colour(aa1,color_scheme)][0][0]+ '">' + aa1 + '</span>')

            #print("============")
        coloured_ref="".join(aa2_list)

        coloured_seq="".join(aa1_list)

        f.write(coloured_ref +"<br>")
        f.write(coloured_seq)
        f.write("</p>")

    f.write("</body></html>")
    f.close()
out="/media/data/AtteR/projects/hiti/AA_3p_mcherry_coloured_alignment.html"


#now you have a df where first col contains the information about the specific amplicon, the cols 1-3 are
#named based on the ref codon frame + the entire translated AA seq

#Go over each alignment col by row, find the longest match and align these next to each other.
#then move onto the next col (ref seq translated using a different frame), find the longest consequtive seq
#in the amplicon and align this to the ref too. repeat with the last col as well.

#when aligning to the ref using e.g. ref1 (in frame) and amplicon 1, cut off the amplicon 1 after
#the longest conseq. seq has been found. Then take the same amplicon but when it has been aligned to the
#ref2 (out of frame), find the longest seq, remove the extra "---" from the beginning relative to the 
#amplicon 1 and then align to the ref. merge the two. 

def visualise_aa_hybrid_alignments(df, output_html): #output as html
    seqinfos=[]
    aligs_merged_seq=[]
    aligs_merged_ref=[]
    frame_info=[]
    for ampl_row in range(round(len(df.index)*0.1)):
        match1 = SequenceMatcher(None, df.columns[1].split(":")[-1], df.iloc[ampl_row,1]).find_longest_match(0, len(df.columns[1].split(":")[-1]), 0, len(df.iloc[ampl_row,1]))
        end_of_match=match1.b+match1.size
        matched_ampl1= len(df.columns[1].split(":")[-1][:match1.a])*"-" + str(df.iloc[ampl_row,1][match1.b:]) + "|"
        seq_inframe1=matched_ampl1[0:end_of_match] + "|"
        ref_inframe1=df.columns[1].split(":")[-1][0:end_of_match] + "|"
        #get the matches with the other reading frames (i.e. the last 2 cols of the df)
        match2_a = SequenceMatcher(None, df.columns[2].split(":")[-1], df.iloc[ampl_row,2]).find_longest_match(0, len(df.columns[2].split(":")[-1]), 0, len(df.iloc[ampl_row,2]))
        match2_b = SequenceMatcher(None, df.columns[3].split(":")[-1], df.iloc[ampl_row,3]).find_longest_match(0, len(df.columns[3].split(":")[-1]), 0, len(df.iloc[ampl_row,3]))
        
        #modify this part so that you include AAs generated by all different frames
        if match2_a.size>match2_b.size:
            matched_ampl2= len(df.columns[2].split(":")[-1][:match2_a.a])*"-" + str(df.iloc[ampl_row,2][match2_a.b:]) + "|"
            ref_outframe2=df.columns[2].split(":")[-1][0:end_of_match] + "|"
            # seq_inframe2=matched_ampl2[0:end_of_match]
            seq_outframe=matched_ampl2[end_of_match:]
            frame_info.append("Frame:+1 " + ref_inframe1.index("|")*" " + "Frame:0")
        else:
            matched_ampl2= len(df.columns[3].split(":")[-1][:match2_b.a])*"-" + str(df.iloc[ampl_row,3][match2_b.b:]) + "|"
            ref_outframe2=df.columns[3].split(":")[-1][end_of_match:] + "|"
            seq_outframe=matched_ampl2[end_of_match:]
            # seq_inframe2=matched_ampl2[0:end_of_match]
            seq_outframe=matched_ampl2[end_of_match:]
            frame_info.append("Frame:+1 " + ref_inframe1.index("|")*" " + "Frame:+2")
            print("Frame:+1 " + ref_inframe1.index("|")*" " + "Frame:+2")
        merged_align_seq=seq_inframe1+seq_outframe

        merged_align_ref=ref_inframe1+ref_outframe2
        header_info=df.iloc[ampl_row,0]
        #add one more line which adds the frames

        seqinfos.append(df.iloc[ampl_row,0])
        #aligs_merged_seq.append(merged_align_seq[0:round(len(aligs_merged_seq)*0.6)])
        aligs_merged_seq.append(merged_align_seq)
        #aligs_merged_ref.append(merged_align_ref[0:round(len(merged_align_ref)*0.6)])
        aligs_merged_ref.append(merged_align_ref)
    aa_coloring(seqinfos,aligs_merged_seq,aligs_merged_ref, frame_info, output_html)


#old approach for removing duplicates
def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return zip(*[iter(iterable)]*n)

def remove_multiples(index_pos, id_keys_dict):
    perc_sum_all=0
    var_sum_all=0
    sd_sum_all=0
    if len(index_pos)%2==0:
        for x,y in grouped(index_pos,2):
            perc_sum=float(id_keys_dict[x].split("_")[0].split(":")[1])+float(id_keys_dict[y].split("_")[0].split(":")[1])
            perc_sum_all+=perc_sum
            var_sum=float(id_keys_dict[x].split("_")[1].split(":")[1])+float(id_keys_dict[y].split("_")[1].split(":")[1])
            var_sum_all+=var_sum
            var_sum=float(id_keys_dict[x].split("_")[2].split(":")[1])+float(id_keys_dict[y].split("_")[2].split(":")[1])
            var_sum_all+=var_sum

    else:
        index_pos=index_pos.append(0)
        for x,y in grouped(index_pos,2):
            perc_sum=float(id_keys_dict[x].split("_")[0].split(":")[1])+float(id_keys_dict[y].split("_")[0].split(":")[1])
            perc_sum_all+=perc_sum
            var_sum=float(id_keys_dict[x].split("_")[1].split(":")[1])+float(id_keys_dict[y].split("_")[1].split(":")[1])
            var_sum_all+=var_sum
            var_sum=float(id_keys_dict[x].split("_")[2].split(":")[1])+float(id_keys_dict[y].split("_")[2].split(":")[1])
            var_sum_all+=var_sum
    
    #fix the sd summing!
    merged_id=">CluSeq:"+str(perc_sum_all)+ "_var:"+str(var_sum_all) +"_sd:"+str(sqrt(var_sum_all))
    return(merged_id)
