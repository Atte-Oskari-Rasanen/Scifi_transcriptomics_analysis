from fnmatch import translate
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

def trimRead_hiti(animal_nr,base_path,transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd):
    animal_nr = str(animal_nr)
    "Filters and trims the reads"
    search_path = base_path+animal_nr+'*'+transgene+'*'+assay_end+'*/'
    
    animal_p5_cat = tempfile.NamedTemporaryFile(suffix = '.fastq.gz').name
    animal_p7_cat = tempfile.NamedTemporaryFile(suffix = '.fastq.gz').name
    test_file_p5_out = tempfile.NamedTemporaryFile(suffix = '.fastq').name
    test_file_p7_out = tempfile.NamedTemporaryFile(suffix = '.fastq').name
    test_file_p5_filter = tempfile.NamedTemporaryFile(suffix = '.fastq').name
    
    if read_fwd:
        animal_p5 = glob.glob(search_path+'*R1*')
        animal_p7 = glob.glob(search_path+'*R2*')
        #display('Forward run Animal: '+animal_nr)
    else:
        animal_p5 = glob.glob(search_path+'*R2*')
        animal_p7 = glob.glob(search_path+'*R1*')
        #display('Reverse run Animal: '+animal_nr)
    

    cat_p5= "cat "+" ".join(animal_p5)+" > "+animal_p5_cat
    call([cat_p5], shell=True)
    cat_p7= "cat "+" ".join(animal_p7)+" > "+animal_p7_cat
    call([cat_p7], shell=True)

    stats_out = export_path+animal_nr+'_'+transgene+'_'+assay_end+'_stats-filter.txt'
    
    kmer = '20'
    hdist = '3'
    param=" k="+kmer+" hdist="+hdist+" rcomp=f skipr2=t threads=32 overwrite=true"
    
    call_sequence = "/media/data/AtteR/Attes_bin/bbmap/bbduk.sh in="+animal_p7_cat+" in2="+animal_p5_cat+" outm1="+test_file_p7_out+" outm2="+test_file_p5_out+" literal="+filterlitteral+" stats="+stats_out + param
    call([call_sequence], shell=True)
    
    call_sequence = "/media/data/AtteR/Attes_bin/bbmap/bbduk.sh in="+test_file_p5_out+" out="+test_file_p5_filter+ " literal=AAAAAAAAA,CCCCCCCCC,GGGGGGGGG,TTTTTTTTT k=9 mm=f overwrite=true minlength=40"
    call([call_sequence], shell=True)
    
    test_file_p5_out_starcode = tempfile.NamedTemporaryFile(suffix = '.tsv').name
    starcode_call= "/media/data/AtteR/Attes_bin/starcode/starcode -i "+test_file_p5_filter+" -t 32 -o "+test_file_p5_out_starcode
    call([starcode_call], shell=True)
    
    df=pd.read_csv(test_file_p5_out_starcode, sep='\t', header=None)
    df = df.rename(columns={0: 'sequence', 1:'count'})
    total_counts = int(df[['count']].sum())
    df = df[df['count'].astype(int)>total_counts/10000]
    total_counts = int(df[['count']].sum())
    df['percent'] = (df['count'] / total_counts)
    df = df.rename(columns={'percent':animal_nr+'_percent','count':animal_nr+'_count',})
    
    return df
def align_to_ref(query_sequence, target_sequence):
    alignment, score, start_end_positions = local_pairwise_align_ssw(
        DNA(target_sequence),
        DNA(query_sequence),
        gap_open_penalty = 3,
        gap_extend_penalty = 1
    )
    out_align = ('-'*start_end_positions[0][0])+str(alignment[1])+('-'*(len(target_sequence)-start_end_positions[0][1]-1))
    
    return out_align

"""
understand how the pairwise alignment was done and put into df
use different ones and compare results
visualise
maybe apply same with trimming methods
translate to AAs
"""
transgene = 'mCherry'
assay_end = '3p'
read_fwd = True
animal_list = [7, 8, 9, 10, 11, 12] 
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT'
lliteral = ' literal=GGCGGCATGGACGAG'
rliteral = ' literal=CATATGACCACCGG'
#base_path = '/home/lcadmin/mnm-lts/SequencingFiles/arc_hiti_asyn/HITI1-8977973/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
#export_path = '/home/lcadmin/mnm-lts/HITI-analysis/'
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path = '/media/data/AtteR/projects/hiti/output/'
target_sequence = "CGGCGGCATGGACGAGCTGTACAAGGTCGGTGCTGCGGCTCCGCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAGGTCTTCTACCGTCTGGAGAGG"

complete_df = pd.DataFrame({'sequence': ['CTGTACAAGGTCGGTGCTGCGGCTCCGCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGAC']})
complete_df

#this line caused the perfect match to occur since we had the target seq included!
#complete_df = pd.DataFrame({'sequence': [target_sequence]})

#trim the reads and make a df of each animal containing the cluster groups, number of seqs they contain and % of the total seqs in
#the cluster and merge everything together based on the seq (if found same clusters, otherwise NA)
for animal in animal_list:
    df_this = trimRead_hiti(animal,base_path,transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd)
    complete_df = pd.merge(complete_df, df_this, on="sequence", how='outer')
complete_df.columns

complete_df = complete_df.fillna(value=0)
perc_cols = [col for col in complete_df.columns if 'percent' in col]
perc_cols #['7_percent', '8_percent', '9_percent', '10_percent', '11_percent', '12_percent']

#sum the percentages of each seq cluster (animals 7-12)
complete_df['percent_sum'] = complete_df[perc_cols].sum(axis=1)

complete_df.sort_values(by=['percent_sum'], ascending=False, inplace=True)

#generate a column with seq alignment wher eyou take the seq cluster from sequence column and map it to the target
complete_df.loc[:,'sequence_align'] = complete_df.loc[:,'sequence'].apply(lambda x: align_to_ref(x, target_sequence))
export_csv = export_path+transgene+'_'+assay_end+'.csv'
complete_df.to_csv(export_csv, index=False)

complete_df.loc[:,'sequence_align']

from Bio.Seq import Seq 

complete_df.iloc[:,0]
s = Seq(complete_df.iloc[1,0])
s

def align_to_ref(query_sequence, target_sequence):
    alignment, score, start_end_positions = local_pairwise_align_ssw(
        DNA(target_sequence),
        DNA(query_sequence),
        gap_open_penalty = 3,
        gap_extend_penalty = 1
    )
    out_align = ('-'*start_end_positions[0][0])+str(alignment[1])+('-'*(len(target_sequence)-start_end_positions[0][1]-1))
    
    return out_align


al = align_to_ref(complete_df.iloc[2,0],target_sequence)
al
alignment, score, start_end_positions = local_pairwise_align_ssw(
    DNA(target_sequence),
    DNA(complete_df.iloc[2,0]),
    gap_open_penalty = 3,
    gap_extend_penalty = 1
)
alignment[0]
out_align = ('-'*start_end_positions[0][0])+str(alignment[1])+('-'*(len(target_sequence)-start_end_positions[0][1]-1))

seqA = Seq("ACCGGT") 
seqB = Seq("ACGT")
alignment, score, start, end = pairwise2.align.globalxx(seqA, seqB, one_alignment_only=True)
al_s = pairwise2.align.globalxx(seqA, seqB, one_alignment_only=True)
al_s[0]
al_s_f =format_alignment(*al_s[0])
al_s_f
#takes in the clustered seqs, makes them into seqrecord objects

#transforms the aligned seqs into seqrecord objects
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#pairwise2.align.globalxx(seq1, seq2)

#extract the aligned seq part, ignoring the scores etc. 
#to process and save the files downstream, must convert the result into seq record
#first take specific alignment from the alignment object via indexing, then format it into a string, split based on line changes
#and take the 3rd element, i.e the NTs matching the template one, i.e. the alignment. Transform this into seq object, save
#objects into a list and then iteratively save into fasta file. visualise them
#try different alignment methods


#we need to know what percentage of the clusters mapped with high accuracy, ideally 100%, to the reference genome. Thus
#we must save the percentages from the starcluster data as these are used as identifiers!!!!
############################
############################

#Need the seq cluster and the % attached to it. maybe instead of aligning all t
complete_df
#Then find the one with the highest match

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

seq_and_perc = {}

for cols in complete_df.columns:
    #print(cols)
    if "percent" in cols and not "sum" in cols:
        print(cols)
        seq_and_perc[cols]=complete_df.loc[:,["sequence", cols]]
seq_and_perc["8_percent"].iloc[:,0]


############################
############################
#if we align using muscle we will need to transfer our seq files into fasta first
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

input = "/media/data/AtteR/projects/hiti/mcherry_p3_seq_clusters.fasta"
output = "/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_clusters_aligned.fasta"
############################
############################


import subprocess
subprocess.call('muscle -align %s -output %s'%(input,output))

muscle_exe = "/home/lcadmin/miniconda3/envs/scifi-analysis/bin/muscle" #specify the location of your muscle exe file



complete_df
seq_and_perc["12_percent"]
alignments = pairwise2.align.globalms(target_sequence, seq_and_perc[group].iloc[i,0]
#returns a dict with percentage value of the certain cluster seq and the aligned seq as the value
def align_global(complete_df,target_sequence):
    seq_and_perc = {}

    #make another dict into which you save the percentage value of the seq into key and the aligned seq as the value
    align_and_perc = {}

    for cols in complete_df.columns:
        #print(cols)
        if "percent" in cols and not "sum" in cols:
            print(cols)
            seq_and_perc[cols]=complete_df.loc[:,["sequence", cols]]

    #you take the group, from the df you take the sequence cluster column and align each sequence 
    for group in seq_and_perc.keys():
        #alignments = pairwise2.align.globalms(target_sequence, seq_and_perc[group][:,0],  2, -1, -.5, -.1)

        for i in range(len(seq_and_perc[group].index)):
            #query=seq_and_perc[group][i,0]
            alignments = pairwise2.align.globalms(target_sequence, seq_and_perc[group].iloc[i,0],  2, -1, -.5, -.1)
            #alignments = pairwise2.align.globalms(target_sequence, seq_and_perc[group][:,0],  2, -1, -.5, -.1)

            alignm=format_alignment(*alignments[0])
            #make string objects, save into a list. then count which one has least ----, cut off the rest of the reference based on this?
            seq_align = alignm.split("\n")[2]
            #nt_count=count_nts(seq_align)
            #seq_and_perc[group]["match"]
            align_and_perc[seq_and_perc[group].iloc[i,1]]=alignm.split("\n")[2]
            #aligned_seqs.append(seq_align)
    return(align_and_perc)

align_and_perc=align_global(complete_df,target_sequence)
align_and_perc

#make a file for especially muscle alignment which you can then visualise. This file contains the 
#percentages of seq clusters

result="/media/data/AtteR/projects/hiti/mcherry_p3_seq_cluster_all.fasta"
id_f = 1
with open(result, "w") as handle:
    seq_obj = SeqRecord(Seq(target_sequence), id=str(0), description="mcherry_p3_seq_ref")
    count = SeqIO.write(seq_obj, handle, "fasta")

    for group in seq_and_perc.keys():
        for seq_i in range(len(seq_and_perc[group])):
            descr="ClusterSeq%: " + str(round((seq_and_perc[group].iloc[seq_i,1]*100),4))
            print(descr)
            seq_obj = SeqRecord(Seq(seq_and_perc[group].iloc[seq_i,0]), id=str(id_f), description=descr)
            count = SeqIO.write(seq_obj, handle, "fasta")
            id_f+=1

#this saves them via seq record formatting. However, when importing this content back for translation,
#there are new line characters in some of the longer seqs so need to create another approach too
id_f=1
result="/media/data/AtteR/projects/hiti/mcherry_p3_seq_cluster_all.fasta"
with open(result, "w") as handle:
    seq_obj = SeqRecord(Seq(target_sequence), id=str(0), description="mcherry_p3_seq_ref")
    header=">0"+" mcherry_p3_seq_ref"
    handle.write(header + "\n" + target_sequence + "\n")

    for group in seq_and_perc.keys():
        for seq_i in range(len(seq_and_perc[group])):
            header=">"+ str(id_f)+" ClusterSeq_%: " + str(round((seq_and_perc[group].iloc[seq_i,1]*100),4))
            handle.write(header + "\n" + seq_and_perc[group].iloc[seq_i,0] + "\n")
            id_f+=1

aligned_seqs = []
#seq_and_perc["8_percent"].iloc[,1]
#extract the df, then align, then save into a file

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

#Translate and save as dictionary containing header as the key and the translated AA seq as value

amplicon=seq_and_perc["8_percent"].iloc[1,0]

#this may be wrong (bottom, 5'-3')
mcherry_full="TACCACTCGTTCCCGCTCCTCCTATTGTACCGGTAGTAGTTCCTCAAGTACGCGAAGTTCCACGTGTACCTCCCGAGGCACTTGCCGGTGCTCAAGCTCTAGCTCCCGCTCCCGCTCCCGGCGGGGATGCTCCCGTGGGTCTGGCGGTTCGACTTCCACTGGTTCCCACCGGGGGACGGGAAGCGGACCCTGTAGGACAGGGGAGTCAAGTACATGCCGAGGTTCCGGATGCACTTCGTGGGGCGGCTGTAGGGGCTGATGAACTTCGACAGGAAGGGGCTCCCGAAGTTCACCCTCGCGCACTACTTGAAGCTCCTGCCGCCGCACCACTGGCACTGGGTCCTGAGGAGGGACGTCCTGCCGCTCAAGTAGATGTTCCACTTCGACGCGCCGTGGTTGAAGGGGAGGCTGCCGGGGCATTACGTCTTCTTCTGGTACCCGACCCTCCGGAGGAGGCTCGCCTACATGGGGCTCCTGCCGCGGGACTTCCCGCTCTAGTTCGTCTCCGACTTCGACTTCCTGCCGCCGGTGATGCTGCGACTCCAGTTCTGGTGGATGTTCCGGTTCTTCGGGCACGTCGACGGGCCGCGGATGTTGCAGTTGTAGTTCAACCTGTAGTGGAGGGTGTTGCTCCTGATGTGGTAGCACCTTGTCATGCTTGCGCGGCTCCCGGCGGTGAGGTGGCCGCCGTACCTGCTCGACATGTTC"

mcherry_full="ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAG"

common=set(mcherry_full).intersection(amplicon)

#make a function that takes in each amplicon along with the full mcherry template, finds the intersection, finds the 
#distance from the starting codon (lies in the start of mcherry) till the start of this intersection, calculate the 
# remainder, if not in frame, start translating the AAs into protein from the Nth position of the amplicon 

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


def find_frame(overlap_bases):
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

def translate_nt_aa(amplicon, mcherry_full):
    overlap=find_overlap(amplicon, mcherry_full)
    frameN=find_frame(overlap)
    seq=Seq(amplicon[frameN:]).translate()
    return(seq)

aa_and_perc={}

#no intersection found for most of the sequences....how to deal with these cases? this would imply that the cluster seq
#exists but downstream
#issues with aligning with muscle downstream
with open(result) as nt_seq:
    header=[]
    for i, line in enumerate(nt_seq):
        if line.startswith(">"):
            header.append(line.strip())
        else:
            print(header)
            amplicon=line.strip()
            aa_and_perc[header[0]]=str(translate_nt_aa(amplicon, mcherry_full))
            # overlap=find_overlap(line, mcherry_full)
            # print(overlap)
            # frameN=find_frame(overlap)
            # aa_and_perc[header[0]]=str(Seq(amplicon[frameN:]).translate())
            header=[]
            #print(str(Seq(next(nt_seq)).translate()))
    print("done")
aa_and_perc

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


#trim the overhangs from muscle file
muscle_alignment="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_clusters_all_aligned_nonrecseq.fasta"
muscle_alignment

from Bio import SeqIO

muscle_seqs = SeqIO.parse(open(muscle_alignment),'fasta')
muscle_seqs

mcherry_full

#find the highest match


#parse the muscle file, put into dict. once we have specific seq, calculate the similarity score, add to the id key for later

#another option is to make a function that takes in all the seqs and finds the highest match, then trims all the seqs based on this.
#function takes in the dict of muscle alignment file
muscle_dict = {}
for seq_record in SeqIO.parse(muscle_alignment, "fasta"):
     print(seq_record.id)
     print(repr(seq_record.seq))

     print(len(seq_record))


with open(output_file) as out_file:
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        new_sequence = some_function(sequence)
        write_fasta(out_file)
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

#you need to find the highest count but you also need to keep the index so save this info into the same dic. Then we define the highest matching seq's len as 
#the cut off for the rest and thus we cut off the extra from the others as well.

#get all the results and save them into the list, find the highest one along with its index, then trim the lengths of all seqs based on this, then calculate the percentage of
#match relative to the reference template, save into list, put into the df

#if several nts with equally good match, get the one which extends longest and trim seqs based on this one
def 
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


