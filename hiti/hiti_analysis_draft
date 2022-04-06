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
template_seq = SeqRecord(Seq(target_sequence), id="TEMPLATE", description="mcherry_p3_reference")
alignments = pairwise2.align.globalms(target_sequence, complete_df.iloc[3,0],  2, -1, -.5, -.1) #localxx earlier
alignments
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
############################
#Need the seq cluster and the % attached to it. maybe instead of aligning all t
complete_df
#Then find the one with the highest match
template_seq = SeqRecord(Seq(target_sequence), id="TEMPLATE", description="mcherry_p3_reference")

#take col names which contain percent in them

#complete_df.loc[:,["sequence","8_percent"]]
seq_and_perc = {}
for cols in complete_df.columns:
    #print(cols)
    if "percent" in cols and not "sum" in cols:
        print(cols)
        seq_and_perc[cols]=complete_df.loc[:,["sequence", cols]]
seq_and_perc["8_percent"].iloc[:,0]

#als = pairwise2.align.globalms(target_sequence, complete_df.iloc[3,0],  2, -1, -.5, -.1) #localxx earlier
#l = format_alignment(*als[0])
#l_s = l.split("\n")[2]


aligned_seqs = []
#seq_and_perc["8_percent"].iloc[,1]
#extract the df, then align, then save into a file

#make another dict into which you save the percentage value of the seq into key and the aligned seq as the value
align_and_perc = {}

def count_nts(seq):
    nts=0
    for base in seq:
        if base.isalpha():
            nts+=1
    return nts
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
df2

#now we have found the longest seq and trimmed based on that and calculated match % to the reference
#next put this into function and compare results. If more than one match of max length, then report how many seqs have this
#we should then take each one of these seqs, make them as seqrecords and save as a fasta file which then can be visualised

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


result="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_aligned_log.txt"
with open(result, "w") as handle:
    count = Bio.AlignIO.write(alignments, handle, format)

with open("example.faa", "w") as handle:
    count = SeqIO.write(alignments, handle, "fasta")

a=str(alignments[0])
b=a.split(",")
b
alignment, score, start, end = alignments[0].split(",")
alignment, score, start, end = str(alignments[0]).split(",")
alignments[0]

p = view_alignment(alignments, plot_width=900)
pn.pane.Bokeh(p)

#figure out a way to put the alignments into a format in which you CAN visualise them!
