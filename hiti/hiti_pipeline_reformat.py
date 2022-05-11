from cmath import nan
from heapq import merge
import pandas as pd
import os 
os.getcwd()
os.chdir("/media/data/AtteR/projects/hiti")
from scripts_hiti import *

'''
Exp1
LP1
#1_GFP       #2_mCherry
1a 5'dna    2a
1b 3'dna    2b
1c 5'aa
1d 3'aa     2d

Exp2

#3_mCherry   #4_mcherry  #5_SATI
3a 5'dna     4a          5a 
3b 3'dna     4b
             4c
3d 3'aa      4d
'''

def calculate_perc_sd2(full_df):
    full_df = full_df.fillna(value=0)
    perc_cols = [col for col in full_df.columns if 'percent' in col]
    count_cols = [col for col in full_df.columns if 'count' in col]

    perc_cols
    #sum the percentages of each seq cluster of each animal together to see how much the certain seq is found 
    full_df['percent_sum_unit'] = full_df[perc_cols].sum(axis=1)  

    total_perc_unit=full_df.iloc[:,-1].sum()
    full_df['percent_sum'] = (full_df['percent_sum_unit'] / total_perc_unit) #so divide the total perc of each seq (summed across animals) with total summed percentage of summed perc column 
    full_df.head()
    count_cols
    #full_df['sd']=full_df[count_cols].std()
    full_df['total_reads_seq'] = full_df[count_cols].sum(axis=1)  
    full_df['variance']=full_df[perc_cols].var(axis=1)
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

    return(full_df_trim.iloc[0:round(len(full_df_trim.index)*0.3),:])

def create_datadict2(base_path, transgene):
    group_folders = [folder for folder in os.listdir(base_path) if transgene in folder]
    animals=[]
    for s in group_folders:
        animals.append("_".join(s.split("_")[:3]))

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


#==EXP1==#
#GFP 3p
#############
transgene = 'GFP'
assay_end = '3p'
read_fwd = True
filterlitteral = 'GGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGG'
lliteral = ' literal=GTGGTCATATGGTCCAGCTCC'
rliteral = ' literal=TCGTCCATGCCGAG'
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path = '/media/data/AtteR/projects/hiti/pipeline_output_reorg/'
target_sequence = "GTGGTCATATGGTCCAGCTCCATCTGGTCGTCGGTGCTGCGGCTCCGCGGAGCCGCAGCACCGACCTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGCTTCTCGTTGGGGTCTTTGCTCAGGGCGGACTGGGTGCTCAGGTAGTGGTTGTCGGGCAGCAGCACGGGGCCGTCGCCGATGGGGGTGTTCTGCTGGTAGTGGTCGGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGATCTTGAAGTTCACCTTGATGCCGTTCTTCTGCTTGTCGGCCATGATATAGACGTTGTGGCTGTTGTAGTTGTACTCCAGCTTGTGCCCCAGGATGTTGCCGTCCTCCTTGAAGTCGATGCCCTTCAGCTCGATGCGGTTCACCAGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGGTCTTGTAGTTGCCGTCGTCC"
data_dict=create_datadict2(base_path, transgene)
full_df=import_reads_process(data_dict, transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd)
full_df_trim_orig=calculate_perc_sd2(full_df)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_3p_GFP_LP1.fasta"
save_fasta(result, full_df_trim_orig, target_sequence)


#########
#NT
####################
output_path="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html"
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/exp1/Exp1_3p_GFP_LP1_local3.fasta"
test_res=aligner(full_df_trim_orig, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/exp1/Exp1_3p_GFP_LP1_local2.fasta"
test_res=aligner(full_df_trim_orig, target_sequence, "align_local2", result, output_path, 3,1)
####################


#AA
####################
corr_frame=1
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_3p_GFP_LP1.fasta"
df_aa=translate_nt_aa(result, 1)
df_aa.columns
df_aa.iloc[1,2]
output_html="/media/data/AtteR/projects/hiti/pipeline_output_reorg/AA_aligned/Exp1_3p_GFP_LP1_AA.html"
visualise_aa_hybrid_alignments(df_aa, output_html)

df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral,read_fwd)
####################

#mCherry 5p
#############
transgene = 'GFP'
assay_end = '5p'
read_fwd = True
filterlitteral = 'GCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTCCTGGACG'
lliteral = ' literal=CCTCAGAGGAGTTCT'
rliteral = ' literal=GGCGAGGAGCTGTT'
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path = '/media/data/AtteR/projects/hiti/pipeline_output_reorg/'
target_sequence = "TTAGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTAACAGGCGCGCCACCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGC"
target_sequence=target_sequence.upper()
#########
data_dict=create_datadict2(base_path, transgene)
full_df=import_reads_process(data_dict, transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd)
full_df_trim_orig=calculate_perc_sd2(full_df)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_5p_GFP_LP1.fasta"
save_fasta(result, full_df_trim_orig, target_sequence)

#NT
####################
output_path="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html"
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_5p_GFP_LP1_local3.fasta"
test_res=aligner(full_df_trim_orig, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_5p_GFP_LP1_local2.fasta"
test_res=aligner(full_df_trim_orig, target_sequence, "align_local2", result, output_path, 3,1)
####################

#AA
####################
corr_frame=1
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_5p_GFP_LP1.fasta"
df_aa=translate_nt_aa(result, 1)
df_aa.columns
df_aa.iloc[1,2]
output_html="/media/data/AtteR/projects/hiti/pipeline_output/AA_aligned/Exp1_5p_GFP_LP1_AA.html"
visualise_aa_hybrid_alignments(df_aa, output_html)
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral,read_fwd)
#############



#mCherry 3p
#############
transgene = 'mCherry'
assay_end = '3p'
read_fwd = True
filterlitteral='CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' 
lliteral=' literal=GGCGGCATGGACGAGC' #added C at the end as its part of tje primer #to check that its on target with mcherry - the primer
rliteral=' literal=CATATGACCACCGG'
#base_path = '/home/lcadmin/mnm-lts/SequencingFiles/arc_hiti_asyn/HITI1-8977973/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
#export_path = '/home/lcadmin/mnm-lts/HITI-analysis/'
base_path='/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path='/media/data/AtteR/projects/hiti/pipeline_output_reorg/'
#with primer
target_sequence="TGTACAAGGtcggtgctgcggctccgCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAGGTCTTCTACCGTCTGGAGAGG"
target_sequence=target_sequence.upper()
data_dict=create_datadict2(base_path, transgene)
full_df=import_reads_process(data_dict, transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd)
full_df_trim_orig=calculate_perc_sd2(full_df)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_3p_mcherry_LP1.fasta"
save_fasta(result, full_df_trim_orig, target_sequence)

#########

#NT
####################
output_path="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html"
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_3p_mcherry_LP1_local3.fasta"
test_res=aligner(full_df_trim_orig, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_3p_mcherry_LP1_local2.fasta"
test_res=aligner(full_df_trim_orig, target_sequence, "align_local2", result, output_path, 3,1)
####################

#AA
####################
corr_frame=1
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_3p_mcherry_LP1.fasta"
df_aa=translate_nt_aa(result, 1)
df_aa.columns
df_aa.iloc[1,2]
output_html="/media/data/AtteR/projects/hiti/pipeline_output/AA_aligned/Exp1_3p_mcherry_LP1_AA.html"
visualise_aa_hybrid_alignments(df_aa, output_html)

df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral,read_fwd)

############
transgene='mCherry'
assay_end='5p'
read_fwd=True
filterlitteral='CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' #region prior to r2 primer
lliteral = ' literal=GGGCGAGGAGGATAACATGG' 
#base_path = '/home/lcadmin/mnm-lts/SequencingFiles/arc_hiti_asyn/HITI1-8977973/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
#export_path = '/home/lcadmin/mnm-lts/HITI-analysis/'
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path = '/media/data/AtteR/projects/hiti/pipeline_output_reorg/'

#target_sequence = "CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTaacaggCGCGCCACCATGGTGAGCAAGGGCGAGGAGGATAACATGG"
target_sequence = "CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTaacaggCGCGCCACCATGGTGAGCAA"
target_sequence = target_sequence.upper()
data_dict=create_datadict2(base_path, transgene)
data_dict
full_df=import_reads_process(data_dict, transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd)
full_df_trim_orig=calculate_perc_sd2(full_df)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_5p_mcherry_LP1.fasta"
save_fasta(result, full_df_trim_orig, target_sequence)

#########
#NT
####################
output_path="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html"
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_5p_mcherry_LP1_local3.fasta"
aligner(full_df_trim_orig, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_5p_mcherry_LP1_local2.fasta"
aligner(full_df_trim_orig, target_sequence, "align_local2", result, output_path, 3,1)
####################

#AA
corr_frame=1
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_5p_mcherry_LP1.fasta"
df_aa=translate_nt_aa(result, 1)
df_aa.columns
df_aa.iloc[1,2]
output_html="/media/data/AtteR/projects/hiti/pipeline_output/AA_aligned/Exp1_5p_mcherry_LP1_AA.html"
visualise_aa_hybrid_alignments(df_aa, output_html)
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral,read_fwd)
#############



#==EXP2==#
#############
#SP4 3' 
trimmed_df="/media/data/AtteR/projects/hiti/dfs/full_df_trim_mcherry_p3_HITI_SP4.csv"
transgene = 'mCherry'
assay_end = '3p'
read_fwd = True
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_3p"
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' #region prior to r2 primer
lliteral=" literal=CGGCGGCATGGACGAGC"
target_sequence="tgtacaagATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAGGTCTTCTACCGTCTGGAGAGG"
target_sequence=target_sequence.upper()
#########
#########

#The counts and percs for each animal regarding cluster seq are identical

#########
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral,read_fwd)
#starcode produces files with exactly same seqs, counts and percs
#files do have different number of reads but for some reason starcode clusters all of them as the same
df_trim_full2=calculate_perc_sd2(df_full)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/unaligned_seqs/Exp2_3p_mcherry_mcherry_SP4_var.fasta"
save_fasta(result, df_trim_full2, target_sequence)


#Checking for duplicates and removing them - a dummy script
#################

test_df=df_full.iloc[0:10,0:]
test_df
test_df2=pd.DataFrame(columns = test_df.columns)
vals1=['CCTCGGCATGGACGAGCTGTACAAGATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCCTCCACGCCTACCC', 1000,1000,1000,1000,1000,1000]
test_df.loc[-1] = vals1
test_df.index = test_df.index + 1
test_df = test_df.sort_index()
test_df.fillna(0)
vals2=['TGTACAAGATGGAGCTGGACCATATGACCCGGTGCCGGCGGCCTCCACGCCTACCCTG', 1000,1000,1000,1000,1000,1000]
test_df.loc[-1] = vals2
test_df.index = test_df.index + 1
test_df = test_df.sort_index()
test_df.fillna(0)
vals3=['CTTCGGCATGGACGAGCTGTACAAGATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCCTCCACGCCTACCC', 1000,1000,1000,1000,1000,1000]
test_df.loc[-1] = vals3
test_df.index = test_df.index + 1
test_df = test_df.sort_index()
test_df.fillna(0)
vals4=['TGTACAAGATGGAGCTGGACCATATGACCCGGTGCTTTGGTCGCCCGGCCTCAGTGAG',1000,1000,1000,1000,1000,1000]
test_df.loc[-1] = vals4
test_df.index = test_df.index + 1
test_df = test_df.sort_index()
test_df.fillna(0)

test_df.columns[1:]
grouped=test_df.groupby("sequence")
grouped.first()
grouped = test_df.groupby(['sequence'], as_index=False)[test_df.columns[1:]].sum()
grouped
#################
def aligner2(full_df, target_sequence, align_method, filename, output_path, gop=3, gep=1):
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
        header=">"+ str(id_f)+"CluSeq:" + str((round(full_df.iloc[seq_i,-4],5))) + "_var:"+str((round(full_df.iloc[seq_i,-2],5))) +"_sd:" + str((round(full_df.iloc[seq_i,-1],5)))
        #seq_obj_align = aligner_init(full_df.iloc[seq_i,0], target_sequence, gop, gep).align()

        print("===SEQ===:" + full_df.iloc[seq_i,0])
        seq_obj_align = aligner_init(full_df.iloc[seq_i,0], target_sequence, gop, gep).align()
        seq_obj_align = re.sub(r'[(\d|\s]', '', seq_obj_align) #remove digits from the string caused by the alignment and empty spaces from the start
        aligned_data[header]=seq_obj_align
        id_f+=1

    aligned_data_trim=align_trimmer(aligned_data, target_sequence)
    write_align(aligned_data_trim, filename, target_sequence)
    #Generate a visual alignment file using mview
    mview_file=output_path + "/" + filename.split("/")[-1].split(".")[-2] + ".html"
    mview_command='/media/data/AtteR/Attes_bin/mview -in fasta -html head -css on -coloring any ' + filename + '>' + mview_file
    call([mview_command], shell=True)
    print("html file created!")
    return(aligned_data_trim)

alignments=aligner2(df_trim_full2, target_sequence, "align_local2", result, output_path, 3,1)
algs=list(alignments.values())

alignments.keys()
len(algs)
#get the prev dict values (consists of the id as the keys and the aligned nt seq as value)
#and make a new dict of it in which the key is the index and the numerical values, perc and sd, are 
#saved as list elements as the value part
id_keys_dict=dict()
for i, id in enumerate(alignments.keys()):
    id_keys_dict[i]=[id.split("_")[0].split(":")[1], id.split("_")[1].split(":")[1]]


id_keys_dict=dict()
for i, id in enumerate(alignments.keys()):
    id_keys_dict[i]=id
# none_match = {0:'>CluSeq:0_var:0_sd:0'} 
# id_keys_dict.update(none_match)
# id_keys_dict[1].split("_")[2].split(":")[1]


perc_sum=float(id_keys_dict[33].split("_")[0].split(":")[1])+float(id_keys_dict[i+1].split("_")[0].split(":")[1])
def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return zip(*[iter(iterable)]*n)
#MITEN ITER TOIMII? NONEE EI HYVÄKSY ELI KEKSI MUU KONSTI...


l=[0, 1, 2, 4, 33]
for x, y in grouped(l, 2):
   print("%d,%d" % (x, y))#the generator function takes in the index_pos list and the dict, then sums them 
import math
#go over the indeces where the duplicates of a seq found, two at a time and sum them all up together
def remove_multiples(index_pos, id_keys_dict):
    perc_sum_all=0
    var_sum_all=0
    sd_sum_all=0
    for i in index_pos:
        perc_sum_all+=float(id_keys_dict[i].split("_")[0].split(":")[1])
        var_sum_all+=float(id_keys_dict[i].split("_")[1].split(":")[1])
        sd_sum_all+=float(id_keys_dict[i].split("_")[2].split(":")[1])
    #fix the sd summing!
    merged_id=">CluSeq:"+str(perc_sum_all)+ "_var:"+str(var_sum_all) +"_sd:"+str(math.sqrt(var_sum_all))
    return(merged_id)
#id_keys_dict[33].split("_")[0].split(":")[1]

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

def remove_key_from_dict(index_pos,id_keys_dict,alignments):
    for pos in index_pos:
        #id_keys_dict[pos] value corresponds to the key of the alignments 
        print(id_keys_dict[pos])
        print(alignments[id_keys_dict[pos]])
        del alignments[id_keys_dict[pos]]
    return(alignments)
t=get_indeces_matching_seqs(algs, 'TGTACAAGATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCCTCCACGCCTACCC-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
t
aa=alignments.copy()
aa
del alignments['>4CluSeq:0.01283_var:0.00049_sd:0.02223']
len(alignments.keys())
import collections
#merged values
dups = [item for item, count in collections.Counter(algs).items() if count >= 2]
merged_multiples=dict()
for dupl_seq in dups:
    index_pos = get_indeces_matching_seqs(algs, dupl_seq)
    #take the seq from any of the matches
    #seq=algs[index_pos[0]]
    print(dupl_seq[:40] + " found at pos " + str(index_pos))
    alignments=remove_key_from_dict(index_pos,id_keys_dict,alignments)
    #integrate the merged values into the dict in the end
    #single_id=remove_multiples(index_pos, id_keys_dict)
    #key includes the merged matches values and the value the seq
    merged_multiples[remove_multiples(index_pos, id_keys_dict)]=algs[index_pos[0]]

#now we have removed the seqs that have duplicates and effectively merged them together into merged_multiples dict. This
#is then added back to the original dict from which the duplicates had been removed.
len(merged_multiples)
len(alignments)
alignments.update(merged_multiples)

#once you get the indices of the duplicates, get the keys based on these indices, merge the key values
#Now we have a script for detecting duplicates from the data dict, removing them via merging the values
#need to remove the seqs from the original dict using the index_pos approach

#NT
####################
output_path="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html"
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_3p_mcherry_mcherry_SP4_local3.fasta"
aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_3p_mcherry_mcherry_SP4_local2.fasta"
aligner(df_trim_full2, target_sequence, "align_local2", result, output_path, 3,1)
####################
from Bio import SeqIO

alig_fasta="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_3p_mcherry_mcherry_SP4_local2.fasta"
fasta_sequences = SeqIO.parse(open(alig_fasta),'fasta')
with open(alig_fasta) as out_file:
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
def find(lst, a, b):
    return [i for i, x in enumerate(lst) if x<a or x>b]
#find indeces of matching seqs, then take their id rows and sum up the percs (and sds?)
#maybe do this right after merging the dfs?
#AA
####################
corr_frame=1
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/unaligned_seqs/Exp2_3p_mcherry_mcherry_SP4.fasta"
output_html="/media/data/AtteR/projects/hiti/pipeline_output_reorg/AA_aligned_html/Exp2_3p_mcherry_mcherry_SP4_AA.html"
translate_nt_aa_hiti2(result, corr_frame, output_html)
####################
####################



####################
#SP1 3' 
trimmed_df="/media/data/AtteR/projects/hiti/dfs/full_df_trim_mcherry_p3_HITI_SP1.csv"
transgene = 'mCherry'
assay_end = '3p'
read_fwd = True
lliteral=" literal=CGGCGGCATGGACGAGC"
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_3p"
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' #region prior to r2 primer

target_sequence="tgtacaagccgaacGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAGGTCTTCTACCGTCTGGAGAGG"
target_sequence=target_sequence.upper()

#########
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral,read_fwd)
df_trim_full2=calculate_perc_sd2(df_full)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/unaligned_seqs/Exp2_3p_mcherry_mcherry_SP1.fasta"
save_fasta(result, df_trim_full2, target_sequence)
####################
#NT
output_path="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html"
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_3p_mcherry_mcherry_SP1_local3.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_3p_mcherry_mcherry_SP1_local2.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local2", result, output_path, 3,1)
####################

#AA
####################
corr_frame=1
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/unaligned_seqs/Exp2_3p_mcherry_mcherry_SP1.fasta"
output_html="/media/data/AtteR/projects/hiti/pipeline_output_reorg/AA_aligned_html/Exp2_3p_mcherry_mcherry_SP1_AA.html"
translate_nt_aa_hiti2(result, corr_frame, output_html)
####################



#5'
#############
#SP4 5' 
transgene = 'mCherry'
assay_end = '5p'
lliteral = ' literal=GGGCGAGGAGGATAACATGG'
filterlitteral = 'CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAG'
target_sequence = "GTGTCTCCGGTCCCCAAAATCCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACGTCATATGGTCCAGCTCGGgtgagcaagggcgaggaggataacatgg"
target_sequence="GTGTCTCCGGTCCCCAAAATCCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACGTCATATGGTCCAGCTCGGgtgagcaa"
target_sequence=target_sequence.upper()
read_fwd = True
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_5p"
#########
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, read_fwd)
df_trim_full2=calculate_perc_sd2(df_full)
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/unaligned_seqs/Exp2_5p_mcherry_mcherry_SP4.fasta"
save_fasta(result, df_trim_full2, target_sequence)
#NT
####################
output_path="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html"
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_5p_mcherry_mcherry_SP4_local3.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_5p_mcherry_mcherry_SP4_local2.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local2", result, output_path, 3,1)
####################

#AA
####################
corr_frame=1
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/unaligned_seqs/Exp2_5p_mcherry_mcherry_SP4.fasta"
output_html="/media/data/AtteR/projects/hiti/pipeline_output_reorg/AA_aligned_html/Exp2_5p_mcherry_mcherry_SP4_AA.html"
translate_nt_aa_hiti2(result, corr_frame, output_html)
####################



####################
#SP1 5' 
assay_end = '5p'
filterlitteral = 'CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAG'
lliteral = ' literal=GTGTCTCCGGTCCCCAAAAT'
target_sequence= "GTGTCTCCGGTCCCCAAAATCCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTaggctaagaactcctccgcgccaccatggtgagca"
target_sequence=target_sequence.upper()
read_fwd = True
lliteral = ' literal=GGGCGAGGAGGATAACATGG'

base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_5p"
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' #region prior to r2 primer
#########
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, read_fwd)
df_trim_full2=calculate_perc_sd2(df_full)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp2_5p_mcherry_mcherry_SP1.fasta"
save_fasta(result, df_trim_full2, target_sequence)


#NT
####################
output_path="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html"
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_5p_mcherry_mcherry_SP1_local3.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_5p_mcherry_mcherry_SP1_local3.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, 3,1)
####################

#AA
####################
corr_frame=1
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp2_5p_mcherry_mcherry_SP1.fasta"
output_html="/media/data/AtteR/projects/hiti/pipeline_output_reorg/AA_aligned_html/Exp2_5p_mcherry_mcherry_SP1_AA.html"
translate_nt_aa_hiti2(result, corr_frame, output_html)
####################



