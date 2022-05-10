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
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral,read_fwd)
#starcode produces files with exactly same seqs, counts and percs
#files do have different number of reads but for some reason starcode clusters all of them as the same

df_trim_full2=calculate_perc_sd2(df_full)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/unaligned_seqs/Exp2_3p_mcherry_mcherry_SP4.fasta"
save_fasta(result, df_trim_full2, target_sequence)

#NT
####################
output_path="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html"
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_3p_mcherry_mcherry_SP4_local3.fasta"
aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_3p_mcherry_mcherry_SP4_local2.fasta"
aligner(df_trim_full2, target_sequence, "align_local2", result, output_path, 3,1)
####################

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



