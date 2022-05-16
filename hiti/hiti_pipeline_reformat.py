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

LP1 mcherry
5p lliteral GGGCGAGGAGGATAACATGG
3p lliteral CGGCGGCATGGACGAG
CCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCAC

LP1 GFP
5p CGAGGTGAAGTTCGAGGGC
3p GGACGACGGCAACTACAAGA


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
    full_df.sort_values(by=['percent_sum'], ascending=False, inplace=True)
    full_df.head()
    full_df.columns
    #discard seqs that contribute less than 0.0001% percentage
    rows_drop=[]
    full_df[count_cols]
    full_df_trim = full_df.drop(full_df[full_df["total_reads_seq"] <= 3].index)

    #return(full_df_trim.iloc[0:round(len(full_df_trim.index)),:])
    return(full_df_trim)

def create_datadict(base_path, transgene, animal_list):

    #hip_folders = [folder for folder in os.listdir(base_path) if "mCherry" in folder and "h_" in folder or "s_" in folder]
    group_folders = [folder for folder in os.listdir(base_path) if transgene in folder]
    #str_folders = [folder for folder in os.listdir(base_path) if "mCherry" in folder and "s_" in folder]
    group_folders


    def animal_names(animal_list):
        for s in animal_list:
            animal_name="_".join(s.split("_")[:3])
            if int(animal_name.split("_")[0]) in animal_list:
                print(animal_name)
                animals.append("_".join(s.split("_")[:3]))
        #animals=list(set(animals))
        return(sorted(list(((set(animals))))))

    animals = animal_names(animal_list)
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
    else:
        animal_p5 = glob.glob(search_path+'*R2*')
        animal_p7 = glob.glob(search_path+'*R1*')
    

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
    test_file_p5_filter2 = tempfile.NamedTemporaryFile(suffix = '.fastq').name #when cutadapt applied

    cutadapt_call="cutadapt -g "+lliteral+" -o " + test_file_p5_filter2 + " " + test_file_p5_filter
    call([cutadapt_call], shell=True)
    cutadapt_call="cutadapt -a "+rliteral+" -o " + test_file_p5_filter2 + " " + test_file_p5_filter2
    call([cutadapt_call], shell=True)


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
def analyze_all():
    complete_df = pd.DataFrame({'sequence': [target_sequence]})
    for animal in animal_list:
        df_this = trimRead_hiti(animal,base_path,transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd)
        complete_df = pd.merge(complete_df, df_this, on="sequence", how='outer')
    
    complete_df = complete_df.fillna(value=0)
    perc_cols = [col for col in complete_df.columns if 'percent' in col]
    #complete_df['percent_sum'] = complete_df[perc_cols].sum(axis=1)
    export_csv = export_path+transgene+'_'+assay_end+'.csv'
    complete_df.to_csv(export_csv, index=False)
    return complete_df

#==EXP1==#
#GFP 3p

#############
transgene = 'GFP'
assay_end = '3p'
animal_list = [20, 21, 22, 23, 24] 

read_fwd = True
filterlitteral = 'GGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGG'
lliteral = ' literal=GGACGACGGCAACTACAAGA'
rliteral = ' literal=GGAGCTGGACCATATGACCAC'
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path = '/media/data/AtteR/projects/hiti/pipeline_output_reorg/'
target_sequence = "CCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGGtcggtgctgcggctccgcggagccgcagcaccgacgaccagAT"
target_sequence=target_sequence.upper()
full_df=analyze_all()

data_dict=create_datadict(base_path,transgene, animal_list)
full_df=import_reads_process(data_dict, transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd)
full_df_trim_orig=calculate_perc_sd2(full_df)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_3p_GFP_LP1_TB.fasta"
save_fasta(result, full_df_trim_orig, target_sequence)

#########
#NT
####################
output_path="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html"
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/exp1/Exp1_3p_GFP_LP1_local3.fasta"
test_res=aligner(full_df_trim_orig, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/exp1/Exp1_3p_GFP_LP1_local2_TB.fasta"
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

#GFP 5p
#############
transgene = 'GFP'
assay_end = '5p'
read_fwd = True
animal_list = [13, 14, 15, 16, 17, 18] 
filterlitteral = 'GCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTCCTGGACG'
lliteral = ' literal=GCTTCTGCCTCAGAGGAGTTCT'
rliteral = ' literal=CGAGGTGAAGTTCGAGGGC'
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path = '/media/data/AtteR/projects/hiti/pipeline_output_reorg/'
target_sequence = "tagcctgttaacaggCGCGCCACCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGC"
target_sequence=target_sequence.upper()
#########
full_df=analyze_all()

data_dict=create_datadict(base_path, transgene, animal_list)
full_df=import_reads_process(data_dict, transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd)
full_df_trim_orig=calculate_perc_sd2(full_df)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_5p_GFP_LP1_TB.fasta"
save_fasta(result, full_df_trim_orig, target_sequence)

#NT
####################
output_path="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html"
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_5p_GFP_LP1_local3.fasta"
test_res=aligner(full_df_trim_orig, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_5p_GFP_LP1_local2_TB.fasta"
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


'''
mCherry LP1
3p lliteral CGGCGGCATGGACGAG rliteral GTCTTCTACCGTCTGGAGAGG
CTGTACAAGGtcggtgctgcggctccgcggagccgcagcaccgacgaccagATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAG

5p lliteral ---- rliteral GGGCGAGGAGGATAACATGG
gcagagctcaagcgagttctcccgcagccgcagtctctgggcctctctagcttcagcggcgacgagcctgccacactcgctaagctcctccggcaccgcacacctgccactgccgctgcagccgccggctctgctcccttccggcttctgcctcagaggagttcttagcctgttaacaggCGCGCCACCATGGTGAGCAA

'''

#mCherry 3p
#############
transgene = 'mCherry'
assay_end = '3p'
animal_list = [7, 8, 9, 10, 11, 12]

read_fwd = True
filterlitteral='CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' 
lliteral=' literal=CGGCGGCATGGACGAG' #added C at the end as its part of tje primer #to check that its on target with mcherry - the primer
rliteral=' literal=GTCTTCTACCGTCTGGAGAGG'
#base_path = '/home/lcadmin/mnm-lts/SequencingFiles/arc_hiti_asyn/HITI1-8977973/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
#export_path = '/home/lcadmin/mnm-lts/HITI-analysis/'
base_path='/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path='/media/data/AtteR/projects/hiti/pipeline_output_reorg/'
#with primer
target_sequence="CTGTACAAGGtcggtgctgcggctccgcggagccgcagcaccgacgaccagATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAG"
target_sequence=target_sequence.upper()
full_df=analyze_all()

data_dict=create_datadict(base_path, transgene, animal_list)
full_df=import_reads_process(data_dict, transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd)
full_df_trim_orig=calculate_perc_sd2(full_df)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_3p_mcherry_LP1_TB.fasta"
save_fasta(result, full_df_trim_orig, target_sequence)

#########

#NT
####################
output_path="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html"
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_3p_mcherry_LP1_local3.fasta"
test_res=aligner(full_df_trim_orig, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_3p_mcherry_LP1_local2_TB.fasta"
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


#mCherry 5p
############
transgene='mCherry'
assay_end='5p'
read_fwd=True
filterlitteral='CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' #region prior to r2 primer
#lliteral = ' literal=GGGCGAGGAGGATAACATGG' 
rliteral = ' literal=GGGCGAGGAGGATAACATGG'
#base_path = '/home/lcadmin/mnm-lts/SequencingFiles/arc_hiti_asyn/HITI1-8977973/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
#export_path = '/home/lcadmin/mnm-lts/HITI-analysis/'
base_path = '/media/data/AtteR/projects/hiti/FASTQ_Generation_2020-03-09_08_30_27Z-13364364/'
export_path = '/media/data/AtteR/projects/hiti/pipeline_output_reorg/'
animal_list = [1, 2, 3, 4, 5, 6] 
#target_sequence = "CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTaacaggCGCGCCACCATGGTGAGCAAGGGCGAGGAGGATAACATGG"
target_sequence = "CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTaacaggCGCGCCACCATGGTGAGCAA"
target_sequence = "gcagagctcaagcgagttctcccgcagccgcagtctctgggcctctctagcttcagcggcgacgagcctgccacactcgctaagctcctccggcaccgcacacctgccactgccgctgcagccgccggctctgctcccttccggcttctgcctcagaggagttcttagcctgttaacaggCGCGCCACCATGGTGAGCAA"
target_sequence = target_sequence.upper()
full_df=analyze_all()
full_df_trim_orig=calculate_perc_sd2(full_df)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_5p_mcherry_LP1_TB.fasta"
save_fasta(result, full_df_trim_orig, target_sequence)
full_df_trim_orig['percent_sum']
#########
#NT
####################
output_path="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html"
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_5p_mcherry_LP1_local3.fasta"
aligner(full_df_trim_orig, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp1_5p_mcherry_LP1_local2_TB.fasta"
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

'''
SP1 ref
ccatcatcaaggagttcatgcgcttcaaggtgcacatggagggctccgtgaacggccacgagttcgagatcgagggcgagggcgagggccgcccctacgagggcacccagaccgccaagctgaaggtgaccaagggtggccccctgcccttcgcctgggacatcctgtcccctcagttcatgtacggctccaaggcctacgtgaagcaccccgccgacatccccgactacttgaagctgtccttccccgagggcttcaagtgggagcgcgtgatgaacttcgaggacggcggcgtggtgaccgtgacccaggactcctccctgcaggacggcgagttcatctacaaggtgaagctgcgcggcaccaacttcccctccgacggccccgtaatgcagaagaagaccatgggctgggaggcctcctccgagcggatgtaccccgaggacggcgccctgaagggcgagatcaagcagaggctgaagctgaaggacggcggccactacgacgctgaggtcaagaccacctacaaggccaagaagcccgtgcagctgcccggcgcctacaacgtcaacatcaagttggacatcacctcccacaacgaggactacaccatcgtggaacagtacgaacgcgccgagggccgccactccac
3p lliteral cggcggcatggacgag
5p lliteral gggcgaggaggataacatgg

SP4 ref
ccatcatcaaggagttcatgcgcttcaaggtgcacatggagggctccgtgaacggccacgagttcgagatcgagggcgagggcgagggccgcccctacgagggcacccagaccgccaagctgaaggtgaccaagggtggccccctgcccttcgcctgggacatcctgtcccctcagttcatgtacggctccaaggcctacgtgaagcaccccgccgacatccccgactacttgaagctgtccttccccgagggcttcaagtgggagcgcgtgatgaacttcgaggacggcggcgtggtgaccgtgacccaggactcctccctgcaggacggcgagttcatctacaaggtgaagctgcgcggcaccaacttcccctccgacggccccgtaatgcagaagaagaccatgggctgggaggcctcctccgagcggatgtaccccgaggacggcgccctgaagggcgagatcaagcagaggctgaagctgaaggacggcggccactacgacgctgaggtcaagaccacctacaaggccaagaagcccgtgcagctgcccggcgcctacaacgtcaacatcaagttggacatcacctcccacaacgaggactacaccatcgtggaacagtacgaacgcgccgagggccgccactccac
3p lliteral cggcggcatggacgag
5p lliteral gggcgaggaggataacatgg
---------------------------------------------------------
GFP
3p lliteral GGACGACGGCAACTACAAGA   rliteral GGAGCTGGACCATATGACCAC
CCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGGtcggtgctgcggctccgcggagccgcagcaccgacgaccagAT

5p lliteral GCTTCTGCCTCAGAGGAGTTCT rliteral CGAGGTGAAGTTCGAGGGC
tagcctgttaacaggCGCGCCACCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGC

mCherry LP1
3p lliteral CGGCGGCATGGACGAG rliteral GTCTTCTACCGTCTGGAGAGG
CTGTACAAGGtcggtgctgcggctccgcggagccgcagcaccgacgaccagATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAG

5p lliteral ---- rliteral GGGCGAGGAGGATAACATGG
gcagagctcaagcgagttctcccgcagccgcagtctctgggcctctctagcttcagcggcgacgagcctgccacactcgctaagctcctccggcaccgcacacctgccactgccgctgcagccgccggctctgctcccttccggcttctgcctcagaggagttcttagcctgttaacaggCGCGCCACCATGGTGAGCAA
-------------------------------
Exp2
SP1
3p
ctgtacaagccgaacgttcggagccgcagcaccgacgaccagATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAG
lliteral cggcggcatggacgag rliteral GTCTTCTACCGTCTGGAGAGG
5p
gcagagctcaagcgagttctcccgcagccgcagtctctgggcctctctagcttcagcggcgacgagcctgccacactcgctaagctcctccggcaccgcacacctgccactgccgctgcagccgccggctctgctcccttccggcttctgcctcagaggagttcttagcctaggctaagaactcctccgcgccaccatggtgagcaa
lliteral gggcgaggaggataacatgg


SP4
5p
gcagagctcaagcgagttctcccgcagccgcagtctctgggcctctctagcttcagcggcgacgagcctgccacactcgctaagctcctccggcaccgcacacctgccactgccgctgcagccgccggctctgctcccttccggcttctgcctcagaggagttcttagcctgttcggagccgcagcaccgacgaccagATGGAGCTGGACCATATGACGTCATATGGTCCAGCTCGGgtgagcaa
lliteral 
rliteral gggcgaggaggataacatgg

3p
ctgtacaagATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCG
lliteral cggcggcatggacgag rliteral CTGGGTCAAGCGTGAGATG

so from amplicons trim both sides lliteral and rliteral
'''
#==EXP2==#
'''
SP4

3p
ctgtacaagATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCG
lliteral cggcggcatggacgag rliteral CTGGGTCAAGCGTGAGATG

'''
#############
#SP4 3' 
trimmed_df="/media/data/AtteR/projects/hiti/dfs/full_df_trim_mcherry_p3_HITI_SP4.csv"
transgene = 'mCherry'
assay_end = '3p'
read_fwd = True
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_3p"
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' #region prior to r2 primer
lliteral=" literal=CGGCGGCATGGACGAG"
rliteral=" literal=CTGGGTCAAGCGTGAGATG"

target_sequence="ctgtacaagATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCG"
target_sequence=target_sequence.upper()
#########
#########

#The counts and percs for each animal regarding cluster seq are identical

#########
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, read_fwd)
#starcode produces files with exactly same seqs, counts and percs
#files do have different number of reads but for some reason starcode clusters all of them as the same
df_trim_full2=calculate_perc_sd2(df_full)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/unaligned_seqs/Exp2_3p_mcherry_mcherry_SP4_var.fasta"
save_fasta(result, df_trim_full2, target_sequence)


#NT
####################
output_path="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html"
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_3p_mcherry_mcherry_SP4_local3.fasta"
aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_3p_mcherry_mcherry_SP4_local2.fasta"
aligner(df_trim_full2, target_sequence, "align_local2", result, output_path, 3,1)
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
'''
SP1
3p
ctgtacaagccgaacgttcggagccgcagcaccgacgaccagATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAG
lliteral cggcggcatggacgag rliteral GTCTTCTACCGTCTGGAGAGG

'''
trimmed_df="/media/data/AtteR/projects/hiti/dfs/full_df_trim_mcherry_p3_HITI_SP1.csv"
transgene = 'mCherry'
assay_end = '3p'
read_fwd = True
lliteral=" literal=CGGCGGCATGGACGAGC"
rliteral=" literal=GTCTTCTACCGTCTGGAGAGG"

base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_3p"
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' #region prior to r2 primer

target_sequence="ctgtacaagccgaacgttcggagccgcagcaccgacgaccagATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAG"
target_sequence=target_sequence.upper()

#########
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral,rliteral,read_fwd)
df_trim_full2=calculate_perc_sd2(df_full)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/unaligned_seqs/Exp2_3p_mcherry_mcherry_SP1.fasta"
save_fasta(result, df_trim_full2, target_sequence)
####################
df_trim_full2['percent_sum']
#NT
output_path="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html"
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_3p_mcherry_mcherry_SP1_local3.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_3p_mcherry_mcherry_SP1_local2.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local2", result, output_path, 3,1)
NT_and_perc={}
for record in SeqIO.parse(result, "fasta"):
    print(record.seq)
    NT_and_perc[">"+str(record.description)]=lliteral.split("=")[1] + "|" + str(Seq(record.seq)) + "|" + rliteral.split("=")[1]
print(NT_and_perc)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_3p_mcherry_mcherry_SP1_local2_prim.fasta"

with open(result, "w") as f: 
    for id, seq_prims in NT_and_perc.items():
        print(id)
        f.write(id + "\n")
        f.write(seq_prims + "\n")

def addprimers(fasta_alig, lliteral,rliteral):
    NT_and_perc={}
    for record in SeqIO.parse(fasta_alig, "fasta"):
        NT_and_perc[">"+str(record.description)]=lliteral.split("=")[1] + "|" + str(Seq(record.seq)) + "|" + rliteral.split("=")[1]
    print(NT_and_perc)
    with open(fasta_alig, "w") as f: 
        for id, seq_prims in NT_and_perc.items():
            f.write(id + "\n")
            f.write(seq_prims + "\n")

#add primers to both ends
addprimers(result, lliteral, rliteral)
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
'''
SP4
5p
gcagagctcaagcgagttctcccgcagccgcagtctctgggcctctctagcttcagcggcgacgagcctgccacactcgctaagctcctccggcaccgcacacctgccactgccgctgcagccgccggctctgctcccttccggcttctgcctcagaggagttcttagcctgttcggagccgcagcaccgacgaccagATGGAGCTGGACCATATGACGTCATATGGTCCAGCTCGGgtgagcaa
lliteral 
rliteral gggcgaggaggataacatgg

'''
transgene = 'mCherry'
assay_end = '5p'
lliteral = ' literal='
rliteral = ' literal=GGGCGAGGAGGATAACATGG'
filterlitteral = 'CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAG'
filterlitteral='AGCCTGTTAACAGGCGCGCCACCATGGTGAGCAAGGGCGAGGAGGATAACATGG'
target_sequence = "GTGTCTCCGGTCCCCAAAATCCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACGTCATATGGTCCAGCTCGGgtgagcaagggcgaggaggataacatgg"
target_sequence="gcagagctcaagcgagttctcccgcagccgcagtctctgggcctctctagcttcagcggcgacgagcctgccacactcgctaagctcctccggcaccgcacacctgccactgccgctgcagccgccggctctgctcccttccggcttctgcctcagaggagttcttagcctgttcggagccgcagcaccgacgaccagATGGAGCTGGACCATATGACGTCATATGGTCCAGCTCGGgtgagcaa"
target_sequence=target_sequence.upper()
read_fwd = True
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_5p"
#########
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, rliteral, read_fwd)
df_trim_full2=calculate_perc_sd2(df_full)
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/unaligned_seqs/Exp2_5p_mcherry_mcherry_SP4.fasta"
save_fasta(result, df_trim_full2, target_sequence)
#NT
####################
output_path="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html"
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_5p_mcherry_mcherry_SP4_local3.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_5p_mcherry_mcherry_SP4_local2.fasta"
#gop and gep got rid of an annoying shift in NTs which made part of the alignments off so changed to 4 and
aligner(df_trim_full2, target_sequence, "align_local2", result, output_path, 4,2)
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
'''
5p
gcagagctcaagcgagttctcccgcagccgcagtctctgggcctctctagcttcagcggcgacgagcctgccacactcgctaagctcctccggcaccgcacacctgccactgccgctgcagccgccggctctgctcccttccggcttctgcctcagaggagttcttagcctaggctaagaactcctccgcgccaccatggtgagcaa
rliteral gggcgaggaggataacatgg
'''
assay_end = '5p'
filterlitteral = 'CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAG'
rliteral = ' literal=GGGCGAGGAGGATAACATGG'
lliteral = ' literal='

target_sequence= "gcagagctcaagcgagttctcccgcagccgcagtctctgggcctctctagcttcagcggcgacgagcctgccacactcgctaagctcctccggcaccgcacacctgccactgccgctgcagccgccggctctgctcccttccggcttctgcctcagaggagttcttagcctaggctaagaactcctccgcgccaccatggtgagcaa"
target_sequence=target_sequence.upper()
read_fwd = True

base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_5p"
#########
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral,rliteral, read_fwd)
df_trim_full2=calculate_perc_sd2(df_full)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp2_5p_mcherry_mcherry_SP1.fasta"
save_fasta(result, df_trim_full2, target_sequence)


#NT
####################
output_path="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html"
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_5p_mcherry_mcherry_SP1_local3.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/2Exp2_5p_mcherry_mcherry_SP1_local2.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local2", result, output_path, 4,2)
####################
fasta_alig="/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_alig_fastas/Exp2_3p_mcherry_mcherry_SP1_local2.fasta"


ref_dict={}
NT_and_perc={}
for record in SeqIO.parse(fasta_alig, "fasta"):
    if record.id=="0":
        ref_dict[">"+str(record.description)]=str(Seq(record.seq))
    else:
        NT_and_perc[">"+str(record.description)]=str(Seq(record.seq))

NT_and_perc['>32CluSeq:0.02716_var:3e-05_sd:0.00548']
f = open("/media/data/AtteR/projects/hiti/pipeline_output_reorg/NT_aligned_html/new_format.html","w")
rliteral = ' literal=GGGCGAGGAGGATAACATGG'
lliteral = ' literal=GGGCGAGGAGGATAACATGG'

#you'll also need to check if the upper and lower alignments have the same NT, if not, give a faulty color. 
color_scheme={"RED": [["#FF0000"], ["F"]], "BLUE":[["#6495ED"],["A"]], "GREEN":[["#9FE2BF"], ["G"]], "ORANGE":[["#FF7F50"],["C"]], "MAGENTA":[["#DE3163"], ["T"]], "BLACK": [["#000000"]]}
for id, alig_nt in list(NT_and_perc.items())[:80]:
    aa2_list=[]
    aa1_list=[]
    alig_nt=lliteral.split("=")[1] + "|" + alig_nt + "|" + rliteral.split("=")[1]
    ref=lliteral.split("=")[1] + "|" + ref_dict['>0 Ref_seq'] + "|" + rliteral.split("=")[1]

    f.write("<html> \n <body> <p>"+ id +"</p> <p>")
    #test_html.append('<!DOCTYPE html><html> <head><link rel="stylesheet" href="format.css"></head> <meta charset="utf-8"> <body> <p>'+ "_".join(info.split("_")[:2]) + "<br>" + frame +'</p> <p>')
#go over the individual AAs of the alignemnts
    #write the html colouring into a string var for each aa and append into a list, then join and make a massive string,
    # then after all AAs have been iterated over for the certain seq, then write into the html 
    for nt1,nt2, in zip(alig_nt,ref):
        if nt1!=nt2 and nt1.isalpha() and nt2.isalpha():
            aa2_list.append('<span style="color:'+ color_scheme[find_colour(nt2,color_scheme)][0][0]+ '">' + nt2 + '</span>')
            aa1_list.append('<span style="color:'+ color_scheme[find_colour("F",color_scheme)][0][0]+ '">' + nt1 + '</span>')
        aa2_list.append('<span style="color:'+ color_scheme[find_colour(nt2,color_scheme)][0][0]+ '">' + nt2 + '</span>')
        aa1_list.append('<span style="color:'+ color_scheme[find_colour(nt1,color_scheme)][0][0]+ '">' + nt1 + '</span>')

        #print("============")
    coloured_ref="".join(aa2_list)

    coloured_seq="".join(aa1_list)

    f.write(coloured_ref +"<br>")
    f.write(coloured_seq)
    f.write("</p>")

f.write("</body></html>")
f.close()

#atm if theres a mismatch with NTs in rows colours the first NT mismatch is coloured correctly as red but the next one is not!


#AA
####################
corr_frame=1
result="/media/data/AtteR/projects/hiti/pipeline_output_reorg/fastas/Exp2_5p_mcherry_mcherry_SP1.fasta"
output_html="/media/data/AtteR/projects/hiti/pipeline_output_reorg/AA_aligned_html/Exp2_5p_mcherry_mcherry_SP1_AA.html"
translate_nt_aa_hiti2(result, corr_frame, output_html)
####################



