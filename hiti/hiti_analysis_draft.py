import pandas as pd
import os 
os.getcwd()
os.chdir("/media/data/AtteR/projects/hiti")
from scripts_hiti import *


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

#############
#SP4 3' 
trimmed_df="/media/data/AtteR/projects/hiti/dfs/full_df_trim_mcherry_p3_HITI_SP4.csv"
transgene = 'mCherry'
assay_end = '3p'
read_fwd = True
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_3p"
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' #region prior to r2 primer
lliteral = ' literal=GGCGGCATGGACGAGC' #added C at the end as its part of tje primer #to check that its on target with mcherry - the primer

target_sequence="tgtacaagATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAGGTCTTCTACCGTCTGGAGAGG"
target_sequence=target_sequence.upper()
#########
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral,read_fwd)
df_animal=[]
complete_df = pd.DataFrame({'sequence': [ref]})

for read in os.listdir(base_path):
    animal_group_name=read.split("_")[3] + "_" + read.split("_")[4]
    if "R1" in read:
        animal_p5_cat = tempfile.NamedTemporaryFile(suffix = '.fastq.gz').name
        animal_p7_cat = tempfile.NamedTemporaryFile(suffix = '.fastq.gz').name
        test_file_p5_out = tempfile.NamedTemporaryFile(suffix = '.fastq').name
        test_file_p7_out = tempfile.NamedTemporaryFile(suffix = '.fastq').name
        test_file_p5_filter = tempfile.NamedTemporaryFile(suffix = '.fastq').name#when bbduk applied


        if read_fwd:
            animal_p5 = glob.glob(base_path+'/*R1*')
            animal_p7 = glob.glob(base_path+'/*R2*')
            #display('Forward run Animal: '+animal_nr)
        else:
            print("===")
            break
            # animal_p5 = glob.glob(search_path+'/*R2*')
            # animal_p7 = glob.glob(search_path+'/*R1*')
            #display('Reverse run Animal: '+animal_nr)

        cat_p5= "cat "+" ".join(animal_p5)+" > "+animal_p5_cat
        print(cat_p5)
        call([cat_p5], shell=True) #call caused the terminal to freeze so switched to os
        cat_p7= "cat "+" ".join(animal_p7)+" > "+animal_p7_cat
        call([cat_p7], shell=True)

        kmer = '20'
        hdist = '3'
        param=" k="+kmer+" hdist="+hdist+" rcomp=f skipr2=t threads=32 overwrite=true"

        #to check if the read is an amplicon
        call_sequence = "/media/data/AtteR/Attes_bin/bbmap/bbduk.sh in="+animal_p7_cat+" in2="+animal_p5_cat+" outm1="+test_file_p7_out+" outm2="+test_file_p5_out+" literal="+filterlitteral+param
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
        df_animal.append(df)
        print(animal_group_name + " done!")
        df["count"]=df.sum(axis=1) #make a column with total count sum of reads and remove the rest. This gives a df that has the seqs and the total counts from all lanes
        df.drop(df.iloc[:, 1:((len(df.columns)-1))], inplace = True, axis = 1)

        #Once you have combined all the lane dfs, then you take the percentage
        total_counts = int(df[['count']].sum())
        df['percent'] = (df['count'] / total_counts)
        df = df.rename(columns={'percent':animal_group_name+'_percent','count':animal_group_name+'_count',})

        complete_df = pd.merge(complete_df, df, on="sequence", how='outer')

#starcode produces files with exactly same seqs, counts and percs
#files do have different number of reads but for some reason starcode clusters all of them as the same

df_trim_full2=calculate_perc_sd2(df_full)
df_full.columns
df_full['20_S16_percent']
df_full['19_S15_percent']

df_trim_full2.head()
df_trim_full2.columns
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p3_seq_SP4_v2.fasta"
save_fasta(result, df_trim_full2, target_sequence)
#NT
output_path="/media/data/AtteR/projects/hiti/pipeline_output/NT_aligned"
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p3_seq_SP4_local3_v2.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p3_seq_SP4_local2_v2.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local2", result, output_path, 3,1)
#AA
corr_frame=1
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p3_seq_SP4_v2.fasta"

df_aa=translate_nt_aa(result, 1)
df_aa2=translate_nt_aa_hiti2(result,1)
df_aa2
df_aa.head()
output_html="/media/data/AtteR/projects/hiti/pipeline_output/AA_aligned/AA_3p_mcherry_coloured_alignment_SP4_v2_fr0.html"
visualise_aa_hybrid_alignments(df_aa, output_html)
#########

#############


#############
#SP1 3' 
trimmed_df="/media/data/AtteR/projects/hiti/dfs/full_df_trim_mcherry_p3_HITI_SP1.csv"
transgene = 'mCherry'
assay_end = '3p'
read_fwd = True
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_3p"
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' #region prior to r2 primer

target_sequence="tgtacaagccgaacGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAGGTCTTCTACCGTCTGGAGAGG"
target_sequence=target_sequence.upper()
#########
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral,read_fwd)
df_trim_full2=calculate_perc_sd2(df_full)
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p3_seq_SP1.fasta"
save_fasta(result, df_trim_full2, target_sequence)
#NT
output_path="/media/data/AtteR/projects/hiti/pipeline_output/NT_aligned"
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p3_seq_SP1_local3.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p3_seq_SP1_local2.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local2", result, output_path, 3,1)
#AA
corr_frame=1
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p3_seq_SP1.fasta"

out_of_frames=[0,1,2]
out_of_frames.remove(corr_frame)
refs_aa_frames={}
aa_and_perc={}
for record in SeqIO.parse(result, "fasta"):
    if record.id=="0":
        refs_aa_frames["Frame:" + str(corr_frame)]=str(Seq(record.seq[corr_frame:]).translate())
    else:
        aa_and_perc[">"+str(record.description) + "_transl.frame:" + str(corr_frame)]=str(Seq(record.seq[corr_frame:]).translate())

aa_and_perc
refs_aa_frames

#overall the code takes the inte ref seq, translates it, then the amplicon and translates it
"VQAERSEPQHRRPDGAGPYDHRRPPRLPCPAGWAGRQTQCDPADW*VPS*DAGTRTEDPPASVDRSVQAGGARAERVAQVGGQAGEQLGRLRAHRRLTALEEVHQGLSLPLPGDHRQPGALGQA*DARVEGGLLPSGE"
"VQAEQFGAAAPTTRWSWTI"

r="TGTACAAGCCGAACGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAGGTCTTCTACCGTCTGGAGAGG"
seq="TGTACAAGCCGAACAGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATA"
aa=Seq(seq[1:]).translate()
aa
#you go over the ref seqs in different frames and align all the amplicons to them. save the alignment into the df's specific column. initially
#save into a list or such
ref_x_alignment={}
refs_aa_frames
alignments_per_ref=[]
refs_aa_frames['Frame:1']
for ampl in aa_and_perc.keys():
    matches=SequenceMatcher(None,refs_aa_frames['Frame:1'],aa_and_perc[ampl])
    seqs=[]
    #you use range_line so that when you fill the remnants from left side of the match, you wont keep adding from
    #beginning since in the end, we merge the seq list elements into the whole alignment of the amplicon against the ref
    range_line=0
    for i in range(len(matches.get_matching_blocks())):
        match=matches.get_matching_blocks()[i]
        seqs.append(len(refs_aa_frames['Frame:1'][range_line:match.a])*"-"+str(aa_and_perc[ampl])[match.b:match.b+match.size])
        range_line=match.a+match.size
    alignments_per_ref.append(''.join(seqs))
ref_x_alignment['Frame:1' + "|Ref:" +refs_aa_frames['Frame:1']]=alignments_per_ref
seq_info={"Seq_info:":aa_and_perc.keys()}
keys=list(aa_and_perc.keys())
ref_x_alignment.keys()
list(aa_and_perc.values())[0]

seq_info=["Seq_info"]+list(aa_and_perc.keys())
ref_x_alig_list=[]
for keys, values in ref_x_alignment.items():
    #make into a list with first being the ref, the rest being the aligned seqs. 
    ref_x_alig_list.append([keys]+list(values))

#so we have a list containing sublist pairs of the ref seq
ref_x_alig_list[0][0]
df = pd.DataFrame(data= {seq_info[0]: seq_info[1:], ref_x_alig_list[0][0]:ref_x_alig_list[0][1:]})
output_html="/media/data/AtteR/projects/hiti/pipeline_output/AA_aligned/AA_3p_mcherry_coloured_alignment_SP1_NOSCAR.html"
df.iloc[:,1]
aa_coloring(list(df['Seq_info']), list(df.iloc[:,1]), df.columns[1].split(":")[2], str(corr_frame), output_html)


#(2910-2210)%3
df_aa=translate_nt_aa(result, 1)
df_aa2=translate_nt_aa_hiti2(result,1)
output_html="/media/data/AtteR/projects/hiti/pipeline_output/AA_aligned/AA_3p_mcherry_coloured_alignment_SP1.html"
visualise_aa_hybrid_alignments(df_aa, output_html)
#########

############
#OLD
############
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
#target_sequence = "GGCGGCATGGACGAGCTGTACAAGGTCGGTGCTGCGGCTCCGCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAGGTCTTCTACCGTCTGGAGAGG"
target_sequence = "TGTACAAGGTCGGTGCTGCGGCTCCGCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAGGTCTTCTACCGTCTGGAGAGG"
data_dict=create_datadict(base_path)
data_dict
df_trim_full2=import_reads_process(data_dict, transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd)
df_trim_full2=calculate_perc_sd2(df_full)
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p3_seq_Exp1.fasta"
save_fasta(result, df_trim_full2, target_sequence)
#NT
output_path="/media/data/AtteR/projects/hiti/pipeline_output/NT_aligned"
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p3_seq_local3_Exp1.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p3_seq_SP1_local2_Exp1.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local2", result, output_path, 3,1)
#AA
corr_frame=1
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p3_seq_Exp1.fasta"
#(2910-2210)%3
df_aa=translate_nt_aa(result, 1)
output_html="/media/data/AtteR/projects/hiti/pipeline_output/AA_aligned/AA_3p_mcherry_coloured_alignment_Exp1.html"
visualise_aa_hybrid_alignments(df_aa, output_html)




#5'
#############
#SP4 5' 
transgene = 'mCherry'
assay_end = '5p'
filterlitteral = 'CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAG'
target_sequence = "CCATGTTATCCTCCTCGCCCTTGCTCACCATGGTGGCGCGCCTGTTAACAGGCTAAGAACTCCTCTGAGGCAGAAGCCGGAAGGGAGCAGAGCCGGCGGCTGCAGCGGCAGTGGCAGGTGTGCGGTGCCGGAGGAGCTTAGCGAGTGTGGCAGGCTCGTCGCCGCTGAAGCTAGAGAGGCCCAGAGACTGCGGCTGCGGGAGAACTCGCTTGAGCTCTGCACCGAAACCGCCACCAGCGGCTATTTATGCTGCGCGGGGCCCGTGGGCGGCAGCTCGTGGAGCCCTGGCTCCCATTGGCGGCGCCCAGGCTCCGCGCCAGAGCCCGCCCGCGCCCTCTCCTGCCGCAGGCCCCGCCTTCCGCGCCTACTCGCTCCCCTCCCGTGGGTGCCCTCAAGGACCCGCCCCTTCACAGCCCCGAGTGACTAATGTGCTCTGCTGCGCGCCTCCCACCGGGAGGGATTTTGGGGACCGGAGACAC"
target_sequence = "GTGTCTCCGGTCCCCAAAATCCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACGTCATATGGTCCAGCTCGGgtgagcaagggcgaggaggataacatgg"
#target_sequence= "GTGTCTCCGGTCCCCAAAATCCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACGTCATATGGTCCAGCTCGGgtgagcaagggcgaggaggataacatgg"
target_sequence=target_sequence.upper()
'''
CCATGTTATCCTCCTCGCCCTTGCTCACCATGGTGGCGCGCCTGTTAACAGGCTAAGAACTCCTCTGAGGCAGAAGCCGGAAGGGAGCAGAGCCGGCGGCTGCAGCGGCAGTGGCAGGTGTGCGGTGCCGGAGGAGCTTAGCGAGTGTGGCAGGCTCGTCGCCGCTGAAGCTAGAGAGGCCCAGAGACTGCGGCTGCGGGAGAACTCGCTTGAGCTCTGCACCGAAACCGCCACCAGCGGCTATTTATGCTGCGCGGGGCCCGTGGGCGGCAGCTCGTGGAGCCCTGGCTCCCATTGGCGGCGCCCAGGCTCCGCGCCAGAGCCCGCCCGCGCCCTCTCCTGCCGCAGGCCCCGCCTTCCGCGCCTACTCGCTCCCCTCCCGTGGGTGCCCTCAAGGACCCGCCCCTTCACAGCCCCGAGTGACTAATGTGCTCTGCTGCGCGCCTCCCACCGGGAGGGATTTTGGGGACCGGAGACAC
TGTACAAGGTCGGTGCTGCGGCTCCGCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAGGTCTTCTACCGTCTGGAGAGG
'''
target_sequence=target_sequence.upper()
read_fwd = True
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP4_5p"
#########
df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral,read_fwd)
df_trim_full2=calculate_perc_sd2(df_full)
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p5_seq_SP4.fasta"
save_fasta(result, df_trim_full2, target_sequence)
#NT
output_path="/media/data/AtteR/projects/hiti/pipeline_output/NT_aligned"
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p3_seq_SP5_local3.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p3_seq_SP5_local2.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local2", result, output_path, 3,1)
#AA
corr_frame=1
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p5_seq_SP4.fasta"

df_aa=translate_nt_aa(result, 1)
output_html="/media/data/AtteR/projects/hiti/pipeline_output/AA_aligned/AA_5p_mcherry_coloured_alignment_SP4.html"
visualise_aa_hybrid_alignments(df_aa, output_html)
#########


#############


#############
#SP1 5' 
assay_end = '5p'
filterlitteral = 'CCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAG'
#target_sequence = "CCATGTTATCCTCCTCGCCCTTGCTCACCATGGTGGCGCGCCTGTTAACAGGCTAAGAACTCCTCTGAGGCAGAAGCCGGAAGGGAGCAGAGCCGGCGGCTGCAGCGGCAGTGGCAGGTGTGCGGTGCCGGAGGAGCTTAGCGAGTGTGGCAGGCTCGTCGCCGCTGAAGCTAGAGAGGCCCAGAGACTGCGGCTGCGGGAGAACTCGCTTGAGCTCTGCACCGAAACCGCCACCAGCGGCTATTTATGCTGCGCGGGGCCCGTGGGCGGCAGCTCGTGGAGCCCTGGCTCCCATTGGCGGCGCCCAGGCTCCGCGCCAGAGCCCGCCCGCGCCCTCTCCTGCCGCAGGCCCCGCCTTCCGCGCCTACTCGCTCCCCTCCCGTGGGTGCCCTCAAGGACCCGCCCCTTCACAGCCCCGAGTGACTAATGTGCTCTGCTGCGCGCCTCCCACCGGGAGGGATTTTGGGGACCGGAGACAC"
target_sequence= "GTGTCTCCGGTCCCCAAAATCCCTCCCGGTGGGAGGCGCGCAGCAGAGCACATTAGTCACTCGGGGCTGTGAAGGGGCGGGTCCTTGAGGGCACCCACGGGAGGGGAGCGAGTAGGCGCGGAAGGCGGGGCCTGCGGCAGGAGAGGGCGCGGGCGGGCTCTGGCGCGGAGCCTGGGCGCCGCCAATGGGAGCCAGGGCTCCACGAGCTGCCGCCCACGGGCCCCGCGCAGCATAAATAGCCGCTGGTGGCGGTTTCGGTGCAGAGCTCAAGCGAGTTCTCCCGCAGCCGCAGTCTCTGGGCCTCTCTAGCTTCAGCGGCGACGAGCCTGCCACACTCGCTAAGCTCCTCCGGCACCGCACACCTGCCACTGCCGCTGCAGCCGCCGGCTCTGCTCCCTTCCGGCTTCTGCCTCAGAGGAGTTCTTAGCCTaggctaagaactcctccgcgccaccatggtgagcaagggcgaggaggataacatgg"
target_sequence=target_sequence.upper()
read_fwd = True
lliteral = ' literal=TTATCCTCCTCGCCC'
base_path="/media/data/AtteR/projects/hiti/220426_NB502004_0185_AHKVHYAFX3_HITI-only/SP1_5p"
filterlitteral = 'CTCCCTCCACACGTGCATCTCACGCTTGACCCAGCGCTCCAGGTTGGCGATGGT' #region prior to r2 primer
#########

df_full=import_reads_process_mini(base_path, target_sequence, filterlitteral, lliteral, read_fwd)
df_trim_full2=calculate_perc_sd2(df_full)
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p5_seq_PS1.fasta"
save_fasta(result, df_trim_full2, target_sequence)
#NT
output_path="/media/data/AtteR/projects/hiti/pipeline_output/NT_aligned"
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p5_seq_SP1_local3.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, -3,-1)
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p5_seq_SP1_local2.fasta"
test_res=aligner(df_trim_full2, target_sequence, "align_local3", result, output_path, 3,1)
#AA
corr_frame=1
result="/media/data/AtteR/projects/hiti/pipeline_output/fastas/mcherry_p5_seq_PS1.fasta"
df_aa=translate_nt_aa(result, 1)
output_html="/media/data/AtteR/projects/hiti/pipeline_output/AA_aligned/AA_5p_mcherry_coloured_alignment_SP1.html"
visualise_aa_hybrid_alignments(df_aa, output_html)
#########


