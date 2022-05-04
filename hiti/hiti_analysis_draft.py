import pandas as pd
import os 
os.getcwd()
os.chdir("/media/data/AtteR/projects/hiti")
from scripts_hiti import *

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

'''
SP1 3' 
lliteral=cggcggcatggacgagc
target_sequence="tgtacaagccgaacGTTCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAGGTCTTCTACCGTCTGGAGAG"
target_sequence=target_sequence.upper()

SP4 3'
target_sequence="tgtacaagATGGAGCTGGACCATATGACCCGGTGCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAGGTCTTCTACCGTCTGGAGAGG"
target_sequence=target_sequence.upper()

'''


def calculate_perc_sd2(full_df):
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
target_sequence = "TGTACAAGGTCGGTGCTGCGGCTCCGCGGAGCCGCAGCACCGACGACCAGATGGAGCTGGACCATATGACCACCGGCGGCCTCCACGCCTACCCTGCCCCGCGGGGTGGGCCGGCCGCCAAACCCAATGTGATCCTGCAGATTGGTAAGTGCCGAGCTGAGATGCTGGAACACGTACGGAGGACCCACCGGCATCTGTTGACCGAAGTGTCCAAGCAGGTGGAGCGAGAGCTGAAAGGGTTGCACAGGTCGGTGGGCAAGCTGGAGAACAACTTGGACGGCTACGTGCCCACCGGCGACTCACAGCGCTGGAAGAAGTCCATCAAGGCCTGTCTTTGCCGCTGCCAGGAGACCATCGCCAACCTGGAGCGCTGGGTCAAGCGTGAGATGCACGTGTGGAGGGAGGTCTTCTACCGTCTGGAGAGG"

#removed the first C from the target seq
"Filters and trims the reads"

trimmed_df="/media/data/AtteR/projects/hiti/dfs/full_df_trim_mcherry_p3_exp2_trimPrimer.csv"
data_dict=create_datadict(base_path)
full_df_test=import_reads_process(data_dict, transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd)
full_df_test

df_trim=calculate_perc_sd2(full_df_test)
df_trim.iloc[0,0]

df_trim.to_csv(trimmed_df)
# tf="/media/data/AtteR/projects/hiti/dfs/full_df_trim_mcherry_p3_exp2_retrimC.csv"
# full_df_trim2=pd.read_csv(trimmed_df)
full_df_trim=pd.read_csv(trimmed_df)
del full_df_trim['Unnamed: 0']
df_trim.columns
full_df_trim.columns
full_df_trim.iloc[2,0]
result="/media/data/AtteR/projects/hiti/mcherry_p3_seq_cluster_all_retrim.fasta"
save_fasta(result, full_df_trim, target_sequence)

#save_fasta(result, df_trim, target_sequence)

###############
from scripts_hiti import *
#Alignments NT
output_path="/media/data/AtteR/projects/hiti/pipeline_output/NT_aligned"

alignment_nt="/media/data/AtteR/projects/hiti/align_output_pipeline/mcherry_p3_seq_aligned_pairwise_local_3_go3_ge1.fasta"
test_res=aligner(full_df_trim, target_sequence, "align_local3", alignment_nt, output_path, -3, -1)

#when running the aligner with biopython's methods, need to add negative penalty scores but when running 
#smith watermans from sk package then need to be pos

#add mview function into the aligner function too.
'''
./mview -in fasta -html head -css on -coloring any /media/data/AtteR/projects/hiti/align_output_pipeline/mcherry_p3_seq_aligned_pairwise_local_3_go3_ge1.fasta > mview_data/aligned_loc_3_refcut.html

'''
alignment_nt="/media/data/AtteR/projects/hiti/align_output_pipeline/mcherry_p3_seq_aligned_pairwise_local_2_go3_ge1.fasta"
test_res=aligner(df_trim, target_sequence, "align_local2", alignment_nt)


corr_frame=1


def translate_nt_aa(result, corr_frame):
    out_of_frames=[0,1,2]
    out_of_frames.remove(corr_frame)
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

    #df = pd.DataFrame([seq_info, ref_x_alig_list[0],ref_x_alig_list[1],ref_x_alig_list[2]], columns =[seq_info[0],ref_x_alig_list[0][0], ref_x_alig_list[1][0], ref_x_alig_list[2][0]])
    df = pd.DataFrame(data= {seq_info[0]: seq_info[1:], ref_x_alig_list[0][0]:ref_x_alig_list[0][1:], ref_x_alig_list[1][0]:ref_x_alig_list[1][1:], ref_x_alig_list[2][0]:ref_x_alig_list[2][1:]})
    return(df)
#to show alignments to ref seq that has been translated in frame
df_aa=translate_nt_aa(result, 1)





#iterate over files
#align
#trim


#in alignment method: apply the iterator and trimmer
#we have deletions for sure, not sure about insertions

'''
RRHGRAVQGRCCGSAEPQHRRPDGAGPYDHRRPPRLPCPAGWAGRQTQCDPADW*VPS*DAGTRTEDPPASVDRSVQAGGARAERVAQVGGQAGEQLGRLRAHRRLTALEEVHQGLSLPLPGDHRQPGALGQA*DARVEGGLLPSGE


'''


from Bio.Data import CodonTable
from Bio.Seq import Seq


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

matched_ampl= len(ref[:match.a])*"-" + str(lines[i+1][match.b:]) + "|"

seqs=[]
range_line=0
for i in range(len(temp.get_matching_blocks())):
    match=temp.get_matching_blocks()[i]
    seqs.append(len(a[range_line:match.a])*"-"+b[match.b:match.b+match.size])
    range_line=match.a+match.size
    #based on the coordinates, extract the start and end, get the seqs between them
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



out="/media/data/AtteR/projects/hiti/AA_3p_mcherry_coloured_alignment2.html"

a=aa_alignments(result,out)
a.align_ampl_x_ref()
from scripts_hiti import *


