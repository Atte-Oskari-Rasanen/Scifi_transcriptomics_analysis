from typing import Iterator
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


data_dict=create_datadict(base_path)
full_df_trim=import_reads_process(data_dict, transgene,assay_end,filterlitteral,lliteral,rliteral,export_path,read_fwd):
full_df_trim.to_csv("/media/data/AtteR/projects/hiti/dfs/full_df_trim_mcherry_p3_exp2_trimPrimer.csv")

result="/media/data/AtteR/projects/hiti/mcherry_p3_seq_cluster_all_retrim.fasta"
save_fasta(result, full_df_trim, target_sequence)


###############
#Alignments
alignment_nt="/media/data/AtteR/projects/hiti/align_output/mcherry_p3_seq_aligned_pairwise_local_3_go3_ge1.fasta"

test_res=aligner(full_df_trim, target_sequence, "align_local2", alignment_nt)




#downstream an issue with visualising the seqs using mview is that all the seqs are not SAME length. thus, needs to fix this


#takes in the trimmed df and aligns each amplicon against the ref based on the aligmment method
#of choice taken from align class. 

def align_and_save(filename, full_df, align_method):
    id_f=1
    aligned_data=dict()
    #align all the data, save into dict, then ensure that all the seqs are same length (take the longest seq). IF not, then add padding!
    for seq_i in range(len(full_df.iloc[:,-1])):
        header=">"+ str(id_f)+"CluSeq:" + str((round(full_df.iloc[seq_i,-3],5))) + "_sd:" + str((round(full_df.iloc[seq_i,-1],5)))
        align_inst=align(full_df.iloc[seq_i,0], target_sequence)
        seq_obj_1=align_inst.align_local3()
        #seq_obj_1= align_local(full_df.iloc[seq_i,0], target_sequence)
        seq_obj_1 = re.sub(r'[(\d|\s]', '', seq_obj_1) #remove digits from the string caused by the alignment and empty spaces from the start
        aligned_data[header]=seq_obj_1
        id_f+=1
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

#iterate over files
#align
#trim


#in alignment method: apply the iterator and trimmer
align_and_save(filename, full_df_trim,target_sequence, align_local3)
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



from Bio.Data import CodonTable
from Bio.Seq import Seq


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


