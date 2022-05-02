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
        full_df = pd.merge(complete_df, df_all_lanes, on="sequence", how='outer')
        print("A full df containing the sum from all lanes of " + animal_group_name + " is done!")
    
    # after generation of the full df containing all animals, the following function calculates the stats (sd, mean,
    # perc. sum per seq cluster and raw counts sum per seq cluster)
        #animal group contains all lanes of the certain data
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
    ref="mcherry_p3_seq_ref"
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

class aa_alignments():

    def __init__(self, NT_seq_fasta, corr_frame, alt_frame,out_html):
        self.NT_seq_fasta=NT_seq_fasta
        self.out_html=out_html
        self.corr_frame=corr_frame
        self.alt_frame=alt_frame

    def find_colour(aa,color_scheme):
        col=[key for key,val in color_scheme.items() if any(aa in s for s in val)]
        out="".join(col)
        if out==None or not col:
            out="BLACK"
            return(out)
        else:
            return(out)

    def align_ampl_x_ref(self):
        refs_aa_frames={}
        aa_and_perc={}
        out_of_frames=[0,1,2]
        out_of_frames.remove(self.corr_frame)

        #Translate the seq into AAs
        for record in SeqIO.parse(self.NT_seq_fasta, "fasta"):
            if record.id=="0":
                refs_aa_frames["Frame:" + str(self.corr_frame)]=str(Seq(record.seq[self.corr_frame:]).translate())
                for alt_frame in out_of_frames:
                    refs_aa_frames["Frame:" + str(self.alt_frame)]=str(Seq(record.seq[self.alt_frame:]).translate())
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

        seq_info=["Seq_info"]+list(aa_and_perc.keys())
        ref_x_alig_list=[]
        for keys, values in ref_x_alignment.items():
            print("".join(list((keys))[:]))
            #ref_x_alig_list.append([("".join(list((keys))))]+list(values))
            ref_x_alig_list.append([keys]+list(values))

        #df = pd.DataFrame([seq_info, ref_x_alig_list[0],ref_x_alig_list[1],ref_x_alig_list[2]], columns =[seq_info[0],ref_x_alig_list[0][0], ref_x_alig_list[1][0], ref_x_alig_list[2][0]])
        df = pd.DataFrame(data= {seq_info[0]: seq_info[1:], ref_x_alig_list[0][0]:ref_x_alig_list[0][1:], ref_x_alig_list[1][0]:ref_x_alig_list[1][1:], ref_x_alig_list[2][0]:ref_x_alig_list[2][1:]})

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
            if match2_a.size>match2_b.size:
                matched_ampl2= len(df.columns[2].split(":")[-1][:match2_a.a])*"-" + str(df.iloc[ampl_row,2][match2_a.b:]) + "|"
                ref_outframe2=df.columns[2].split(":")[-1][0:end_of_match] + "|"
                # seq_inframe2=matched_ampl2[0:end_of_match]
                seq_outframe=matched_ampl2[end_of_match:]
            else:
                matched_ampl2= len(df.columns[3].split(":")[-1][:match2_b.a])*"-" + str(df.iloc[ampl_row,3][match2_b.b:]) + "|"
                ref_outframe2=df.columns[3].split(":")[-1][end_of_match:] + "|"
                seq_outframe=matched_ampl2[end_of_match:]
                # seq_inframe2=matched_ampl2[0:end_of_match]
                seq_outframe=matched_ampl2[end_of_match:]
                frame_info.append("Frame:+1 " + ref_inframe1.index("|")*" " + "Frame:+2")
            merged_align_seq=seq_inframe1+seq_outframe

            merged_align_ref=ref_inframe1+ref_outframe2

            seqinfos.append(df.iloc[ampl_row,0])
            #aligs_merged_seq.append(merged_align_seq[0:round(len(aligs_merged_seq)*0.6)])
            aligs_merged_seq.append(merged_align_seq)
            #aligs_merged_ref.append(merged_align_ref[0:round(len(merged_align_ref)*0.6)])
            aligs_merged_ref.append(merged_align_ref)


            #using lesk color code scheme:
            #color_scheme={"RED": [["#FF0000"], ["D", "E"]], "BLUE":[["#6495ED"],["K", "R"]], "GREEN":[["#9FE2BF"], ["C", "V", "I", "L", "P", "F", "Y", "M", "W"]], "ORANGE":[["#FF7F50"],["G", "A", "S", "T"]], "MAGENTA":[["#DE3163"], ["N", "Q", "H"]], "BLACK": [["#000000"],["-","|", "*"]]}
            color_scheme={"RED": [["#FF0000"], ["D", "E"]], "BLUE":[["#6495ED"],["K", "R"]], "GREEN":[["#9FE2BF"], ["C", "V", "I", "L", "P", "F", "Y", "M", "W"]], "ORANGE":[["#FF7F50"],["G", "A", "S", "T"]], "MAGENTA":[["#DE3163"], ["N", "Q", "H"]], "BLACK": [["#000000"]]}
            def find_colour(aa,color_scheme):
                col=[key for key,val in color_scheme.items() if any(aa in s for s in val)]
                out="".join(col)
                if out==None or not col:
                    out="BLACK"
                    return(out)
                else:
                    return(out)


            f = open(self.out_html,"w")
            for info, a1,a2, frame in zip(seqinfos, aligs_merged_seq, aligs_merged_ref,frame_info):
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
                    try:
                        aa2_list.append('<span style="color:'+ color_scheme[find_colour(aa2,color_scheme)][0][0]+ '">' + aa2 + '</span>')
                    except KeyError:
                        print("Keyerror with aa2 being: " + aa2)
                    try:
                        aa1_list.append('<span style="color:'+ color_scheme[find_colour(aa1,color_scheme)][0][0]+ '">' + aa1 + '</span>')
                    except KeyError:
                        print("Keyerror with aa1 being: " + aa1)
                    
                    # aa2_list.append('<span style="color:'+ color_scheme[find_colour(aa2,color_scheme)][0][0]+ '">' + aa2 + '</span>')
                    # aa1_list.append('<span style="color:'+ color_scheme[find_colour(aa1,color_scheme)][0][0]+ '">' + aa1 + '</span>')

                    #print("============")
                coloured_ref="".join(aa2_list)

                coloured_seq="".join(aa1_list)

                f.write(coloured_ref +"<br>")
                f.write(coloured_seq)
                f.write("</p>")

            f.write("</body></html>")
            f.close()
        print("alignment file saved into " + "/".join(self.out_html.split("/")[:-1]) + " as " + self.out_html.split("/")[-1])

