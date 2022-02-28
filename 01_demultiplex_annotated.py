import os
import warnings
from subprocess import call
import subprocess
import glob
import pandas as pd
import tempfile
import ntpath
import numpy as np
import matplotlib.pyplot as plt
import re
from datetime import date

run_date = date.today()



in_path = "'/media/data/AtteR/210517_A00681_0370_AHF2JHDRXY/demultiplexed'"
out_path = "'/media/data/AtteR/scifi-analysis/pilot1_output'"
base_path = '/media/data/AtteR/210517_A00681_0370_AHF2JHDRXY/demultiplexed'
export_path = '/media/data/AtteR/scifi-analysis/export_data'
output_data2 = '/media/data/AtteR/scifi-analysis/' 
output_tables = '/media/data/AtteR/scifi-analysis/pilot1_output/tables'

#the symlinks are generated to ensure path without spaces which bbmap cannot handle
""" os.system("mkdir /media/data/AtteR/scifi-analysis/pilot1_output")
!rm $base_path
!ln -s $in_path $base_path
!rm $export_path
!ln -s $out_path """ #$export_path
base_path = base_path +'/'
export_path = export_path + '/'


### Fuse Scifi and Bead barcode reads
def fuseBCs(sample_id) :
    print(sample_id)
    '''
    Takes in the sample_id i.e. the read name. In our case read 1 (UMI+sample bc) and fuses it with read 2
    '''
    revCompR1 = in_path + '/' + os.path.basename(sample_id).replace('_R1_001.fastq.gz','_RC_R1_001.fastq.gz')
    rc_call= "reformat.sh overwrite=t rcomp=t in='"+ sample_id +"' out='"+revCompR1+"'"
    #print(rc_call)
  #  !$rc_call
    
    inPathR2 = base_path + ntpath.basename(sample_id).replace('_R1_001.fastq.gz','_R2_001.fastq.gz')
    outPath = base_path + ntpath.basename(sample_id).replace('_R1_001.fastq.gz','_R21_001.fastq.gz')
    fuse_call= "/usr/local/bin/bbmap/fuse.sh in1='"+ inPathR2 +"' in2='"+revCompR1+"' out='"+outPath+"' overwrite=t interleaved=f fusepairs=t pad=0"
   # !$fuse_call 
    #os.remove(revCompR1)
    print(fuse_call)

#get the sample 
samples_list = sorted(glob.glob(base_path+'*R1*.fastq.gz'))
#samples_list= [item for item in samples_list if 'Undetermined' not in item]
#samples_list= [item for item in samples_list if 'BC' in item]
print(samples_list)


#Then fuse the read1 with the read 2 (bead bc) OR with read 3 (cdna).
#When you fuse read1 (contains the umi and BC dictating which well it came from) with read 2 (bead BC) then you receive
#the unique identity of the transcriptome 
[fuseBCs(sample_id) for sample_id in samples_list]


#Trim polyA
def runTrimPA(sample_id) :
    print(sample_id)
    outPath = export_path + ntpath.basename(sample_id).replace('R21_001.fastq.gz','R21_trimPA.fastq.gz')
    samplePair = base_path + ntpath.basename(sample_id).replace('R21_001.fastq.gz','R3_001.fastq.gz')
    outPathR3 = export_path + ntpath.basename(samplePair).replace('R3_001.fastq.gz','R3_trimPA.fastq.gz')
    cutAdapt_call="cutadapt --cores=0 --nextseq-trim=25 --minimum-length 25 -g GCTGCACCGTCAGA -O 12 -a 'A{30}N{104}' -o "+ outPathR3 + " -p " + outPath + " " + samplePair + " " + sample_id
    print(cutAdapt_call)
   # !$cutAdapt_call 
    
samples_list = sorted(glob.glob(base_path+'*R21_001.fastq.gz'))
samples_list
[runTrimPA(sample_id) for sample_id in samples_list[2:]]



#find the nucleus identity
#Scifi demultiplex
samples_list = sorted(glob.glob(export_path+'*R21_trimPA.fastq.gz'))
#samples_list= [item for item in samples_list if 'Undetermined' not in item]
#samples_list= [item for item in samples_list if 'BC' in item]
samples_list
def deMultiFuse(sample_id) :
    print(sample_id)
    #sample_id = samples_list[6]
    inPathR3 = export_path + ntpath.basename(sample_id).replace('_R21_trimPA.fastq.gz','_R3_trimPA.fastq.gz')
    outPathR1 = export_path + ntpath.basename(sample_id).replace('_R21_trimPA.fastq.gz','_trimmed-{name}_R21.fastq.gz')
    outPathR3 = export_path + ntpath.basename(sample_id).replace('_R21_trimPA.fastq.gz','_trimmed-{name}_R3.fastq.gz')
    trimLog = export_path + ntpath.basename(sample_id).replace('_R21_trimPA.fastq.gz','_trimLog.tsv')
    param = "-e 0.15 --cores=0 --action none -a file:\'" + export_path + "../ArrayIndex_211217.fasta\' --report=minimal"
    call_sequence = "cutadapt "+param+" -o \'"+outPathR1+"\' -p \'" + outPathR3 + "\' \'"+ sample_id + "\' \'" + inPathR3 +"\' > \'" + trimLog +"\'"
    print(call_sequence)
    #!($call_sequence)      
[deMultiFuse(sample_id) for sample_id in samples_list]



###Generate a QC table

qc_list = sorted(glob.glob(export_path+'*_trimLog.tsv'))
allSampleTable = pd.DataFrame()

for qc_file in qc_list :
    df=pd.read_csv(qc_file, sep='\t') # , header=0, index_col=0,skiprows=None
    df = df.set_index(pd.Series(ntpath.basename(qc_file).replace('_trimLog.tsv','')))
    allSampleTable = pd.concat([allSampleTable, df], axis=0)

print(allSampleTable)

with open(output_tables + "/QC" +str(run_date)+".csv", "w") as out: #write the dataframe into an output file
#write the dataframe into an output file
    allSampleTable.to_csv(out, sep='\t')
    # df.to_string(out, index=None)
    print('output info file saved!')


###Remove false cDNA indeces
samples_list = sorted(glob.glob(export_path+'*trimmed*.fastq.gz'))
for sample_file in samples_list :
    print(sample_file)
    result = subprocess.run(['pigz -dc '+ sample_file+ ' | wc -l'], stdout=subprocess.PIPE, shell=True)
    sampleReads = int(int(result.stdout)/4)
    if sampleReads < 4000 :
        os.remove(sample_file)
    print(sampleReads)
    
###List selected samples
samples_list = sorted(glob.glob(export_path+'*trimmed*R21.fastq.gz'))
for sample_file in samples_list :
    print(sample_file)
    result = subprocess.run(['pigz -dc '+ sample_file+ ' | wc -l'], stdout=subprocess.PIPE, shell=True)
    sampleReads = int(int(result.stdout)/4)
    print(sampleReads)
    
    
#Merge the individual files from the same cDNA index: 
def catSameID(cDNA_id) :
    print(cDNA_id)
    catSeq = "cat "+export_path+"*oDT_"+cDNA_id+"_R21.fastq.gz > /tmp/outdata/SciFi_oDT_" +cDNA_id+"_R21.fastq.gz"
    catSeqRev = "cat "+export_path+"*oDT_"+cDNA_id+"_R3.fastq.gz > /tmp/outdata/SciFi_oDT_" +cDNA_id+"_R3.fastq.gz"
    print(catSeq)
    #!($catSeq)
    #!($catSeqRev) 

print("Running weird cats:")
#!cat /media/data/AtteR/scifi-analysis/export_data/*R21_trimPA.fastq.gz > /media/data/AtteR/scifi-analysis/export_data/SciFi_oDT_Complete_R21_trimPA.fastq.gz
#!cat /media/data/AtteR/scifi-analysis/export_data/*R3_trimPA.fastq.gz > /media/data/AtteR/scifi-analysis/export_data/SciFi_oDT_Complete_R3_trimPA.fastq.gz

