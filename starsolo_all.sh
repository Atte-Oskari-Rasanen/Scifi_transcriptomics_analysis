WHITELIST=/home/data/737K-cratac-v1.txt
WHITELIST=/media/data/AtteR/scifi-analysis/737K-cratac-v1_RC.txt


solo_exec () {
filename=${EXP}${ID} 
folder_name=$(echo $filename | cut -f1-6 -d'_') #get the range from first to 6th occurence based on delim _

echo "Folder name: $folder_name"

DIR=$DATAOUT/$folder_name             
echo "output starsolo dir: $DIR"
/media/data/AtteR/Attes_bin/STAR-2.7.10a/bin/Linux_x86_64/STAR --genomeDir $MAP --outTmpDir solo_tmp --runThreadN 24 \
    --soloType CB_UMI_Simple \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloCBwhitelist /media/data/AtteR/scifi-analysis/737K-cratac-v1_RC.txt --soloBarcodeMate 0 --soloBarcodeReadLength 0 \
    --soloFeatures GeneFull --readFilesCommand zcat --soloMultiMappers Uniform EM PropUnique Unique \
    --limitBAMsortRAM 220000000000 --soloUMIdedup 1MM_All --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --genomeLoad LoadAndRemove --soloStrand Unstranded \
    --clip3pAdapterSeq polyA --alignIntronMin 20 --alignIntronMax 1000000 \
    --readFilesIn $RV, $FW  
    #--soloOutFileNames $DIR/ features.tsv barcodes.tsv matrix.mtx \
    

mkdir $DIR
chmod -R 755 $DIR

gzip Solo.out/GeneFull/filtered/*

mv Solo.out $DIR


mv Log.final.out Log.out Log.progress.out SJ.out.tab Aligned.sortedByCoord.out.bam Aligned.out.sam $DIR/
rm -rf solo_tmp
#if test -f "Solo.out"; then
#    echo "Solo.out exists!"
#fi
--soloCellFilter TopCells 1 0.99 10 \

}

# cell_filter ()  {
# STAR --runMode soloCellFiltering path_raw/ path_feat/ --soloCellFilter TopCells 1 0.99 10

# }


#DATAPTH=/media/data/AtteR/scifi-analysis
DATAPTH=/media/data/AtteR/scifi-analysis/scifi6/2nd_try/output_enriched
DATAOUT=/media/data/AtteR/scifi-analysis/scifi6/2nd_try/starsolo_outputs_enriched



#DATAPTH=/media/data/AtteR/scifi-analysis/scifi6/output
#DATAOUT=/media/data/AtteR/scifi-analysis/scifi6/scifi_output

mkdir $DATAOUT
cd $DATAOUT


MAP=/media/data/AtteR/scifi-analysis/ref_genome/rat_index_aav
echo "GENOME DIR: $MAP"

#POOL1
EXP="SciFi_new_cDNA_Pool1_trimmed-"
ID="oDT_A9_50"
#we take the files with the above expressions, then run the solo_exec function
#A9
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz

RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
echo "RV strand: $RV"

solo_exec 

ID="WP_A9_16"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
echo "RV strand: $RV"

solo_exec 


#F8
ID="oDT_F8_255"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 
echo "RV strand: $RV"

ID="WP_F8_237"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 
echo "RV strand: $RV"


echo "Done F8_Asyn, pool1"

#G8
ID="oDT_G8_287"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
echo "RV strand: $RV"

solo_exec 
ID="WP_G8_295"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 


#H8
ID="oDT_H8_351"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

ID="WP_H8_350"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

#old 5
EXP="SciFi_new_cDNA_Pool1_trimmed-old"
ID="oDT_5"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

ID="WP_5"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

#old 6 
EXP="SciFi_new_cDNA_Pool1_trimmed-old"
ID="oDT_6"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

ID="WP_6"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

echo "Pool1 done!"

########################################
########################################

#Pool2
EXP="SciFi_old_cDNA_Pool2_trimmed"
ID="-oDT_A9_50"
#we take the files with the above expressions, then run the solo_exec function

#A9
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz

RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
echo "RV strand: $RV"

solo_exec 

ID="-WP_A9_16"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
echo "RV strand: $RV"

solo_exec 


#F8
ID="-oDT_F8_255"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 
echo "RV strand: $RV"

ID="-WP_F8_237"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 
echo "RV strand: $RV"


echo "Done F8_Asyn, pool1"

#G8
ID="-oDT_G8_287"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
echo "RV strand: $RV"

solo_exec 
ID="-WP_G8_295"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 


#H8
ID="-oDT_H8_351"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

ID="-WP_H8_350"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

#old 5
ID="-oldoDT_5"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

ID="-oldWP_5"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

#old 6 
ID="-oldoDT_6"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

ID="-oldWP_6"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 
echo "Pool2 done!"

########################################
#Pool 1-2 AAV
#F8
EXP="SciFi_Pool1-2_EnrichedBCs_All_AAV-aSyn_trimmed-"
ID="WP_F8_237"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec
########################################

########################################
#Pool 1-2 tagBFP

#WP5
EXP="SciFi_Pool1-2_EnrichedBCs_All_AAV-tagBFP_trimmed"
ID="-oldWP_5"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec

#WP6
ID="-oldWP_6"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec

#WP_A9
ID="-WP_A9_16"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec

ID="-WP_G8_295"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec

ID="-WP_H8_350"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec


echo "BCs done!"
########################################




STAR --genomeDir $MAP --genomeLoad Remove
