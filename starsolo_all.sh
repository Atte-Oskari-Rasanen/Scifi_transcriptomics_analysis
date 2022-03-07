WHITELIST=/home/data/737K-cratac-v1.txt
WHITELIST=/media/data/AtteR/scifi-analysis/737K-cratac-v1_RC.txt


solo_exec () {
DIR=$DATAOUT/${EXP}${ID}             
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
DATAPTH=/media/data/AtteR/scifi-analysis/scifi6/2nd_try/output
DATAOUT=/media/data/AtteR/scifi-analysis/scifi6/2nd_try/starsolo_outputs



#DATAPTH=/media/data/AtteR/scifi-analysis/scifi6/output
#DATAOUT=/media/data/AtteR/scifi-analysis/scifi6/scifi_output

mkdir $DATAOUT
cd $DATAOUT


MAP=/media/data/AtteR/scifi-analysis/ref_genome/rat_index_aav
echo "GENOME DIR: $MAP"
# #### A9, G8, H8 F8_Asyn, 5, 6

#POOL1
EXP="SciFi_new_cDNA_Pool1__trimmed-Pool1"
ID="_oDT_A9"
#we take the files with the above expressions, then run the solo_exec function

#A9
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz

RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
echo "RV strand: $RV"

solo_exec 
#give the starsolo the path where to find the raw files and then the output path
# path_raw = $DATAOUT$EXP$ID/Solo.out/GeneFull/raw
# path_feat = $DATAOUT$EXP$ID/Solo.out/GeneFull/filtered
# cell_filter
####################
ID="_WP_A9"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
echo "RV strand: $RV"

solo_exec 
# ####################

#F8_Asyn
ID="_oDT_F8_Asyn"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 
echo "RV strand: $RV"

ID="_WP_F8_Asyn"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 
echo "RV strand: $RV"

echo "Done F8_Asyn, pool1"

#G8
ID="_oDT_G8"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
echo "RV strand: $RV"

solo_exec 
ID="_WP_G8"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

#H8
ID="_oDT_H8"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

ID="_WP_H8"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

# #5
# ID="_oDT_5"
# FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
# solo_exec 

# ID="_WP_5"
# FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
# solo_exec 

# #6
# ID="_oDT_6"
# FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
# solo_exec 

# ID="_WP_6"
# FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
# solo_exec 


echo "Pool1 done!"

#Pool2
EXP="SciFi_old_cDNA_Pool2__trimmed-Pool2"
ID="_oDT_A9"
#we take the files with the above expressions, then run the solo_exec function

# #A9
# FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
# solo_exec 

# ID="_WP_A9"
# FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
# solo_exec 

# echo "Done A9, pool1"
# #F8_Asyn
# ID="_oDT_F8_Asyn"
# FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
# solo_exec 

# ID="_WP_F8_Asyn"
# FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
# solo_exec 


# #G8
# ID="_oDT_G8"
# FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
# solo_exec 
# ID="_WP_G8"
# FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
# solo_exec 

# #H8
# ID="_oDT_H8"
# FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
# solo_exec 

# ID="_WP_H8"
# FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
# RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
# solo_exec 

#5
ID="_oDT_5"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

ID="_WP_5"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

#6
ID="_oDT_6"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

ID="_WP_6"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 


echo "Pool2 done!"




STAR --genomeDir $MAP --genomeLoad Remove
