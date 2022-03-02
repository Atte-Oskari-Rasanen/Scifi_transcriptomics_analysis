WHITELIST=/home/data/737K-cratac-v1.txt
WHITELIST=/media/data/AtteR/scifi-analysis/737K-cratac-v1_RC.txt


solo_exec () {
DIR=$DATAOUT/${EXP}_${ID}             
echo "output starsolo dir: $DIR"
/media/data/AtteR/Attes_bin/STAR-2.7.10a/bin/Linux_x86_64/STAR --genomeDir $MAP --runThreadN 24 \
    --soloType CB_UMI_Simple \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloCBwhitelist $WHITELIST --soloBarcodeMate 0 --soloBarcodeReadLength 0 \
    --soloFeatures GeneFull --readFilesCommand zcat --soloMultiMappers Uniform EM PropUnique Unique \
    --limitBAMsortRAM 220000000000 --soloUMIdedup 1MM_All --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --genomeLoad LoadAndKeep --soloStrand Unstranded \
    --clip3pAdapterSeq polyA --alignIntronMin 20 --alignIntronMax 1000000 \
    --readFilesIn $RV, $FW 
mkdir $DIR
mv Solo.out $DIR
mv Log.final.out Log.out Log.progress.out SJ.out.tab Aligned.sortedByCoord.out.bam Aligned.out.sam $DIR/
chmod -R 755 $DIR
rm -rf _STARtmp
}


#DATAPTH=/media/data/AtteR/scifi-analysis
DATAPTH=/media/data/AtteR/scifi-analysis/scifi6/2nd_try/output
DATAOUT=/media/data/AtteR/scifi-analysis/scifi6/2nd_try/starsolo_outputs

#DATAPTH=/media/data/AtteR/scifi-analysis/scifi6/output
#DATAOUT=/media/data/AtteR/scifi-analysis/scifi6/scifi_output

mkdir $DATAOUT
cd $DATAOUT


MAP=/media/data/AtteR/scifi-analysis/ref_genome/rat_index_aav

# #### A9, G8, H8 F8_Asyn, 5, 6

#POOL1
EXP="SciFi_new_cDNA_Pool1__trimmed-Pool1"
ID="_oDT_A9"
#we take the files with the above expressions, then run the solo_exec function

#A9
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz

echo "FW strand: $FW"
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

ID="_WP_A9"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
echo $FW
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

echo "Done A9, pool1"
#F8_Asyn
ID="_oDT_F8_Asyn"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

ID="_WP_F8_Asyn"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

echo "Done F8_Asyn, pool1"

#G8
ID="_oDT_G8"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
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


echo "Pool1 done!"

#Pool2
EXP="SciFi_old_cDNA_Pool2__trimmed-Pool2"
ID="_oDT_A9"
#we take the files with the above expressions, then run the solo_exec function

#A9
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

ID="_WP_A9"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

echo "Done A9, pool1"
#F8_Asyn
ID="_oDT_F8_Asyn"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 

ID="_WP_F8_Asyn"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
solo_exec 


#G8
ID="_oDT_G8"
FW=$DATAPTH/${EXP}${ID}_R21.fastq.gz
RV=$DATAPTH/${EXP}${ID}_R3.fastq.gz
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
