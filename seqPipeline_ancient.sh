#!/bin/bash -x
####################################################################################################
#The objective of this bash script is to call the pipeline.sub script n times, one for each sample.#
####################################################################################################

# usage > loop.sh fname
#Needed folders are created to store outputs.

#*********************************************************************************#
#*****************************Things to consider**********************************#
#*********************************************************************************#
#*****--export flag exports the variable $LINE in $SAMPLE in the .sub script.*****#


#$1 receives the file name passed when calling this script.
#This file has to be a list of sample names, one per line, without the file extension.

FILE=$1
## USAGE sh seqPipeline_ancient.sh /projects/125-emmer/scripts/ancient_exome_pipeline/sampleList.txt

#A counter is set.
counter=0

cd ~/ancient_exome_pipeline
OUTPATH=/scratch/125-emmer/exome_capture
INPATH=/projects/125-emmer/berlin-ancient/exomes # where the fastqs are
REFERENCE=/projects/125-emmer/durum_ref/Triticum_turgidum.Svevo.v1.dna_rm.toplevel.fa

#FOR PICARD
RGLB=$(echo $FILE | cut -f6 -d'/' | cut -f1 -d'.')
echo $RGLB

#For each line of that file, which is open,
for LINE in $(cat < $FILE)
do
        FNAME=$(echo $LINE | cut -f1-4 -d'_')
        indID=$(echo $FNAME | cut -f1 -d'_')_exome # This is the name that the g.vcf file will have
        echo $indID
        RGSM=$(echo $indID | cut -f1 -d'_' | cut -f1 -d'-')
        echo $RGSM
#        jid1=$(sbatch --parsable --export=FNAME=$FNAME,OUTPATH=$OUTPATH,INPATH=$INPATH adapRem.sub)
#        echo $jid1
#        jid2=$(sbatch --parsable --dependency=afterok:$jid1 --export=INDID=$indID,OUTPATH=$OUTPATH,REFERENCE=$REFERENCE,FQFILE=${indID}_paired.collapsed.gz bwaAlignSE.sub)
#        echo $jid2
#        jid3=$(sbatch --parsable --dependency=afterok:$jid2 --export=INDID=$indID,OUTPATH=$OUTPATH,REFERENCE=$REFERENCE,BAMFILE=${indID}_sorted.bam mapDamage.sub)        
#        jid4=$(sbatch --parsable --dependency=afterok:$jid2 --export=INDID=$indID,OUTPATH=$OUTPATH,BAMFILE=${indID}_sorted.bam,RGLB=$RGLB,RGSM=$RGSM processBam.sub)
         jid4=$(sbatch --parsable --export=INDID=$indID,OUTPATH=$OUTPATH,BAMFILE=${indID}_sorted.bam,RGLB=$RGLB,RGSM=$RGSM processBam.sub)
#        jid5=$(sbatch --dependency=afterok:$jid4 --export=INDID=$indID,OUTPATH=$OUTPATH,REFERENCE=$REFERENCE,BAMFILE=${indID}_sorted_dupMarked_RG.bam gatk_HaploCall_loop.sub)
        counter=$((counter+1))	
	sleep 2
done

echo $counter " jobs have been sent."
