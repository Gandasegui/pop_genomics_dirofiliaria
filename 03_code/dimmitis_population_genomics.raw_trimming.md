#!/bin/bash

# run_trimgalore.sh

mkdir RAW

while read LANE NAME; do

mkdir ${NAME}_${LANE}_trimdata

trim_galore --paired --fastqc --length 50 --output_dir ${NAME}_${LANE}_trimdata ${LANE}_1.fastq.gz ${LANE}_2.fastq.gz

mv ${NAME}_${LANE}_trimdata/${LANE}_1_val_1.fq.gz ${NAME}_1.fastq.gz
mv ${NAME}_${LANE}_trimdata/${LANE}_2_val_2.fq.gz ${NAME}_2.fastq.gz

mv ${LANE}_1.fastq.gz RAW
mv ${LANE}_2.fastq.gz RAW

done < lanes_samples.list



#bsub.py --queue long 10 trim "bash ./run_trimgalore.sh"



#!/bin/bash

# run_trimmomatic.sh

module load trimmomatic/0.39--1

mkdir RAW

while read LANE NAME; do

mkdir ${NAME}_${LANE}_trimdata

trimmomatic PE \
-threads 10 -phred33 \
${LANE}_1.fastq.gz ${LANE}_2.fastq.gz \
${NAME}_1.fastq.gz ${NAME}_${LANE}_trimdata/${NAME}.unpaired_1.fastq.gz \
${NAME}_2.fastq.gz ${NAME}_${LANE}_trimdata/${NAME}.unpaired_2.fastq.gz \
ILLUMINACLIP:/nfs/users/nfs_s/sd21/lustre_link/databases/trimmomatic_Illumina-adapters.fa:2:30:10 \
SLIDINGWINDOW:10:20 MINLEN:50

mv ${LANE}_1.fastq.gz RAW
mv ${LANE}_2.fastq.gz RAW

done < lanes_samples.list


# bsub.py --queue long --threads 10 20 trim "bash ./run_trimmomatic.sh"