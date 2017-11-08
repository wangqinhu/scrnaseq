#!/bin/bash

#$ -t 1-1920
#$ -cwd
#$ -N sc_rnaseq.014
#$ -j y

data="data/clean/*/*"
fastqs=(`ls $data`)
fastq=${fastqs[$(expr $SGE_TASK_ID - 1)]}

cell=`basename $fastq`
tiss=`dirname $fastq`
tiss=`basename $tiss`
prefix=${cell%.fq}
mkdir -p data/sam/$tiss
mkdir -p data/bam/$tiss

echo $tiss/$prefix
hisat2 -x data/refseq/mm10_ercc92.genome -U $fastq -S data/sam/$tiss/$prefix.sam
samtools view -bS data/sam/$tiss/$prefix.sam > data/bam/$tiss/$prefix.bam
samtools sort data/bam/$tiss/$prefix.bam data/bam/$tiss/$prefix.sorted
echo "Done"

