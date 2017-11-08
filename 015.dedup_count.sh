#!/bin/bash

#$ -t 1-1920
#$ -cwd
#$ -N sc_rnaseq.015
#$ -j y

data="data/bam/*/*.sorted.bam "
bams=(`ls $data`)
bam=${bams[$(expr $SGE_TASK_ID - 1)]}

cell=`basename $bam`
tiss=`dirname $bam`
tiss=`basename $tiss`
well=${cell%.sorted.bam}

mkdir -p data/count/$tiss

echo $tiss/$well

cd data/bam/$tiss

samtools index $well.sorted.bam 
umi_tools dedup -I $well.sorted.bam -S $well.dedup.bam
samtools index $well.dedup.bam
featureCounts -a ../../refseq/mm10_ercc92.gtf -o $well.gene_assigned -R BAM $well.dedup.bam
samtools sort $well.dedup.bam.featureCounts.bam $well.gene_assigned
samtools index $well.gene_assigned.bam
umi_tools count --per-gene --gene-tag=XT --per-cell -I $well.gene_assigned.bam -S ../../count/$tiss/$well.count.tsv

echo "Done"
