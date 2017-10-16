#!/bin/bash

#$ -t 1-1920
#$ -cwd
#$ -N sc_rnaseq.012
#$ -j y

data="data/plate2well/*/*"
fastqs=(`ls $data`)
fastq=${fastqs[$(expr $SGE_TASK_ID - 1)]}

cell=`basename $fastq`

tiss=`dirname $fastq`
tiss=`basename $tiss`
outdir="data/clean/$tiss"
mkdir -p $outdir

perl 012.filter_fastq.pl $fastq $outdir/$cell

