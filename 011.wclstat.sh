#!/bin/bash

#$ -t 1-5
#$ -cwd
#$ -N sc_rnaseq.011-stat
#$ -j y

data="data/plate2well"

subdirs=(`ls $data`)
subdir=${subdirs[$(expr $SGE_TASK_ID - 1)]}

echo "Counting $subdir ..."
wc -l $data/$subdir/* && echo Done

