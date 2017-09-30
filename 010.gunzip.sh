#!/bin/bash

#$ -t 1-4
#$ -cwd
#$ -N sc_rnaseq.010
#$ -j y

rawDir="data/raw"

files=(`ls $rawDir`)
file=${files[$(expr $SGE_TASK_ID - 1)]}

echo "Unzip $file ..."
gunzip $rawDir/$file && echo Done

