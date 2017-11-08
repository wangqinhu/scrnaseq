#!/bin/bash

#$ -cwd
#$ -N sc_rnaseq.013
#$ -j y

cd data/refseq
hisat2-build mm10_ercc92.genome.fa mm10_ercc92.genome
