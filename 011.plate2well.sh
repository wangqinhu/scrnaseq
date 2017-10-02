#!/bin/bash

#$ -cwd
#$ -N sc_rnaseq.011
#$ -j y
#$ -q great.q

perl 011.plate2well.pl
