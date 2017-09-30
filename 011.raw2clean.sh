#!/bin/bash

#$ -cwd
#$ -N sc_rnaseq.011
#$ -j y
#$ -q great.q

perl 011.raw2clean.pl
