#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -notify
#$ -P gentech-rqc.p
#$ -N rqc_metag_bam
#$ -l h_rt=12:00:00
#$ -l ram.c=100G
#$ -t 1-2
#$ -M wandreopoulos@lbl.gov
#$ -m abe
./run.sh
