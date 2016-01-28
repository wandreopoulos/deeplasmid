#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -notify
#$ -P gentech-rqc.p
#$ -N rqc_metag_bin
#$ -l h_rt=12:00:00
#$ -l ram.c=100G
#$ -M wandreopoulos@lbl.gov
#$ -m abe
./run_bin.sh
