#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -notify
#$ -P gentech-rqc.p
#$ -N rqc_NOSETESTS_assemblyqc
#$ -l h_rt=48:00:00
#$ -l exclusive.c
#$ -M wandreopoulos@lbl.gov
#$ -m abe
./run.sh

