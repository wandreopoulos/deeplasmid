#!/usr/bin/env python
""" Formatting of training data: read raw input data
both sequences and features (yaml files)
write out to a yaml file the normalized, split into segments dataset.
"""
__author__ = "Bill Andreopoulos, Jan Balewski"
__email__ = "wandreopoulos@lbl.gov, janstar1122@gmail.com"

#exit(1)  # block 3a
from Plotter_Plasmid import Plotter_Plasmid
from DL_Model import DL_Model
from Util_Plasmid import get_glob_info_files

import argparse

import Constants


def get_parser():
    parser = argparse.ArgumentParser(
        description='Format raw scafolding for Plasmid Oracle training',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v","--verbosity",type=int,choices=[0, 1, 2],
                        help="increase output verbosity", default=1, dest='verb')
    parser.add_argument("--dataPath",
                        default='data',help="path to ML input/output")
    parser.add_argument("--outPath",
                        default='outPR',help="output path for plots/output")
    parser.add_argument( "-X","--no-Xterm", dest='noXterm',
                         action='store_true', default=False,
                         help="disable X-term for batch mode")

    parser.add_argument("-i","--inputfasta", dest='inputfasta',
                         help="The first class input fasta file", required=False)
    parser.add_argument("-y","--inputyml", dest='inputyml',
                         help="The first class input yml directory", required=False)
    parser.add_argument("-i2","--inputfasta2", dest='inputfasta2',
                         help="The second class input fasta file", required=False)
    parser.add_argument("-y2","--inputyml2", dest='inputyml2',
                         help="The second class input yml directory", required=False)

    # These are unused in this script for now
    parser.add_argument("-a", "--arrIdx", type=int, default=1,
                        help="slurm array index")
    parser.add_argument("-k", "--kfoldOffset", type=int, default=0,
                        help="decides which segments merge for training")

    args = parser.parse_args()
    for arg in vars(args):  print( 'myArg:',arg, getattr(args, arg))
    return args


# plasmids/main
#commPath='/global/projectb/sandbox/rqc/andreopo/Jan_genome_plasmid/v1_apr62018/'
#plasmFN=commPath+'plasmid/aclame_plasmid_sequences.fasta'
#plasmGFDir=commPath+'plasmid/yml/'
#mainFN=commPath+'genomic/refseq.bacteria.nonplasmid.nonmito.fasta.subsam40kreads.fasta'
#mainGFDir=commPath+'genomic/yml/'

#The plasmids
commPath='/global/projectb/sandbox/rqc/andreopo/src/bitbucket/jgi-ml/classifier/DATA/ACLAME.REFSEQMICROB/'
plasmFN=commPath+'aclame_plasmid_sequences.fasta.MIN1kMAX330k.fasta'
plasmGFDir=commPath+'aclame_plasmid_sequences.fasta.MIN1kMAX330k.fasta.yml/2019-02-27_08:55:26/yml/'

#The main (chromosome) files
mainFN=commPath+'refseq.bacteria.nonplasmid.nonmito.fasta.subsam40kreads.fasta.MIN1kMAX330k.fasta'
mainGFDir=commPath+'refseq.bacteria.nonplasmid.nonmito.fasta.subsam40kreads.fasta.MIN1kMAX330k.fasta.yml/2019-02-27_09:45:48/yml/'

#=================================
#=================================
#  M A I N 
#=================================
#=================================
args=get_parser()
#Get the plasmid class input
if args.inputfasta:
    plasmFN=args.inputfasta
if args.inputyml:
    plasmGFDir=args.inputyml
#Get the main class input
if args.inputfasta2:
    mainFN=args.inputfasta2
if args.inputyml2:
    mainGFDir=args.inputyml2


ppp=Plotter_Plasmid(args )
deep=DL_Model(args)

# Get the files (yaml) with the stats features
plasmGFD=get_glob_info_files(plasmGFDir)
mainGFD=get_glob_info_files(mainGFDir)

#######################
#Format the plasmid class
#max number of plasmids to use in the training
mxScaf=Constants.mxScaf_train_class1
#read the fasta and the features
recP=deep.read_scaffolds_fasta(plasmFN,plasmGFD,mxScaf)
#split the seqs into folds
recPP=deep.split_into_folds(recP,'plasm',info=[plasmFN,mxScaf])

if not args.noXterm:
    ppp.plot_float_featuresA(deep,recPP,'plasm',16)
#ppp.pause(args,'format')

#######################
#Format the chromosome class
#max number of chromosomal sequences to use in training
mxScaf=Constants.mxScaf_train_class2
#read the fasta and the features
recM=deep.read_scaffolds_fasta(mainFN,mainGFD,mxScaf)
#split the seqs into folds
recMM=deep.split_into_folds(recM,'main',info=[mainFN,mxScaf])
#print('skip plots, ok'); exit(1)  

#No need for plotting the stats
if not args.noXterm:
    ppp.plot_scaff_len(deep,recPP,'plasm',6)
    ppp.plot_scaff_len(deep,recMM,'main',7)
    ppp.plot_float_featuresA(deep,recMM,'main',17)
    ppp.pause(args,'format')

