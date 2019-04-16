#!/usr/bin/env python
""" Formatting of testing data: read raw input data
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
                         help="The input fasta file", required=False)
    parser.add_argument("-y","--inputyml", dest='inputyml',
                         help="The input yml directory", required=False)

    # These are unused in prediction for now
    parser.add_argument("-a", "--arrIdx", type=int, default=1,
                        help="slurm array index")
    parser.add_argument("-k", "--kfoldOffset", type=int, default=0,
                        help="decides which segments merge for training")

    args = parser.parse_args()
    for arg in vars(args):  print( 'myArg:',arg, getattr(args, arg))
    return args


# plasmids/main
commPath='/global/projectb/scratch/andreopo/AsafPlasmids/'  ###'/global/projectb/scratch/andreopo/AsafPlasmids/IMG_plasmids/'  ###651053060/'  ###20180425.651053060.fna/global/projectb/sandbox/rqc/andreopo/Jan_genome_plasmid/v1_apr62018/'
plasmFN=commPath+'genome_list.fasta' ###genome_listb.fasta' ###'651053060.genes.fna'
plasmGFDir=commPath+'yml/'  ###'yml2/' ###'yml.Asaf_Plasmid.fasta/' ###CP002622.yml' ###features.txt'  ###plasmid/yml/'

#=================================
#=================================
#  M A I N 
#=================================
#=================================
args=get_parser()
if args.inputfasta:
    plasmFN=args.inputfasta
if args.inputyml:
    plasmGFDir=args.inputyml

#not needed unless we plan to plot the stats
ppp=Plotter_Plasmid(args )

#just used for reading the information for formatting
deep=DL_Model(args)

#get all files with statistics
plasmGFD=get_glob_info_files(plasmGFDir)

#read the files : fasta sequences and statistics files
recP=deep.read_scaffolds_fasta(plasmFN,plasmGFD,Constants.mxScaf_pred)

#split the input into folds
recPP=deep.split_into_folds(recP,'test',info=[plasmFN,Constants.mxScaf_pred])

if not args.noXterm:
    ###Plot the feaatures and the scaffold length
    ppp.plot_float_featuresA(deep,recPP,'test',16)
    #ppp.pause(args,'format')

    #print('skip plots, ok'); exit(1)  
    ppp.plot_scaff_len(deep,recPP,'test',6)
    ppp.pause(args,'format')

