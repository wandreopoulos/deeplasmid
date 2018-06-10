#!/usr/bin/env python
"""  read raw input data
sanitize, split, ballance correction
write data as tensors as 10 segments
"""
__author__ = "Jan Balewski"
__email__ = "janstar1122@gmail.com"

#exit(1)  # block 3a
from Plotter_Plasmid import Plotter_Plasmid
from Deep_Plasmid import Deep_Plasmid
from Util_Plasmid import get_glob_info_files

import argparse
def get_parser():
    parser = argparse.ArgumentParser(
        description='Format raw scafolding for Plasmid Oracle training',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v","--verbosity",type=int,choices=[0, 1, 2],
                        help="increase output verbosity", default=1, dest='verb')
    parser.add_argument("--project",
                        default='assayer4',dest='prjName',
                        help="core name used everywhere")
    parser.add_argument("--dataPath",
                        default='data',help="path to ML input/output")
    parser.add_argument("--outPath",
                        default='out',help="output path for plots")
    parser.add_argument( "-X","--no-Xterm", dest='noXterm',
                         action='store_true', default=False,
                         help="disable X-term for batch mode")
    parser.add_argument("-i","--inputfasta", dest='inputfasta',
                         help="The input fasta file", required=False)
    parser.add_argument("-y","--inputyml", dest='inputyml',
                         help="The input yml directory", required=False)
    args = parser.parse_args()
    for arg in vars(args):  print( 'myArg:',arg, getattr(args, arg))
    args.arrIdx=0 # for plotter, not needed here
    args.kfoldOffset=0 # for training, not needed here
    args.events=0 # for training, not needed here
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

ppp=Plotter_Plasmid(args )
deep=Deep_Plasmid(args)

plasmGFD=get_glob_info_files(plasmGFDir)

mxScaf=122000 # set it below 1k for testing

recP=deep.read_scaffolds_fasta(plasmFN,plasmGFD,mxScaf)

recPP=deep.split_species(recP,'test',info=[plasmFN,mxScaf])

#ppp.plot_float_featuresA(deep,recPP,'test',16)
#ppp.pause(args,'format')

#print('skip plots, ok'); exit(1)  
#ppp.plot_scaff_len(deep,recPP,'test',6)
#ppp.pause(args,'format')

# - -- - - ver 4f was the last which used it
# plasmids/main
plasmFN='/global/projectb/scratch/andreopo/GAA-1290_plasmids/Dingl/plasmids/aclame_plasmid_sequences.fasta'
plasmGFDir='/global/projectb/scratch/andreopo/GAA-1290_plasmids/Dingl/plasmids/20180327.aclame_plasmid_sequences.fasta/yml/'

mainFN='/global/projectb/scratch/andreopo/GAA-1290_plasmids/Dingl/microbial/NEW3_ALL_FEATURES/refseq.bacteria.NOplasmids_NOmito.fasta.1'
mainGFDir='/global/projectb/scratch/andreopo/GAA-1290_plasmids/Dingl/microbial/NEW3_ALL_FEATURES/20180327.refseq.bacteria.NOplasmids_NOmito.fasta.1/yml/'
