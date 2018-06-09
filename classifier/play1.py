#!/usr/bin/env python
""" read input 
train net
write net + weights as HD5

"""
__author__ = "Jan Balewski"
__email__ = "janstar1122@gmail.com"

from Plotter_Plasmid import Plotter_Plasmid
from Util_Plasmid import write_yaml, read_yaml

import argparse
def get_parser():
    parser = argparse.ArgumentParser(
        description='Perform Plasmid Oracle training',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--arrIdx", type=int, default=1,
                        help="slurm array index")
    parser.add_argument("--project",
                        default='assayer4',dest='prjName',
                        help="core name used to store outputs")

    parser.add_argument( "-X","--noXterm", dest='noXterm',
                         action='store_true', default=False,
                         help="disable X-term for batch mode")

    parser.add_argument("--outPath",
                        default='out/',help="output path for plots and tables")

    args = parser.parse_args()
    for arg in vars(args):  print( 'myArg:',arg, getattr(args, arg))
    return args

#=================================
#=================================
#  M A I N 
#=================================
#=================================
args=get_parser()

inpF='/global/cscratch1/sd/balewski/plasmid_sum/4g/plasmid4g-12/assayer4.history.yml'
inpD=read_yaml(inpF)
#print(inpD)
args.train_hirD=inpD
args.name='abc'
args.train_sec=123
ppp=Plotter_Plasmid(args )
ppp.plot_train_history(args,args,10)

ppp.pause(args,'play') 




