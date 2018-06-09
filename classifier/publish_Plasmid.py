#!/usr/bin/env python
""" 
read yaml predictions for a species-singles and publish it on the web
"""
__author__ = "Jan Balewski"
__email__ = "janstar1122@gmail.com"


from PubAll_Plasmid import PubAll_Plasmid

import argparse
def get_parser():
    parser = argparse.ArgumentParser(
        description='publish mito/main classification on the web',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-v","--verbosity",type=int,choices=[0, 1, 2],
                        help="increase output verbosity", default=1, dest='verb') 
    parser.add_argument("--project",
                        default='assayer3',dest='prjName',
                        help="core name used to store outputs")

    parser.add_argument("--dataPath",
                        default='out/main/',help="path to input")
    parser.add_argument("--webPath",
                        default='/project/projectdirs/mpccc/www/balewski/tryAny/assayer3/',help="output path for plots")
    

    parser.add_argument('-X', "--no-Xterm", dest='noXterm',
                         action='store_true', default=False,
                         help="disable X-term for batch mode")

    args = parser.parse_args()
    args.deepName='Genom Assayer v1.0'
    for arg in vars(args):  print( 'myArg:',arg, getattr(args, arg))
    return args

#=================================
#=================================
#  M A I N 
#=================================
#=================================
args=get_parser()
#args.webPath='web/'

specNameL=['NC_010119.1','NC_017232.1'] # hard
specNameL=['NZ_AVBO01000047.1','NZ_FIXO01000005.1','NZ_JUHW01000036.1'] #main

#pub=PubAll_Plasmid(args ,specNameL=specNameL[:])
#pub=PubAll_Plasmid(args , maxSpec=20)
pub=PubAll_Plasmid(args )

pub.species_loop()
pub.summary()

pub.coverHTML()

