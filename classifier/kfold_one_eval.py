#!/usr/bin/env python
""" 
 evaluates one kfold for  HK-RR 
 output: yaml
 uses pre trained net
"""
__author__ = "Jan Balewski"
__email__ = "janstar1122@gmail.com"

from Deep_Plasmid import Deep_Plasmid
from Util_Plasmid import write_yaml, read_yaml 

import time
import argparse
def get_parser():
    parser = argparse.ArgumentParser(
        description='recommend best HK-RR pairs for given species',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-v","--verbosity",type=int,choices=[0, 1, 2],
                        help="increase output verbosity", default=1, dest='verb') 
    parser.add_argument("--project",
                        default='assayer4',dest='prjName',
                        help="core name used to store outputs")
    parser.add_argument("-i", "--arrIdx", type=int, default=1,
                        help="slurm array index")

    parser.add_argument("--dataPath",
                        default='data',help="path to input")
    parser.add_argument("--outPath",
                        default='outKF',help="output path for plots/output")

    parser.add_argument("-k", "--kfoldOffset", type=int, default=-1,
                        help="decides which segments merge for training")

    parser.add_argument("-n", "--events", type=int, default=8000*8,
                        help="num HK-RR pairs for training, use 0 for all")    

    parser.add_argument("--seedModel",
                        default='out/',
                        help="seed model and weights")

    parser.add_argument('-X', "--no-Xterm", dest='noXterm',
                         action='store_true', default=False,
                         help="disable X-term for batch mode")

    args = parser.parse_args()
    for arg in vars(args):  print( 'myArg:',arg, getattr(args, arg))
    args.noiseFract=0  # only for training
    args.tf_timeline=False
    return args


#=================================
#=================================
#  M A I N 
#=================================
#=================================
args=get_parser()
#args.seedModel='/global/cscratch1/sd/balewski/'

deep=Deep_Plasmid(args)

for role in deep.mr12:
    deep.prep_labeled_input(role,'val')  
    deep.prep_labeled_input(role,'test')
deep.build_training_data()

deep.load_Kmodels(path=args.seedModel)

domL=['val','test']
out={}
for dom in domL:
    start = time.time()
    print('\n------ process dom:',dom)
    sumD=deep.model_predict_domain(dom)
    sumD['kfold_offset']=args.kfoldOffset
    sumD['idx']=args.arrIdx
    out[dom]=sumD

write_yaml(out,deep.outPath+'/kfold_idx%d.yml'%args.arrIdx)

