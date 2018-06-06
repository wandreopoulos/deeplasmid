#!/usr/bin/env python
""" read input 
train net
write net + weights as HD5

"""
__author__ = "Jan Balewski"
__email__ = "janstar1122@gmail.com"

from Plotter_Plasmid import Plotter_Plasmid
from Deep_Plasmid import Deep_Plasmid

import argparse
def get_parser():
    parser = argparse.ArgumentParser(
        description='Perform Plasmid Oracle training',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--arrIdx", type=int, default=1,
                        help="slurm array index")
    parser.add_argument("-v","--verbosity",type=int,choices=[0, 1, 2],
                        help="increase output verbosity", default=1, dest='verb')
 
    parser.add_argument("--project",
                        default='assayer4',dest='prjName',
                        help="core name used to store outputs")

    parser.add_argument("--dataPath",
                        default='data',help="path to input/output")
    parser.add_argument("--outPath",
                        default='out/',help="output path for plots and tables")

    parser.add_argument("--seedModel",
                        default=None,
                        help="seed model and weights")
 
    parser.add_argument("--seedWeights",
                        default=None,
                        help="seed weights only, after model is created, for re-trainig only ")

    parser.add_argument("-k", "--kfoldOffset", type=int, default=0,
                        help="decides which segments merge for training")

    parser.add_argument( "-X","--noXterm", dest='noXterm',
                         action='store_true', default=False,
                         help="disable X-term for batch mode")
    parser.add_argument("-e", "--epochs", type=int, default=2,
                        help="num epochs")

    parser.add_argument("-b", "--batch_size", type=int, default=200,
                        help="fit batch_size")
    parser.add_argument("-n", "--events", type=int, default=3000,
                        help="num seqences for training")
    parser.add_argument("--dropFrac", type=float, default=0.1,
                        help="drop fraction at all layers")

    parser.add_argument( "-s","--earlyStop", type=int,
                         dest='earlyStopPatience', default=10,
                         help="early stop:  epochs w/o improvement (aka patience), 0=off")
    parser.add_argument( "--checkPt", dest='checkPtOn',
                         action='store_true',default=False,help="enable check points for weights")

    parser.add_argument( "--reduceLr", dest='reduceLearn',
                         action='store_true',default=False,help="reduce learning at plateau")


    args = parser.parse_args()
    for arg in vars(args):  print( 'myArg:',arg, getattr(args, arg))
    return args

#=================================
#=================================
#  M A I N 
#=================================
#=================================
args=get_parser()

ppp=Plotter_Plasmid(args )
deep=Deep_Plasmid(args)

for role in deep.mr12:
    deep.prep_labeled_input(role,'val')  
    deep.prep_labeled_input(role,'train')
    
deep.build_training_data()

'''
ppp.plot_float_featuresB(deep,'val',0,20)
ppp.plot_float_featuresB(deep,'val',1,21)
ppp.plot_float_featuresB(deep,'train',0,22)
ppp.plot_float_featuresB(deep,'train',1,23)
'''

if args.seedModel==None:
    print('start fresh model')
    deep.build_compile_model(args)
else:
    deep.load_compile_1model(path=args.seedModel)

if args.seedWeights:
    deep.load_weights_4_retrain(path=args.seedWeights)


ppp.plot_model(deep,1) # depth:0,1,2

if args.epochs >3:  deep.save_model() 

if args.epochs >0:
    deep.train_model(args) 
    deep.save_model() 
    deep.save_training_history() 
    ppp.plot_train_history(deep,args,10)
    ppp.plot_AUC('val',deep,args,10)
else:
    y_true, y_pred=ppp.plot_AUC('val',deep,args)
    ppp.plot_labeled_scores(y_true, y_pred,'seg=val',score_thr=0.65)

ppp.pause(args,'train') 




