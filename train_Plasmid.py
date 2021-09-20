#!/usr/bin/env python
""" read input 
train net
write net + weights as HD5
Run as:
   /usr/bin/time -v python -u ./train_Plasmid.py    --epochs 20 --events 601000 --arrIdx $arrIdx  --checkPt --kfoldOffset $kfoldOff --verbosity 0  --noXterm  --reduceLr  --outPath ./   >&fit_${arrIdx}.log
Note: see batchTrain.slr . This assumes pwd is the run directory for the specific fold, kFoldOffset is the index of the fold, dataPath (data subdirectory) has all the data.
"""
__author__ = "Bill Andrepoulos, Jan Balewski"
__email__ = "wandreopoulos@lbl.gov, janstar1122@gmail.com"

from Plotter_Plasmid import Plotter_Plasmid
from DL_Model import DL_Model

import argparse

import Constants


def get_parser():
    parser = argparse.ArgumentParser(
        description='Perform Plasmid Oracle training',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-v","--verbosity",type=int,choices=[0, 1, 2],
                        help="increase output verbosity", default=1, dest='verb')

    parser.add_argument("--dataPath",
                        default='data',help="path to input/output")

    parser.add_argument("--outPath",
                        default='outPR',help="output path for plots/output")

    parser.add_argument("-a", "--arrIdx", type=int, default=1,
                        help="slurm array index")

    parser.add_argument("-k", "--kfoldOffset", type=int, default=0,
                        help="decides which segments merge for training")

    parser.add_argument( "-X","--noXterm", dest='noXterm',
                         action='store_true', default=False,
                         help="disable X-term for batch mode")


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
deep=DL_Model(args)

# For each fold k determines which segments are train and which are val. Do this for each class of "main" and "plasmid"
for classif in Constants.classes:
    deep.prep_labeled_input(classif,'val')
    deep.prep_labeled_input(classif,'train')

deep.build_training_data()

'''
ppp.plot_float_featuresB(deep,'val',0,20)
ppp.plot_float_featuresB(deep,'val',1,21)
ppp.plot_float_featuresB(deep,'train',0,22)
ppp.plot_float_featuresB(deep,'train',1,23)
'''

print('start fresh model')
deep.build_compile_model(args)


ppp.plot_model(deep,1) # depth:0,1,2

deep.save_model() 
deep.train_model(args) 
deep.save_model() 
deep.save_training_history() 


#####Plots##############
###Plot various statistics about the training, such as descrease in loss over the epochs, and AUC.
ppp.plot_train_history(deep,args,10)
ppp.plot_AUC('val',deep,args,10)
y_true, y_pred=ppp.plot_AUC('val',deep,args)
ppp.plot_labeled_scores(y_true, y_pred,'seg=val',score_thr=0.65)

ppp.pause(args,'train') 


