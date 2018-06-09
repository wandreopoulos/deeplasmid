#!/usr/bin/env python
""" 
lenear single scaffold classifier
 uses a series of pre trained models (kModels), 
 run prediction with each model, then average the score
"""
__author__ = "Jan Balewski"
__email__ = "janstar1122@gmail.com"

from Oracle_Plasmid import Oracle_Plasmid
from Util_Plasmid import write_yaml, read_yaml 

import yaml
import math
from pprint import pprint

import argparse
def get_parser():
    parser = argparse.ArgumentParser(
        description='classify Plasmids based on training',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-v","--verbosity",type=int,choices=[0, 1, 2],
                        help="increase output verbosity", default=1, dest='verb') 
    parser.add_argument("--project",
                        default='assayer3',dest='prjName',
                        help="core name used to store outputs")

    parser.add_argument("--dataPath",
                        default='data',help="path to input/output")

    parser.add_argument('-g',"--given", choices=['main','plasm'],
                        default='plasm',help="fixed protein role")

    parser.add_argument('-s',"--dataSegment", type=int,
                        default=-1, help="which data segment to process")

    parser.add_argument("-k","--kModelList",nargs="+",
                        default=['12'],
                        help=" blank separated list of kfold IDs, takes n1-n2 n3 n5")

    parser.add_argument("--outPath",
                        default='out',help="output path for plots/output")
 
    parser.add_argument("-n", "--events", type=int, default=100,
                        help="num seqences for classification")
 
    parser.add_argument("--seedModel",
                        default='/global/cscratch1/sd/balewski/plasmid_sum/3a/plasmid3a-',
                        help="trained model and weights")

    parser.add_argument("-i", "--arrIdx", type=int, default=1,
                        help="slurm array index")

    parser.add_argument('-X', "--no-Xterm", dest='noXterm',
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
#args.dataPath='dataPlasmSmall/'

role=args.given
inpF=args.dataPath+'/'+args.prjName+'.%s-scaff-split.yml'%role
bulk=read_yaml(inpF)


# select how many scaffolds to be processed
if args.dataSegment>=0:
    scaffD=bulk[args.dataSegment]
else:
    scaffD={}
    for seg in range(10): 
        #print('ww',seg,type(bulk[seg]),len(scaffD))
        scaffD.update(bulk[seg])


print('M: seg:',args.dataSegment,' num scaffolds:', len(scaffD))

ora=Oracle_Plasmid(args)
ora.load_models(path=args.seedModel,kL=args.kModelList)

#print('work on ', scaffD.keys())

cnt={'inp':0,'plasmid':0,'main':0,'NNjunk':0,'ambig':0}
for scaffN in scaffD:
    cnt['inp']+=1
    seqStr=scaffD[scaffN]['seq']
    seqLen=len(seqStr)
    if seqLen<5e4 : continue
    #if seqLen>3000 : continue
    print( cnt['inp'],'M:scaff',role,scaffN,seqLen)
    #if scaffN!='NC_009431.1': continue #ambigous
    if scaffN!='NC_005307.1': continue #ambigous+double peak
    #if cnt['inp']<750: continue

    sampL=ora.sample_scaffold(seqStr,args.events,isRandom=0)
    if len(sampL)< args.events/3. :
        cnt['NNjunk']+=1
        continue
    rec={'data_info':{'size':len(seqStr),'scaffName':scaffN,'given':role}}

    Xhot=ora.build_data_one(sampL)
    rec['ora_inp']={'nSample':Xhot.shape[0],'seqLen':Xhot.shape[1]}
    Yscore,avr_str,rec1=ora.classify_one_scaffold(Xhot,isLinear=1,verb=0)

    if rec1['avr']>0.5+rec1['err']*2:
        decision='PLASM'
        cnt['plasmid']+=1
    elif rec1['avr']<0.5-rec1['err']*2:
        decision= 'MAIN'
        cnt['main']+=1
    else:
        decision= 'AMBIG'
        cnt['ambig']+=1

    print(cnt['inp'],scaffN,'decision=%s, avr score:'%decision, avr_str,' len=%.1fk sampl=%d'%(len(seqStr)/1000.,len(sampL)))
    rec['score']=rec1
    rec['model_info']=ora.info
    #print('out rec='); pprint(rec)

    write_yaml(rec,args.outPath+'/%s/%s.assayer.yml'%(role,scaffN),0)
    #exit(4)
nClass=cnt['plasmid']+cnt['main']+cnt['ambig']
print('M:%s endCnt:'%role,cnt,'  fraction: Plasm=%.3f Ambig=%.3f  Main=%.3f'%(cnt['plasmid']/nClass,cnt['ambig']/nClass,cnt['main']/nClass))
