#!/usr/bin/env python
""" 
 classify a scaffold 
 uses a series of pre trained models (kfolds), 
 run prediction with each model, then average the score
"""
__author__ = "Jan Balewski"
__email__ = "janstar1122@gmail.com"

from Oracle_Plasmid import Oracle_Plasmid
from Util_Plasmid import write_yaml, read_yaml , dump_ana_details

import yaml
import math
from pprint import pprint
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import numpy as np


import argparse
def get_parser():
    parser = argparse.ArgumentParser(
        description='classify Plasmids based on training',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-v","--verbosity",type=int,choices=[0, 1, 2],
                        help="increase output verbosity", default=1, dest='verb') 
    parser.add_argument("--project",
                        default='assayer4',dest='prjName',
                        help="core name used to store outputs")

    parser.add_argument("--dataPath",
                        default='dataBig',help="path to input")

    parser.add_argument('-g',"--given", choices=['main','plasm','test'],
                        default='plasm',help="fixed protein role")

    parser.add_argument('-s',"--dataSegment", type=int,
                        default=5,help="which data segment to process")

    parser.add_argument("-k","--kModelList",nargs="+",
                        default=['10'],
                        help=" blank separated list of kfold IDs, takes n1-n2")

    parser.add_argument("--outPath",
                        default='outPR',help="output path for plots/output")
 
    parser.add_argument("-n", "--events", type=int, default=100,
                        help="approx num samples per scaffolds")
 
    parser.add_argument("--seedModel",
                        default='/global/cscratch1/sd/andreopo/plasmid_sum/4g/plasmid4g-',
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
    for seg in range(6):
        if seg in bulk:
            #print('ww',seg,type(bulk[seg]),len(scaffD))
            scaffD.update(bulk[seg])
            ###it_will_crash_make_5

score_thr=0.50
max_scaff=100000
print('M: seg:',args.dataSegment,' num scaffolds:', len(scaffD),' score_thr=',score_thr,' max_scaff=',max_scaff)

ora=Oracle_Plasmid(args)
ora.load_Kmodels(path=args.seedModel,kL=args.kModelList)

#print('work on ', scaffD.keys())

Yscore_sum=[]
Yscore_sum1=[]  # per sample
Yscore_sum2=[]  # per scaffold
Yclass = []  #Known class based on header text

#import numpy as np
cnt={'inp':0,'plasmid':0,'main':0,'NNjunk':0,'ambig':0}

rootdir = os.getcwd()
output_predix = os.path.join( rootdir , "predictions.txt" )
f_predix = open(output_predix, 'a')
#f_predix.write("name,pred,conf\n")

for scaffN in scaffD:
    if cnt['inp'] >=max_scaff: break
    cnt['inp']+=1
    seqStr=scaffD[scaffN]['seq']
    txt = scaffD[scaffN]['text'].lower()
    
    sampFact= 0.5 + 0.5*math.sqrt(len(seqStr)/1e4)
    sampL=ora.sample_scaffold(seqStr,args.events*sampFact)
    print("len(sampL) %s" % ( len(sampL) ))
    assert len(sampL) >10 # must ahve few samples to compute the average

    floatD=scaffD[scaffN]['features']
    floatL=[ floatD[x] for x in ora.globFeatureL ]
    Xfloat=[floatL]*len(sampL) # the same scaffold features for all samples
    
    if len(sampL)< args.events/3. :
        cnt['NNjunk']+=1
        f_predix.write("%s,%s,%s\n" %(txt,"NNJUNK",1))
        continue
    rec={'data_info':{'size':len(seqStr),'scaffName':scaffN,'given':role}}
    XhotA,XfloatA=ora.build_data_one(sampL,Xfloat)
    #myNoise=0.02
    #XfloatA+=np.random.uniform(-myNoise,myNoise,size=(XfloatA.shape))

    #dump_ana_details(scaffN,sampL,Xfloat,XhotA,XfloatA)

    rec['ora_inp']={'nSample':XhotA.shape[0],'seqLen':XhotA.shape[1],'floatLen':XfloatA.shape[1]}
    Yscore,avr_str,rec1=ora.classify_one_scaffold(XhotA,XfloatA,verb=0)
    Yscore_sum+=Yscore.tolist()
    Yscore_sum1+=Yscore.tolist()
    Yscore_sum2.append(rec1['avr'])
    if txt.find('plasmid') > -1:
         Yclass.append(1)
    else:
         Yclass.append(0)

    if rec1['avr']>score_thr+rec1['err']*2:
        decision='PLASM'
        cnt['plasmid']+=1
        f_predix.write("%s,%s,%s\n" %(txt,"PLASMID",avr_str))
    elif rec1['avr']<score_thr-rec1['err']*2:
        decision= 'MAIN'
        cnt['main']+=1
        f_predix.write("%s,%s,%s\n" %(txt,"GENOME",avr_str))
    else:
        decision= 'AMBIG'
        cnt['ambig']+=1
        f_predix.write("%s,%s,%s\n" %(txt,"AMBIGUOUS",avr_str))

    print(cnt['inp'],scaffN,'decision=%s, avr score:'%decision, avr_str,' len=%.1fk sampl=%d'%(len(seqStr)/1000.,len(sampL)))
    rec['score']=rec1
    rec['model_info']=ora.info
    #print('out rec='); pprint(rec)

    write_yaml(rec,args.outPath+'/%s.assayer.yml'%(scaffN),0)
    
nClass=cnt['plasmid']+cnt['main']+cnt['ambig']
#print('M:%s endCnt:'%role,cnt,'  fraction: Plasm=%.3f Ambig=%.3f  Main=%.3f'%(cnt['plasmid']/nClass,cnt['ambig']/nClass,cnt['main']/nClass))
print('Counts: Plasm=%s  Ambig=%s  Main=%s  nCount=%s'%(cnt['plasmid'],cnt['ambig'],cnt['main'],nClass))

f_predix.close()

# make plot of all scores
#ROC curve
#print("Yclass %s " % (Yclass))
#print("Yscore_sum2 %s " % (Yscore_sum2))

'''
fpr, tpr, _ = roc_curve(np.array(Yclass).reshape(-1,1), np.array(Yscore_sum2).reshape(-1,1))
roc_auc = auc(fpr, tpr)
plt.figure()
lw = 2
plt.plot(fpr, tpr, color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Deep Learning on test set')
plt.legend(loc="lower right")
#plt.show()
plt.savefig("roc_curve.png")
'''

plt.figure()
#histogram
bins = np.linspace(0, 1, 100)
y_score_plas = np.array(Yscore_sum2)[np.array(Yclass) == 1]
y_score_gen = np.array(Yscore_sum2)[np.array(Yclass) == 0]
plt.hist(y_score_plas, bins, alpha=0.5, label='plasmid')
plt.hist(y_score_gen, bins, alpha=0.5, label='genomic')
plt.legend(loc="upper left", prop={'size': 11})
plt.xlabel('Score')
plt.ylabel('Counts')
plt.title('Scores Deep Learning on test set')
plt.show()
plt.savefig("hist.png")


from Plotter_Plasmid import Plotter_Plasmid
ppp=Plotter_Plasmid(args )
y1=0
if args.given=='plasm': y1=1

#Ytrue1=Yclass ###[y1]*len(Yscore_sum1)
#ppp.plot_labeled_scores(Ytrue1, Yscore_sum1,'samples:'+args.given,score_thr=score_thr,figId=21)
Ytrue2=Yclass ###[y1]*len(Yscore_sum2)
ppp.plot_labeled_scores(Ytrue2, Yscore_sum2,'scaffolds:'+args.given,score_thr=score_thr,figId=22)
ppp.pause(args,'predict') 
#Ytrue=[y1]*len(Yscore_sum)
#ppp.plot_labeled_scores(Ytrue, Yscore_sum,'given='+args.given,score_thr=score_thr)
#ppp.pause(args,'predict') 

