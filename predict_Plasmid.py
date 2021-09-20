#!/usr/bin/env python
""" 
 classify sequences as Plasmids or Chromosomes based on a training model.
 classify a set of scaffold from an input fasta file.
 uses a series of pre trained models (kfolds), usually 10 or 12.
 runs prediction with each model, then average the score over the models.
 outputs a text file named predictions.txt.
 How to run:
     python3 $PARENT/predict_Plasmid.py --dataPath $OUT/dlDataPath  --given test  --kModelList  18-22  --dataSegment -1  --seedModel   /global/cscratch1/sd/andreopo/plasmid4o-  
"""
__author__ = "Bill Andreopoulos, Jan Balewski"
__email__ = "wandreopoulos@lbl.gov, janstar1122@gmail.com"

from DL_Model import DL_Model
from Util_Plasmid import write_yaml, read_yaml

import yaml
import math
from pprint import pprint
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse

import Constants

def get_parser():
    parser = argparse.ArgumentParser(
        description='classify sequences as Plasmids or Chromosomes based on a training model',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-v","--verbosity",type=int,choices=[0, 1, 2],
                        help="increase output verbosity", default=0, dest='verb') 

    parser.add_argument("--dataPath",
                        default='dataBig',help="path to input")

    parser.add_argument('-g',"--given", choices=['main','plasm','test'],
                        default='test',help="whether the file contains known classified items (main or plasm) or test (unknown) items")

    parser.add_argument('-s',"--dataSegment", type=int,
                        default=-1,help="which data segment to process, -1 for process all")

    parser.add_argument("-m","--kModelList",nargs="+",
                        default=['10'],
                        help=" blank separated list of kfold IDs, takes n1-n2")

    parser.add_argument("--outPath",
                        default='outPR',help="output path for plots/output")
 
    parser.add_argument("--seedModel",
                        default='/global/cscratch1/sd/andreopo/plasmid_sum/4g/plasmid4g-',
                        help="trained model and weights")

    parser.add_argument("-a", "--arrIdx", type=int, default=1,
                        help="slurm array index (unused in prediction at the moment)")

    parser.add_argument("-k", "--kfoldOffset", type=int, default=0,
                        help="decides which segments merge for training (unused in prediction at the moment)")

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
given_class=args.given
#This file was written in the function DL_Model::split_into_folds
inpF=args.dataPath+'/'+ Constants.project +'.%s-scaff-split.yml'%given_class
assert(os.path.exists(inpF))
bulk=read_yaml(inpF)


# select how many scaffolds to be processed
if args.dataSegment>=0:
    all_scaffolds=bulk[args.dataSegment]
else:
    all_scaffolds={}
    for seg in range(Constants.numSegm_total):
        if seg in bulk:
            #print('ww',seg,type(bulk[seg]),len(all_scaffolds))
            all_scaffolds.update(bulk[seg])
            ###it_will_crash_make_5

score_thr=0.50
max_scaff=Constants.mxScaf_pred 
print('M: seg:',args.dataSegment,' num scaffolds:', len(all_scaffolds),' score_thr=',score_thr,' max_scaff=',max_scaff)

deep=DL_Model(args)
deep.load_Kmodels(path=args.seedModel,kL=args.kModelList)

#print('work on ', all_scaffolds.keys())

scores_samples_list_class1=[]  # per sample score
scores_samples_list_class0=[]  # per sample score
scores_scaffs_list=[]  # per scaffold score
Yclass = []  #Known class based on header text

#import numpy as np
#The cnt dictionary counts how many sequences are in each class
cnt={'inp':0,'plasmid':0,'main':0,'NNjunk':0,'ambig':0}

#create the predictions output file
rootdir = os.getcwd()
output_predix = os.path.join( args.outPath , "predictions.txt" )
f_predix = open(output_predix, 'a')
#f_predix.write("name,pred,conf\n")

#iterate through all sequences
for scaffN in all_scaffolds:
    if cnt['inp'] >=max_scaff: break
    cnt['inp']+=1
    sequenceString=all_scaffolds[scaffN]['seq']
    txt = all_scaffolds[scaffN]['text'].lower()
    
    ############################
    ###Subsample!
    ###The number of sequences to sample depends on sampling_rate (length-dependent)
    ###times the target_samples_per_contig_pred or target_samples_per_contig_train
    ############################
    sampling_rate= 0.5 + 0.5*math.sqrt(len(sequenceString)/1e4)
    sampList=deep.sample_scaffold(sequenceString,Constants.target_samples_per_contig_pred*sampling_rate)
    print("len(sampList) %s" % ( len(sampList) ))
    if args.verb > 0: print("sample: %s" % (sampList), Constants.target_samples_per_contig_pred, sampling_rate)
    assert len(sampList) >1 #0 # must ahve few samples to compute the average

    floatD=all_scaffolds[scaffN]['features']
    is_digit = lambda str: str.lstrip('-').replace('.', '').isdigit()
    flatten = lambda l: [float(item) for sublist in l for item in sublist  if is_digit(item)]
    #flatten = lambda l: [float(item) for sublist in l for item in sublist]
    #flatten = lambda l: [float(item) for sublist in l for item in sublist]
    b=flatten( [ list(floatD[x]) for x in Constants.globFeatureL if type(floatD[x]) is not float ] )
    a=[ floatD[x] for x in Constants.globFeatureL if type(floatD[x]) is float ]
    featureList=a+b #[ floatD[x] for x in Constants.globFeatureL if type(floatD[x]) is float else list(floatD[x]) ]
    print("featureList %s" % (featureList))
    # replicate the same scaffold features for all samples (for sequence scaffN)
    featureListXsamples=[featureList]*len(sampList)
    
    #If too few samples from scaffN, skip
    if len(sampList)< Constants.target_samples_per_contig_pred/3. :
        cnt['NNjunk']+=1
        f_predix.write("%s,%s,%s\n" %(txt,"NNJUNK",1))
        continue
    classif_details_yaml={'data_info':{'size':len(sequenceString),'scaffName':scaffN,'given':given_class}}
    #Get samples one hot encoded, and feature list for each sample
    sampListHotEncodedA,featureListXsamplesA=deep.build_data_classify_one(sampList,featureListXsamples)
    #myNoise=0.02
    #featureListXsamplesA+=np.random.uniform(-myNoise,myNoise,size=(featureListXsamplesA.shape))

    #dump_ana_details(scaffN,sampList,featureListXsamples,sampListHotEncodedA,featureListXsamplesA)

    classif_details_yaml['ora_inp']={'nSample':sampListHotEncodedA.shape[0],'seqLen':sampListHotEncodedA.shape[1],'featureListen':featureListXsamplesA.shape[1]}
    Yscores_samples,classif_avrg_score_str,classif_score=deep.classify_one_scaffold(sampListHotEncodedA,featureListXsamplesA,verb=args.verb)

    scores_scaffs_list.append(classif_score['avr'])
    if txt.find('plasmid') > -1:
         Yclass.append(1)
         scores_samples_list_class1+=Yscores_samples.tolist()
    else:
         Yclass.append(0)
         scores_samples_list_class0+=Yscores_samples.tolist()

    if classif_score['avr']>score_thr+classif_score['err']*2:
        decision='PLASM'
        cnt['plasmid']+=1
        f_predix.write("%s,%s,%s\n" %(txt,"PLASMID",classif_avrg_score_str))
    elif classif_score['avr']<score_thr-classif_score['err']*2:
        decision= 'MAIN'
        cnt['main']+=1
        f_predix.write("%s,%s,%s\n" %(txt,"GENOME",classif_avrg_score_str))
    else:
        decision= 'AMBIG'
        cnt['ambig']+=1
        f_predix.write("%s,%s,%s\n" %(txt,"AMBIGUOUS",classif_avrg_score_str))

    print(cnt['inp'],scaffN,'decision=%s, avr score:'%decision, classif_avrg_score_str,' len=%.1fk samples=%d'%(len(sequenceString)/1000.,len(sampList)))
    if args.verb > 0: print(Yscores_samples,classif_avrg_score_str,classif_score  , scores_samples_list_class1, scores_samples_list_class0, scores_scaffs_list, Yclass , sampListHotEncodedA,featureListXsamplesA)

    classif_details_yaml['score']=classif_score
    classif_details_yaml['model_info']=deep.info
    #print('out classif_details_yaml='); pprint(classif_details_yaml)

    #write under outPR the info for this scaffold prediction in a yaml file
    write_yaml(classif_details_yaml,args.outPath+'/%s.assayer.yml'%(scaffN),0)
    
nClass=cnt['plasmid']+cnt['main']+cnt['ambig']
#print('M:%s endCnt:'%given_class,cnt,'  fraction: Plasm=%.3f Ambig=%.3f  Main=%.3f'%(cnt['plasmid']/nClass,cnt['ambig']/nClass,cnt['main']/nClass))
print('Counts: Plasm=%s  Ambig=%s  Main=%s  nCount=%s'%(cnt['plasmid'],cnt['ambig'],cnt['main'],nClass))

f_predix.close()

# make plot of all scores
#ROC curve
#print("Yclass %s " % (Yclass))
#print("scores_scaffs_list %s " % (scores_scaffs_list))
####Print the AUC-ROC and FP-TP rate list to a file
if 1 in Yclass and 0 in Yclass:  ###Need both classes for this to be meaningful
    fpr, tpr, _ = roc_curve(np.array(Yclass).reshape(-1,1), np.array(scores_scaffs_list).reshape(-1,1))
    roc_auc = auc(fpr, tpr)
    fptpauc = open(os.path.join( args.outPath , "predictions_fp_tp_auc.txt" ), 'w' )
    fptpauc.write("#"+str(roc_auc)+"\n" + "\n".join(map(str, zip(fpr, tpr)))   )
    fptpauc.close()


'''
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

#Only save plots if the user requested them.
if not args.noXterm:
    plt.figure()
    #Histogram plot. This hist plots the per-scaffold scores for each class (plasmid vs. chrom) with different colors.
    bins = np.linspace(0, 1, 100)
    y_score_plas = np.array(scores_scaffs_list)[np.array(Yclass) == 1]
    y_score_gen = np.array(scores_scaffs_list)[np.array(Yclass) == 0]
    plt.hist(y_score_plas, bins, alpha=0.5, label='plasmid', density=True)
    plt.hist(y_score_gen, bins, alpha=0.5, label='genomic', density=True)
    plt.legend(loc="upper left", prop={'size': 11})
    plt.xlabel('Score')
    plt.ylabel('Counts')
    plt.title('Per-scaffold average Scores Deep Learning on test set')
    plt.show()
    plt.savefig(os.path.join( args.outPath , "scaffscore_hist.png"))
    
    plt.figure()
    #Hist plot for per-sample score for each class. scores_samples_list_class0 and scores_samples_list_class1.
    bins = np.linspace(0, 1, 100)
    y_score_plas = np.array(scores_samples_list_class1)
    y_score_gen = np.array(scores_samples_list_class0)
    plt.hist(y_score_plas, bins, alpha=0.5, label='plasmid', density=True)
    plt.hist(y_score_gen, bins, alpha=0.5, label='genomic', density=True)
    plt.legend(loc="upper left", prop={'size': 11})
    plt.xlabel('Score')
    plt.ylabel('Counts')
    plt.title('Per-sample Scores Deep Learning on test set')
    plt.show()
    plt.savefig(os.path.join( args.outPath , "samplescore_hist.png"))
    
    
    from Plotter_Plasmid import Plotter_Plasmid
    ppp=Plotter_Plasmid(args)
    y1=0
    if args.given=='plasm': y1=1
    
    #Ytrue1=Yclass ###[y1]*len(scores_samples_list)
    #ppp.plot_labeled_scores(Ytrue1, scores_samples_list,'samples:'+args.given,score_thr=score_thr,figId=21)
    #Ytrue2=Yclass ###[y1]*len(scores_scaffs_list)
    ppp.plot_labeled_scores(Yclass, scores_scaffs_list,'scaffolds:'+args.given,score_thr=score_thr,figId=22)
    ppp.pause(args,'predict') 
    #Ytrue=[y1]*len(Yscore_sum)
    #ppp.plot_labeled_scores(Ytrue, Yscore_sum,'given='+args.given,score_thr=score_thr)
    #ppp.pause(args,'predict') 
    
