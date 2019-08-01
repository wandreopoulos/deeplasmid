# small operations w/o dependencies
# keep double indent - in case you want to make a class of it

import os, time
import shutil
import yaml

import os, time, sys
(va,vb,vc,vd,ve)=sys.version_info ; assert(va==3)  # needes Python3
import warnings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' #Hide messy TensorFlow warnings
warnings.filterwarnings("ignore") #Hide messy Numpy warnings

start = time.time()

from keras.models import Model, load_model
from sklearn.metrics import roc_curve, auc,confusion_matrix
from scipy.stats import skew

import datetime
from pprint import pprint
import numpy as np
import yaml

print('deep-libs1 imported elaT=%.1f sec'%(time.time() - start))

__author__ = "Jan Balewski"
__email__ = "janstar1122@gmail.com"




    #............................
    #write a yaml file
    #............................
def write_yaml(rec,ymlFn,verb=1):
        start = time.time()
        ymlFd = open(ymlFn, 'w')
        yaml.dump(rec, ymlFd, Dumper=yaml.CDumper)
        ymlFd.close()
        xx=os.path.getsize(ymlFn)/1048576
        if verb:
                print('  closed  yaml:',ymlFn,' size=%.2f MB'%xx,'  elaT=%.1f sec'%(time.time() - start))

    #............................
    #read a yaml file
    #............................
def read_yaml(ymlFn):
        start = time.time()
        ymlFd = open(ymlFn, 'r')
        bulk=yaml.load( ymlFd, Loader=yaml.CLoader)
        ymlFd.close()

        print('  read  yaml:',ymlFn,' size=%d'%len(bulk),'  elaT=%.1f sec'%(time.time() - start))
        return bulk


    #............................
    #collect all yaml files under a directory
    #............................
def get_glob_info_files(dir0):
        allL=os.listdir(dir0)
        print('get_glob_info from',dir0)
        outD={}
        for x in allL:
            if 'yml' not in x: continue
            #print('aa',x)
            core=x[:-4]
            fn=dir0+'/'+x
            #print('cc',core,fn)
            outD[core]=fn

        print('found %d glob-features files '%len(outD), 'in ',dir0)
        assert len(outD)>0
        return outD


    #............................
    #normalize the features
    #TODO Check for correctness
    #............................
def normalize_features(rawD):
        outD={}
        for xN in rawD:
            val=rawD[xN]
            if xN=='len_sequence':
                xx=min(1.,val/3e5) #30k is the max length
                val=max(0.007,xx) #0.007 is the min accepted value for this feature
            elif xN=='genecount':
                xx=min(1.,val/1e3) #1k is the max genecount
                val=max(0.007,xx) #0.007 is the min accepted value for this feature
            elif xN=='aalenavg':
                xx=min(1.,val/1e3) #1k is the max aalenavg
                val=max(0.007,xx) #0.007 is the min accepted value for this feature
            elif 'longestHomopol' in xN:
                val=min(1,val/15)  #we don't expect homopols of length >15 and 1 is the max accepted value for this feature
            elif 'totalLongHomopol'  in xN:
                val=min(1,val/1000.)   #we don't expect homopols of length >1000 and 1 is the max accepted value for this feature
            else: 
                #val-=0.5
                outD[xN]=val #other vals such as GC plassketch pentamer are normalized from [0,1] to [-0.5,0.5]
                continue
            val-=0.5
            outD[xN]=val  # now val is in range [-0.5,0.5]
        return outD


