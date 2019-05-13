'''
########################
###This is the file of constants that the user does not need to provide every time.
##All constants used throughout the software 
########################
'''
import numpy as np

__author__ = "Bill Andreopoulos"
__email__ = "wandreopoulos@lbl.gov"


#Number of processes to use.
PROCESSES = 32 ### multiprocessing.cpu_count() ###32



#core name used to store outputs
project="assayer4"


#maximum scaffold to predict on (for small datasets doesn't matter)
mxScaf_pred=122000 # set it below 1k for testing

#maximum scaffold to train on
mxScaf_train_class1=1100 # set it below 1k for testing
mxScaf_train_class2=mxScaf_train_class1*50 # set it below 1k for testing




#trained model and weights. Usually unused. Only used if retraining is to be done.
#Note the prediction script has its own seedModel (pretrained) for predictions.
seedModel_trained="/global/cscratch1/sd/andreopo/plasmid_sum/4g/plasmid4g-"



#Sample size in base-pairs
seqLenCut=300


#The minimum and maximum sequence lengths
minSeqLen = 1e3
maxSeqLen = 330e3


#Names of bases
padChar='N'
basesL=list('ACTG')
bases2L=list('NYKMVSHRWBD')# other partial decoddings, all will map to 'N'
basesSet=set(basesL)



# prepare 1-hot base for basess:
myDim=len(basesL)
oneHotBase={} # len=1+myDim, pad=00000000...
i=0
for x in basesL:
    oneHotBase[x]=np.eye(myDim)[i]
    i+=1
oneHotBase[padChar]=np.zeros(myDim)
for x in bases2L:
    oneHotBase[x]=np.zeros(myDim)
print('oneHot base, size=',len(oneHotBase), ', sample:')
for x in basesL:
    print('base:',x,'1-hot:',oneHotBase[x].tolist())
print('all bases :', list(oneHotBase.keys()))
print('use seqLenCut=',seqLenCut)

#Feature names
globFeatureL=['gc_content','len_sequence']
for x in basesL: #Note dont change this to the basesSet, which is unordered! Keep basesL
    globFeatureL.append(x+'_longestHomopol')
    globFeatureL.append(x+'_totalLongHomopol')



#The class names we are looking for
classes=['plasm','main']

#Segments used in training
#Total segments (training and validation)
numSegm_total=6

#Training segments 
# For example, 5 training+1 validation K-fold partition of scaffolds
numSegm_train=numSegm_total-1




#DL variables
# DL - epochs
epochs = 30

# DL - drop fraction
dropFrac = 0.1

# Early stop:  epochs w/o improvement (aka patience), 0=off
earlyStopPatience=10
                         
# whether to enable check points for weights
checkPtOn = False

#whether to reduce learning at plateau
reduceLearn = False

# Fit batch_size
batch_size = 200


####Sampling constants
#approx num samples per contig for training
target_samples_per_contig_training=601000

#approx num samples per contig for prediction
target_samples_per_contig_pred=100

