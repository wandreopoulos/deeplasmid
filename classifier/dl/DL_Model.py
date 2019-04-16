import os, time
import warnings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' #Hide messy TensorFlow warnings
warnings.filterwarnings("ignore") #Hide messy Numpy warnings
from sklearn.metrics import roc_curve, auc, confusion_matrix

start = time.time()

from keras.callbacks import EarlyStopping, ModelCheckpoint ,ReduceLROnPlateau, Callback 
from keras.layers import Dense, Dropout,  LSTM, Input, concatenate
from keras.models import Model, load_model
import keras.backend as K
#import tensorflow as tf

import random
import numpy as np
import yaml
from scipy.stats import truncnorm, skew

from Util_Plasmid import write_yaml, read_yaml ,normalize_features

import datetime
from pprint import pprint

import Constants

#print('deep-libs2 imported, TF.ver=%s, elaT=%.1f sec'%(tf.__version__,(time.time() - start)))

__author__ = "Bill Andreopoulos, Jan Balewski"
__email__ = "wandreopoulos@lbl.gov, janstar1122@gmail.com"

#.............................
# This callback is used by the DL model for what happens after each epoch.
# For now it stores the learning rate for plotting. Callbacks are automatically called when an epoch ends.
##############################
class MyLearningTracker(Callback):
    def __init__(self):
        self.hir=[]   
    def on_epoch_end(self, epoch, logs={}):
        optimizer = self.model.optimizer
        lr = K.eval(optimizer.lr)
        self.hir.append(float(lr))

#............................
#............................
#............................
class DL_Model(object):

#............................
#............................
###The next code originated from Oracle class. #class Oracle_Plasmid(object):
#............................

    def __init__(self,args):
        self.runDate=datetime.datetime.now()

        self.name=Constants.project
        self.verb=args.verb
        print(self.__class__.__name__,', prj:',self.name)


        # - - - - - - - - - data containers - - - - - - - -
        # data divided by train/val , all species chopped to NNNbp samples
        self.trainvalid_data={} # [dom][classif][ seqL,globL,nameL]
        # 1-hot encoded data ready for training
        self.trainvalid_1hot_data={} # ???[X[],Y[]]*[seg]


        self.outPath=args.outPath
        self.arrIdx=args.arrIdx


        self.kfoldOffset=args.kfoldOffset
        self.trainvalid_1hot_dataPath=args.dataPath
        self.trainNoise=0.05 # added truncated gAuss to only to train-data gloft
        self.train_hirD={'acc': [],'loss': [],'lr': [],'val_acc': [],'val_loss': []}

        print('globFeatureL %d fixed order:'%len(Constants.globFeatureL),Constants.globFeatureL)

        assert os.path.isdir(self.outPath) # check it now to avoid lost time

        '''
        moved to Constants
        self.padChar='N'
        self.basesL=list('ACTG')
        Constants.bases2L=list('NYKMVSHRWBD')# other partial decoddings, all will map to 'N'
        Constants.basesSet=set(self.basesL)
        Constants.seqLenCut=300 # sample size in base-pairs
        Constants.globFeatureL=['gc_content','len_sequence']
        for x in 'ACTG' :
            Constants.globFeatureL.append(x+'_longestHomopol')
            Constants.globFeatureL.append(x+'_totalLongHomopol')

        self.mr12=['plasm','main']

        # prepare 1-hot base for basess:
        myDim=len(self.basesL)
        self.oneHotBase={} # len=1+myDim, pad=00000000...
        i=0
        for x in self.basesL:
            self.oneHotBase[x]=np.eye(myDim)[i]
            i+=1
        self.oneHotBase[self.padChar]=np.zeros(myDim)
        for x in Constants.bases2L:
            self.oneHotBase[x]=np.zeros(myDim)
        print('oneHot base, size=',len(self.oneHotBase), ', sample:')
        for x in self.basesL:
            print('base:',x,'1-hot:',self.oneHotBase[x].tolist())
        print('all bases :', list( self.oneHotBase.keys()))
        print('use seqLenCut=',Constants.seqLenCut)
        '''




    #............................
    # A list with samples from the the scaffold that is input
    # sourceSequence - the source sequence for the sampling
    # numSamplesF - the number of samples to collect as a float
    # returns seqSampleList sampling of seqs
    # Called from: training and prediction.
    #............................
    def sample_scaffold(self,sourceSequence,numSamplesF):
        numSamples=int(numSamplesF)
        seqLen=Constants.seqLenCut
        #the max possible number of chunks is the length of the source scaffold - sampled seq length
        maxChunks=len(sourceSequence) - seqLen
        #maxSampples is a bit less than maxChunks
        maxSamples=int(maxChunks*0.9)
        if numSamples>maxSamples: numSamples=maxSamples
        #print('qqq',maxChunks,numSamples)
        # pick numSamples random numbers of where to start the samples in the original scaffold from 0..maxChunks
        idxL=np.random.choice(maxChunks, numSamples, replace=False)
        #print('sample: idxL',idxL)
        seqSampleList = []
        for i0 in idxL:
            chunk=sourceSequence[i0:i0+seqLen]
            if 'NNNN' in chunk: continue  # not sure if it works - need larger input
            #print(i0,chunk, len(chunk))
            seqSampleList.append(chunk)
        #print('sss',len(seqL))
        return seqSampleList

    #............................
    # convert sequences to 1-hot 2-D arrays
    # called from add_scaff_to_1hot, which is called from build_data_classify_one, which is called from predict script.
    # Called from: training and prediction.
    #............................
    def encode1hotBase(self,seqPad): 
        hot2D=[]
        for y in seqPad:
            hot2D.append(Constants.oneHotBase[y])
        return np.array(hot2D).astype(np.float32)

    #............................
    #called from build_data_classify_one, which is called from predict script.
    #The purpose is just to copy the data from inSeq,inFeatures to outHotSeq,outFeatures and to convert seqs to one-hot encodings.
    #Returns outHotSeq,outFeatures where outHotSeq contains 1-hot encoded sequence.
    # Called from: training and prediction.
    #............................
    def list_to_nparray(self,outHotSeq,outFeatures,inSeq,inFeatures, offset=0):
        offs = offset
        #print('wq1',outHotSeq.shape,outFeatures.shape,len(inSeq),len(inFeatures))
        assert len(inSeq)==len(inFeatures)
        for seq,featureV in zip(inSeq,inFeatures):
            XhotOne=self.encode1hotBase(seq)
            outHotSeq[offs,:]=XhotOne[:]
            outFeatures[offs,:]=featureV[:]
            offs+=1
        #print(' add  shape:',len(Xseq),' ieve=',ieve)

    #............................
    #called from predict script.
    #Turn sampList into one how encoded vector and numpy array
    #Turn featureListXsamples into numpy array
    # Called from: prediction.
    #............................
    def build_data_classify_one(self, sampList,featureListXsamples):  ###seqL,Xfloat):
        #These are used for resizing the lists into numpy arrays
        numSamples=len(sampList)
        sampleSeqLen=Constants.seqLenCut
        num_bases=len(Constants.basesSet)
        num_features=len(Constants.globFeatureL)

        # clever list-->numpy conversion, Thorsten's idea
        # First turn into np.zeros of desired size, then fill with data
        XhotAll=np.zeros([numSamples,sampleSeqLen,num_bases],dtype=np.float32)
        featureListXsamplesAll=np.zeros([numSamples,num_features],dtype=np.float32)
        self.list_to_nparray(XhotAll,featureListXsamplesAll,sampList,featureListXsamples)
        #print('build_data_classify_one:, X:',XhotAll.shape,'done')
        return XhotAll,featureListXsamplesAll



    #.............................
    # This is called from predict script. It takes as input the resulting vectors from build_data_classify_one
    # Returns the results from the classification: all scores and the average string.
    # Called from: prediction.
    #.............................
    def classify_one_scaffold(self,Xs,Xf,verb=1):

        nK=len(self.model)
        #print('clsf  inputs shapes:',Xs.shape,Xf.shape,' using kModel=%d models'%(nK))
        start = time.time()

        Yscore=np.array(0)  # accumulator
        for k in self.model:
            Yscore=np.add(Yscore,self.model[k].predict([Xs,Xf]).flatten())
            if verb>1: print(k,' kModel done, score[0].sum=',Yscore[0])
        #print('   prediction  done,  elaT=%.1f sec'%((time.time() - start)))
        # renormalize score if multi-model predictionw as made
        if nK>1:
            Yscore=np.multiply(Yscore, 1./nK)

        avr=float(Yscore.mean())
        std=float(Yscore.std())
        err=float(std/np.sqrt(len(Yscore)-1))
        skew1=float(skew(Yscore))
        avr_str="%.3f +/- %.3f"%(avr,err)

        bins = np.linspace(0.0, 1., 21)
        hist_np,_= np.histogram(Yscore, bins)
        hist=[ float(x) for x in hist_np]

        rec={'avr':avr,'err':err,'std':std,'skew':skew1,'hist':hist}
        if verb >0: print(rec)

        return Yscore,avr_str,rec




    #............................
    #kL is the array of h5 model indices. For example it may contain 18-22,25.
    #It is called from the prediction script. It will parse the kL model indices and 
    # store the indices in kL.
    #load all the h5 models into self.model
    # modifies self.model
    ############################
    def load_Kmodels(self,path='.',kL=['']):
        # expand list if '-' are present
        kkL=[]
        for x in kL:
            if '-' not in x:
                kkL.append(x) ; continue
            xL=x.split('-')
            for i in range(int(xL[0]),int(xL[1])+1):
                kkL.append(i)
        #print(kL,'  to ',kkL)
        #store the indices in kL.
        kL=kkL

        nK=len(kL)
        assert nK>0
        self.model={}
        runDateStr=self.runDate.strftime("%Y-%m-%d_%H.%M")
        self.info={'modelPath':path,'kModelList':kL,'oracle':self.name,'runDate':runDateStr}

        print('load_kModels =',nK,kL)
        start = time.time()
        for k in kL:
            inpF5m=path+'%s/'%k+self.name+'.model.h5'
            print('load %d model and weights  from'%len(self.model),inpF5m,'  ... ')
            self.model[k]=load_model(inpF5m) # creates mode from HDF5
            if len(self.model)==1 :
                self.model[k].summary()
            else:
                assert self.model[k].count_params() == self.model[kL[0]].count_params()
        print('%d models loaded, elaT=%.1f sec'%(nK,(time.time() - start)))



    ##################################################
    ######End of the code that originated from Oracle class.
    ##################################################



    #............................
    #Read the scaffolds fasta file
    # returns scafInfoList - list of scaffolds, excluding the long and the short ones.
    # Called from: training and prediction format scripts.
    #
    # Parameters:
    # ------------
    # nPresc : If set to > 1 it sets that the nPresc'th sequences will be processed only. It will only process seqs such that the mod nPresc == 0, which means the nPresc'th sequence.
    #################################
    def read_scaffolds_fasta(self,fName,globFD,mxRec=0, nPresc=1):
        #minNucleoLen= 1000 #500
        #maxNucleoLen= 330e3 #2100e3   
        #set minimum and maximum length low enough and high enough respectively,
        #since the fasta file contains all sequences we want (filtered)
        minNucleoLen = 500
        maxNucleoLen = 2100e3
        minNucleoLen = Constants.minSeqLen #2e3
        maxNucleoLen = Constants.maxSeqLen #300e3



    
        #timestamp = strftime("%m%d%Y-%H%M%S")
        rootdir = os.getcwd()
        output_predix = os.path.join( self.outPath , "predictions.txt" )
        f_predix = open(output_predix, 'w')
        f_predix.write("name,pred,conf\n")
        

        # just used to validate the chars in a line.
        bSet=Constants.basesSet | set(Constants.bases2L)
        scafInfoList=[]     
        scafName=''
        sequenceLengthList=[] # only for accepted
        cnt={'any':0,'short':0,'long':0,'presc':0,'noFeatures':0,'accept':0}
        print('read %d from %s'%(mxRec,fName))   
        fp = open(fName, 'r')
        for line in fp:
            line = line.strip()
            if len(line)<1: continue  # skip empty lines
            if line.startswith('>'):#  header: Name+ Species+soem info
                if mxRec and len(scafInfoList)>=mxRec: break
                #print('raw',len(scafInfoList),mxRec)
                if len(scafName)>0:
                    #archive last protein
                    cnt['any']+=1
                    seqLen=len(seq)

                    if seqLen<= minNucleoLen:
                        cnt['short']+=1 
                        f_predix.write("%s,%s,%s\n" %(scafName,"SHORTER_%s"%(minNucleoLen),1))
                    elif seqLen >= maxNucleoLen:
                        cnt['long']+=1
                        f_predix.write("%s,%s,%s\n" %(scafName,"LONGER_%s"%(maxNucleoLen),1))
                    elif  cnt['any']% nPresc!=0:
                        cnt['presc']+=1
                    elif  scafName not in globFD:
                        cnt['noFeatures']+=1 
                        f_predix.write("%s,%s,%s\n" %(scafName,"NOFEATURE",1))
                    else:
                        gloftRaw=self.read_global_features(scafName,globFD) 
                        eps=np.abs(gloftRaw['len_sequence']-seqLen)
                        if eps>5:
                            print(scafName,len(scafInfoList),'bad eps=%.1g'% eps,gloftRaw['len_sequence'],seqLen)                       
                        # accpeted scaffold
                        cnt['accept']+=1
                        featureV=normalize_features(gloftRaw)
                        sequenceLengthList.append(seqLen)
                        scafInfoList.append( [scafName,text,featureV,seq])  
 
                        if cnt['accept']%1000==0:
                            print(len(scafInfoList),seqLen,scafName,featureV,seq[:100],'...')
                                
                      
                # init new scaffold
                seq=''
                # header is different for mito/main vs. plasmid/main
                lineL=line[1:].split()
                #print('r2',lineL)
                scafName = lineL[0].replace("/", "_")
                #assert len(lineL)>=2
                text=' '.join(lineL)
                
            else:    # - - - - it is the sequence
                line=line.upper()
                if not  bSet.issuperset(set(line)):
                        print( 'new base:',set(line) - bSet)
                
                assert bSet.issuperset(set(line))
                assert len(scafName)>0
                seq+=  line


        if len(scafName)>0:
            #archive last protein
            cnt['any']+=1
            seqLen=len(seq)

            if seqLen<= minNucleoLen:
                cnt['short']+=1
                f_predix.write("%s,%s,%s\n" %(scafName,"SHORTER_%s"%(minNucleoLen),1))
            elif seqLen >= maxNucleoLen:
                cnt['long']+=1
                f_predix.write("%s,%s,%s\n" %(scafName,"LONGER_%s"%(maxNucleoLen),1))
            elif  cnt['any']% nPresc!=0:
                cnt['presc']+=1
            elif  scafName not in globFD:
                cnt['noFeatures']+=1
                f_predix.write("%s,%s,%s\n" %(scafName,"NOFEATURE",1))
            else:
                gloftRaw=self.read_global_features(scafName,globFD)
                eps=np.abs(gloftRaw['len_sequence']-seqLen)
                if eps>5:
                    print(scafName,len(scafInfoList),'bad eps=%.1g'% eps,gloftRaw['len_sequence'],seqLen)
                # accpeted scaffold
                cnt['accept']+=1
                featureV=normalize_features(gloftRaw)
                sequenceLengthList.append(seqLen)
                scafInfoList.append( [scafName,text,featureV,seq])

                if cnt['accept']%1000==0:
                    print(len(scafInfoList),seqLen,scafName,featureV,seq[:100],'...')

        f_predix.close()
        fp.close()
        
        sequenceLengthListnp=np.array(sequenceLengthList)
        xavr=np.mean(sequenceLengthListnp)/1000.
        xstd=np.std(sequenceLengthListnp)/1000.
        print('  read proteins inp list len=',len(scafInfoList),', avr seq len(k)=%.1f +/- %.1f, '%(xavr,xstd),cnt)
        
        return scafInfoList


    #............................
    #Read the features
    # Called from: training and prediction (read_scaffolds_fasta).
    ##############################
    def read_global_features(self, scafName, globFD ):  #return dictionary
        ymlF=globFD[scafName]
        with open(ymlF,'r') as fd:
            print('  fd %s' % fd)
            bulk=yaml.load(fd)
        inpD=bulk['sequence']
        #print(scafName, 'bb',inpD.keys())
        outD={}
        for xN in Constants.globFeatureL:
            outD[xN]=inpD[xN]
        return outD

  
    #............................
    # Called from the formatting script. Split the sequences into folds randomly
    # Called from: training and prediction.
    # 
    # Parameteres:
    # -------------
    #      inpL : list 
    #           The items to split into segments
    #      classif : string
    #           The class of the items to be used as dict key. 
    #
    # Output:
    # ---------
    # This function writes the yaml file with the split of scaffolds
    #............................
    def split_into_folds(self, inpL, classif,save=True, excludeL=None,info=''):
        
        out_allinfo={}
        out_lengths={}
        numSegm=Constants.numSegm_total # one of the segments will be validation and the rest used in training
        
        for seg in range(numSegm):  
            out_allinfo[seg]={}
            out_lengths[seg]=0

        for scafName,text,fetureV,seq in inpL:
            randomSegment=np.random.randint(numSegm)
            out_allinfo[randomSegment][scafName]={'text':text,'seq':seq,'len':len(seq),'features':fetureV} #
            out_lengths[randomSegment]+=len(seq)

        if self.verb>0:
            print('  achieved split for  %d scaffolds to segments:'%len(inpL),[ (x,len(out_allinfo[x]),out_lengths[x]) for x in out_allinfo ])

        for seg in range(numSegm):
            sum=out_lengths[seg]/1.e6
            print('seg:',seg,' size=%.2f (MB) numScaf=%d'%(sum,len(out_allinfo[seg])))

        out_allinfo['info']=info

        # save segments
        if save==True: 
            out_file=self.trainvalid_1hot_dataPath+'/'+self.name+'.%s-scaff-split.yml'%classif
            write_yaml(out_allinfo,out_file)
        return out_allinfo



    #............................
    ''' Final formula for required number of samples per species.
    It holds both for majority and minority label.
    given: totSamp, M=numScaff, {scaffLen}_i
    sum_species[samplesPerSpe]=1.0

    Formula:
    0.2/numScaff + 0.8 * sqrt(len)/SUM[ sqrt(len)]

    Note, replace values in place in specL[.]
    returns sample_rates, which is the sample rate for each scaffold name.
    '''
    def compute_sampling_rate(self, specL):
        sqrt_lens={}
        sum_sqrt_lens=0.
        # find sum of weights

        for specN in specL:
            #print('qqq',specN, specL[specN]['text'])
            sqrtlen=pow(specL[specN]['len'],0.5)
            sqrt_lens[specN] = sqrtlen
            sum_sqrt_lens += sqrtlen
            #print('www1',specL[specN],w,sum,specN)

        # apply the formula
        sample_rates={}
        M=len(specL)
        sum2=0
        for specN in specL:
            frac=0.2/M + 0.8*sqrt_lens[specN]/ sum_sqrt_lens
            sample_rates[specN]=frac
            sum2+=frac
            #print('www2',specL[specN],frac*1000,specN)

        if self.verb>1:
            print('cross check, sum2=',sum2)

        # testing formula
        #sum3=0
        #for specN in specL:
        #    sum3+=int(1e5 * sample_rates[specN] +0.5)
        #print('sum3=',sum3,' numScaff=',len(specL))

        return sample_rates




    #............................
    # For each fold k determines which segments are train and which are val
    # dataset_trainval means val train test
    # classif means the class: main, plasm....
    # Called from: training.
    # Tasks: determine if train or validation data segment.
    # Then subsample. Save in dictionary self.trainvalid_data.
    ###############################
    def prep_labeled_input(self,classif,dataset_trainval):  #dom==domain
        # The input file as yaml
        inpF=self.trainvalid_1hot_dataPath+'/'+self.name+'.%s-scaff-split.yml'%classif
        print('\nprep_labeled scaffolds from:',inpF)
        start=time.time()

        #read it in
        inpD=read_yaml(inpF)
        numSegments=len(inpD)
        if self.verb>1:
            print(' found %d split-segments:'%numSegments,list(inpD.keys()),', kfoldOffset=', self.kfoldOffset)
            print('load_labeled_input_yaml dataset_trainval=',dataset_trainval,'kfoldOffset=', self.kfoldOffset)
        assert numSegments>0

        numTrainingSegments=Constants.numSegm_train
        numTotalSegments=Constants.numSegm_total

        assert numTrainingSegments>0  # makes no sense to split to train/eval/test
        kFold=self.kfoldOffset
        currentWrkSegment = {}

        #As example we are at kFold=1. Assuming we have 6 segments and 5 training segments:
        #This will set the train segment for this fold to 1..5 and the val to 0.
        #If we are at kFold=0, train will be 0..4 and val will be 5.
        if dataset_trainval=='train':
            jL=[ (kFold+j)%numTotalSegments for j in range(numTrainingSegments)]
            currentWrkSegment={}
            for k in jL:  currentWrkSegment.update(inpD[k])

        if dataset_trainval=='val':
            jL=(kFold+numTrainingSegments)%numTotalSegments # one element
            currentWrkSegment=inpD[jL]

        #We dont use the test segment at this point
        #if dataset_trainval=='test':
        #    jL=(kFold+numTrainingSegments+1)%self.numSegm # one element
        #    wrkL=inpD[jL]

        print('  class:',classif,'segL:',jL,' dataset_trainval:', dataset_trainval,', numSpec:',len(currentWrkSegment),', kfoldOffset=', self.kfoldOffset)

        assert len(currentWrkSegment)>0
        #print('  dump segg:',segg,wrkL)

        ############################
        ###The number of sequences to sample depends on sampling_rate (length-dependent)
        ###times the target_samples_per_contig_pred or target_samples_per_contig_train
        ############################
        sampling_rate=self.compute_sampling_rate(currentWrkSegment)

        if dataset_trainval not in self.trainvalid_data:
            self.trainvalid_data[dataset_trainval]={}

        #For validation and testing it is ok to use fewer samples, makes no difference.
        target_samples=int(Constants.target_samples_per_contig_training /2)
        if dataset_trainval!='train' : target_samples=int(target_samples/8)

        print('   request  target_samples=', target_samples, dataset_trainval,classif)
        # use all scaffolds, sample them to desired quota
        self.trainvalid_data[dataset_trainval][classif]=self.sample_labeled_scaffolds(currentWrkSegment,sampling_rate,classif,target_samples)
        print('prep:',dataset_trainval,classif,' completed, elaT=%.1f sec'%(time.time() - start),', gotSamples=',len(self.trainvalid_data[dataset_trainval][classif][0]))



    #............................
    #Called from prep_labeled_input to do subsampling.
    #Performs a subsampling of what is in currentWrkSegment.
    # Called from: training.
    # Parameters:
    # -----------
    #    currentWrkSegment : dict
    #       The current segment (0...5) to do subsampling on.
    #        
    #    sampling_rate : float 
    #       The rate to subsample
    #    
    #    classif : string
    #       main or plasmid
    # 
    #    target_samples : int
    #       How many samples in total could be subsampled.
    ##############################
    def sample_labeled_scaffolds(self, currentWrkSegment,sampling_rate,classif,target_samples):

        if target_samples < len(currentWrkSegment): # at least 1 sample per scaff on average
            print('Abort: partition_labeled_scaffold: samples:',target_samples ,', num scafold:',len(currentWrkSegment),', Hint: increase number of events\n')
            assert 1==2
        outSequences=[]; outFeatures=[] #1  outLname=[]
        j=0
        for specN in currentWrkSegment:
            j+=1
            seqStr=currentWrkSegment[specN]['seq']
            floatD=currentWrkSegment[specN]['features']
            floatL=[ floatD[x] for x in Constants.globFeatureL ]

            ############################
            ###The number of sequences to sample depends on sampling_rate (length-dependent)
            ###times the target_samples_per_contig_pred or target_samples_per_contig_train
            ############################
            numSamples=int(target_samples * sampling_rate[specN] +0.5)
            sequence_samples=self.sample_scaffold(seqStr,numSamples)

            feature_vector=[floatL]*len(sequence_samples) # the same scaffold features for all samples
            #tmpL3=[specN]*len(sequence_samples)

            outSequences[0:0]=sequence_samples
            outFeatures[0:0]=feature_vector
            #1 outLname[0:0]=tmpL3


            if j<5 or j%500==0:
                print(j,'scaff: %s %s size/B=%d numSampl=%d got=%d=%d'%(specN,classif, len(seqStr),numSamples,len(sequence_samples),len(feature_vector)))


            #exit(1)
        print('   completed class=',classif,',got numSamples=',len(outSequences))
        return [outSequences,outFeatures] #,outLname


 
        
    #...................................
    #This is called to build the training data. Called from the training script.
    # Assumes prep_labeled_input has been called.
    # dataset_trainval means val train test
    # classif means the class: main, plasm....
    #returns XhotAll,XfloatAll,YAll. It builds these from the data that is in self.trainvalid_data, which is created in prep_labeled_input
    # Called from: training.
    #
    # Tasks: call function self.list_to_nparray to convert seqs to one-hot encodings.
    #  Then sabe the 1-hot encoded data in self.trainvalid_1hot_data dictionary, for both classes main and plasm.
    #...................................
    def build_training_data(self):
        for dataset_trainval in self.trainvalid_data:
            num_main_seqs=len(self.trainvalid_data[dataset_trainval]['main'][0])
            num_plasm_seqs=len(self.trainvalid_data[dataset_trainval]['plasm'][0])
            seq_len=Constants.seqLenCut
            num_bases=len(Constants.basesSet)
            num_features=len(Constants.globFeatureL)

            # clever list-->numpy conversion, Thorsten's idea
            XhotAll=np.zeros([num_main_seqs+num_plasm_seqs,seq_len,num_bases],dtype=np.float32)
            XfloatAll=np.zeros([num_main_seqs+num_plasm_seqs,num_features],dtype=np.float32)
            #XnameAll=[]

            #The classes are a vector of 0s and 1s.
            YAll=np.zeros([num_main_seqs+num_plasm_seqs],dtype=np.float32)
            YAll[num_main_seqs:]=np.ones(num_plasm_seqs) # before randomization all plasmids are in 2nd half of array

            #print(dom,'wq0',XhotAll.shape,XfloatAll.shape)
            off={'main':0,'plasm':num_main_seqs}
            for classif in Constants.classes: 
                seqL,gloftL=self.trainvalid_data[dataset_trainval][classif]
                numSeq=len(seqL)
                print('   build_training_data:',dataset_trainval,classif,' numSeq=%d ...'%numSeq)
                self.list_to_nparray(XhotAll,XfloatAll,seqL,gloftL,offset=off[classif])
                #print('tt dom=',dom,' classif=',classif,' len=',len(seqL))

            print('build_training_data for dataset_trainval:',dataset_trainval,'Xs,Xf,Y:',XhotAll.shape,XfloatAll.shape,YAll.shape,'SNR=%.3f'%(num_plasm_seqs/num_main_seqs),'done')

            if self.trainNoise>0:
                sigNoise= self.trainNoise
                if dataset_trainval=='val' :  sigNoise=0.01
                print(' add noise ',sigNoise,' to gloft in',dataset_trainval)
                XfloatAll+=truncnorm.rvs(-2, 2, loc=0, scale=sigNoise, size=XfloatAll.shape)

            self.trainvalid_1hot_data[dataset_trainval]=[XhotAll,XfloatAll,YAll]


    #............................
    # Build the DL model. The dict trainvalid_1hot_data was created in build_training_data
    # to contain the [XhotAll,XfloatAll,YAll] for each training or validation segment.
    # Called from: training.
    #............................
    def build_compile_model(self,args):
        X_1hot_seqsamples,X_features_floats,Y=self.trainvalid_1hot_data['train']
        
        shA=X_1hot_seqsamples.shape
        shB=X_features_floats.shape
        inputA = Input(shape=(shA[1],shA[2]), name='seq_%d_x_%d'%(shA[1],shA[2]))
        inputB = Input(shape=(shB[1],), name='globF_%d'%(shB[1]))
        print('build_model inpA:',inputA.get_shape(),'  inpB:',inputB.get_shape())
        
        lstm_dim=40
        dens_dim=20
        densAct='relu'

        layerDropFrac=Constants.dropFrac
        recDropFrac=layerDropFrac/2.

        print('Dens act=',densAct,' recurDropFrac=',recDropFrac,' layerDropFrac=',layerDropFrac,' idx=',self.arrIdx)
        
        #inputA is the 1-hot encoded sequences and it goes into an LSTM.
        net= LSTM(lstm_dim, activation='tanh',recurrent_dropout=recDropFrac,dropout=layerDropFrac,name='A_%d'%lstm_dim,return_sequences=True) (inputA)
        net= LSTM(lstm_dim, activation='tanh',recurrent_dropout=recDropFrac,dropout=layerDropFrac,name='B_%d'%lstm_dim) (net)

        #inputB is the statistical features and it goes into a Dense (CNN) network.
        net2=Dense(4, activation='tanh', name='G_4')(inputB)
        net2=Dense(4, activation='tanh', name='H_4')(net2)
 
        #the two above get concatenated.
        net = concatenate([net,net2],name='seq_glob')

        net=Dense(dens_dim*2, activation=densAct, name='C_%d'%(dens_dim*2))(net)
        net=Dropout(layerDropFrac,name='fr_%.1f'%layerDropFrac)(net)
        net=Dense(dens_dim, activation=densAct, name='D_%d'%dens_dim)(net)
        net=Dropout(layerDropFrac,name='fr_same')(net)
        outputs=Dense(1, activation='sigmoid', name='score')(net) 
        model = Model(inputs=[inputA,inputB], outputs=outputs)
        self.model=model
        self.compile_model()


    #............................
    # This is called from build_compile_model
    #.........................
    #...........................
    def compile_model(self):
        """ https://machinelearningmastery.com/save-load-keras-deep-learning-models/
        It is important to compile the loaded model before it is used.
        This is so that predictions made using the model can use
        the appropriate efficient computation from the Keras backend.
        """
        print(' (re)Compile model')
        start = time.time()
        self.model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
        self.model.summary() # will print
        print('model (re)compiled elaT=%.1f sec'%(time.time() - start))



    #............................
    # This is called to train. Called from the training script.
    # Called from: training.
    #............................
    def train_model(self,args):
        X_1hot_seqsamples,X_features_floats,Y=self.trainvalid_1hot_data['train']
        X_1hot_seqsamples_val,X_features_floats_val,Y_val=self.trainvalid_1hot_data['val']

        epochs=Constants.epochs
 
        if args.verb==0:
            print('train for epochs:', Constants.epochs)
     
        callbacks_list = []
        lrCb=MyLearningTracker()
        callbacks_list.append(lrCb)

        if Constants.earlyStopPatience>0:
            earlyStop=EarlyStopping(monitor='val_loss', patience=Constants.earlyStopPatience, verbose=1, min_delta=1.e-5, mode='auto')
            callbacks_list.append(earlyStop)
            print('enabled EarlyStopping, patience=',Constants.earlyStopPatience)

        if Constants.checkPtOn:
            outF5w=self.outPath+'/'+self.name+'.weights_best.h5'
            chkPer=1
            ckpt=ModelCheckpoint(outF5w, monitor='val_loss', save_best_only=True, save_weights_only=True, verbose=1,period=chkPer)
            callbacks_list.append(ckpt)
            print('enabled ModelCheckpoint, period=',chkPer)

        if Constants.reduceLearn:
            redu_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.3, patience=4, min_lr=0.0, verbose=1,epsilon=0.003)
            callbacks_list.append(redu_lr)
            print('enabled ReduceLROnPlateau')
        
        # final check of positive/negative label balance
        sumT=sum(Y); snrT= sumT/(Y.shape[0]-sumT)
        sumV=sum(Y_val); snrV= sumV/(Y_val.shape[0]-sumV)
        print('\nTrain_model X:',X_1hot_seqsamples.shape, ' earlyStop=', Constants.earlyStopPatience,' epochs=',epochs,' batch=', Constants.batch_size,', data SNR(train)=%.4f, SNR(valid)=%.4f'%(snrT,snrV))

        fitVerb=1   # prints live:  [=========>....] - ETA: xx s 
        if args.verb==2: fitVerb=1 
        if args.verb==0: fitVerb=2  # prints 1-line summary at epoch end

        startTm = time.time()
        hir=self.model.fit([X_1hot_seqsamples,X_features_floats],Y, callbacks=callbacks_list,
                 validation_data=([X_1hot_seqsamples_val,X_features_floats_val],Y_val),  shuffle=True,
                 batch_size=Constants.batch_size, nb_epoch=epochs, 
                 verbose=fitVerb)
        fitTime=time.time() - start

        ###Save the training history!
        hir=hir.history
        for obs in hir:
            rec=[ float(x) for x in hir[obs] ]
            self.train_hirD[obs].extend(rec)
        
        # this is a hack, 'lr' is returned by fit only when --reduceLr is used
        if 'lr' not in  hir: 
            self.train_hirD['lr'].extend(lrCb.hir)

        loss=self.train_hirD['val_loss'][-1]
        acc=self.train_hirD['val_acc'][-1]
        acc0=self.train_hirD['val_acc'][0]
        nEpoch=len(self.train_hirD['val_loss'])
      
        print('\n Validation Acc:%.3f -->%.3f'%(acc0,acc), ', end-loss:%.3f'%loss,' %d epochs, idx=%d'%(nEpoch, self.arrIdx),', fit time=%.1f min'%(fitTime/60.))
        self.train_sec=fitTime
        self.acc=acc


    #............................
    # Save the model used.
    # Called from: training.
    #............................
    def save_model(self):
        outF=self.outPath+'/'+self.name+'.model.h5'
        print('save model full to',outF)
        self.model.save(outF)
        xx=os.path.getsize(outF)/1048576
        print('  closed  hdf5:',outF,' size=%.2f MB'%xx)
 
    #............................
    # Save the training history.
    # Called from: training.
    #............................
    def save_training_history(self) :
        outD=self.train_hirD
        outF=self.outPath+'/'+self.name+'.history.yml'
        write_yaml(outD,outF)

