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
#............................
#............................
class Oracle_Plasmid(object):

    def __init__(self,args):
        self.name=args.prjName
        self.verb=args.verb
        print(self.__class__.__name__,', prj:',self.name)
        self.outPath=args.outPath     
        self.arrIdx=args.arrIdx
        self.padChar='N'
        self.basesL=list('ACTG')  
        self.bases2L=list('NYKMVSHRWBD')# other partial decoddings, all will map to 'N'
        self.basesSet=set(self.basesL)
        self.seqLenCut=300 # sample size in base-pairs
        self.globFeatureL=['gc_content','len_sequence']
        for x in 'ACTG' :
            self.globFeatureL.append(x+'_longestHomopol')
            self.globFeatureL.append(x+'_totalLongHomopol')

        self.mr12=['plasm','main']
        self.runDate=datetime.datetime.now()

        # - - - - - - - - - data containers - - - - - - - -
                        
        # data divided by train/val , all species chopped to NNNbp samples
        self.bases_data={} # [dom][role][ seqL,globL,nameL]

        # 1-hot encoded data ready for training
        self.data={} # ???[X[],Y[]]*[seg]
        
        # prepare 1-hot base for basess:
        myDim=len(self.basesL)
        self.oneHotBase={} # len=1+myDim, pad=00000000...
        i=0
        for x in self.basesL:
            self.oneHotBase[x]=np.eye(myDim)[i]
            i+=1
        self.oneHotBase[self.padChar]=np.zeros(myDim)
        for x in self.bases2L:
            self.oneHotBase[x]=np.zeros(myDim)
        print('oneHot base, size=',len(self.oneHotBase), ', sample:')
        for x in self.basesL:
            print('base:',x,'1-hot:',self.oneHotBase[x].tolist())
        print('all bases :', list( self.oneHotBase.keys()))
        print('use seqLenCut=',self.seqLenCut)

        #print('Cnstr ready:%s\n'%self.name)
 
   #............................
    def load_one_scaffold_fasta(self,fName):
        #print('read ',fName)
        fp = open(fName, 'r')
        out=[]
        for line in fp:
            if '>'==line[0] : continue
            #print('ll=%s='%line)
            line=line[:-1]
            setACTG=set(line)
            #print(fName,setACTG)
            assert self.basesSet.issuperset(setACTG)
            out.append(line)
        fp.close()
        return  ''.join(out)

    #............................
    def sample_scaffold(self,scaff,numSamplesF, isRandom=1):
        numSamples=int(numSamplesF)
        seqLen=self.seqLenCut
        maxChunks=len(scaff) - seqLen
        maxSamples=int(maxChunks*0.9)
        if numSamples>maxSamples: numSamples=maxSamples 
        #print('qqq',maxChunks,numSamples)
        if isRandom:
            idxL=np.random.choice(maxChunks, numSamples, replace=False)
        else:
            stepLen=int(maxChunks/numSamplesF)
            idxL=np.arange(0,maxChunks,stepLen)
        #print('sample: idxL',idxL)
        seqL=[]
        for i0 in idxL:
            chunk=scaff[i0:i0+seqLen]
            if 'NNNN' in chunk: continue  # not sure if it works - need larger input
            #print(i0,chunk, len(chunk))
            seqL.append(chunk)
        #print('sss',len(seqL))
        return seqL

    #............................
    def encode1hotBase(self,seqPad): # convert sequences to 1-hot 2-D arrays
        hot2D=[]
        for y in seqPad:
            hot2D.append(self.oneHotBase[y])
        return np.array(hot2D).astype(np.float32)

    #............................
    def add_scaff_to_1hot(self,XhotAll,XfloatAll,off,Xseq,Xfloat):
        ieve=off
        #print('wq1',XhotAll.shape,XfloatAll.shape,len(Xseq),len(Xfloat))
        assert len(Xseq)==len(Xfloat)
        for seqStr,floatV in zip(Xseq,Xfloat):
            XhotOne=self.encode1hotBase(seqStr)
            XhotAll[ieve,:]=XhotOne[:]
            XfloatAll[ieve,:]=floatV[:]
            ieve+=1
        #print(' add  shape:',len(Xseq),' ieve=',ieve)

    #............................
    def build_data_one(self,seqL,Xfloat):
        num0=len(seqL)
        seq_len=self.seqLenCut        
        num_bases=len(self.basesSet)
        num_features=len(self.globFeatureL)

        # clever list-->numpy conversion, Thorsten's idea
        XhotAll=np.zeros([num0,seq_len,num_bases],dtype=np.float32)
        XfloatAll=np.zeros([num0,num_features],dtype=np.float32)
        self.add_scaff_to_1hot(XhotAll,XfloatAll,0,seqL,Xfloat)
        #print('build_data_one:, X:',XhotAll.shape,'done')
        return XhotAll,XfloatAll

    #............................
    def classify_one_scaffold(self,Xs,Xf,verb=1,isLinear=0):
        
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
        if isLinear==0:
            return Yscore,avr_str,rec

        print('clsf: do linear classification along the scaffold ')
        for i in range(0,Yscore.shape[0]):
            print('i=',i,Yscore[i],Xf[i])
        window=50
        cnt={'any':0,'plasm':0,'main':0}
        for i in range(0,Yscore.shape[0]-window,window):
            vals=Yscore[i:i+window]
            avr=float(vals.mean())
            err=float(vals.std()/np.sqrt(window-1))
            pred=''
            cnt['any']+=1
            if avr<0.5-2*err: 
                pred='MAIN';  cnt['main']+=1
            if avr>0.5+2*err: 
                pred='PLASM' ;  cnt['plasm']+=1
            print('i=%d  score=%.2f  +/- %.2f'%(i,avr,err),pred)
        print('linear cnt',cnt)

    #............................
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
