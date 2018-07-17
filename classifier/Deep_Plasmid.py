import os, time
import warnings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' #Hide messy TensorFlow warnings
warnings.filterwarnings("ignore") #Hide messy Numpy warnings
from sklearn.metrics import roc_curve, auc, confusion_matrix

start = time.time()
from Oracle_Plasmid import Oracle_Plasmid

from keras.callbacks import EarlyStopping, ModelCheckpoint ,ReduceLROnPlateau, Callback 
from keras.layers import Dense, Dropout,  LSTM, Input, concatenate
from keras.models import Model, load_model
import keras.backend as K
#import tensorflow as tf
from time import strftime

import random
import numpy as np
import yaml
from scipy.stats import truncnorm

from Util_Plasmid import write_yaml, read_yaml ,normalize_features

#print('deep-libs2 imported, TF.ver=%s, elaT=%.1f sec'%(tf.__version__,(time.time() - start)))

__author__ = "Jan Balewski"
__email__ = "janstar1122@gmail.com"

#............................
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
class Deep_Plasmid(Oracle_Plasmid):

    def __init__(self,args):
        Oracle_Plasmid.__init__(self,args)
        self.name=args.prjName
        print(self.__class__.__name__,', prj:',self.name)

        self.kfoldOffset=args.kfoldOffset
        self.events=args.events
        self.numSegm=5  # 5+1 K-fold partition of scaffolds
        self.dataPath=args.dataPath
        self.trainNoise=0.05 # added truncated gAuss to only to train-data gloft
        self.train_hirD={'acc': [],'loss': [],'lr': [],'val_acc': [],'val_loss': []}
        
        print('globFeatureL %d fixed order:'%len(self.globFeatureL),self.globFeatureL)
        
        assert os.path.isdir(self.outPath) # check it now to avoid lost time

    #............................
    def read_scaffolds_fasta(self,fName,globFD,mxRec=0, nPresc=1):
        minNucleoLen=2e3    
        maxNucleoLen= 300e3   
    
        #timestamp = strftime("%m%d%Y-%H%M%S")
        rootdir = os.getcwd()
        output_predix = os.path.join( rootdir , "predictions.txt" )
        f_predix = open(output_predix, 'a')
        f_predix.write("name,pred,conf\n")
        

        bSet=self.basesSet | set(self.bases2L)
        scafDBL=[]     
        scafN=''
        seqLenL=[] # only for accepted
        cnt={'any':0,'short':0,'long':0,'presc':0,'noFeatures':0,'accept':0}
        print('read %d from %s'%(mxRec,fName))   
        fp = open(fName, 'r')
        for line in fp:
            line = line.rstrip()
            if len(line)<1: continue  # skip empty lines
            if line[0]=='>':#  header: Name+ Species+soem info
                if mxRec and len(scafDBL)>=mxRec: break
                #print('raw',len(scafDBL),mxRec)
                if len(scafN)>0:
                    #archive last protein
                    cnt['any']+=1
                    seqLen=len(seq)

                    if seqLen<= minNucleoLen:
                        cnt['short']+=1 
                        f_predix.write("%s,%s,%s\n" %(scafN,"TOOSHORT",1))
                    elif seqLen >= maxNucleoLen:
                        cnt['long']+=1
                        f_predix.write("%s,%s,%s\n" %(scafN,"TOOLONG",1))
                    elif  cnt['any']% nPresc!=0:
                        cnt['presc']+=1
                    elif  scafN not in globFD:
                        cnt['noFeatures']+=1 
                        f_predix.write("%s,%s,%s\n" %(scafN,"NOFEATURE",1))
                    else:
                        gloftRaw=self.read_global_features(scafN,globFD) 
                        eps=np.abs(gloftRaw['len_sequence']-seqLen)
                        if eps>5:
                            print(scafN,len(scafDBL),'bad eps=%.1g'% eps,gloftRaw['len_sequence'],seqLen)                       
                        # accpeted scaffold
                        cnt['accept']+=1
                        featureV=normalize_features(gloftRaw)
                        seqLenL.append(seqLen)
                        scafDBL.append( [scafN,text,featureV,seq])  
 
                        if cnt['accept']%1000==0:
                            print(len(scafDBL),seqLen,scafN,featureV,seq[:100],'...')
                                
                      
                # init new scaffold
                seq=''
                # header is different for mito/main vs. plasmid/main
                lineL=line[1:].split()
                #print('r2',lineL)
                scafN = lineL[0]
                #assert len(lineL)>=2
                text=' '.join(lineL)
                
            else:    # - - - - it is the sequence
                line=line.upper()
                if not  bSet.issuperset(set(line)):
                        print( 'new base:',set(line) - bSet)
                
                assert bSet.issuperset(set(line))
                assert len(scafN)>0
                seq+=  line


        if len(scafN)>0:
            #archive last protein
            cnt['any']+=1
            seqLen=len(seq)

            if seqLen<= minNucleoLen:
                cnt['short']+=1
                f_predix.write("%s,%s,%s\n" %(scafN,"TOOSHORT",1))
            elif seqLen >= maxNucleoLen:
                cnt['long']+=1
                f_predix.write("%s,%s,%s\n" %(scafN,"TOOLONG",1))
            elif  cnt['any']% nPresc!=0:
                cnt['presc']+=1
            elif  scafN not in globFD:
                cnt['noFeatures']+=1
                f_predix.write("%s,%s,%s\n" %(scafN,"NOFEATURE",1))
            else:
                gloftRaw=self.read_global_features(scafN,globFD)
                eps=np.abs(gloftRaw['len_sequence']-seqLen)
                if eps>5:
                    print(scafN,len(scafDBL),'bad eps=%.1g'% eps,gloftRaw['len_sequence'],seqLen)
                # accpeted scaffold
                cnt['accept']+=1
                featureV=normalize_features(gloftRaw)
                seqLenL.append(seqLen)
                scafDBL.append( [scafN,text,featureV,seq])

                if cnt['accept']%1000==0:
                    print(len(scafDBL),seqLen,scafN,featureV,seq[:100],'...')

        f_predix.close()
        fp.close()
        
        xar=np.array(seqLenL)
        xavr=np.mean(xar)/1000.
        xstd=np.std(xar)/1000.
        print('  read proteins inp list len=',len(scafDBL),', avr seq len(k)=%.1f +/- %.1f, '%(xavr,xstd),cnt)
        
        return scafDBL


  #............................
    def read_global_features(self, scafN, globFD ):  #return dictionary
        ymlF=globFD[scafN]
        with open(ymlF,'r') as fd:
            bulk=yaml.load(fd)
        inpD=bulk['sequence']
        #print(scafN, 'bb',inpD.keys())
        outD={}
        for xN in self.globFeatureL:
            outD[xN]=inpD[xN]
        return outD

 #............................
    def split_species(self, inpL, role,save=True, excludeL=None,info=''):
        
        out0={}
        out1={}
        numSegm=self.numSegm+1 # the last one is not used in any training
        
        for seg in range(numSegm):  
            out0[seg]={}
            out1[seg]=0

        for scafN,text,fetureV,seq in inpL:
            ix=np.random.randint(numSegm)
            out0[ix][scafN]={'text':text,'seq':seq,'len':len(seq),'features':fetureV} #
            out1[ix]+=len(seq)
        if self.verb>0:
            print('  achieved split for  %d scaffolds to segments:'%len(inpL),[ (x,len(out0[x]),out1[x]) for x in out0 ])


        for seg in range(numSegm):
            sum=out1[seg]/1.e6
            print('seg:',seg,' size=%.2f (MB) numScaf=%d'%(sum,len(out0[seg])))

        out0['info']=info

        # save segments
        if save==True: 
            outF=self.dataPath+'/'+self.name+'.%s-scaff-split.yml'%role
            write_yaml(out0,outF)
        return out0

 
    #............................
    def prep_labeled_input(self,role,dom):  #dom==domain
        inpF=self.dataPath+'/'+self.name+'.%s-scaff-split.yml'%role
        print('\nprep_labeled scaffolds from:',inpF)
        start=time.time()

        inpD=read_yaml(inpF)
        nSeg=len(inpD) 
        if self.verb>1:
            print(' found %d split-segments:'%nSeg,list(inpD.keys()),', kfoldOffset=', self.kfoldOffset)
            print('load_labeled_input_yaml dom=',dom,'kfoldOffset=', self.kfoldOffset)
        assert nSeg>0

        numTrainSeg=self.numSegm-1

        assert numTrainSeg>0  # makes no sense to split to train/eval/test
        n0=self.kfoldOffset

        if dom=='val':
            jL=(n0+numTrainSeg)%self.numSegm # one element                
            wrkL=inpD[jL]

        if dom=='test':
            jL=(n0+numTrainSeg+1)%self.numSegm # one element                
            wrkL=inpD[jL]

        if dom=='train':
            jL=[ (n0+j)%self.numSegm for j in range(numTrainSeg)]
            wrkL={}            
            for k in jL:  wrkL.update(inpD[k])

        print('  role:',role,'segL:',jL,' dom:',dom,', numSpec:',len(wrkL),', kfoldOffset=', self.kfoldOffset)
            
        assert len(wrkL)>0
        #print('  dump segg:',segg,wrkL)

        dynL=self.compute_sampling_rate(wrkL)

        if dom not in self.bases_data:
            self.bases_data[dom]={}

        numSamples=int(self.events/2)
        if dom!='train' : numSamples=int(numSamples/8)

        print('   request  numSamples=',numSamples, dom,role)
        # use all scaffolds, sample them to desired quota
        self.bases_data[dom][role]=self.partition_labeled_scaffolds(wrkL,dynL,role,numSamples)
        print('prep:',dom,role,' completed, elaT=%.1f sec'%(time.time() - start),', gotSamples=',len(self.bases_data[dom][role][0]))

        
        
    #............................
    def partition_labeled_scaffolds(self, scaffL,dynL,role,totSamples): 
        
        if totSamples < len(scaffL): # at least 1 sample per scaff on average
            print('Abort: partition_labeled_scaffold: samples:',totSamples ,', num scafold:',len(scaffL),', Hint: increase number of events\n')
            assert 1==2
        outLseq=[]; outLgloft=[] #1  outLname=[]
        j=0        
        for specN in scaffL:
            j+=1            
            seqStr=scaffL[specN]['seq']
            floatD=scaffL[specN]['features']
            floatL=[ floatD[x] for x in self.globFeatureL ]

            numSamples=int(totSamples * dynL[specN] +0.5) 
            tmpL1=self.sample_scaffold(seqStr,numSamples)

            tmpL2=[floatL]*len(tmpL1) # the same scaffold features for all samples
            tmpL3=[specN]*len(tmpL1)

            outLseq[0:0]=tmpL1
            outLgloft[0:0]=tmpL2
            #1 outLname[0:0]=tmpL3


            if j<5 or j%500==0:
                print(j,'scaff: %s %s size/B=%d numSampl=%d got=%d=%d'%(specN,role, len(seqStr),numSamples,len(tmpL1),len(tmpL2)))


            #exit(1)
        print('   completed role=',role,',got numSamples=',len(outLseq))
        return [outLseq,outLgloft] #,outLname  
   

  
    #............................
    def build_training_data(self):
        for dom in self.bases_data:
            num0=len(self.bases_data[dom]['main'][0])
            num1=len(self.bases_data[dom]['plasm'][0])
            seq_len=self.seqLenCut
            num_bases=len(self.basesSet)
            num_features=len(self.globFeatureL)

            # clever list-->numpy conversion, Thorsten's idea
            XhotAll=np.zeros([num0+num1,seq_len,num_bases],dtype=np.float32)
            XfloatAll=np.zeros([num0+num1,num_features],dtype=np.float32)
            #XnameAll=[]

            YAll=np.zeros([num0+num1],dtype=np.float32)
            YAll[num0:]=np.ones(num1) # before randomization all plasmids are in 2nd half of array

            #print(dom,'wq0',XhotAll.shape,XfloatAll.shape)
            off={'main':0,'plasm':num0}
            for role in self.mr12: 
                seqL,gloftL=self.bases_data[dom][role]
                numSeq=len(seqL)
                print('   build_training_data:',dom,role,' numSeq=%d ...'%numSeq)
                self.add_scaff_to_1hot(XhotAll,XfloatAll,off[role],seqL,gloftL)
                #print('tt dom=',dom,' role=',role,' len=',len(seqL))
            print('build_training_data for dom:',dom,'Xs,Xf,Y:',XhotAll.shape,XfloatAll.shape,YAll.shape,'SNR=%.3f'%(num1/num0),'done')
            if self.trainNoise>0:
                sigNoise= self.trainNoise
                if dom=='val' :  sigNoise=0.01
                print(' add noise ',sigNoise,' to gloft in',dom)
                XfloatAll+=truncnorm.rvs(-2, 2, loc=0, scale=sigNoise, size=XfloatAll.shape)
            self.data[dom]=[XhotAll,XfloatAll,YAll]

    #............................
    def compute_sampling_rate(self, specL):
        ''' Final formula for required number of samples per species.
        It holds both for majority and minority label.
        given: totSamp, M=numScaff, {scaffLen}_i
        sum_species[samplesPerSpe]=1.0
        
        Formula:
        0.2/numScaff + 0.8 * sqrt(len)/SUM[ sqrt(len)]

        Note, replace values in place in specL[.]
        '''

        wD={}
        sum=0.
        # find sum of weights
        
        for specN in specL:
            #print('qqq',specN, specL[specN]['text'])
            w=pow(specL[specN]['len'],0.5)
            wD[specN]=w
            sum+=w
            #print('www1',specL[specN],w,sum,specN)

        # applay the formula
        dynL={}
        M=len(specL)
        sum2=0
        for specN in specL:
            frac=0.2/M + 0.8*wD[specN]/ sum
            dynL[specN]=frac
            sum2+=frac
            #print('www2',specL[specN],frac*1000,specN)
        if self.verb>1:
            print('cross check, sum2=',sum2)

        # testing formula
        sum3=0
        for specN in specL:
            sum3+=int(1e5 * dynL[specN] +0.5) 
        print('sum3=',sum3,' numScaff=',len(specL))    
            
        return dynL

    #............................
    def build_compile_model(self,args):
        XA,XB,Y=self.data['train']
        
        shA=XA.shape
        shB=XB.shape
        inputA = Input(shape=(shA[1],shA[2]), name='seq_%d_x_%d'%(shA[1],shA[2]))
        inputB = Input(shape=(shB[1],), name='globF_%d'%(shB[1]))
        print('build_model inpA:',inputA.get_shape(),'  inpB:',inputB.get_shape())
        
        lstm_dim=40
        dens_dim=20
        densAct='relu'

        layerDropFrac=args.dropFrac
        recDropFrac=layerDropFrac/2.

        print('Dens act=',densAct,' recurDropFrac=',recDropFrac,' layerDropFrac=',layerDropFrac,' idx=',args.arrIdx)
        
        net= LSTM(lstm_dim, activation='tanh',recurrent_dropout=recDropFrac,dropout=layerDropFrac,name='A_%d'%lstm_dim,return_sequences=True) (inputA)

        net= LSTM(lstm_dim, activation='tanh',recurrent_dropout=recDropFrac,dropout=layerDropFrac,name='B_%d'%lstm_dim) (net)

        net2=Dense(4, activation='tanh', name='G_4')(inputB)
        net2=Dense(4, activation='tanh', name='H_4')(net2)
 
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
    def load_compile_1model(self,path='.'):  # for re-training
        try:
            del self.model
            print('delet old model')
        except:
            a=1
        start = time.time()
        inpF5m=path+'/'+self.name+'.model.h5'
        print('load model and weights  from',inpF5m,'  ... ')
        self.model=load_model(inpF5m) # creates mode from HDF5
        #self.model.summary()
        print(' model loaded, elaT=%.1f sec'%(time.time() - start))
        self.compile_model()



    #............................
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
    def train_model(self,args):
        XA,XB,Y=self.data['train']
        XA_val,XB_val,Y_val=self.data['val']

        epochs=args.epochs
 
        if args.verb==0:
            print('train for epochs:',args.epochs)
     
        callbacks_list = []
        lrCb=MyLearningTracker()
        callbacks_list.append(lrCb)

        if args.earlyStopPatience>0:
            earlyStop=EarlyStopping(monitor='val_loss', patience=args.earlyStopPatience, verbose=1, min_delta=1.e-5, mode='auto')
            callbacks_list.append(earlyStop)
            print('enabled EarlyStopping, patience=',args.earlyStopPatience)

        if args.checkPtOn:
            outF5w=self.outPath+'/'+self.name+'.weights_best.h5'
            chkPer=1
            ckpt=ModelCheckpoint(outF5w, monitor='val_loss', save_best_only=True, save_weights_only=True, verbose=1,period=chkPer)
            callbacks_list.append(ckpt)
            print('enabled ModelCheckpoint, period=',chkPer)

        if args.reduceLearn:
            redu_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.3, patience=4, min_lr=0.0, verbose=1,epsilon=0.003)
            callbacks_list.append(redu_lr)
            print('enabled ReduceLROnPlateau')
        
        # final check of positive/negative label balance
        sumT=sum(Y); snrT= sumT/(Y.shape[0]-sumT)
        sumV=sum(Y_val); snrV= sumV/(Y_val.shape[0]-sumV)
        print('\nTrain_model X:',XA.shape, ' earlyStop=',args.earlyStopPatience,' epochs=',epochs,' batch=',args.batch_size,', data SNR(train)=%.4f, SNR(valid)=%.4f'%(snrT,snrV))

        fitVerb=1   # prints live:  [=========>....] - ETA: xx s 
        if args.verb==2: fitVerb=1 
        if args.verb==0: fitVerb=2  # prints 1-line summary at epoch end

        startTm = time.time()
        hir=self.model.fit([XA,XB],Y, callbacks=callbacks_list,
                 validation_data=([XA_val,XB_val],Y_val),  shuffle=True,
                 batch_size=args.batch_size, nb_epoch=epochs, 
                 verbose=fitVerb)
        fitTime=time.time() - start

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
      
        print('\n Validation Acc:%.3f -->%.3f'%(acc0,acc), ', end-loss:%.3f'%loss,' %d epochs, idx=%d'%(nEpoch, args.arrIdx),', fit time=%.1f min'%(fitTime/60.))
        self.train_sec=fitTime
        self.acc=acc


    #............................
    def save_model(self):
        outF=self.outPath+'/'+self.name+'.model.h5'
        print('save model full to',outF)
        self.model.save(outF)
        xx=os.path.getsize(outF)/1048576
        print('  closed  hdf5:',outF,' size=%.2f MB'%xx)
 
    #............................
    def load_weights_4_retrain(self,path='.'):
        start = time.time()
        inpF5m=path+'/'+self.name+'.weights_best.h5'
        print('load  weights  from',inpF5m,end='... ')
        self.model.load_weights(inpF5m) # creates mode from HDF5
        print('loaded, elaT=%.2f sec'%(time.time() - start))


    #............................
    def save_training_history(self) :
        outD=self.train_hirD
        outF=self.outPath+'/'+self.name+'.history.yml'
        write_yaml(outD,outF)


   #............................
    def model_predict_domain(self,dom):
        start = time.time()
        kf='' # here only one kflod is expected to be loaded
        (XA,XB,y_true)=self.data[dom]
        kf='' # here only one kflod is expected to be loaded
        print('model_predict_domain :%s   Y-size=%d'%(dom,len(y_true)))
        y_score = self.model[kf].predict([XA,XB]) 

        fpr, tpr, _ = roc_curve(y_true, y_score)
        roc_auc = float(auc(fpr, tpr))

        # this is wastefull - 2nd time computes score inside 'evaluate(..)'
        score = self.model[kf].evaluate([XA,XB],y_true, verbose=0)
        #print(' loss:', score[0])
        #print('accuracy:', score[1])
        loss= float(score[0])
        accuracy=float(score[1])

        nAll=y_true.shape[0]
        L1=int(y_true.sum())
        L0=nAll - L1

        print('  segment done, loss=%.3f, acc=%.3f, AUC=%.3f   elaT=%.1f sec'%(loss,accuracy,roc_auc,time.time() - start))

        sumD={'accuracy':accuracy,'AUC':roc_auc,'loss':loss,'nNeg':L0,'nPos':L1}

        return sumD


