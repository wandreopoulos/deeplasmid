from sklearn.metrics import roc_curve, auc, confusion_matrix
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib import cm as CM
import socket  # for hostname
from keras.utils import plot_model
import yaml
import Constants


__author__ = "Jan Balewski"
__email__ = "janstar1122@gmail.com"

#............................
#............................
#............................
class Plotter_Plasmid(object):
    """Graphic ops related to formatin,training, eval of deep net """

    def __init__(self, args):  #,xinch=8,yinch=6
        if args.noXterm:
            print('disable Xterm')
            import matplotlib as mpl
            mpl.use('Agg')  # to plot w/o X-server
        import matplotlib.pyplot as plt
        print(self.__class__.__name__,':','Graphics started')
        self.plt=plt
        self.nr_nc=(3,3)
        self.figL=[]

#............................
    def pause(self,args,ext,pdf=1):
        if len(self.figL)<=0: return
        self.plt.tight_layout()
        if pdf:
            for fid in self.figL:
                self.plt.figure(fid)
                self.plt.tight_layout()
                figName='%s/idx%d_%s_%s_f%d'%(args.outPath,args.arrIdx, Constants.project ,ext,fid)
                print('Graphics saving to %s PDF ...'%figName)
                self.plt.savefig(figName+'.pdf')
        self.plt.show()

#............................
    def plot_model(self,deep,flag=0):
        #print('plot_model broken in TF=1.4, skip it'); return

        if 'cori' not in socket.gethostname(): return  # software not installed
        fname=deep.outPath+'/'+deep.name+'_graph.svg'
        plot_model(deep.model, to_file=fname, show_layer_names=flag>0, show_shapes=flag>1)
        print('Graph saved as ',fname,' flag=',flag)



#............................
    def plot_scores(self,Yscore,score_thr,figId,title,xtxt='',ytxt=''):
        self.figL.append(figId)
        fig=self.plt.figure(figId,facecolor='white', figsize=(6,4))
        ax=self.plt.subplot(1,1, 1)
        
        m=0
        for y in Yscore:
            if y >score_thr: m+=1
        
        bins = np.linspace(0.0, 1., 20)

        ax.hist(Yscore, bins, alpha=0.5,label='n=%d --> %d'%(len(Yscore),m))
     
        ax.axvline(x=score_thr,linewidth=2, color='red')

        ax.set(xlabel='predicted score , avr:'+xtxt, ylabel='num samples '+ytxt,title=title)
        #ax.set_yscale('log')
        ax.grid(True)
        ax.legend(loc='upper center', title='score thr=%.2f'%score_thr)


#............................
    def plot_confusion_matrix(self,Ytrue,Yscore,score_thr,figId,title):
        self.figL.append(figId)
        fig=self.plt.figure(figId,facecolor='white', figsize=(5,4))
        ax=self.plt.subplot(1,1, 1)

        # must convert continuus score to 0/1 decision
        Ypred=Yscore > score_thr
        cm = confusion_matrix(Ytrue,Ypred)
        print('Confusion matrix for thr=%.3f\n'%score_thr,confusion_matrix(Ytrue,Ypred))
        ax.imshow(cm, interpolation='nearest', cmap=self.plt.cm.Blues)

        catL=('main DNA','mito DNA')
        tickV=np.arange(len(catL))

        ax.set_xticks( tickV+0.25)
        ax.set_xticklabels(catL)
        ax.set_yticks(tickV)
        ax.set_yticklabels(catL, rotation=90)

        ax.set( xlabel='True label', ylabel='Predicted label')
        ax.set_title(title , size=9)
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

        thresh = cm.max() / 2.
        for i in range(cm.shape[0]):
            for j in  range(cm.shape[1]):
                ax.text(j, i, '%.3f'%cm[i, j], horizontalalignment="center", size=12,color="red" if cm[i, j] > thresh else "black")
                

#............................
    def plot_float_featuresA(self,deep,recHH,roleN,figId):
        self.figL.append(figId)
        fig=self.plt.figure(figId,facecolor='white', figsize=(12,6))

        num_gloft=len(Constants.globFeatureL)#+19 #TODO change 19!
        nrow,ncol=3,4
        assert nrow*ncol >= num_gloft
        #  grid is (yN,xN) - y=0 is at the top,  so dumm
        k=0
        for gloftN in Constants.globFeatureL:
            
            k+=1
            ax=self.plt.subplot(nrow,ncol,k)
            uV=[]
            for seg in range(Constants.numSegm_total):
                rec=recHH[seg]
                for scafN in rec:
                    
                    uV.append( rec[scafN]['features'][gloftN])

            u1=min(uV); u2=max(uV)
            ax.hist(uV, 50,alpha=0.9)
            ax.set(xlabel= gloftN ,ylabel='%d scaffolds'%len(uV),title='range: %.1g, %.1g'%(u1,u2))
            if gloftN=='gc_content':                
                ax.text(0.05, 0.5,roleN,rotation='90', transform=ax.transAxes)
            #break

#............................
    def plot_float_featuresB(self,deep,dom,y_true,figId):
        self.figL.append(figId)
        fig=self.plt.figure(figId,facecolor='white', figsize=(12,6))
        print('plot_float_featuresB y_true=',y_true,dom)
        
        XA,XB,Y=deep.trainvalid_1hot_data[dom]
        num_gloft=len(Constants.globFeatureL)#+19 #TODO change 19!
        nrow,ncol=3,4
        assert nrow*ncol >= num_gloft
        #  grid is (yN,xN) - y=0 is at the top,  so dumm
        k=0
        for gloftN in Constants.globFeatureL:
            ax=self.plt.subplot(nrow,ncol,k+1)
            uV=[]
            for i in range(Y.shape[0]):
                if Y[i]!=y_true: continue
                uV.append( XB[i][k])

            u1=min(uV); u2=max(uV)
            ax.hist(uV, 50,alpha=0.9)
            ax.set(xlabel= gloftN ,ylabel='%d scaffolds'%len(uV),title='range: %.1g, %.1g'%(u1,u2))
            if k==1: ax.text(0.05, 0.8,'label=%d'%y_true, transform=ax.transAxes)
            if k==0: ax.text(0.05, 0.8,'dom=%s'%dom, transform=ax.transAxes)
            #break
            k+=1
#............................
    def plot_scaff_len(self,deep,recHH,text,figId):
        self.figL.append(figId)
        fig=self.plt.figure(figId,facecolor='white', figsize=(6,4))

        u1=[]
        for seg in range(Constants.numSegm_total):
            recH=recHH[seg]
            for scafN in recH:
                u1.append( recH[scafN]['len']/1000.)

        ax=self.plt.subplot(1,1, 1)
        ax.hist(u1, 50,alpha=0.5,label=' nScaf=%d '%(len(u1)))
        ax.set(xlabel= 'seq len (kb)' ,ylabel='scaffolds',title=text+', scaffolds :%d'%len(u1))
        ax.grid(color='brown', linestyle='--')
        #ax.set_yscale('log')
 


#............................
    def plot_train_history(self,deep,args,figId):
        self.figL.append(figId)
        self.plt.figure(figId,facecolor='white', figsize=(13,6))
        nrow,ncol=self.nr_nc
        #  grid is (yN,xN) - y=0 is at the top,  so dumm
        ax1 = self.plt.subplot2grid((nrow,ncol), (1,0), colspan=2 )

        DL=deep.train_hirD
        nEpochs=len(DL['val_loss'])
        epoV=[i for i in range(1,nEpochs+1)]
        print('sss',nEpochs,len(epoV),len(DL['loss']))
        loss=DL['val_loss'][-1]

        tit1='%s, train %.2f h, end-val-loss=%.3f'%(deep.name,deep.train_sec/3600.,loss)

        try:
            nTrain=deep.trainvalid_1hot_data['train'][2].shape[0]
            nValid=deep.trainvalid_1hot_data['val'][2].shape[0]
        except:
            nTrain=nValid=-1

        ax1.set(ylabel='loss',title=tit1)
        ax1.plot(epoV,DL['loss'],'.-.',label='train, n=%d'%nTrain)
        ax1.plot(epoV,DL['val_loss'],'.-',label='valid, n=%d'%nValid)
        ax1.legend(loc='best')
        ax1.grid(color='brown', linestyle='--',which='both')
        ax1.set_xlim(0.5,nEpochs+0.5)
        if 'acc' not in DL :
            ax1.set_yscale('log')
            return

        # this part is only for classifier
        acc=DL['val_acc'][-1]

        tit2='arrIdx=%d, end-val-acc=%.3f'%(args.arrIdx,acc)
        if 'use_encoder'  in dir(deep):
            tit2='encoder=%d, '%deep.use_encoder+tit2
        ax2 = self.plt.subplot2grid((nrow,ncol), (2,0), colspan=2,sharex=ax1)
        ax2.set(xlabel='epochs',ylabel='accuracy',title=tit2)
        ax2.plot(epoV,DL['acc'],'.-',label='train')
        ax2.plot(epoV,DL['val_acc'],'.-',label='valid')
        ax2.legend(loc='lower right')
        ax2.grid(color='brown', linestyle='--',which='both')

        if 'lr' not in DL: return
        print('sss2',len(DL['lr']))
        #if 2*nEpochs==len(DL['lr']): DL['lr']=DL['lr'][:nEpochs] #tmp hack

        ax3 = self.plt.subplot2grid((nrow,ncol), (0,0), colspan=2,sharex=ax1 )
        ax3.plot(epoV,DL['lr'],'.-',label='learn rate')
        ax3.legend(loc='best')
        ax3.grid(color='brown', linestyle='--',which='both')
        ax3.set_yscale('log')
        ax3.set(ylabel='learning rate')

#............................
    def plot_AUC(self,dom,deep,args,figId=20):         
        name=dom
        fpr_cut=0.10
        if figId==10 : 
            self.plt.figure(figId)
        else:
            self.plt.figure(figId,facecolor='white', figsize=(12,6))
            self.figL.append(figId)
        nrow,ncol=self.nr_nc
        (Xs,Xf,y_true)=deep.trainvalid_1hot_data[dom]
    
        print('Produce AUC of ROC, domain=',name,'Y shape',y_true.shape,y_true[:6])
        m=len(y_true)
        y_pred = deep.model.predict([Xs,Xf]).flatten()

        fpr, tpr, _ = roc_curve(y_true, y_pred)
        roc_auc = auc(fpr, tpr)
        outF=args.outPath+'/'+name+'.AUC.csv'
        print('AUC: %f' % roc_auc,', save:',outF)
        
        with open(outF,'w') as file:
            file.write('# FPR, TPR \n')
            for x,y in zip(fpr,tpr):
                file.write('%.5f,%.5f\n'%(x,y))
            file.close()
 
        plr=np.divide(tpr,fpr)
        for x,y,z in zip(fpr,tpr,plr):            
            if x <fpr_cut :continue
            print('found fpr=',x, ', tpr=',y,', LR+=',z)
            break

        #  grid is (yN,xN) - y=0 is at the top,  so dumm
        ax1 = self.plt.subplot2grid((nrow,ncol), (1,2), rowspan=2)
        ax1.plot(fpr, tpr, label='ROC n=%d'%len(y_pred),color='seagreen' )

        ax1.plot([0, 1], [0, 1], 'k--', label='coin flip')
        ax1.axvline(x=x,linewidth=1, color='blue')
        ax1.set(xlabel='False Positive Rate',ylabel='True Positive Rate',title='ROC , area = %.3f' % roc_auc)
        ax1.plot([0], [1], linestyle='None', marker='+', color='magenta', markersize=15,markeredgewidth=2, label='Singularity')
        ax1.legend(loc='center right', title='input set:'+name)
        ax1.grid(color='brown', linestyle='--',which='both')
        ax1.text(x*1.1,y*0.9,"%.3f @ %.3f"%(y,x),color='blue')
        #ax1.set_ylim(0.4,1.05) ;   ax1.set_xlim(-0.05,0.6)


        ax2 = self.plt.subplot(nrow,ncol, 3)
        ax2.plot(fpr,plr, label='ROC', color='teal')
        ax2.plot([0, 1], [1, 1], 'k--',label='coin flip')
        ax2.set(ylabel='Pos. Likelih. Ratio',xlabel='False Positive Rate',title='LR+(%.3f)=%.2f, %s'%(x,z,name))
        ax2.set_xlim([0.,fpr_cut+0.05])
        ax2.set_ylim([0.,3*z])

        ax2.axvline(x=x,linewidth=1, color='blue')
        ax2.legend(loc='bottom middle')
        ax2.grid(color='brown', linestyle='--',which='both')

        return y_true, y_pred

#............................
    def plot_labeled_scores(self,Ytrue,Yscore,segName,score_thr=0.5,figId=21):
        self.figL.append(figId)
        fig=self.plt.figure(figId,facecolor='white', figsize=(7,4))
        ax=self.plt.subplot(1,1, 1)
        #print('nn',Yscore.ndim)
        #assert Yscore.ndim==1
        u={0:[],1:[]}
        #mAcc={0:0,1:0}
        for ygt,ysc in  zip(Ytrue,Yscore):
            u[ygt].append(ysc)
            #if ysc> score_thr : 
            #mAcc[ygt]+=1

        mInp={ 0:len(u[0]), 1:len(u[1]) }

        print('Labeled scores found mInp',mInp, ' thr=',score_thr)

        bins = np.linspace(0.0, 1., 50)
        txt='plasmids %d' % (mInp[1])
        #txt='TPR=%.2f, '%(mAcc[1]/mInp[1])
        ax.hist(u[1], bins, alpha=0.6,label=txt) ###+'%d plasmids out of %d'%(mAcc[1],mInp[1]))
        txt='genome %d' % (mInp[0])
        #txt='FPR=%.2f, '%(mAcc[0]/mInp[0])
        ax.hist(u[0], bins, alpha=0.5,label=txt) ###+'%d main out of %d'%(mAcc[0],mInp[0]))

        ax.axvline(x=score_thr,linewidth=2, color='blue', linestyle='--')

        ax.set(xlabel='predicted score', ylabel='num samples')
        #ax.set_yscale('log')
        ax.grid(True)
        ax.set_title('Labeled scores for %s'%(segName))
        ax.legend(loc='upper right', title='score thr > %.2f'%score_thr)

