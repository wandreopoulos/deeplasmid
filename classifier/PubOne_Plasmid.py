import numpy as np
import os
import shutil
import datetime
import yaml
from matplotlib.colors import LogNorm
from pprint import pprint


__author__ = "Jan Balewski"
__email__ = "janstar1122@gmail.com"


#............................
#............................
#............................
class PubOne_Plasmid(object):

    def __init__(self, args,coreName,pubPath):  
        self.name=coreName
        self.deepName=args.deepName
        self.pubPath=pubPath
        relDataPath='../data/'
        self.data={}
        self.figL=[]

        inpF=coreName+'.assayer.yml'
        inpFF=args.dataPath+'/'+inpF 
        with open(inpFF, 'r') as fd:
            rec=yaml.load( fd)
            rec.update({'yaml':relDataPath+inpF})
            print ('loaded  predictions for: %s'%(inpFF))
            outFF=self.pubPath+'/'+relDataPath+inpF
            shutil.copy(inpFF,outFF)
            os.chmod(outFF,0o764) # a+r 
            self.data=rec
            
        if args.noXterm:
            print('disable Xterm')
            import matplotlib as mpl
            mpl.use('Agg')  # to plot w/o X-server

        import matplotlib.pyplot as plt
        self.plt=plt

    #......................
    def save_plot(self, coreFig ):
        figName='%s/%s.png'%(self.pubPath,coreFig)
        print('Saving  %s ...'%figName)
        self.plt.savefig(figName)
        self.plt.close()
        os.chmod(figName,0o764) # a+r 


    #----------------
    def doQA(self): #,qaTh
        #print('do QA', self.name)
      
        data=self.data['data_info']
        role=data['given']
        rec=self.data['score']
        avr=rec['avr']
        err=rec['err'] 
            
        if role=='plasm':  #flip
            avr=1-avr
        
        if avr < 0.5- 2*err : 
            qa='good'
        elif  avr < 0.5+ 2*err : 
            qa='ambig'
        else:
            qa='bad'
    
        rec['qa']=qa

        

    #----------------
    def coverHTML(self):
        outF=self.pubPath+'/index.html'

        dateStop=datetime.datetime.now()
        dateNowStr=dateStop.strftime("%Y-%m-%d_%H.%M")
        
        #print(self.name,self.outRec)
        spac3=3*'&nbsp;'
        with open(outF,'w') as fd:
            fd.write("<html>\n<head></head>\n<body>\n")
            fd.write('\n<p>%s , produced: '%self.deepName+dateNowStr)
            
            fd.write('\n<br>species: <b> %s </b> '%self.name+spac3)

            fd.write('\n <a href="%s">  YAML</a>  '%(self.data['yaml'])+spac3)

            #common figures
            for id,name,caption in self.figL:
                #print('HH1', id,name,caption)
                txtHtml='\n  <p> <img src="%s.png" /> <br>Fig.%d.  %s<hr>  '%( name,id,caption)
                fd.write(txtHtml)


            fd.write("<p>Scaffold  analysis details<pre>\n")
            pprint(self.data,fd)
            fd.write("</pre>\n")


            fd.write("</body> \n</html>\n")
            
        os.chmod(outF,0o764) # a+r 


#............................
    def plot_labeled_scores(self,figId):
        
        fig=self.plt.figure(figId,facecolor='white', figsize=(7,4))
        ax=self.plt.subplot(1,1, 1)
        col12={'plasm':'darkorange','main':'blue'}
        bins=None
        
        nSample=self.data['ora_inp']['nSample']
        seqLen=self.data['ora_inp']['seqLen']
        rec=self.data['score']
        avr=rec['avr']
        err=rec['err']
        qa=rec['qa']
        hist=rec['hist']

        nbin=len(hist)
        if bins==None:  bins = np.linspace(0.0, 1., nbin)
        thBin=int(0.5*nbin)

        role=self.data['data_info']['given']
        if role=='main':
            sum2=sum(hist[:thBin])
            yh=10
        else:
            thBin=nbin-thBin
            sum2=sum(hist[thBin:])
            yh=20
 
        thX=bins[thBin]
        width = 1/nbin
        nh2=int(nbin/2)

        #print('plot_scores found:',sum2,thBin,thX)
            
        ax.bar(bins,hist,width, alpha=0.5,color=col12[role],label=' %d out of %s'%(sum2,nSample))
        ax.axvline(x=thX,linewidth=1,linestyle='--')
        ax.text(avr,yh+7,'%s=%s'%(role,qa))

        ax.errorbar([avr], [yh], markersize=14,linewidth=3, xerr=err, fmt='o',color=col12[role])

        ax.grid(True)
        tt1='samples above thr'
        if role=='main' :         tt1='samples below thr'
        ax.legend(loc='upper center',title=tt1)
        ax.set(xlabel='(main<--)     score      (-->mito)', ylabel='num samples, seq=%db'%seqLen, title=self.name+', prior='+role)

        self.plt.tight_layout()
        coreFig='fig%d'%(figId)
        self.save_plot( coreFig )
        self.figL.append((figId,coreFig,'Distribution of scores.'))

        return

