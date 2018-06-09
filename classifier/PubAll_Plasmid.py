import numpy as np
import os, time
import datetime
from pprint import pprint

from matplotlib.colors import LogNorm
import matplotlib.patches as patches


__author__ = "Jan Balewski"
__email__ = "janstar1122@gmail.com"

from PubOne_Plasmid import PubOne_Plasmid
from Util_Plasmid import  createPubDir

#............................
#............................
#............................
class PubAll_Plasmid(object):

    def __init__(self, args,specNameL=None, maxSpec=9999): 
        self.deepName=args.deepName
        self.verb=args.verb
        self.mr12=['main','plasm']
        tmpPath='tmp1'      
        self.pubPath=createPubDir(args.webPath,tmpPath)
        self.dataPath=createPubDir(self.pubPath,'data')
        
        if specNameL==None:
            allL=os.listdir(args.dataPath)
            specNameL=[]
            for x in allL:
                if 'assayer.yml' not in x: continue
                specNameL.append(x[:-12])
            print('found %d yaml-species'%len(specNameL))
            assert len(specNameL) >0

        if maxSpec < len(specNameL): 
            specNameL=specNameL[:maxSpec]
            print('reduced initial spec list to ',len(specNameL))

        specNameL=sorted(specNameL)
        self.specPub={}
        for specN in specNameL:            
            specPubPath=createPubDir(self.pubPath,specN)
            self.specPub[specN]=PubOne_Plasmid( args,specN,specPubPath)
        print(' process %d species ...'%len(specNameL))
            
        if args.noXterm:
            print('disable Xterm')
            import matplotlib as mpl
            mpl.use('Agg')  # to plot w/o X-server
        import matplotlib.pyplot as plt
        print(self.__class__.__name__,':','Graphics started')

        self.plt=plt
        self.figL=[]
       
    #............................
    def species_loop(self,mxPrt=9999):
        k=0
        start = time.time()
        for specN in self.specPub:
            specPub=self.specPub[specN]
            specPub.doQA() #self.qaTh)
            specPub.plot_labeled_scores(10)
            #specPub.plot_boxplot(11)
            specPub.coverHTML()
            k+=1
            #print('done with %d of %d : %s\n'%(k,len(self.specPub),specN))
        print('species loop end  elaT=%.1f sec :'%(time.time() - start))


    #----------------
    def summary(self):
        self.plot_2Dscore(12)

    #----------------
    def coverHTML(self):
        dateStop=datetime.datetime.now()
        dateNowStr=dateStop.strftime("%Y-%m-%d_%H.%M")
        
        outF=self.pubPath+'/index.html'
        
        with open(outF,'w') as fd:
            fd.write("<html>\n<head></head>\n<body>\n")
            fd.write('\n<p>%s  &nbsp; &nbsp; &nbsp; generated:%s'%(self.deepName,dateNowStr))

            #common figures
            fd.write('\n<table> <tr>')
            for id,name,caption in self.figL:
                #print('HH2', id,name,caption)
                txtHtml='\n  <td> <img src="%s.png" /> <br>Fig.%d.  %s<hr>  '%( name,id,caption)
                fd.write(txtHtml)
            fd.write('\n</table> ')

            self.printSpeciesTable(fd)

            fd.write("</body> \n</html>\n")
            fd.close()

        os.chmod(outF,0o764) # a+r 

    #............................
    def plot_2Dscore(self,figId):
        
        fig=self.plt.figure(figId,facecolor='white', figsize=(5,3.5))
        ax=self.plt.subplot(1,1, 1)
        
        xV=[]; yV=[]
        for specN in self.specPub:
            data=self.specPub[specN].data
            
            xV.append(data['score']['avr'])
            yV.append(data['score']['err'])

        if len(xV) <2 :
            print(len(xV),' is too few data for plot_2Dscore, skip')
            return

        assert len(xV) >1
        hb = ax.hist2d(xV, yV, bins=20,range=[[0.,1.],[0.0,0.01]],cmin=1, cmap=self.plt.cm.coolwarm)
        
        cb = fig.colorbar(hb[3])
        cb.set_label('scaffold cnt')
        ax.set_xlabel ('score average ')
        ax.set_ylabel ('score avr err')
        ax.set_title('Model score for %d scaffolds'%len(xV))
        #ax.set_xlim(0.0,1.)
        #ax.set_ylim(0.0,0.45)

        #ax.axhline(y=1.0, linewidth=2, color='r', linestyle='-.')
        ax.text(0.3,0.008,'ambigous scaffolds')
        ax.text(0.02,0.005,'best MAIN',rotation='70')
        ax.text(0.87,0.005,'best PLASM',rotation='-70')
        ax.axvline(0.5,linewidth=1,linestyle='--')
        ax.grid()
        self.plt.tight_layout()
        coreFig='fig%d'%(figId)
        self.save_plot( coreFig )
        self.figL.append((figId,coreFig,'Average scores for %d scaffolds'%len(xV)))

  

   #............................
    def printSpeciesTable(self,fd):

        k=0
        fd.write('\n <p>  <table border="0" > <col width="30"><col width="200">  \n')
        recT=['prior designation','size KB','score avr','score err','QA','']
        colL={'good':'yellow','ambig':'silver','bad':'tomato'}


        fd.write('\n <tr> <th> idx<th>   Scaffold ')
        cnt={}

        for role in self.mr12:
            cnt[role]={'good':0,'bad':0,'ambig':0,'any':0}
        for x in recT:
            fd.write(' <th  width="40" > %s'%(x))
            
        for specN in self.specPub:
            k+=1
            fd.write('\n<tr  align="center"><td> %d <td>      <a href="%s/index.html "> %s </a>  '%(k,specN,specN))
            data=self.specPub[specN].data
            role=data['data_info']['given']

            xrole=role

            fd.write(' <td > %s '%(xrole))
            fd.write(' <td > %.0f '%(data['data_info']['size']/1e3))
            valL=[]
            valL.append(data['score']['avr'])
            valL.append(data['score']['err'])
            qa=data['score']['qa']
            cnt[role][qa]+=1;  cnt[role]['any']+=1
            col=colL[qa]

            for x in valL:
                fd.write(' <td  bgcolor="%s"> %.2f '%(col,x)) 
            fd.write(' <td bgcolor="%s">  %s <td>'%(col,qa)) 
        fd.write('\n </table> ')
        fd.write("<hr> Oracle QA stats :<pre>\n")
        pprint(cnt,fd)
        fd.write("</pre>\n")


 
    #......................
    def save_plot(self, coreFig ):
        figName='%s/%s.png'%(self.pubPath,coreFig)
        print('Saving  %s ...'%figName)
        self.plt.savefig(figName)
        self.plt.close()
        os.chmod(figName,0o764) # a+r 

