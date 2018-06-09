#!/usr/bin/env python
""" 
agregate kfold analysis for kinase Oracle
read yaml kfold_one summaries,
"""
import numpy as np
import os, time
import shutil
import datetime
import yaml

#from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker


__author__ = "Jan Balewski"
__email__ = "janstar1122@gmail.com"

import argparse
def get_parser():
    parser = argparse.ArgumentParser(
        description='summary of AUC & Acc from kfold method for Kinase Oracle',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument("--inpPath",
                        default='outKF',help="path to input")

    parser.add_argument("--plotPath",
                        default='plotKF',help="path to out plots")

    parser.add_argument('-X', "--no-Xterm", dest='noXterm',
                         action='store_true', default=False,
                         help="disable X-term for batch mode")

    parser.add_argument('--venue', dest='venue', choices=['web','icml'],
                         default='web',
                         help="change fromat of output plots to match publisher requirements")

    args = parser.parse_args()
    for arg in vars(args):  print( 'myArg:',arg, getattr(args, arg))
    return args


#............................
#............................
#............................
class Kfold_average(object):

    def __init__(self, args): 
        self.kfold={}
        self.venue=args.venue

        allL=os.listdir(args.inpPath)
        yamlL=[]
        for x in allL:
            if 'yml' not in x: continue
            if 'kfold' not in x: continue
            yamlL.append(x)
            #print('found %d yaml-kfolds'%len(yamlL))
            #if  len(yamlL) >5:  break
        assert len(yamlL) >0


        for x in sorted(yamlL):
            self.read_one(args.inpPath+'/'+x)
        #print(self.kfold['val_natur'])

        if args.noXterm:
            print('disable Xterm')
            import matplotlib as mpl
            mpl.use('Agg')  # to plot w/o X-server
        import matplotlib.pyplot as plt
        print(self.__class__.__name__,':','Graphics started')

        self.plt=plt
        self.figL=[]
        self.pltName=[]

    #............................
    def pause(self,args,ext,pdf=1):
        if len(self.figL)<=0: return
        self.plt.tight_layout()
        if pdf:
            for fid in self.figL:
                self.plt.figure(fid)
                self.plt.tight_layout()
                figName='%s/%s_f%d'%(args.plotPath,ext,fid)
                print('Graphics saving to %s PDF ...'%figName)
                self.plt.savefig(figName+'.pdf')
        self.plt.show()


    #............................
    def read_one(self,inpF):
        print('load kfold from',inpF)
        ymlf = open(inpF, 'r')
        inpD=yaml.load( ymlf)
        ymlf.close()
        #print('   found %d records'%len(inpD));
        assert len(inpD) >0
        for x in inpD:
            if x not in self.kfold: self.kfold[x]=[]
            self.kfold[x].append(inpD[x])

    #............................
    def plot_one(self,segName,obsN,symCol):
        figId=10+len(self.figL)
        self.figL.append(figId)
        fig=self.plt.figure(figId,facecolor='white', figsize=(4.3,2.8))
        nrow,ncol=1,4
        sym,col=symCol
        
        inp=self.kfold[segName]
        vK=[];  vY=[] 
        for rec in inp:
            k=rec['kfold_offset']
            vK.append(k)
            if obsN=='balance neg/pos':
                vY.append(rec['nNeg']/rec['nPos'])
            else:
                vY.append(rec[obsN])
        
        ax=self.plt.subplot(nrow,ncol,1)
        bp=ax.boxplot(vY, 0,'')
        #self.plt.setp(ax.get_xticklabels() , visible=False)
        #ax.set_xticks([])
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
        mm=bp['medians'][0].get_ydata()[0]

        if self.venue!='web': # correct ticks for boxplot
            tick_spacing=0.02
            if obsN=='AUC':
                tick_spacing=0.01
                #obsN='AUC test data
            ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
        else:
            ax.set_title('%s=%.3f'%(obsN,mm))
        #ax.locator_params(axis='y', nticks=4)

        #print('bp',bp,type(bp))
 
        #print('mm',mm)
        ax.set(ylabel=obsN, xlabel=segName)
        ax.yaxis.tick_right()
        ax.set_xticks([])
        

        #  grid is (yN,xN) - y=0 is at the top,  so dumm
        ax = self.plt.subplot2grid((nrow,ncol), (0,1), colspan=ncol-1, sharey=ax )
        ax.plot(vK, vY,symCol, markersize=10);
        ax.set(xlabel='kModel index')

        if self.venue!='web': # correct ticks for boxplot
            a=1
        else:
            ax.set_title('%s for segm=%s'%(obsN, segName))
            ax.grid(axis='y')

        #ax.set_yticklabels([])
        self.plt.setp(ax.get_yticklabels() , visible=False)
        #ax.set_ylim(0.95,1.0)  # modify Y-range


#=================================
#=================================
#  M A I N 
#=================================
#=================================
args=get_parser()

ppp=Kfold_average(args)
obsL={'AUC':'*b','accuracy':'8r'}
obsL={'AUC':'*b'}
#obsL={'balance neg/pos':'db','nPos':'vm','loss':'>g'}

ppp.plot_one('val','AUC','*g')
ppp.plot_one('test','AUC','*b')
ppp.plot_one('val','accuracy','8r')
ppp.plot_one('test','accuracy','8m')

ppp.pause(args,'eval')
