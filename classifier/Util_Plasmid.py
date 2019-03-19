# small operations w/o dependencies
# keep double indent - in case you want to make a class of it

import os, time
import shutil
import yaml

    #............................
def dump_ana_details(scaffN,sampL,floatD,Xhot,Xfloat):
        print('\n======== dump_ana_details ',scaffN,' num sampl:',len(sampL),len(floatD),' inp shapes:',Xhot.shape,Xfloat.shape)
        print('\nfloadD:',floatD[:2])

        print('\n sample of samples:\n',sampL[0])
        for j in [0,5,10,50]:
           print( 'j=',j,sampL[j][:20],Xfloat[j],Xhot[j])
        

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
def read_yaml(ymlFn):
        start = time.time()
        ymlFd = open(ymlFn, 'r')
        bulk=yaml.load( ymlFd, Loader=yaml.CLoader)
        ymlFd.close()

        print('  read  yaml:',ymlFn,' size=%d'%len(bulk),'  elaT=%.1f sec'%(time.time() - start))
        return bulk



    #............................
def createPubDir(path0,child):
        print(' createPubDir path0=',path0)
        assert(os.path.isdir(path0))
        assert(len(child)>1)
        path1=path0+'/'+child
        #print('create dir=',path1)
        if os.path.isdir(path1):
                path2=path1+'_Old'
                if os.path.isdir(path2):
                        shutil.rmtree(path2)
                os.rename(path1,path2)
        try:
            os.makedirs(path1)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        os.chmod(path1,0o765) # a+r a+x
        print(" created dir "+path1)
        return path1


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
def normalize_features(rawD):
        outD={}
        for xN in rawD:
            val=rawD[xN]
            if xN=='len_sequence':
                xx=min(1.,val/3e5)
                val=max(0.007,xx)
            elif 'longestHomopol' in xN:
                    val=min(1,val/15)
            elif 'totalLongHomopol'  in xN:
                    val=min(1,val/1000.)
            val-=0.5
            outD[xN]=val  # now val is in range [-0.5,0.5]
        return outD

