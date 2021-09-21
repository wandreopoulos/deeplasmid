'''
The JGI_Pipeline class is used in multiprocessing to run one process per core on a node.
'''
__author__ = "Bill Andreopoulos"
__email__ = "wandreopoulos@lbl.gov"

class JGI_Pipeline:
    
    def __init__(self, cmd = '', pars = [], paramvals = [], seqin = "", headerin = ""):
        self.command = cmd
        self.params = pars
        self.paramvalues = paramvals
        self.sequence = seqin
        self.header = headerin

    def getCommand(self):
        return self.command
        
