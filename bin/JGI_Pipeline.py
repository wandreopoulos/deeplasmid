class JGI_Pipeline:
    
    def __init__(self, cmd = '', pars = [], paramvals = [], seqin = "", headerin = ""):
        self.command = cmd
        self.params = pars
        self.paramvalues = paramvals
        self.sequence = seqin
        self.header = headerin

    def getCommand(self):
        return self.command
        