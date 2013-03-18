import util
"""
Global variables for the CARMA data reduction pipeline
Current Version (this file): 1.1
Released with pipeline version: 1.1
Author: D. N. Friedel
"""

oFreq = None            # observing frequency
oDate = None            # observation date decimal year form (e.g. 2008.4556)
aBaseline = None        # average baseline length in klambda
aShortBaseline = None   # average baseline for short baselines (sci2 only)
aLongBaseline = None    # average baseline for long baselines (sci2 only)
ref = None              # reference antenna
hybridMode = None       # is the data set taken in hybrid mode
antpos = None           # list of antenna positions in ns
email = None            # email of PI
header = None           # header for images

counters = {"SOURCE":0,"GAINCAL":0,"FLUXCAL":0,"PASSCAL":0,"POLCAL":0,"FLUX":0}

GAINCALEND = "gc"
PASSCALEND = "pc"
FLUXCALEND = "flx"
MINCONTINUUMBW = 400

isSci2 = False

STARTWINDOW = 1
ENDWINDOW = 1

scriptVarList = dict()

# class to hold variable data for command line execution
class Variable :
    def __init__(self,option="",arg="",postfix="",prefix="") :
        self._option = option      # left hand side of argument
        self._arg = arg            # right hand side of argument
        self._postfix = postfix    # anything that needs to be appended to the RHS
        self._prefix = prefix      # anything that needs to be prepended to the RHS

    def getOption(self) :
        """ Method to return the option
            input :
                none
            returns :
                the value of the option
        """
        return self._option

    def getArg(self) :
        """ Method to return the argument
            input :
                none
            returns :
                the value of the argument
        """
        return self._arg

    def setArg(self,arg) :
        """ Method to set the value of the argument
            input :
                arg - the value of the argument
            returns :
                none
        """
        self._arg = arg

    def getPostfix(self) :
        """ Method to return the post-fix value
            input :
                none
            returns :
                the value of the post-fix
        """
        return self._postfix

    def prependPostfix(self,pre) :
        """ Method to prepend a string to the postfix and return the result
            intput :
                pre - the value to prefix
            returns :
                the result of the concatenation
        """
        self._postfix = pre + self._postfix

    def getPrefix(self) :
        """ Method to get the prefix for the argument
            input :
                none
            returns :
                the prefix
        """
        return self._prefix

    def getAll(self) :
        """ Method to get all argument components
            input :
                none
            returns :
                a tuple containing the option, argument, post-fix and pre-fix
        """
        return self._option,self._arg,self._postfix,self._prefix

# class for scaling the amount of cpu and memory the pipeline uses
# this is currently not implemented, but is under development
class SystemInfo :
    def __init__(self) :
        self._numberOfCPUs = 0
        self._systemMemory = 0
        self._availableMemory = 0

    def setCPU(self,cpu) :
        self._numberOfCPUs = cpu

    def setMemory(self,mem,avail) :
        self._systemMemory = mem
        self._availableMemory = avail

    def getCPUs(self) :
        return self._numberOfCPUs

    def getTotalMemory(self) :
        return self._systemMemory

    def getAvailableMemory(self) :
        return self._availableMemory

    def getMemory(self) :
        return self._systemMemory,self._availableMemory

    def getSystemInfo(self) :
        self._numberOfCPUs = util.determineNumberOfCPUs()
        self._systemMemory,temp,self._availableMemory = util.getMemoryInfo()

sysinfo = SystemInfo()
bandpassBreakLow = 0.95 #MHz
bandpassBreakHi = 10.0 #MHz

def setObsInfo(freq,  date,  baseline, pos, ant, short, long) :
    """ Method to set the variables from listobs
        input :
            freq - the observing frequency in GHz
            date - the observation date
            baseline - average baseline in lambda
            ant - the reference antenna
            short - the short baseline average in lambda (sci2 only)
            long - the long baseline average in lambda (sci2 only)
    """
    global oFreq
    oFreq = freq
    global oDate
    oDate = date
    global aBaseline
    aBaseline = baseline
    global ref
    ref = ant
    global antpos
    antpos = pos
    global aShortBaseline
    aShortBaseline = short
    global aLongBaseline
    aLongBaseline = long

def setHeader(head) :
    """ Method for setting the header
        input :
            head - the header
        returns :
            none
    """
    global header
    header = head

def setHybrid(value) :
    """ Method to set the value of the hybridMode variable
        input :
            value - True/False
        returns :
            none
    """
    global hybridMode
    hybridMode = value

def setEmail(address) :
    """Method to set the value of the email
        input :
            address - the email address(es)
        returns :
            none
    """
    global email
    email = address

def obsFreq() :
    """ Method to return the observing frequency
        input :
            none
        returns :
            the value of the observing frequency
    """
    return oFreq

def obsDate() :
    """ Method to return the observation date
        input :
            none
        returns :
            the value of the observation date
    """
    return oDate

def avgBaseline():
    """ Method to return the average baseline
        input :
            none
        returns :
            the value of the average baseline
    """
    return aBaseline

def avgLongBaseline():
    """ Method to return the average long baseline
        input :
            none
        returns :
            the value of the average long baseline
    """
    return aLongBaseline

def avgShortBaseline():
    """ Method to return the average short baseline
        input :
            none
        returns :
            the value of the average short baseline
    """
    return aShortBaseline

def refant() :
    """ Method to return the reference antenna
        input :
            none
        returns :
            the value of the reference antenna
    """
    return ref

def hybrid() :
    """ Method to return the hybrid mode
        input :
            none
        returns :
            the value of the hybrid mode
    """
    return hybridMode

def antennaPositions() :
    """ Method to return the antenna positions list
        input :
            none
        returns :
            the antenna positions list
    """
    return antpos

def piEmail() :
    """ Method to return the PI email
        input :
            none
        returns :
            the value of the PI email
    """
    return email

def getHeader() :
    """ Method to return the header
        input :
            none
        returns :
            the header
    """
    return header

def setScriptVar(var,type) :
    """ Method to check and append to the global script variables
        input :
            var - the script variable value
            type - what type the var is (SOURCE, GAINCAL, ...)
        returns :
            none
    """
    global counters
    global scriptVarList
    # if the variable is already in the list then just return
    if(var in scriptVarList) :
        return
    # otherwise add it to the list and increment the index if necessary
    if(type in counters) :
        counters[type] += 1
        scriptVarList[var] = type + str(counters[type])
    else :
        scriptVarList[var] = type
