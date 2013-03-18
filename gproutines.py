import calculations
import math
import random
from pipeline_miriadwrap import *

"""
Module for reading in gains from a miriad data set
behaves like gplist with options=amp
"""

class Gains:
    """
    Class to hold the gains for all antennas and a time stamp
    """
    def __init__(self, label):
        self.label = label
        self.gains = dict()

    def addGain(self,ant,gain):
        self.gains[ant] = gain

    def getFreq(self) :
        return self.label

def gplist(file):
    """Method to list the amplitude gains for a uv data set
        input :
            file - the name of the file to use
        returns :
            The gains, mean gain, median gain, and gain rms as a tuple
    """
    # open the file
    handle, iostat = hopen(file, "old")
    if(iostat != 0):
        raise Exception, "File %s not found" % (file)
    if(not hexists(handle, "gains")):
        raise Exception,  "No gains present in %s" % (file)
    ngains = rdhdi(handle, "ngains")
    nfeeds = rdhdi(handle, "nfeeds", 1)
    ntau = rdhdi(handle, "ntau")
    if(nfeeds <= 0 or nfeeds > 2 or ngains%(nfeeds+ntau)!= 0 or ntau > 1 or ntau < 0) :
        raise Exception, "Bad number of gains or feeds in %s" % (file)
    nants = ngains / (nfeeds + ntau)
    nantsd = nants
    antsd = [0] * nants
    for i in range(0,nants) :
        antsd[i] = i + 1
    # do the work
    tGains, iostat = haccess(handle,"gains","read")
    if(iostat != 0) :
        raise Exception, "Error opening gains table in %s" % (file)
    nsols = (hsize(tGains)-8)/(8*nants*(nfeeds+ntau)+8)
    ngains = nants*(nfeeds+ntau)
    offset = 8
    pnt = 0
    time = [0.0] * nsols
    gains = [complex()] * 100000
    for i in range(0,nsols) :
        temp,iostat = hreadd(tGains,offset,8)
        time[i] = temp[0]
        offset += 8
        if(iostat == 0) :
            temp,iostat = hreadr(tGains,offset,8*ngains)
        for j in range(0,ngains*2,2) :
            gains[pnt + (j/2)] = complex(temp[j],temp[j+1])
        pnt += ngains
        offset = offset + 8*ngains
        if(iostat != 0) :
            raise Exception, "Error reading gain table of %s" % (file)
    radtodeg=180.0/3.14159
    iostat = hdaccess(tGains)
    if(iostat != 0) :
        raise Exception, "Error closing gain table in %s" % (file)
    MeanGain = [0.0] * (nants)
    MednGain = [0.0] * (nants)
    GainRms = [0.0] * (nants)
    MedArr = [0.0] * (nsols)
    jant = [0] * (nants)
    GainArr = []
    for i in range(0,nants) :
        GainArr.append([0.0]* nsols)
    gainList = []
    for i in range(0,nsols) :
        calday = julday(time[i],'H').strip()
        ctime = calday[8:18]
        Gain = Gains(ctime)
        for m in range(0,nantsd) :
            Gain.addGain(antsd[m],calculations.cabs(gains[i*nants+antsd[m] - 1]))
        for j in range(0,nants) :
            if (calculations.cabs(gains[i*nants+j]) > 0.0) :
                MeanGain[j] += calculations.cabs(gains[i*nants+j])
                GainRms[j] += calculations.cabs(gains[i*nants+j]**2)
                GainArr[j][jant[j]]= calculations.cabs(gains[i*nants+j])
                jant[j] += 1
        gainList.append(Gain)
    doMed = False
    for j in range(0,nants) :
        if(jant[j] > 0) :
            MeanGain[j] = MeanGain[j]/jant[j]
        if(jant[j] > 2) :
            doMed = True
            for k in range(0,jant[j]) :
                MedArr[k] = GainArr[j][k]
            jind = calculations.sortidx(MedArr)
            MednGain[j]=MedArr[jind[int(jant[j]/2)]]
            GainRms[j]=math.sqrt((GainRms[j]- jant[j]*MeanGain[j]*MeanGain[j])/(jant[j]-1))
        else :
            GainRms[j]=0.0
    hclose(handle)
    # convert outputs to dictionary
    mean = dict()
    median = dict()
    rms = dict()
    for i in range(0,nants) :
        mean[i+1] = MeanGain[i]
        median[i+1] = MednGain[i]
        rms[i+1] = GainRms[i]
    return gainList,mean,median,rms

def gpplt(file) :
    """ python version of miriad's gpplt
        input :
            file -  the file to work on
        returns :
            a list of the gains
    """
    handle, iostat = hopen(file, "old")
    if(iostat != 0):
        raise Exception, "File %s not found" % (file)
    if(not hdprsnt(handle, "bandpass")):
        raise Exception,  "No bandpass present in %s" % (file)
    nfeeds = rdhdi(handle,"nfeeds",1)
    ngains = rdhdi(handle,"ngains",1)
    ntau = rdhdi(handle,"ntau")
    nchan = rdhdi(handle,"nchan0")
    nspect = rdhdi(handle,"nspect0")
    if(nfeeds <= 0 or ngains <= 0) :
        raise Exception, "Bad gain table size information"
    nants = ngains / (nfeeds+ntau)
    if(nants*(nfeeds+ntau) != ngains) :
        raise Exception, "Number of gains does equal nants*nfeeds"
    if(nchan <= 0) :
        raise Exception, "Bad number of frequencies"
    if(nspect <= 0 or nspect > nchan) :
        raise Exception, "Bad number of frequency spectral windows"
    item, iostat = haccess(handle,"freqs","read")
    if(iostat != 0) :
        raise Exception, "Error accessing the bandpass frequency table"
    n = 0
    off = 8
    freq = [0.0] * nspect*nchan
    select = [True] * nspect
    for i in range(0,nspect) :
        nschan,iostat = hreadi(item,off,4)
        off = off + 8
        if(iostat == 0) :
            freqs,iostat = hreadd(item,off,2*8)
            print freqs,iostat
        off = off + 2*8
        if(iostat != 0) :
            raise Exception, "Error reading bandpass frequency table"
        for j in range(1,nschan[0] + 1) :
            freq[n] = freqs[0] + (j-1)*freqs[1]
            n += 1
    iostat = hdaccess(item)
    if(iostat != 0) :
        raise Exception, "" + iostat
    item,iostat = haccess(handle,'bandpass','read')
    if(iostat != 0) :
        raise Exception, "Error accessing the bandpass table"

    off = 8
    temp,iostat = hreadr(item,off,8*nants*nfeeds*nchan)
    pnt = 0
    Gains2 = [complex(0.0,0.0)] * (nants*nfeeds*nchan)
    for j in range(0,nants*nfeeds*nchan*2,2) :
        Gains2[pnt + (j/2)] = complex(temp[j],temp[j+1])
    pnt += nants*nfeeds*nchan
    if(iostat != 0) :
        raise Exception, "Error reading the bandpass table"
    iostat = hdaccess(item)
    if(iostat != 0) :
        raise Exception, "" + iostat

    offi = 0
    for k in range(0,nants) :
        for j in range(0,nfeeds) :
            for i in range(0,nchan) :
                if(abs(Gains2[offi].real) + abs(Gains2[offi].imag) > 0.0) :
                    Gains2[offi] = 1.0/Gains2[offi]
                offi += 1
    offi -= 1
    Value = [0] * (nants*nfeeds)

    gainList = []
    y = [0.0] * (nants*nfeeds*nchan)
    for ichan in range(0,nchan) :
        offset = 0
        for iant in range(0,nants) :
            for ifeed in range(0,nfeeds) :
                Gn = Gains2[ichan+nchan*offset]
                y[offset],Value[offset] = GetPhase(Gn,Value[offset])
                offset = offset + 1
        temp = Gains(freq[ichan])
        for ant in range(0,nfeeds*nants) :
            temp.addGain(ant+1,y[ant])
        gainList.append(temp)
    return gainList

def GetPhase(G,pPhase) :
    if(abs(G.real) + abs(G.imag) == 0) :
        return 0.0,pPhase
    else :
        phase = 180/math.pi * math.atan2(G.imag,G.real)
        phase = phase - 360 * int((phase-pPhase)/360.)
        pPhase = 0.5*(phase + pPhase)
        return phase,pPhase

