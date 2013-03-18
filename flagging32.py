import calculations
import math
import logger as log
import random
import gproutines
import globals
from threading import Thread
try :
    import preferences as p
except ImportError:
    import defaultPreferences as p
"""
Module for flagging miriad data
Part of the CARMA data reduction pipeline
Author: D. N. Friedel
"""

GOOD = "good"
BAD = "bad"
breakPoints = []

antflags = dict()

class FlagBadAmpsThread(Thread) :
    def __init__(self,object,objects,onlyBirdies=False) :
        Thread.__init__(self)
        self.object = object
        self.objects = objects
        self.onlyBirdies = onlyBirdies
    def run(self) :
        flagBadAmps(self.object,self.objects,self.onlyBirdies)

def applyFlagging(file, flagString,fileEnd) :
    """ Method to apply flagging to the given miriad data set
        input :
            file - the name of the miriad data set
            flagString - the string used for the select keyword in uvflag
        Note: the string is parsed based on '|' to separate individual
        selections
        returns :
            none
    """
    if(flagString == "" or flagString == GOOD) :
        return
    flagSplit = flagString.split("|")
    for i in flagSplit :
        args = []
        args.append(globals.Variable("vis",file,fileEnd))
        args.append(globals.Variable("select",i))
        args.append(globals.Variable("flagval","flag"))
        log.run("uvflag",args)

def flagShadowing(file) :
    """ Method to flag antenna shadowing
        input :
            file - the name of the file to flag
        returns :
            none
    """
    log.writeComment("Flagging shadowed data with the following options: OVRO shadow fraction %f, BIMA shadow fraction %f, SZA shadow fraction %f" % (p.preferences.get("OVROShadowFraction"),p.preferences.get("BIMAShadowFraction"),p.preferences.get("SZAShadowFraction")))
    args = []
    args.append(globals.Variable("vis",file))
    args.append(globals.Variable("carma","true"))
    args.append(globals.Variable("cfraction",str(p.preferences.get("OVROShadowFraction")) + "," + str(p.preferences.get("BIMAShadowFraction")) + "," + str(p.preferences.get("SZAShadowFraction"))))
    args.append(globals.Variable("sarray","0"))
    log.run("csflag",args)

# class for gains flagging
class GainsFlag :
    def __init__(self) :
        self.flags = dict()            # for the gains flags
        self.allGood = True            # are all data points good
        self.allBad = True             # are all data points bad

    def setFlag(self, time, flagval) :
        """ Method to set the flag value for a given time
            input :
                time - the time for the flag
                flagval - the flag value
            returns :
                none
        """
        self.flags[time] = flagval
        if(flagval) :
            self.allBad = False
        else :
            self.allGood = False

    def getFlags(self) :
        """ Method to return a string for the select= keyword in uvflag
            inputs :
                none
            returns :
                a string containing the select= value
        """
        global breakPoints
        if(self.allGood) :
            return GOOD
        if(self.allBad) :
            return BAD
        first = True
        start = 0.0
        end = 0.0
        flagString = ""
        times = self.flags.keys()
        times.sort()
        # go through all times in increasing order
        for time in times :
            if(not self.flags.get(time)) :
                if(first) :
                    start = time - (p.preferences.get("selfcalInterval") / 2.0)/60.0
                    first = False
                    breakPoints.append(start)
                if(time == times[-1]) :
                    end = time + (p.preferences.get("selfcalInterval") / 2.0)/60.0
                    flagString = flagString + ",time'('%s,%s')'" % (calculations.unconvertTime(start), calculations.unconvertTime(end))
                    breakPoints.append(end)
            else :
                if(start != 0.0) :
                    end = time + (p.preferences.get("selfcalInterval") / 2.0)/60.0
                    flagString = flagString + ",time'('%s,%s')'" % (calculations.unconvertTime(start), calculations.unconvertTime(end))
                    breakPoints.append(end)
                    start = 0.0
                    end = 0.0
                    first = True
        return flagString[1:]

def flagByGains(file, window) :
    """ Method to flag visibilities if the gains are outside of the range specified bu the preferences file
        input :
            file - the file to scan
            window - the window to scan
        returns :
            True/False - whether any flagging was done
    """
    # get the gains info into a text file
    global breakPoints
    anyFlagged = False
    logFile = ""
    if(window == "LSB" or window == "USB") :
        log.writeLog("%s: Flagging on gains based on thresholds given in defaultPreferences.py" % (window))
        log.run("gplist vis=%s.%s options=amp > gp.%s.log" % (file, window,window),[], logit=False)
        logFile = "gp.%s.log" % (window)
    else :
        log.writeLog("Window %i: Flagging on gains based on thresholds given in defaultPreferences.py" % (window))
        log.run("gplist vis=%s.w%i options=amp > gp.%i.log" % (file, window,window),[], logit=False)
        logFile = "gp.%i.log" % (window)
    input = open(logFile)
    fileList = input.readlines()
    input.close()
    fileList.reverse()

    gainFlags = dict()
    means = dict()
    rms = dict()
    fullMean = 0.0
    numAnts = 0
    numGoodEntries = 0
    numTimes = 0
    # read the file and record the gains
    while(len(fileList) > 0) :
        line = fileList.pop()
        if("Found" in line) :
            splitLine = line.split()
            numAnts = int(splitLine[4])
            for i in range (1, numAnts + 1) :
                gainFlags[i] = GainsFlag()
        elif("Time" in line) :
            line=fileList.pop()
            while(len(fileList) > 0 and not("---------" in line)) :
                splitLine = line.split()
                tempSplit = splitLine[0].replace(":", "")
                if(not "." in tempSplit) :
                    tempSplit = tempSplit + ".0"
                numTimes += 1
                time = calculations.convertTime(tempSplit)
                for i in range(1, len(splitLine)) :
                    if("*" in splitLine[i]) :
                        gain = 100.0
                    else :
                        gain = float(splitLine[i])
                    if(gain > 0.0) :
                        numGoodEntries += 1
                    gainFlags[i].setFlag(time, ((gain == 0.0) or (gain >= p.preferences.get("amplitudeGainRange")[0]) and (gain <= p.preferences.get("amplitudeGainRange")[1])))
                line = fileList.pop()
        elif("Means:" in line) :
            numGoodEntries /= numTimes
            splitLine = line.split()
            for ant in range(1,  numAnts + 1) :
                means[ant] = float(splitLine[ant])
                fullMean += means.get(ant)
            fullMean /= float(numGoodEntries)
        elif("Rms:" in line) :
            splitLine = line.split()
            for ant in range(1,  numAnts + 1) :
                rms[ant] = float(splitLine[ant])

    # check for bad gains and flag appropriately
    for ant in range(1, numAnts + 1) :
        select = gainFlags.get(ant).getFlags()
        if(select == BAD) :
            if(window == "LSB" or window == "USB") :
                applyFlagging(file,"antenna'('%i')'" % (ant),".%s" % (window))
            else :
                applyFlagging(file,"antenna'('%i')'" % (ant),".w%i" % (window))
            anyFlagged = True
        elif(select != GOOD) :
            if(window == "LSB" or window == "USB") :
                applyFlagging(file,"antenna'('%i')',%s" % (ant, select), ".%s" % (window))
            else :
                applyFlagging(file,"antenna'('%i')',%s" % (ant, select),".w%i" % (window))
            anyFlagged = True
            bp = ""
            for time in breakPoints :
                bp = bp + ",%s" % calculations.unconvertTime(time)
            if(bp != "") :
                args = []
                if(window == "LSB" or window == "USB") :
                    args.append(globals.Variable("vis",file,".%s" % (window)))
                else:
                    args.append(globals.Variable("vis",file,".w%i" % (window)))
                args.append(globals.Variable("break",bp[1:]))
                args.append(globals.Variable("ants",str(ant)))
                log.run("gpbreak",args)

            del breakPoints[:]
    if(anyFlagged) :
        return True
    for ant in range(1, numAnts + 1) :
        if(rms.get(ant) > p.preferences.get("maxGainRms")) :
            if(window == "LSB" or window == "USB") :
                applyFlagging(file, "antenna'('%i')'" % (ant),".%s" % (window))
            else :
                applyFlagging(file, "antenna'('%i')'" % (ant),".w%i" % (window))
            anyFlagged = True
    if(anyFlagged) :
        return True
    for ant in range(1, numAnts + 1) :
        if(means.get(ant) != 0.0 and (means.get(ant)/fullMean > p.preferences.get("maxAmplitudeGainFactor") or (means.get(ant)/fullMean < 1.0/p.preferences.get("maxAmplitudeGainFactor")))) :
            if(window == "LSB" or window == "USB") :
                applyFlagging(file, "antenna'('%i')'" % (ant),".%s" % (window))
            else :
                applyFlagging(file, "antenna'('%i')'" % (ant),".w%i" % (window))
            anyFlagged = True
    if(not anyFlagged) :
        if(window == "LSB" or window == "USB") :
            log.writeComment("%s: All gains look ok, no flaging done" % (window))
        else :
            log.writeComment("Window %i: All gains look ok, no flaging done" % (window))
    return anyFlagged

def flagByBandpass(file,fEnd,chanList) :
    """ Method to flag visibilities based on badly behaving bandpass (phase only)
        the method looks to see if the rms of the phase gains is too large and then flags the bad window(s)
        input :
            file - the file to scan
            bandpassCal - the bandpassCal object
        returns :
            True/False - whether any flagging was done
    """
    log.writeComment("Flagging on bandpass (very high phase scatter)")
    # check for overlaps first
    gain = []
    rms = []
    flagAnts = []
    anyFlagged = False
    # then read in the bandpass gains and flag as necessary
    sys = log.run("gpplt vis=%s yaxis=phase options=bandpass log=bp.log" % (file + fEnd),[],logit=False)
    if(sys != 0) :
        log.writeComment("Cannot flag based on bandpass solution")
        return anyFlagged
    input = open("bp.log")
    gainsList = input.readlines()
    input.close()
    gainsList.reverse()
    numAnts = 0
    numPoints = 0
    anyFlagged = False
    startChan = 0
    # go over each window
    tx = []
    ty = []
    while(len(gainsList) > 0) :
        line = gainsList.pop()
        if("Number of" in line) :
            splitLine = line.split()
            numAnts = int(splitLine[4])
            for j in range(0, numAnts) :
                gain.append(0.0)
                rms.append(0.0)
                ty.append([])
            continue
        elif("#" in line) :
            continue
        numPoints += 1
        splitLine = line.split()
        tx.append(float(splitLine[0]))
        splitLine = splitLine[1:]
        k = 0
        first = True
        while(k < numAnts) :
            if(not first) :
                line = gainsList.pop()
                splitLine = line.split()
            for j in splitLine :
                gain[k] += float(j)
                ty[k].append(float(j))
                rms[k] += float(j)**2
                k += 1
            first = False

        # see if we have a good solution for each antenna
    for c in range(0,len(chanList)) :
        x = None
        y = []
        if(len(chanList) == 1) :
            x = tx
            y = ty
        else :
            numPoints = chanList[c]
            x = tx[startChan:startChan+numPoints]
            for k in range(0,numAnts) :
                y.append(ty[k][startChan:startChan+numPoints])
            startChan += numPoints
        for k in range(0,numAnts) :
            if(gain[k] != 0.0) :
                y[k] = calculations.unwrap(y[k])
                # try to fit polynomials of increasing order to the
                # gains (this is necessary for windows that do not have
                # phase flattening)
                for order in range(0,4) :
                    z = calculations.fitPoly(x,y[k],order)
                    Nrms = calculations.getRms(x,y[k],z)
                    if(Nrms > 50.0 and order == 3) :
                        flagAnts.append(k+1)
                    else :
                        break
    flagAnts = list(set(flagAnts))
    # apply flagging, but only if there are not too many bad antennas
    if(len(flagAnts) > 0 and float(len(flagAnts))/float(numAnts) < 0.5) :
        for k in flagAnts :
            applyFlagging(file, "antenna'('%i')'" % (k),fEnd)
            anyFlagged = True
    return anyFlagged


""" Class for flagging bad amplitudes """
class Amplitude :
    def __init__(self,amp,rms) :
        self._amp = amp
        self._rms = rms

""" Class for flagging birdies """
class Birdie :
    def __init__(self,numChans) :
        self._highChans = []         # list of peak channels for each integration
        self._intTimes = []          # integration times
        self._numInt = 0             # number of integrations
        for i in range(0,numChans) :
            self._highChans.append(0)

    def addChans(self,chanList, time) :
        """ Method to track peak channels
            input :
                chanList - the list of peak channels from the last integration
                time - the time of the last integration
            returns :
                none
        """
        for chan in chanList:
            self._highChans[chan - 1] += 1
        if(not(time in self._intTimes)) :
            self._numInt += 1
            self._intTimes.append(time)

""" Methods to convert/unconvert the baseline numbers """
def blconvert(ant1,ant2) :
    """ Method to convert from antenna numbers to baseline number, based on MIRIAD formulation
        input :
            ant1 - number of first antenna
            ant2 - number of second antenna
        returns :
            the baseline number
    """
    return ant1*256+ant2

def blunconvert(blNumber) :
    """ Method to convert from baseline number to antenna pair, based on MIRIAD formulation
        input :
            blNumber - the baseline number
        returns :
            the antenna pair as a list [antenna1,antenna2]
    """
    ant1 = blNumber / 256
    ant2 = blNumber % 256
    return [ant1,ant2]

# ONLY FLAG ON SOURCE AND GAINCAL - they are the only ones observed long enough for this to be reasonable
#   apply their solutions to the other calibrators based on bandwidth
def flagBadAmps(object,objects,onlyBirdies=False) :
    """ Method to flag bad amplitudes
        Should only be done on source(s) and gaincal(s) as other objects will not have the s/n
        inputs :
            object - the object to work on
            objects - the full objects list
            onlyBirdies - True/False only flag the birdies
    """
    fileName = str(random.randint(1,100000))
    log.writeComment("Flagging bad amplitudes (including birdies) in %s" % (object._name))
    birdies = dict()
    numChans = 0
    LSBDone = False
    USBDone = False
    for chans in object._numChans :
        numChans += chans
    for window in range(globals.STARTWINDOW, globals.ENDWINDOW + 1) :
        vis = object._file
        fEnd = ".w%i" % (window)
        if(object.haveSuper() and object.isSuper(window)) :
            if(not LSBDone) :
                fEnd = ".LSB"
                LSBDone = True
            elif(not USBDone) :
                fEnd = ".USB"
                USBDone = True
            else :
                continue
        amps = dict()
        numPoints = 0
        avgRms = 0.0
        avgAmp = 0.0
        # get all of the visibility averages and rms values
        log.run("uvlist vis=%s options=stat log=%s recnum=0" % (vis,"uvlist." + fileName),[],logit=False)
        input = open("uvlist." + fileName)
        dataList = input.readlines()
        input.close()
        dataList.reverse()
        # accumulate the data
        while(len(dataList) > 0) :
            line = dataList.pop()
            splitLine = line.split()
            try :
               visnum = int(splitLine[0])
            except:
                pass # do nothing
            else :
                ant1 = int(splitLine[2][:-1])
                ant2 = int(splitLine[3])
                blNumber = blconvert(ant1,ant2)
                time = calculations.convertTime(splitLine[1].replace(":", ""))
                temp = None
                if(blNumber in amps) :
                    temp = amps.get(blNumber)
                else :
                    temp = dict()
                avgRms += float(splitLine[7])
                avgAmp += float(splitLine[6])
                numPoints += 1
                temp[time] = Amplitude(float(splitLine[6]),float(splitLine[7]))
                amps[blNumber] = temp
                # deal with the birdies, gather all high channels
                if(int(splitLine[9]) != 0) :
                    line = dataList.pop()
                    splitLine = line.split()
                    hiChans = []
                    for i in range(1,len(splitLine)) :
                        hiChans.append(int(splitLine[i]))
                    temp = None
                    if(blNumber in birdies) :
                        temp = birdies.get(blNumber)
                    else :
                        temp = Birdie(numChans)
                    temp.addChans(hiChans,time)
                    birdies[blNumber] = temp
        if(not onlyBirdies) :
            avgRms /= float(numPoints)
            avgAmp /= float(numPoints)
            for bl in amps :
                detectBadAmps(bl,amps.get(bl),avgAmp,avgRms,vis,fEnd)
    log.run("rm -rf uvlist." + fileName, [],logit=False)

    # do birdie detection
    for bl in birdies :
        birdieList = birdieDetect(birdies.get(bl))
        if(len(birdieList) == 0) :
            log.writeComment("No birdies found in %s" % (object._name))
            return
        for window in range(globals.STARTWINDOW, globals.ENDWINDOW + 1) :
            flagString = ""
            bw = object._bandwidths[window - 1]
            nchan = object._numChans[window - 1]
            isSuper = False
            for chan in birdieList :
                if(chan >= object._channels[window - 1] and chan < (object._channels[window - 1] + nchan)) :
                    if(object.isSuper(window)) :
                        allSuper = True
                        nonList = []
                        chanShift = 0
                        for w in range(globals.STARTWINDOW,window) :
                            allSuper = allSuper and object.isSuper(w)
                        if(not allSuper) :
                            for w in nonList :
                                chanShift += object._channels[w - 1]
                        flagstring = flagString + "|channel,1,%i,1,1" % (chan - chanShift + 1)
                    else :
                        flagString = flagString + "|channel,1,%i,1,1" % (chan - object._channels[window - 1] + 1)
            if(flagString == "") :
                continue
            for source in objects._sources :
                if(source._bandwidths[window - 1] == bw and source._numChans[window - 1] == nchan) :
                    flagBirdie(source._file,flagString,".w%i" % (window))
            for fcal in objects._fluxcals :
                if(fcal._bandwidths[window - 1] == bw and fcal._numChans[window - 1] == nchan) :
                    flagBirdie(fcal._file,flagString,".w%i" % (window))
            for pcal in objects._passcals :
                if(pcal._hybrid) :
                    if(pcal._hybridConf.get(bw)) :
                        flagBirdie(pcal._file,flagString,".%i.w%i" % (bw,window))
                else :
                    flagBirdie(pcal._file,flagString,".w%i" % (window))

def flagBirdie(file,flagString,endString) :
    """ Method to flag birides
        input :
            file - the file name
            flagString - the string containing the flagging info
        returns :
            none
    """
    if(flagString == "") :
        return
    splitLine = flagString.split("|")
    for i in splitLine :
        args = []
        args.append(globals.Variable("vis",file,endString))
        args.append(globals.Variable("line",str(i)))
        args.append(globals.Variable("flagval","flag"))
        log.run("uvflag",args)

def birdieDetect(baseline) :
    """ Method to detect birdies
        if a channel appears in both USB and LSB as the high channel 90% of the time then the channel is flagged
        input :
            baseline - the baseline number
        returns :
            list of birdies
    """
    halfChan = len(baseline._highChans) / 2
    cutoff = float(baseline._numInt) * 0.9 # birdies must appear in 90% of the data to be considered a birdie
    birdieList = []
    for chan in range(1, halfChan + 1) :
        if(baseline._highChans[chan - 1] >= cutoff and baseline._highChans[chan - 1 + halfChan] >= cutoff) :
            birdieList.append(chan)
    return birdieList

def detectBadAmps(baseline,data,avgAmp,avgRms,visFile,fEnd) :
    """ Method to detect bad amplitides by time
        any amplitude which is > 2*rms from the average amplitude gets flagged
        input :
            baseline - the baseline number
            data - the list of amplitudes
            avgAmp - average amplitude
            avgRms - average rms
            visFile - the visibility file
        returns :
            none
    """
    first = True
    start = 0.0
    end = 0.0
    flagString = ""
    times = data.keys()
    times.sort()
    for time in times :
        datum = data.get(time)
        if(abs(datum._amp - avgAmp) > 2.0 * avgRms) :
            if(first) :
                start = time - 1.0/3600.0
                first = False
            if(time == times[-1]) :
                end = time + 1.0/3600.0
                flagString = flagString + ",time'('%s,%s')'" % (calculations.unconvertTime(start),calculations.unconvertTime(end))
        else :
            if(start != 0.0) :
                end = time - 1.0/3600.0
                flagString = flagString + ",time'('%s,%s')'" % (calculations.unconvertTime(start),calculations.unconvertTime(end))
                start = 0.0
                end = 0.0
                first = True
    if(flagString != "") :
        bl = blunconvert(baseline)
        applyFlagging(visFile,"antennae'('%i,%i')',%s" % (bl[0],bl[1],flagString[1:]),fEnd)

def flagTsys(file) :
    """ Method to flag the system temperatures of uv data
        input :
            file - the file to flag
        returns :
            none
    """
    log.writeComment("Flagging high system temperatures. Cutoff is %i K" % (p.preferences.get("tsysThreshold")))
    args = []
    args.append(globals.Variable("vis",file))
    args.append(globals.Variable("tsys",str(p.preferences.get("tsysThreshold"))))
    args.append(globals.Variable("flagval","flag"))
    log.run("uvflag",args)

def flagElevation(file) :
    """ Method to flag high elevations in data
        input :
            file - the file to flag
        returns :
            none
    """
    log.writeComment("Flagging data taken at high elevations (between 85 and 90 deg)")
    args = []
    args.append(globals.Variable("vis",file))
    args.append(globals.Variable("select","'elev(85,90)'"))
    args.append(globals.Variable("flagval","flag"))
    log.run("uvflag",args)
