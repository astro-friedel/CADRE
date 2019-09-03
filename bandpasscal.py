import sources
import logger as log
import globals
import random
import flagging
from bootflux import getCatalogFlux
from threading import Thread
import calculations


"""
Module to calculate and apply bandpass solutions
Part of the CARMA data reduction pipeline
Author: D. N. Friedel
"""

class bandpassThread(Thread) :
    def __init__(self,objects,refant,window) :
        Thread.__init__(self)
        self.objects = objects
        self.objects = objects
        self.refant = refant
        self.window = window
    def run(self) :
        bandpass(self.objects,self.refant,self.window)

def bandpass(objects,refant,window) :
    """ Method to apply bandpasses, works for both wideband and narrowband windows
        input :
            objects - the observed objects
            refant - the reference antenna
    """
    fileEnd = str(random.randint(1,100000))
    # find an appropriate passband calibrator (quasars are preferred but a planet will do if we have no choice)
    passcal = ""
    passcalFile = ""      # in case there are multiple passband calls
    passcalEnd = ""
    selectedPcal = None
    passcalFlux = 0.0
    found = False
    numChans = 1
    startChan = 1
    chanList = [-1]
    # do superwide windows
    if(objects._passcals[0].isSuper(window)) :
        tempPcal = None
        for pcal in objects._passcals :
            if(pcal._bandwidths[window - 1] == objects._sources[0]._bandwidths[window - 1] or (globals.hybrid() and pcal._hybridConf.get(objects._sources[0]._bandwidths[window - 1]))) :
                if(pcal._type == sources.QUASAR) :
                    flux = getCatalogFlux(pcal._name,calculations.gregorianToNumeric(globals.obsDate()), globals.obsFreq())
                    if(flux < passcalFlux or pcal._type == sources.PLANET) :
                        continue
                    found = True
                    tempPcal = pcal
                    if(globals.hybrid() and pcal._hybridConf.get(objects._sources[0]._bandwidths[window - 1])) :
                        passcal = pcal._file + ".%i.w%i" % (objects._sources[0]._bandwidths[window - 1],window)
                        passcalFile = pcal._file
                        passcalEnd = ".%i.w%i" %(objects._sources[0]._bandwidths[window - 1],window)
                        selectedPcal = pcal
                        numChans = pcal._hybridChans.get(objects._sources[0]._bandwidths[window - 1])
                    else :
                        chanList = []
                        passcal = pcal._file + ".w%i" %(window)
                        if(window <= len(pcal._bandwidths)/2) :
                            passcalFile = pcal._file
                            passcalEnd = ".LSB"
                            selectedPcal = pcal
                            for i in range(0,len(pcal._superwindow)/2) :
                                win = pcal._superwindow[i]
                                chanList.append(pcal.getNumChan(win))
                        else :
                            passcalFile = pcal._file
                            passcalEnd = ".USB"
                            selectedPcal = pcal
                            for i in range(len(pcal._superwindow)/2,len(pcal._superwindow)) :
                                win = pcal._superwindow[i]
                                chanList.append(pcal.getNumChan(win))
                        numChans = pcal.getSuperNumChans()
        if(not found) :
            log.writeComment("Did not find a quasar to use as a passband calibrator, now searching for a planet.")
            for pcal in objects._passcals :
                if(found) :
                    break
                if(pcal._bandwidths[window - 1] == objects._sources[0]._bandwidths[window - 1] or (globals.hybrid() and pcal._hybridConf.get(objects._sources[0]._bandwidths[window - 1]))) :
                    if(pcal._type == sources.PLANET) :
                        found = True
                        # if the data are taken in hybrid mode then make sure we get the right files
                        if(globals.hybrid() and pcal._hybridConf.get(objects._sources[0]._bandwidths[window - 1])) :
                            passcal = + pcal._file + ".%i.w%i" % (objects._sources[0]._bandwidths[window - 1],window)
                            passcalFile = pcal._file
                            passcalEnd = ".%i.w%i" %(objects._sources[0]._bandwidths[window - 1],window)
                            selectedPcal = pcal
                            numChans = pcal._hybridChans.get(objects._sources[0]._bandwidths[window - 1])
                        else :
                            chanList = []
                            passcal = pcal._file + ".w%i" %(window)
                            if(window <= len(pcal._bandwidths)/2) :
                                passcalFile = pcal._file
                                passcalEnd = ".LSB"
                                selectedPcal = pcal
                                for i in range(0,len(pcal._superwindow)/2) :
                                    win = pcal._superwindow[i]
                                    chanList.append(pcal.getNumChan(win))
                            else :
                                passcalFile = pcal._file
                                passcalEnd = ".USB"
                                selectedPcal = pcal
                                for i in range(len(pcal._superwindow)/2,len(pcal._superwindow)) :
                                    win = pcal._superwindow[i]
                                    chanList.append(pcal.getNumChan(win))
                            numChans = pcal.getSuperNumChans()
                        log.writeComment("Found a planet to use as passband.")
        if(tempPcal != None) :
            pcal = tempPcal
    # noise source is no good for 500 MHz windows
    elif(objects._sources[0].getBandwidth(window) >= 300) :
        tempPcal = None
        for pcal in objects._passcals :
            if(pcal._bandwidths[window - 1] == objects._sources[0]._bandwidths[window - 1] or (globals.hybrid() and pcal._hybridConf.get(objects._sources[0]._bandwidths[window - 1]))) :
                if(pcal._type == sources.QUASAR) :
                    flux = getCatalogFlux(pcal._name,calculations.gregorianToNumeric(globals.obsDate()), globals.obsFreq())
                    if(flux < passcalFlux or pcal._type == sources.PLANET) :
                        continue
                    found = True
                    tempPcal = pcal
                    if(globals.hybrid() and pcal._hybridConf.get(objects._sources[0]._bandwidths[window - 1])) :
                        passcal = pcal._file + ".%i.w%i" % (objects._sources[0]._bandwidths[window - 1],window)
                        passcalFile = pcal._file
                        passcalEnd = ".%i.w%i" %(objects._sources[0]._bandwidths[window - 1],window)
                        selectedPcal = pcal
                        numChans = pcal._hybridChans.get(objects._sources[0]._bandwidths[window - 1])
                    else :
                        passcal = pcal._file + ".w%i" %(window)
                        passcalFile = pcal._file
                        passcalEnd = ".w%i" %(window)
                        selectedPcal = pcal
                        numChans = pcal._numChans[window - 1]
        if(not found) :
            log.writeComment("Did not find a quasar to use as a passband calibrator, now searching for a planet.")
            for pcal in objects._passcals :
                if(found) :
                    break
                if(pcal._bandwidths[window - 1] == objects._sources[0]._bandwidths[window - 1] or (globals.hybrid() and pcal._hybridConf.get(objects._sources[0]._bandwidths[window - 1]))) :
                    if(pcal._type == sources.PLANET) :
                        found = True
                        # if the data are taken in hybrid mode then make sure we get the right files
                        if(globals.hybrid() and pcal._hybridConf.get(objects._sources[0]._bandwidths[window - 1])) :
                            passcal = pcal._file + ".%i.w%i" % (objects._sources[0]._bandwidths[window - 1],window)
                            passcalFile = pcal._file
                            passcalEnd = ".%i.w%i" %(objects._sources[0]._bandwidths[window - 1],window)
                            selectedPcal = pcal
                            numChans = pcal._hybridChans.get(objects._sources[0]._bandwidths[window - 1])
                        else :
                            passcal = pcal._file + ".w%i" %(window)
                            passcalFile = pcal._file
                            passcalEnd = ".w%i" %(window)
                            selectedPcal = pcal
                            numChans = pcal._numChans[window - 1]
                        log.writeComment("Found a planet to use as passband.")
        if(tempPcal != None) :
            pcal = tempPcal
    # astronomical source is preferred for 62 MHz, but noise source will do if need be
    elif(objects._sources[0].getBandwidth(window) >= 50) :
        tempPcal = None
        for pcal in objects._passcals :
            if(pcal._bandwidths[window - 1] == objects._sources[0]._bandwidths[window - 1] or (globals.hybrid() and pcal._hybridConf.get(objects._sources[0]._bandwidths[window - 1]))) :
                if(pcal._type == sources.QUASAR) :
                    flux = getCatalogFlux(pcal._name,calculations.gregorianToNumeric(globals.obsDate()), globals.obsFreq())
                    if(flux < passcalFlux or pcal._type == sources.PLANET) :
                        continue
                    found = True
                    tempPcal = pcal
                    startChan = 4
                    # if the data are taken in hybrid mode then make sure we get the right files
                    if(globals.hybrid() and pcal._hybridConf.get(objects._sources[0]._bandwidths[window - 1])) :
                        passcal = pcal._file + ".%i.w%i" % (objects._sources[0]._bandwidths[window - 1],window)
                        passcalFile = pcal._file
                        passcalEnd = ".%i.w%i" %(objects._sources[0]._bandwidths[window - 1],window)
                        selectedPcal = pcal
                        numChans = pcal._hybridChans.get(objects._sources[0]._bandwidths[window - 1])
                    else :
                        if(pcal.getBandwidth(window) == 62) :
                            numChans = pcal._numChans[window - 1] - 6
                        else :
                            numChans = pcal._numChans[window - 1]
                            starChan = 1
                        passcal = pcal._file + ".w%i" %(window)
                        passcalFile = pcal._file
                        passcalEnd = ".w%i" %(window)
                        selectedPcal = pcal
                        numChans = pcal._numChans[window - 1] - 6
        if(not found) :
            log.writeComment("Did not find a quasar to use as a passband calibrator, locating noise source.")
            for pcal in objects._passcals :
                if(found) :
                    break
                if(pcal._bandwidths[window - 1] == objects._sources[0]._bandwidths[window - 1] or (globals.hybrid() and pcal._hybridConf.get(objects._sources[0]._bandwidths[window - 1]))) :
                    if(pcal._type == sources.NOISE) :
                        found = True
                        # if the data are taken in hybrid mode then make sure we get the right files
                        if(globals.hybrid() and pcal._hybridConf.get(objects._sources[0]._bandwidths[window - 1])) :
                            passcal = pcal._file + ".%i.w%i" % (objects._sources[0]._bandwidths[window - 1],window)
                            passcalFile = pcal._file
                            passcalEnd = ".%i.w%i" %(objects._sources[0]._bandwidths[window - 1],window)
                            selectedPcal = pcal
                            numChans = pcal._hybridChans.get(objects._sources[0]._bandwidths[window - 1])
                        else :
                            passcal = pcal._file + ".w%i" %(window)
                            passcalFile = pcal._file
                            passcalEnd = ".w%i" %(window)
                            selectedPcal = pcal
                            numChans = pcal._numChans[window - 1] - 6
                        startChan = 4
                        log.writeComment("Found the noise source to use as passband.")
        if(not found) :
            log.writeComment("Did not find a quasar or noise source to use as a passband calibrator, now searching for a planet.")
            for pcal in objects._passcals :
                if(found) :
                    break
                if(pcal._bandwidths[window - 1] == objects._sources[0]._bandwidths[window - 1] or (globals.hybrid() and pcal._hybridConf.get(objects._sources[0]._bandwidths[window - 1]))) :
                    if(pcal._type == sources.PLANET) :
                        found = True
                        startChan = 4
                        # if the data are taken in hybrid mode then make sure we get the right files
                        if(globals.hybrid() and pcal._hybridConf.get(objects._sources[0]._bandwidths[window - 1])) :
                            passcal = pcal._file + ".%i.w%i" % (objects._sources[0]._bandwidths[window - 1],window)
                            passcalFile = pcal._file
                            passcalEnd = ".%i.w%i" %(objects._sources[0]._bandwidths[window - 1],window)
                            selectedPcal = pcal
                            numChans = pcal._hybridChans.get(objects._sources[0]._bandwidths[window - 1])
                        else :
                            if(pcal.getBandwidth(window) == 62) :
                                numChans = pcal._numChans[window - 1] - 6
                            else :
                                numChans = pcal._numChans[window - 1]
                                startChan = 1
                            passcal = pcal._file + ".w%i" %(window)
                            passcalFile = pcal._file
                            passcalEnd = ".w%i" %(window)
                            selectedPcal = pcal
                            numChans = pcal._numChans[window - 1] - 6
                        log.writeComment("Found a planet to use as passband.")
    # only the noise source is useable for the narrow windows
        if(tempPcal != None) :
            pcal = tempPcal
    else :
        for pcal in objects._passcals :
            if(found) :
                break
            if(pcal._bandwidths[window - 1] == objects._sources[0]._bandwidths[window - 1] or (globals.hybrid() and pcal._hybridConf.get(objects._sources[0]._bandwidths[window - 1]))) :
                if(pcal._type == sources.NOISE) :
                    found = True
                    # if the data are taken in hybrid mode then make sure we get the right files
                    if(globals.hybrid() and pcal._hybridConf.get(objects._sources[0]._bandwidths[window - 1])) :
                        passcal = pcal._file + ".%i.w%i" % (objects._sources[0]._bandwidths[window - 1], window)
                        passcalFile = pcal._file
                        passcalEnd = ".%i.w%i" % (objects._sources[0]._bandwidths[window - 1], window)
                        selectedPcal = pcal
                        numChans = pcal._hybridChans.get(objects._sources[0]._bandwidths[window - 1])
                    else :
                        passcal = pcal._file + ".w%i" %(window)
                        passcalFile = pcal._file
                        passcalEnd = ".w%i" %(window)
                        selectedPcal = pcal
                        numChans = pcal._numChans[window - 1]

    if(not found) :
        #cannot do passband
        log.writeComment("No appropriate sources found for passband. No passband correction done.")
        return

    log.writeComment("Using %s as passband calibrator" % (passcal))
    # calculate the passband
    log.writeComment("Calculate the bandpass.")
    anyFlag = True
    pcfile = ""
    pcend = ""
    pcal = selectedPcal
    if("passcal" in passcalFile) :
        pcfile = "passcal"
        pcend = passcalFile[7:]
    else :
        pcfile = passcalFile
        pcend = passcalEnd

    # if the channel width is too small, average more channels together to get good S/N
    avgChan = 1

    if(pcal._type != sources.NOISE and pcal.getChannelWidth(window) < 0.75) :
        temp = 0.75/pcal.getChannelWidth(window)
        temp2 = int(temp/2.0) + 1
        avgChan = int(temp2*2)
        numChans /= int(temp2*2)
        objects.updateIndividualNumChans(window,numChans * avgChan + ((startChan - 1) * 2))

    args = []
    while(anyFlag) :
        log.run("rm -rf %s/gains %s/bandpass" % (passcalFile,passcalFile),[],logit=False)
        args = []
        args.append(globals.Variable("vis",pcfile,pcend))
        args.append(globals.Variable("interval","1"))
        args.append(globals.Variable("line","channel,%i,%i,%i,%i" % (numChans,startChan,avgChan,avgChan)))
        args.append(globals.Variable("refant",str(refant)))
        args.append(globals.Variable("tol","0.001"))
        sys = log.run("mfcal",args,logit=False)

        if(sys != 0) :
            return
        # flag bad bandpass solutions
        #if(not globals.isSci2 or not objects._passcals[0].isSuper(window)) :

        anyFlag = flagging.flagByBandpass(pcfile,pcend,chanList)

    sys = log.run("mfcal",args,logit=True,execute=False)

# apply the passband to all sources, gaincals, and fluxcals
    fEnd = ""
    if(pcal.isSuper(window)) :
        if(window <= len(pcal._bandwidths)/2) :
            fEnd = ".LSB"
        else :
            fEnd = ".USB"
    else :
        fEnd = ".w%i" % (window)
    for source in objects._sources :
        log.writeAll("\n")
        log.writeComment("Copying bandpass to source file %s" % (source._file))
        args = []
        args.append(globals.Variable("vis",pcfile,pcend))
        args.append(globals.Variable("out",source._file,fEnd))
        args.append(globals.Variable("options","nocal"))
        log.run("gpcopy",args)
        args = []
        args.append(globals.Variable("vis",source._file,fEnd))
        args.append(globals.Variable("out",source._file,"." + fileEnd + "; sleep 3; rm -r"))
        args.append(globals.Variable(None,source._file,"%s/*; rm -r " % (fEnd)))
        args.append(globals.Variable(None,source._file,"%s; sleep 3; mv" % (fEnd)))
        args.append(globals.Variable(None,source._file,"." + fileEnd))
        args.append(globals.Variable(None,source._file,fEnd))
        log.run("uvcat",args)
        log.writeAll("\n")

    for gcal in objects._gaincals :
        if(gcal._file != passcal) :
            log.writeAll("\n")
            log.writeComment("Copying bandpass to gaincal file %s" % (gcal._file))
            args = []
            args.append(globals.Variable("vis",pcfile,pcend))
            args.append(globals.Variable("out",gcal._file,fEnd))
            args.append(globals.Variable("options","nocal"))
            log.run("gpcopy",args)
            args = []
            args.append(globals.Variable("vis",gcal._file,fEnd))
            args.append(globals.Variable("out",gcal._file,"." + fileEnd + "; sleep 3; rm -r"))
            args.append(globals.Variable(None,gcal._file,"%s/*; rm -r " % (fEnd)))
            args.append(globals.Variable(None,gcal._file,"%s; sleep 3; mv" % (fEnd)))
            args.append(globals.Variable(None,gcal._file,"." + fileEnd))
            args.append(globals.Variable(None,gcal._file,fEnd))
            log.run("uvcat",args)
            log.writeAll("\n")

    for fcal in objects._fluxcals :
        if(fcal._file != passcal) :
            log.writeAll("\n")
            log.writeComment("Copying bandpass to flux file %s" % (fcal._file))
            args = []
            args.append(globals.Variable("vis",pcfile,pcend))
            args.append(globals.Variable("out",fcal._file,fEnd))
            args.append(globals.Variable("options","nocal"))
            log.run("gpcopy",args)
            args = []
            args.append(globals.Variable("vis",fcal._file,fEnd))
            args.append(globals.Variable("out",fcal._file,"." + fileEnd + "; sleep 3; rm -r"))
            args.append(globals.Variable(None,fcal._file,"%s/*; rm -r " % (fEnd)))
            args.append(globals.Variable(None,fcal._file,"%s; sleep 3; mv" % (fEnd)))
            args.append(globals.Variable(None,fcal._file,"." + fileEnd))
            args.append(globals.Variable(None,fcal._file,fEnd))
            log.run("uvcat",args)
            log.writeAll("\n")
