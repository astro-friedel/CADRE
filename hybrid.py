import sources
import globals
import logger as log
import random
try :
    import preferences as p
except ImportError:
    import defaultPreferences as p

"""
Module for handling data taken in a hybrid configuration (i.e. source correlator bandwidths are different from gain calibrator correlator bandwidths
Part of the CARMA data reduction pipeline
Current Version (this file): 1.1
Released with pipeline version: 1.1
Author: D. N. Friedel
"""

def splitHybrid(passbandcal, numberOfWins) :
    """ Method to split hybrid data into seperate files for each bandwidth and window
        inputs :
            passbandcal - a passband calibrator object
            numberOfWins - the total number of windows
        returns :
            none
    """
    log.writeComment("Splitting %s into separate hybrid components" % (passbandcal._name))
    array = []
    halfWin = numberOfWins/2
    if(globals.isSci2) :
        halfWin = globals.ENDWINDOW
    for i in range(globals.STARTWINDOW, halfWin + 1) :
        array.append(0)
    for window in range(1, halfWin + 1) :
        # select out the individual windows
        fileEnd = str(random.randint(1,100000))
        sys = log.run("rm -rf temp.w%i.%s; sleep 3",[],fatal=True)
        args = []
        args.append(globals.Variable("vis",passbandcal._file))
        args.append(globals.Variable("select","win'('%i')'" % (window)))
        args.append(globals.Variable("out","temp.w%i.%s" % (window,fileEnd)))
        sys = log.run("uvcat",args,fatal=True)

        if(sys != 0) :
            log.writeLog("FAILED: uvcat vis=%s select=win'('%i')' out=temp.w%i%s" % (passbandcal._file, window, window,fileEnd))
            exit(0)
        # list all spectral data for the window
        args = []
        args.append(globals.Variable("vis","temp.w%i.%s" % (window,fileEnd)))
        args.append(globals.Variable("options","spec"))
        args.append(globals.Variable("recnum","0"))
        args.append(globals.Variable("log","pcal.log"))
        log.run("uvlist",args,fatal=True)

        input = open("pcal.log",'r')
        uvList = input.readlines()
        input.close()
        uvList.reverse()
        bandwidths = []
        nchan = 0
        # gather the different bandwidth data
        while(len(uvList) > 0) :
            line = uvList.pop()
            if("number of channels" in line) :
                splitLine = line.split()
                nchan = int(splitLine[4])
            elif("frequency interval" in line) :
                splitLine = line.split()
                bw = sources.findBandwidth(float(splitLine[3])*1000.0*nchan)
                bandwidths.append(bw)
                passbandcal.setHybridChans(bw, nchan)
        first = True
        # for each bandwidth detected separate out the data for both USB and LSB windows
        for bw in bandwidths :
            if(bw != bandwidths[0]) :
                if(first) :
                    passbandcal.setHybridConf(bandwidths[0])
                    first = False
                passbandcal.setHybridConf(bw)
        for config in passbandcal._hybridConf :
            if(passbandcal._hybridConf.get(config)) :
                fileEnd = str(random.randint(1,100000))
                array[window - 1] = config
                arrayString = "%i" % (array[0])
                for m in range(1, len(array)) :
                    arrayString = arrayString + ",%i" % (array[m])
                args = []
                args.append(globals.Variable("vis",passbandcal._file))
                args.append(globals.Variable("bw",arrayString))
                args.append(globals.Variable("out","temp.pbcal." + fileEnd))
                log.run("bwsel",args)

                args = []
                args.append(globals.Variable("vis","temp.pbcal." + fileEnd))
                args.append(globals.Variable("select","win'('%i')'" % (window)))
                args.append(globals.Variable("out",passbandcal._file,".%i.w%i" % (config, window)))
                log.run("uvcat",args)

                args = []
                args.append(globals.Variable("vis","temp.pbcal." + fileEnd))
                args.append(globals.Variable("select","win'('%i')'" % (window + halfWin)))
                args.append(globals.Variable("out",passbandcal._file,".%i.w%i" % (config, window + halfWin)))
                log.run("uvcat",args)

                array[window - 1] = 0
                log.run("rm -rf temp.*",[])

def calculateOffsets(passcals, refant) :
    """ Method to calculate the phase offsets between wideband and narrow band observations
        input :
            passcals - the passband calibrator objects
            refant - the reference antenna
        returns :
            True/False if the offset was found
    """
    log.writeComment("Calculating phase offsets between windows")
    found = False
    for pcal in passcals :
        # we do not want to use the noise source
        if(pcal._type != sources.NOISE) :
            found = True
            # go through each window
            for window in range(globals.STARTWINDOW, globals.ENDWINDOW + 1) :
                # find the widest band size for the window
                largest = 0
                startChan = 1
                numChans = 1
                if(pcal._hybridConf.get(500)) :
                    largest = 500
                    numChans = pcal._numChans[window - 1]
                elif(pcal._hybridConf.get(62)) :
                    largest = 62
                    startChan = 3
                    numChans = pcal._numChans[window - 1] - 4
                elif(pcal._hybridConf.get(31)) :
                    largest = 31
                    numChans = pcal._numChans[window - 1]
                else :
                    found = False
                # selfcal the widest band to get rid of atmospheric variations
                if(found) :
                    args = []
                    args.append(globals.Variable("vis",pcal._file,".%i.w%i" % (largest,window)))
                    args.append(globals.Variable("interval",str(p.preferences.get("selfcalInterval"))))
                    args.append(globals.Variable("options","phase,noscale,apriori"))
                    args.append(globals.Variable("refant",str(refant)))
                    args.append(globals.Variable("line","line=chan,1,%i,%i" % (startChan, numChans)))
                    log.run("mselfcal",args)

                    for config in pcal._hybridConf :
                        fileEnd = str(random.randint(1,100000))
                        # copy the gains to get rid of the atmospheric variations
                        if(pcal._hybridConf.get(config) and config != largest) :
                            args = []
                            args.append(globals.Variable("vis",pcal._file,".%i.w%i" %(largest,window)))
                            args.append(globals.Variable("mode","merge"))
                            args.append(globals.Variable("out",pcal._file,".%i.w%i" %(config,window)))
                            log.run("gpcopy",args,fatal=True)
                            args = []
                            args.append(globals.Variable("in",pcal._file,".%i.w%i/interval" %(config,window)))
                            args.append(globals.Variable("value","1"))
                            log.run("puthd",args)

                            args = []
                            args.append(globals.Variable("vis",pcal._file,".%i.w%i" % (config,window)))
                            args.append(globals.Variable("out","temp." + fileEnd + "; sleep 3; rm -rf "))
                            args.append(globals.Variable(None,pcal._file,".%i.w%i/*; rm -rf " % (config,window)))
                            args.append(globals.Variable(None,pcal._file,".%i.w%i; sleep 3; mv " % (config,window)))
                            args.append(globals.Variable(None,"temp." + fileEnd))
                            args.append(globals.Variable(None,pcal._file,".%i.w%i" % (config,window)))
                            log.run("uvcat",args)

                            # calculate the phase offsets
                            args = []
                            args.append(globals.Variable("vis",pcal._file,".%i.w%i" % (largest,window)))
                            args.append(globals.Variable("interval","99999"))
                            args.append(globals.Variable("options","phase,noscale,apriori"))
                            args.append(globals.Variable("refant",str(refant)))
                            args.append(globals.Variable("line","line=chan,1,%i,%i" % (startChan, numChans)))
                            log.run("mselfcal",args)

    return found

def applyOffsets(passcals, insources) :
    """ Method to apply the band phase offsets to the source data
        input :
            passcals - the passband calibrators
            source - the source(s)
        returns :
            none
    """
    # correct the offset for each source
    for source in insources :
        log.writeComment("Applying phase offset corrections to %s" % (source._name))
        applied = [False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False]   # so we do not apply more than one offset solution to each window
        # loop over all passband cals
        for pcal in passcals :
            # ignore the noise source
            if(pcal._type != sources.NOISE and len(pcal._bandwidths) == len(source._bandwidths)) :
                # go through each window
                for window in range(globals.STARTWINDOW, globals.ENDWINDOW + 1) :
                    fileEnd = str(random.randint(1,100000))
                    # make sure only the matching bandwidth is selected
                    if(source._bandwidths[window - 1] == pcal._bandwidths[window - 1] and not(applied[window - 1])) :
                        args = []
                        args.append(globals.Variable("vis",pcal._file,".%i.w%i" %(source._bandwidths[window - 1],window)))
                        args.append(globals.Variable("mode","merge"))
                        args.append(globals.Variable("out",source._file,".%i.w%i" %(source._bandwidths[window - 1],window)))
                        log.run("gpcopy",args,fatal=True)
                        args = []
                        args.append(globals.Variable("in",source._file,".%i.w%i/interval" %(source._bandwidths[window - 1],window)))
                        args.append(globals.Variable("value","1"))
                        log.run("puthd",args)

                        args = []
                        args.append(globals.Variable("vis",source._file,".%i.w%i" % (source._bandwidths[window - 1],window)))
                        args.append(globals.Variable("out","temp." + fileEnd + "; sleep 3; rm -rf "))
                        args.append(globals.Variable(None,source._file,".%i.w%i/*; rm -rf " % (source._bandwidths[window - 1],window)))
                        args.append(globals.Variable(None,source._file,".%i.w%i; sleep 3; mv " % (source._bandwidths[window - 1],window)))
                        args.append(globals.Variable(None,"temp." + fileEnd))
                        args.append(globals.Variable(None,source._file,".%i.w%i" % (source._bandwidths[window - 1],window)))
                        log.run("uvcat",args)

                        applied[window - 1] = True
        for window in range(globals.STARTWINDOW, globals.ENDWINDOW + 1) :
            if(not applied[window - 1]) :
                log.writeComment("Could not get band offsets for %s window %i" % (source._name, window))
