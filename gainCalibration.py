import logger as log
import flagging
import calculations
from threading import Thread
import random
import globals
try :
    import preferences as p
except ImportError:
    import defaultPreferences as p

"""
Module to do gain calibration on both wide band and narrow band windows
Part of the CARMA data reduction pipeline
Author: D. N. Friedel
"""

class continuumThread(Thread) :
    def __init__(self,objects,refant,window) :
        Thread.__init__(self)
        self.objects = objects
        self.refant = refant
        self.window = window
    def run(self) :
        continuum(self.objects, self.refant, self.window)

class bootstrapThread(Thread) :
    def __init__(self,objects,refant,refWindow,window) :
        Thread.__init__(self)
        self.objects = objects
        self.refant = refant
        self.refWindow = refWindow
        self.window = window
    def run(self) :
        bootstrap(self.objects, self.refant, self.refWindow, self.window)

class bootcalThread(Thread) :
    def __init__(self,objects,gcal,refant,refWindow,window) :
        Thread.__init__(self)
        self.gcal = gcal
        self.refant = refant
        self.refWindow = refWindow
        self.window = window
        self.objects = objects
    def run(self) :
        runboot(self.objects,self.gcal,self.refant,self.refWindow,self.window)

def runboot(objects,gcal,refant,refWindow,window) :
    """ Method for bootstrapping the gain solution
        input :
            objects - the objects
            gcal - the gain calibrator
            refant - the reference antenna
            refWindow - the reference window
            window - the window being bootstrapped
        returns :
            none
    """
    fileEnd = str(random.randint(1,100000))
    if(gcal.isSuper(window)) :
        sideband = "LSB"
        if(window > len(gcal._bandwidths)/2) :
            sideband = "USB"

        log.writeComment("Calibrating window %i of %s by bootstrapping from window %s" % (window, gcal._name,sideband))
        startChan = 1
        numChans = gcal.getSuperNumChans()
        # copy gains from wideband window
        args = []
        args.append(globals.Variable("vis",gcal._file,".%s" % (sideband)))
        args.append(globals.Variable("mode","merge"))
        args.append(globals.Variable("out",gcal._file,".w%i" %(window)))
        log.run("gpcopy",args,fatal=True)
        args = []
        args.append(globals.Variable("in",gcal._file,".w%i/interval" %(window)))
        args.append(globals.Variable("value","1"))
        log.run("puthd",args)

        args = []
        args.append(globals.Variable("vis",gcal._file,".w%i" % (window)))
        args.append(globals.Variable("out",gcal._file,"." + fileEnd + "; sleep 3; rm -rf "))
        args.append(globals.Variable(None,gcal._file,".w%i/*; rm -rf " % (window)))
        args.append(globals.Variable(None,gcal._file,".w%i; sleep 3; mv " % (window)))
        args.append(globals.Variable(None,gcal._file,"." + fileEnd))
        args.append(globals.Variable(None,gcal._file,".w%i" % (window)))
        log.run("uvcat",args)

    # take care of phase offset
        args = []
        args.append(globals.Variable("vis",gcal._file,".w%i" % (window)))
        args.append(globals.Variable("interval","99999"))
        args.append(globals.Variable("options","amplitude,noscale,apriori"))
        args.append(globals.Variable("refant",str(refant)))
        args.append(globals.Variable("line","chan,1,%i,%i" % (startChan, numChans)))
        log.run("mselfcal",args)

    # do the same for each source
        for source in objects._sources :
            args = []
            args.append(globals.Variable("vis",gcal._file,".%s" %(sideband)))
            args.append(globals.Variable("mode","merge"))
            args.append(globals.Variable("out",source._file,".w%i" %(window)))
            log.run("gpcopy",args,fatal=True)
            args = []
            args.append(globals.Variable("in",source._file,".w%i/interval" %(window)))
            args.append(globals.Variable("value","1"))
            log.run("puthd",args)

            args = []
            args.append(globals.Variable("vis",source._file,".w%i" % (window)))
            args.append(globals.Variable("out",source._file,"." + fileEnd + "; sleep 3; rm -rf "))
            args.append(globals.Variable(None,source._file,".w%i/*; rm -rf " % (window)))
            args.append(globals.Variable(None,source._file,".w%i; sleep 3; mv " % (window)))
            args.append(globals.Variable(None,source._file,"." + fileEnd))
            args.append(globals.Variable(None,source._file,".w%i" % (window)))
            log.run("uvcat",args)

            args = []
            args.append(globals.Variable("vis",gcal._file,".w%i" %(window)))
            args.append(globals.Variable("mode","merge"))
            args.append(globals.Variable("out",source._file,".w%i" %(window)))
            log.run("gpcopy",args,fatal=True)
            args = []
            args.append(globals.Variable("in",source._file,".w%i/interval" %(window)))
            args.append(globals.Variable("value","1"))
            log.run("puthd",args)
    else :
        if("LSB" in refWindow or "USB" in refWindow) :
            log.writeComment("Calibrating window %i of %s by bootstrapping from window %s" % (window, gcal._name,refWindow))
        else :
            log.writeComment("Calibrating window %i of %s by bootstrapping from window %i" % (window, gcal._name,refWindow))
        if(gcal._bandwidths[window - 1] == 62) :
            startChan = 4
            numChans = gcal._numChans[window - 1] - 6
        else:
            startChan = 1
            numChans = gcal._numChans[window - 1]
        # copy gains from wideband window
        args = []
        if("LSB" in refWindow or "USB" in refWindow) :
            args.append(globals.Variable("vis",gcal._file,".%s" %(refWindow)))
        else :
            args.append(globals.Variable("vis",gcal._file,".w%i" %(refWindow)))
        args.append(globals.Variable("mode","merge"))
        args.append(globals.Variable("out",gcal._file,".w%i" %(window)))
        log.run("gpcopy",args,fatal=True)
        args = []
        args.append(globals.Variable("in",gcal._file,".w%i/interval" %(window)))
        args.append(globals.Variable("value","1"))
        log.run("puthd",args)

        args = []
        args.append(globals.Variable("vis",gcal._file,".w%i" % (window)))
        args.append(globals.Variable("out",gcal._file,"." + fileEnd + "; sleep 3; rm -rf "))
        args.append(globals.Variable(None,gcal._file,".w%i/*; rm -rf " % (window)))
        args.append(globals.Variable(None,gcal._file,".w%i; sleep 3; mv " % (window)))
        args.append(globals.Variable(None,gcal._file,"." + fileEnd))
        args.append(globals.Variable(None,gcal._file,".w%i" % (window)))
        log.run("uvcat",args)

    # take care of phase offset
        args = []
        args.append(globals.Variable("vis",gcal._file,".w%i" % (window)))
        args.append(globals.Variable("interval","99999"))
        args.append(globals.Variable("options","amplitude,noscale,apriori"))
        args.append(globals.Variable("refant",str(refant)))
        args.append(globals.Variable("line","chan,1,%i,%i" % (startChan, numChans)))
        log.run("mselfcal",args)

    # do the same for each source
        for source in objects._sources :
            args = []
            if("LSB" in refWindow or "USB" in refWindow) :
                args.append(globals.Variable("vis",gcal._file,".%s" %(refWindow)))
            else :
                args.append(globals.Variable("vis",gcal._file,".w%i" %(refWindow)))
            args.append(globals.Variable("mode","merge"))
            args.append(globals.Variable("out",source._file,".w%i" %(window)))
            log.run("gpcopy",args,fatal=True)
            args = []
            args.append(globals.Variable("in",source._file,".w%i/interval" %(window)))
            args.append(globals.Variable("value","1"))
            log.run("puthd",args)

            args = []
            args.append(globals.Variable("vis",source._file,".w%i" % (window)))
            args.append(globals.Variable("out",source._file,"." + fileEnd + "; sleep 3; rm -rf "))
            args.append(globals.Variable(None,source._file,".w%i/*; rm -rf " % (window)))
            args.append(globals.Variable(None,source._file,".w%i; sleep 3; mv " % (window)))
            args.append(globals.Variable(None,source._file,"." + fileEnd))
            args.append(globals.Variable(None,source._file,".w%i" % (window)))
            log.run("uvcat",args)

            args = []
            args.append(globals.Variable("vis",gcal._file,".w%i" %(window)))
            args.append(globals.Variable("mode","merge"))
            args.append(globals.Variable("out",source._file,".w%i" %(window)))
            log.run("gpcopy",args,fatal=True)
            args = []
            args.append(globals.Variable("in",source._file,".w%i/interval" %(window)))
            args.append(globals.Variable("value","1"))
            log.run("puthd",args)

def continuum(objects, refant, window) :
    """ method to do the gain calibration for a set of windows
        input :
            objects - the objects to be calibrated
            refant - reference antenna
            window - list of the widows to be calibrated
        returns :
            none
    """
    # loop over the gaincals
    for gcal in objects._gaincals :
        if(gcal.isSuper(window)) :
            sideband = "LSB"
            if(window > len(gcal._bandwidths)/2) :
                sideband = "USB"
            log.writeComment("Self calibrating the %s of %s" %(sideband, gcal._name))
            # if we have no flux info
            if(gcal._flux == 0.0) :
                # make sure to ignore the edge channels for the 62 MHz band
                startChan = 1
                numChans = gcal.getSuperNumChans()
                allGood = False
                # while the gains do not all meet the preferences criteria, alternate between selfcal and flagging
                sys = 0
                count = 0
                log.writeComment("Flagging any bad gains")
                while(not allGood and count < 15) :
                    sys = log.run("mselfcal vis=%s interval=%f options=amplitude,apriori,noscale refant=%i line=chan,1,%i,%i" % (gcal._file + "." + sideband,p.preferences.get("selfcalInterval"), refant, startChan,numChans),[],logit=False)
                    allGood = (not flagging.flagByGains(gcal._file, sideband)) or (sys != 0)
                    count += 1
                if(sys == 0) :
                    args = []
                    args.append(globals.Variable("vis",gcal._file,".%s" % (sideband)))
                    args.append(globals.Variable("interval",str(p.preferences.get("selfcalInterval"))))
                    args.append(globals.Variable("options","amplitude,apriori,noscale"))
                    args.append(globals.Variable("refant",str(refant)))
                    args.append(globals.Variable("line","chan,1,%i,%i" % (startChan, numChans)))
                    log.run("mselfcal",args,execute=False)

                    times = calculations.unconvertTime(gcal._start[0] - 1/3600.0) + "," + calculations.unconvertTime(gcal._end[-1] + 1/3600.0)
                    if(times != "") :
                        args = []
                        args.append(globals.Variable("vis",gcal._file,".%s" % (sideband)))
                        args.append(globals.Variable("break",times))
                        log.run("gpbreak",args)

            # copy the gain info to the source files
                for source in objects._sources :
                    args = []
                    args.append(globals.Variable("vis",gcal._file,".%s" % (sideband)))
                    args.append(globals.Variable("mode","merge"))
                    args.append(globals.Variable("out",source._file,".%s" % (sideband)))
                    log.run("gpcopy",args,fatal=True)
                    args = []
                    args.append(globals.Variable("in",source._file,".%s/interval" %(sideband)))
                    args.append(globals.Variable("value","1"))
                    log.run("puthd",args,fatal=True)

        # if we have the flux info
            else :
            # make sure to ignore the edge channels of the 62 MHz bands
                startChan = 1
                numChans = gcal.getSuperNumChans()
                allGood = False
                # while the gains do not all meet the preferences criteria, alternate between selfcal and flagging
                sys = 0
                count = 0
                log.writeComment("Flagging any bad gains")
                while(not allGood and count < 15) :
                    sys = log.run("mselfcal vis=%s interval=%f options=amplitude,noscale flux=%f refant=%i line=chan,1,%i,%i" % (gcal._file + ".%s" %(sideband), p.preferences.get("selfcalInterval"),gcal._flux, refant, startChan, numChans),[],logit=False)
                    allGood = (not flagging.flagByGains(gcal._file, sideband)) or (sys != 0)
                    count += 1
                if(sys == 0) :
                    args = []
                    args.append(globals.Variable("vis",gcal._file,".%s" %(sideband)))
                    args.append(globals.Variable("interval",str(p.preferences.get("selfcalInterval"))))
                    args.append(globals.Variable("options","amplitude,noscale"))
                    args.append(globals.Variable("flux",str(gcal._flux)))
                    args.append(globals.Variable("refant",str(refant)))
                    args.append(globals.Variable("line","chan,1,%i,%i" % (startChan, numChans)))
                    log.run("mselfcal",args,execute=False)

                # add breakpoints to the begining and end of the data
                    times = calculations.unconvertTime(gcal._start[0] - 1/3600.0) + "," + calculations.unconvertTime(gcal._end[-1] + 1/3600.0)
                    if(times != "") :
                        args = []
                        args.append(globals.Variable("vis",gcal._file,".%s" % (sideband)))
                        args.append(globals.Variable("break",times))
                        log.run("gpbreak",args)

            # copy gains to source
                for source in objects._sources :
                    args = []
                    args.append(globals.Variable("vis",gcal._file,".%s" %(sideband)))
                    args.append(globals.Variable("mode","merge"))
                    args.append(globals.Variable("out",source._file,".%s" %(sideband)))
                    log.run("gpcopy",args,fatal=True)
                    args = []
                    args.append(globals.Variable("in",source._file,".%s/interval" %(sideband)))
                    args.append(globals.Variable("value","1"))
                    log.run("puthd",args,fatal=True)
        else :
            log.writeComment("Self calibrating window %i of %s" %(window, gcal._name))
            # if we have no flux info
            if(gcal._flux == 0.0) :
                # make sure to ignore the edge channels for the 62 MHz band
                if(gcal._bandwidths[window - 1] == 62) :
                    startChan = 4
                    numChans = gcal._numChans[window - 1] - 6
                else:
                    startChan = 1
                    numChans = gcal._numChans[window - 1]
                allGood = False
                # while the gains do not all meet the preferences criteria, alternate between selfcal and flagging
                sys = 0
                count = 0
                log.writeComment("Flagging any bad gains")
                while(not allGood and count < 15) :
                    sys = log.run("mselfcal vis=%s interval=%f options=amplitude,apriori,noscale refant=%i line=chan,1,%i,%i" % (gcal._file + ".w%i" %(window),p.preferences.get("selfcalInterval"), refant, startChan,numChans),[],logit=False)
                    allGood = (not flagging.flagByGains(gcal._file, window)) or (sys != 0)
                    count += 1
                if(sys == 0) :
                    args = []
                    args.append(globals.Variable("vis",gcal._file,".w%i" %(window)))
                    args.append(globals.Variable("interval",str(p.preferences.get("selfcalInterval"))))
                    args.append(globals.Variable("options","amplitude,apriori,noscale"))
                    args.append(globals.Variable("refant",str(refant)))
                    args.append(globals.Variable("line","chan,1,%i,%i" % (startChan, numChans)))
                    log.run("mselfcal",args,execute=False)

                    times = calculations.unconvertTime(gcal._start[0] - 1/3600.0) + "," + calculations.unconvertTime(gcal._end[-1] + 1/3600.0)
                    if(times != "") :
                        args = []
                        args.append(globals.Variable("vis",gcal._file,".w%i" % (window)))
                        args.append(globals.Variable("break",times))
                        log.run("gpbreak",args)

                # copy the gain info to the source files
                for source in objects._sources :
                    args = []
                    args.append(globals.Variable("vis",gcal._file,".w%i" %(window)))
                    args.append(globals.Variable("mode","merge"))
                    args.append(globals.Variable("out",source._file,".w%i" %(window)))
                    log.run("gpcopy",args,fatal=True)
                    args = []
                    args.append(globals.Variable("in",source._file,".w%i/interval" %(window)))
                    args.append(globals.Variable("value","1"))
                    log.run("puthd",args,fatal=True)

            # if we have the flux info
            else :
                # make sure to ignore the edge channels of the 62 MHz bands
                if(gcal._bandwidths[window - 1] == 62) :
                    startChan = 4
                    numChans = gcal._numChans[window - 1] - 6
                else:
                    startChan = 1
                    numChans = gcal._numChans[window - 1]
                allGood = False
                # while the gains do not all meet the preferences criteria, alternate between selfcal and flagging
                sys = 0
                count = 0
                log.writeComment("Flagging any bad gains")
                while(not allGood and count < 15) :
                    sys = log.run("mselfcal vis=%s interval=%f options=amplitude,noscale flux=%f refant=%i line=chan,1,%i,%i" % (gcal._file + ".w%i" %(window), p.preferences.get("selfcalInterval"),gcal._flux, refant, startChan, numChans),[],logit=False)
                    allGood = (not flagging.flagByGains(gcal._file, window)) or (sys != 0)
                    count += 1
                if(sys == 0) :
                    args = []
                    args.append(globals.Variable("vis",gcal._file,".w%i" %(window)))
                    args.append(globals.Variable("interval",str(p.preferences.get("selfcalInterval"))))
                    args.append(globals.Variable("options","amplitude,noscale"))
                    args.append(globals.Variable("flux",str(gcal._flux)))
                    args.append(globals.Variable("refant",str(refant)))
                    args.append(globals.Variable("line","chan,1,%i,%i" % (startChan, numChans)))
                    log.run("mselfcal",args,execute=False)

                    # add breakpoints to the begining and end of the data
                    times = calculations.unconvertTime(gcal._start[0] - 1/3600.0) + "," + calculations.unconvertTime(gcal._end[-1] + 1/3600.0)
                    if(times != "") :
                        args = []
                        args.append(globals.Variable("vis",gcal._file,".w%i" % (window)))
                        args.append(globals.Variable("break",times))
                        log.run("gpbreak",args)

                # copy gains to source
                for source in objects._sources :
                    args = []
                    args.append(globals.Variable("vis",gcal._file,".w%i" %(window)))
                    args.append(globals.Variable("mode","merge"))
                    args.append(globals.Variable("out",source._file,".w%i" %(window)))
                    log.run("gpcopy",args,fatal=True)
                    args = []
                    args.append(globals.Variable("in",source._file,".w%i/interval" %(window)))
                    args.append(globals.Variable("value","1"))
                    log.run("puthd",args,fatal=True)

def bootstrap(objects, refant, refWindow, window) :
    """ Method to bootstrap gain solutions from a wideband window to a narrowband window
        input :
            objects - the objects
            refant - reference antenna
            refWindow - the wideband reference window to be bootstrapped from
            window - the list of windows to be bootstrapped
        returns :
            none
    """
    threadList = []
    for gcal in objects._gaincals :
        current = bootcalThread(objects,gcal,refant,refWindow,window)
        threadList.append(current)
        current.start()
    for thread in threadList :
        thread.join()
