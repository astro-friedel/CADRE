#!/usr/bin/python

"""
Copyright (C) 2011 Board of Trustees of the University of Illinois

Current Version (this file): 1.1
Released with pipeline version: 1.1
Author: D. N. Friedel
"""
import sys
config = True
for i in sys.argv :
    if(i == "--no-config") :
        config = False

if(config) :
    import configure

import sources
import math
import os
import bandpasscal
import bootflux
import calculations
import baselines
import flagging
import logger as log
import listobs
import globals
import gainCalibration
import invert
import startupTeardown
import version
import traceback
import time
import random
import hybrid
from threading import Thread

prjData = dict()

"""
Main module for the CARMA data reduction pipeline
Use: python pipeline.py <vis_file> [type=c/s]
    <vis_file> is the miriad file to be reduced
    type is what type of data reduction to do
        c = continuum only
        s = spectral line and continuum (if there are wider bands used)
Author: D. N. Friedel
"""
objects = sources.Objects()
try :
    sys.path.append(".")

    # import the preferences file
    try :
        import preferences as p
    except ImportError:
        import defaultPreferences as p

    try:
        import numpy
        numpy_loaded = True
    except ImportError:
        numpy_loaded = False
        print "The NumPy module was not found on this machine. The pipeline will continue, however continuum subtraction and smart clean region detection will not be available. The NumPy module can be downloaded from http://www.scipy.org/Download"
        print "Continuing in..."
        for i in range(15,0,-1) :
            sys.stdout.write("%i%s\r" % (i," "*2))
            sys.stdout.flush()
            time.sleep(1)


    globals.setScriptVar(str(p.preferences.get("selfcalInterval")),"SEFLCALINTERVAL")

    if(len(sys.argv) < 2 or sys.argv[1] == "-h" or sys.argv[1] == "--h" or sys.argv[1] == "help" or sys.argv[1] == "-help" or sys.argv[1] == "--help") :
        version.help()
        os._exit(0)
    elif(sys.argv[1] == "v" or sys.argv[1] == "-v" or sys.argv[1] == "--v") :
        version.version()
        os._exit(0)

    visFile = sys.argv[1]
    CONTINUUM = False
    SPECTRAL = True
    type = "spectral line"
    if(len(sys.argv) > 2) :
        temp= sys.argv[2]
        mode = temp.split("=")
        if(mode[1] == "C" or mode[1] == "c") :
            CONTINUUM = True
            SPECTRAL = False
            type = "continuum"

    args = []
    if(".tar.gz" in visFile) :
        args.append(globals.Variable(None,"zxvf"))
        visFile = visFile[:visFile.index(".tar.gz")]
        args.append(globals.Variable(None,visFile,".tar.gz"))
    elif(".tar" in visFile) :
        args.append(globals.Variable(None,"xvf"))
        visFile = visFile[:visFile.index(".tar")]
        args.append(globals.Variable(None,visFile,".tar"))

    globals.setScriptVar(visFile,"MAINFILE")

    fluxTrack = False
    prjData = calculations.splitObsblockId(visFile)
    project = prjData["project"]
    if(project == "flux") :
        fluxTrack = True
        SPECTRAL = False
        CONTINUUM = True
    if("S" in prjData["configuration"]) :
        globals.isSci2 = True

    # initialize what needs to be initialized
    proj = ""
    startupTeardown.startup(visFile,args)

    log.writeLog("Pipeline is doing %s data reduction" % (type))

    # run listobs
    listobs.runListobs(visFile, objects)

    # get global variables
    obsFreq = globals.obsFreq()
    obsDate = globals.obsDate()
    avgBaseline = globals.avgBaseline()
    refant = globals.refant()

    # determine which windows to process
    if(len(objects._sources) > 0) :
        numberOfWins = len(objects._sources[0]._bandwidths)
    elif(len(objects._gaincals) > 0) :
        numberOfWins = len(objects._gaincals[0]._bandwidths)

    globals.ENDWINDOW = numberOfWins
    if(globals.isSci2) :
        if(obsFreq < 50.0) :
            globals.ENDWINDOW = numberOfWins/2
        else :
            globals.STARTWINDOW = numberOfWins/2 + 1

    # flag based on system temperature and shadowing
    log.writeHeader(["Flagging bad data"])
    flagging.flagTsys(visFile)
    log.writeAll("\n")
    flagging.flagShadowing(visFile)
    log.writeAll("\n")
    flagging.flagElevation(visFile)
    # apply the correct baseline solution
    if(p.preferences.get("doBaselines")) :
        log.writeHeader(["Determine if a new baseline solution is available,","   and apply as necessary"])
        baselines.applyBaselines(visFile,calculations.gregorianToNumeric(obsDate))

    # separate the objects into individual files and apply linelength calibration
    # some uvcat processes may fail due to the source being in several lists
    sources.sort(visFile, objects,refant)

    # look for bad amplitudes and  birdies in source(s) and gaincalibrator(s)
    # flagged birides are applied to all matching (bandwidth and number of channels) windows on all objects
    threadList = []
    for object in (objects._sources + objects._gaincals) :
        if(object in objects._sources) :
            current = flagging.FlagBadAmpsThread(object,objects,onlyBirdies=True)
        else :
            current = flagging.FlagBadAmpsThread(object,objects)
        threadList.append(current)
        current.start()
    for thread in threadList :
        thread.join()

    # bandpass calibration
    log.writeHeader(["Bandpass calibration"])
    threadList = []
    LSBdone = False
    USBdone = False
    for window in range(globals.STARTWINDOW, globals.ENDWINDOW + 1) :
        if(objects._passcals[0].isSuper(window)) :
            if(window <= len(objects._passcals[0]._bandwidths)/2 and not LSBdone) :
                LSBdone = True
            elif(window > len(objects._passcals[0]._bandwidths)/2 and not USBdone) :
                USBdone = True
            else :
                continue
        current = bandpasscal.bandpassThread(objects,refant,window)
        threadList.append(current)
        current.start()
    for thread in threadList :
        thread.join()

    if(globals.hybrid() and hybrid.calculateOffsets(objects._passcals, refant)) :
        hybrid.applyOffsets(objects._passcals,objects._sources)
    # now do bootflux
    # determine if there is a planet for the flux calibration

    # choose the best first and work toward the least desirable
    log.writeHeader(["Determining flux of gain calibrators"])
    if(len(objects._fluxcals) != 0) :
        primaryFluxcal = sources.getPrimaryFluxCal(objects)
        for gcal in objects._gaincals :
            flux = bootflux.runBootflux(prjData["project"]+"."+prjData["obsblock"]+"."+prjData["subObsblock"]+"."+str(prjData["trial"]),gcal,primaryFluxcal, calculations.gregorianToNumeric(obsDate), obsFreq)
            gcal.setFlux(flux)
            globals.setScriptVar(str(flux),"FLUX")
            log.writeLog("Found a flux of %f for %s" % (gcal._flux, gcal._name))
    else : # need to find flux of gain calibrator(s)
        log.writeComment("No primary flux calibrator found. Attempting to obtain fluxes of gaincals from catalog")
        for gcal in objects._gaincals :
            flux = bootflux.getCatalogFlux(gcal._name,calculations.gregorianToNumeric(obsDate),obsFreq)
            gcal.setFlux(flux)
            globals.setScriptVar(str(flux),"FLUX")
            log.writeComment("Found a flux of %f for %s" %(gcal._flux, gcal._name))

    # need to determine primary gaincal(s)
    if(len(objects._gaincals) > 1) :
        log.writeComment("There are multiple gain calibrators, sorting them by observation time")
        sources.calSort(objects)

    # gain calibration
    lsbRefWindows = []
    usbRefWindows = []
    specWindows = []
    log.writeHeader(["Gain calibration","","Note that this section of the reduction is threaded","and the commands for each window will not necessarily appear together"])
    if(SPECTRAL) :
        maxBw = objects._gaincals[0].getMaxBw()
        minBw = objects._gaincals[0].getMinBw()
        midpoint = len(objects._gaincals[0]._bandwidths)/2
        LSBdone = False
        USBdone = False
        if(maxBw != minBw) :
            threadList = []
            for window in range(globals.STARTWINDOW, globals.ENDWINDOW + 1) :
                if(objects._gaincals[0]._bandwidths[window - 1] == maxBw) :
                    if(window <= midpoint):
                        lsbRefWindows.append(window)
                        if(objects._gaincals[0].isSuper(window)) :
                            if(not LSBdone) :
                                LSBdone = True
                            else :
                                continue
                    else:
                        usbRefWindows.append(window)
                        if(objects._gaincals[0].isSuper(window)) :
                            if(not USBdone) :
                                USBdone = True
                            else :
                                continue
                    current = gainCalibration.continuumThread(objects, refant, window)
                    threadList.append(current)
                    current.start()
                else :
                    specWindows.append(window)
            for thread in threadList :
                thread.join()
            threadList = []
            for window in specWindows :
                refWindow = None
                if(window <= midpoint) :
                    refWindow = lsbRefWindows[0]
                    if(objects._sources[0].haveSuper()) :
                        refWindow = "LSB"
                else:
                    refWindow = usbRefWindows[0]
                    if(objects._sources[0].haveSuper()) :
                        refWindow = "USB"

                current = gainCalibration.bootstrapThread(objects, refant, refWindow, window)
                threadList.append(current)
                current.start()
            for thread in threadList :
                thread.join()
        else :
            threadList = []
            for window in range(globals.STARTWINDOW, globals.ENDWINDOW + 1) :
                if(window <= midpoint):
                    lsbRefWindows.append(window)
                    if(objects._gaincals[0].isSuper(window)) :
                        if(not LSBdone) :
                            LSBdone = True
                        else :
                            continue
                else:
                    usbRefWindows.append(window)
                    if(objects._gaincals[0].isSuper(window)) :
                        if(not USBdone) :
                            USBdone = True
                        else :
                            continue

                current = gainCalibration.continuumThread(objects, refant, window)
                threadList.append(current)
                current.start()
            for thread in threadList :
                thread.join()
    else :
        LSBdone = False
        USBdone = False
        threadList = []
        for window in range(globals.STARTWINDOW, globals.ENDWINDOW + 1) :
            if(window <= midpoint):
                if(objects._gaincals[0].isSuper(window)) :
                    if(not LSBdone) :
                        LSBdone = True
                    else :
                        continue
            else :
                if(objects._gaincals[0].isSuper(window)) :
                    if(not USBdone) :
                        USBdone = True
                    else :
                        continue

            current = gainCalibration.continuumThread(objects, refant, window)
            threadList.append(current)
            current.start()
        for thread in threadList :
            thread.join()

    if(not fluxTrack) :
        # determine if any sources are secondary calibrators and do a selfcal on them and apply to other sources
        secondaries = []
        srcs = []
        calDone = False
        uvcatDone = False
        for source in objects._sources :
            if(("3C" == source._name[0:2] or ( len(source._name) == 8 and ("+" == source._name[4] or "-" == source._name[4]))) and len(objects._sources) > 1) :
                secondaries.append(source)
            else :
                srcs.append(source)
        if(len(secondaries) > 0) :
            log.writeComment("Found secondary calibrators, selfcaling and applying their solution to all sources")
            if(len(secondaries) > 1) :
                log.writeComment("Multiple secondaries found, sorting by distance and incrementally applying gain solutions")
                # sort the secondary calibrators in order of distance from the source so that the solutions get incrementally better as we approach
                sources.sortDistances(secondaries,srcs)
            # Cannot thread this part as it must be done incrementally
            for i in range(0,len(secondaries)) :
                source = secondaries[i]
                list = []
                wideList = []
                bootstrap = False
                LSBdone = False
                USBdone = False
                sideband = "LSB"
                for window in range(globals.STARTWINDOW, globals.ENDWINDOW + 1) :
                    if(source.isSuper(window)) :
                        if(window <= len(source._bandwidths)/2 and not LSBdone) :
                            LSBdone = True
                        elif(window > len(source._bandwidths)/2 and not USBdone):
                            sideband = "USB"
                            USBdone = True
                        else :
                            continue
                        fileEnd = str(random.randint(1,100000))
                        args = []
                        args.append(globals.Variable("vis",source._file,".%s" % (sideband)))
                        args.append(globals.Variable("options","unflagged"))
                        args.append(globals.Variable("out","temp." + fileEnd, "; sleep 3; rm -rf "))
                        args.append(globals.Variable(None,source._file,".%s/*; rm -rf " % (sideband)))
                        args.append(globals.Variable(None,source._file,".%s; sleep 3; mv temp.%s " % (sideband,fileEnd)))
                        args.append(globals.Variable(None,source._file,".%s" % (sideband)))
                        log.run("uvcat",args,fatal=True)

                        startChan = 1
                        numChans = source.getSuperNumChans()
                        wideList.append(window)
                        bootstrap = True
                        args = []
                        args.append(globals.Variable("vis",source._file,".%s" % (sideband)))
                        args.append(globals.Variable("interval",str(p.preferences.get("selfcalInterval"))))
                        args.append(globals.Variable("options","phase,apriori,noscale"))
                        args.append(globals.Variable("refant",str(refant)))
                        args.append(globals.Variable("line","chan,1,%i,%i" % (startChan,numChans)))
                        log.run("mselfcal",args,fatal=True)

                    else :
                        fileEnd = str(random.randint(1,100000))
                        args = []
                        args.append(globals.Variable("vis",source._file,".w%i" % (window)))
                        args.append(globals.Variable("options","unflagged"))
                        args.append(globals.Variable("out","temp." + fileEnd, "; sleep 3; rm -rf "))
                        args.append(globals.Variable(None,source._file,".w%i/*; rm -rf " % (window)))
                        args.append(globals.Variable(None,source._file,".w%i; sleep 3; mv temp.%s " % (window,fileEnd)))
                        args.append(globals.Variable(None,source._file,".w%i" % (window)))

                        log.run("uvcat",args,fatal=True)

                        # see if we need to bootstrap
                        if((source.getMaxBw != source.getMinBw) and source._bandwidths[window - 1] < source.getMaxBw) :
                            list.append(window)
                            continue
                        if(source._bandwidths[window - 1] == 62) :
                            startChan = 4
                            numChans = source._numChans[window - 1] - 6
                        else:
                            startChan = 1
                            numChans = source._numChans[window - 1]
                        wideList.append(window)
                        bootstrap = True
                        args = []
                        args.append(globals.Variable("vis",source._file,".w%i" % (window)))
                        args.append(globals.Variable("interval",str(p.preferences.get("selfcalInterval"))))
                        args.append(globals.Variable("options","phase,apriori,noscale"))
                        args.append(globals.Variable("refant",str(refant)))
                        args.append(globals.Variable("line","chan,1,%i,%i" % (startChan,numChans)))
                        log.run("mselfcal",args,fatal=True)

                # bootstrap the secondary if necessary
                if(SPECTRAL and bootstrap) :
                    LSBdone = False
                    USBdone = False
                    for window in list :
                        if(window <= len(source._bandwidths)/2) :
                            ending = ".w%i" % (lsbRefWindows[0])
                            if(source.haveSuper()) :
                                if(not LSBdone) :
                                    ending = "LSB"
                                    LSBdone = True
                                else :
                                    continue
                            args = []
                            args.append(globals.Variable("vis",source._file,".%s" % (ending)))
                            args.append(globals.Variable("mode","merge"))
                            args.append(globals.Variable("out",source._file,".w%i" % (window)))
                            log.run("gpcopy",args,fatal=True)
                            args = []
                            args.append(globals.Variable("in",source._file,".w%i/interval" % (window)))
                            args.append(globals.Variable("value","1"))
                            log.run("puthd",args,fatal=True)

                        else :
                            ending = ".w%i" % (usbRefWindows[0])
                            if(source.haveSuper()) :
                                if(not USBdone) :
                                    ending = "USB"
                                    USBdone = True
                                else :
                                    continue
                            args = []
                            args.append(globals.Variable("vis",source._file,".%s" % (ending)))
                            args.append(globals.Variable("mode","merge"))
                            args.append(globals.Variable("out",source._file,".w%i" % (window)))
                            log.run("gpcopy",args,fatal=True)
                            args = []
                            args.append(globals.Variable("in",source._file,".w%i/interval" % (window)))
                            args.append(globals.Variable("value","1"))
                            log.run("puthd",args,fatal=True)

                list += wideList
                # copy the gains from this secondary to the others that have not been done yet
                for j in range(i+1,len(secondaries)):
                    bootsource = secondaries[j]
                    LSBdone = False
                    USBdone = False
                    for window in list :
                        # apply any old gains tables so there is no confusion
                        fileEnd = str(random.randint(1,100000))
                        ending = ".w%i" % (lsbRefWindows[0])
                        if(source.haveSuper()) :
                            if(window <= len(source._bandwidths)/2 and not LSBdone) :
                                ending = "LSB"
                                LSBdone = True
                            elif(window > len(source._bandwidths)/2 and not USBdone) :
                                ending = "USB"
                                USBdone = True
                            else :
                                continue
                        args = []
                        args.append(globals.Variable("vis",bootsource._file,".%s" % (ending)))
                        args.append(globals.Variable("out","temp." + fileEnd + "; rm -rf"))
                        args.append(globals.Variable(None,bootsource._file,".%s/*; rm -rf " % (ending)))
                        args.append(globals.Variable(None,bootsource._file,".%s; sleep 3; mv temp.%s " % (ending,fileEnd)))
                        args.append(globals.Variable(None,bootsource._file,".%s" % (ending)))
                        log.run("uvcat",args,fatal=True)

                        args = []
                        args.append(globals.Variable("vis",source._file,".%s" % (ending)))
                        args.append(globals.Variable("mode","merge"))
                        args.append(globals.Variable("out",bootsource._file,".%s" % (ending)))
                        log.run("gpcopy",args,fatal=True)
                        args = []
                        args.append(globals.Variable("in",bootsource._file,".%s/interval" % (ending)))
                        args.append(globals.Variable("value","1"))
                        log.run("puthd",args,fatal=True)

                        uvcatDone = True
        # copy and apply all gain corrections to the source(s)
        if(len(secondaries) > 0 and len(srcs) > 0) :
            for source in srcs :
                LSBdone = False
                USBdone = False
                for window in range(globals.STARTWINDOW, globals.ENDWINDOW + 1) :
                    fend = "w%i" % (window)
                    if(source.isSuper(window)) :
                        if(window <= len(source._bandwidths)/2 and not LSBdone) :
                            LSBdone = True
                            fend = "LSB"
                        elif(window > len(source._bandwidths)/2 and not USBdone) :
                            USBdone = True
                            fend = "USB"
                        else :
                            continue
                    for sec in secondaries :
                        fileEnd = str(random.randint(1,100000))
                        args = []
                        args.append(globals.Variable("vis",sec._file,".%s" % (fend)))
                        args.append(globals.Variable("mode","merge"))
                        args.append(globals.Variable("out",source._file,".%s" % (fend)))
                        log.run("gpcopy",args,fatal=True)
                        args = []
                        args.append(globals.Variable("in",source._file,".%s/interval" % (fend)))
                        args.append(globals.Variable("value","1"))
                        log.run("puthd",args,fatal=True)

                        args = []
                        args.append(globals.Variable("vis",source._file,".%s" % (fend)))
                        args.append(globals.Variable("out","temp." + fileEnd + "; rm -rf"))
                        args.append(globals.Variable(None,source._file,".%s/*; rm -rf " % (fend)))
                        args.append(globals.Variable(None,source._file,".%s; sleep 3; mv temp.%s " % (fend,fileEnd)))
                        args.append(globals.Variable(None,source._file,".%s" % (fend)))
                        log.run("uvcat",args,fatal=True)

        # we should now be calibrated
        log.writeHeader(["Generate continuum and specrtal line maps of source(s) and","   continuum map of gain calibrator(s)","","Note that this section of the reduction is threaded","and the commands for each window will not necessarily appear together"])
        if(SPECTRAL) :
            for source in objects._sources :
                if(source.haveSuper()) :
                    midpoint = len(source._superwindow)/2
                    for i in range(0,midpoint) :
                        if(not globals.isSci2 or (globals.isSci2 and obsFreq < 50.0)) :
                            args = []
                            args.append(globals.Variable("vis",source._file,".LSB"))
                            args.append(globals.Variable("select","win'('%i')'" % (i+1)))
                            args.append(globals.Variable("out",source._file,".w%i" % (source._superwindow[i])))
                            log.run("uvcat",args,fatal=True)
                        if(not globals.isSci2 or (globals.isSci2 and obsFreq > 50.0)) :
                            args = []
                            args.append(globals.Variable("vis",source._file,".USB"))
                            args.append(globals.Variable("select","win'('%i')'" % (i+1)))
                            args.append(globals.Variable("out",source._file,".w%i" % (source._superwindow[i+midpoint])))
                            log.run("uvcat",args,fatal=True)

            # only do the continuum if there are some windows that are at least 62 MHz wide
            if(maxBw > 60) :
                invert.invertContinuum(objects, obsFreq, avgBaseline, lsbRefWindows + usbRefWindows)
                if(globals.isSci2) :
                    threadList = []
                    for window in lsbRefWindows + usbRefWindows :
                        current = invert.invertContinuumThread(objects, obsFreq, avgBaseline ,window)
                        threadList.append(current)
                        current.start()
                    for thread in threadList :
                        thread.join()

            specWindows = specWindows + lsbRefWindows + usbRefWindows
            invert.invertSpectra(objects, obsFreq, avgBaseline, specWindows)
        else :
            invert.invertContinuum(objects, obsFreq, avgBaseline, [])

    # clean up temporary files and finish up logs
    startupTeardown.cleanup(objects)
except :
    # if anything went wrong
    message  = "%s : %s" % (sys.exc_type,sys.exc_value)
    trace = traceback.print_exc()
    try :
        trace = open("trace.message","w")
        traceback.print_exc(file=trace)
        trace.close()
        traceback.print_exc()
    except :
        print "An exception occurred, however trace.message could not be opened, traceback is :\n"
        traceback.print_exc()
        print "\n\nOriginal exception :\n"
        print message
        if(log.messageLog != None) :
            log.writeComment("An exception occurred, however trace.message could not be opened, traceback is :")
            traceback.print_exc(file=log.messageLog)
            log.writeComment("\n\nOriginal exception :")
            log.writeComment(message)
            startupTeardown.cleanup(True)
        os._exit(0)
    print "An exception occurred, traceback is in trace.message\n"
    print "Exception :",message
    if(log.messageLog != None) :
        log.writeComment("An exception occurred, traceback is in trace.message")
        log.writeComment(message)
        startupTeardown.cleanup(objects,True)
    os._exit(0)
