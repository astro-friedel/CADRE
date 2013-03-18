import sources
import calculations
import globals
import logger as log
from pipeline_miriadwrap import *
import calflux
import math
try :
    import preferences as p
except ImportError:
    import defaultPreferences as p
haveLog = False
try :
    import fluxlog
    haveLog = True
except :
    pass
"""
Module to calculate the flux of the gain calibrator(s) by bootstrapping from flux calibrator(s) or by looking up the flux in FluxSource.cat as a fallback
Part of the CARMA data reduction pipeline
Author: D. N. Friedel
"""

averagingTime = p.preferences.get("bootfluxInterval")

def runBootflux(project,gaincal, fluxcal, obsDate, obsFreq=95.0) :
    """ Method to calculate the flux of the gain calibrator(s)
        by bootstrapping from flux calibrator(s) or by looking up
        the flux in FluxSource.cat as a fallback
        input :
            gaincal - the gain calibrator
            fluxcal -the flux calibrator
            obsDate - observation date
            obsFreq - observation frequency in GHz
        returns :
            the flux of gaincal in Jy
    """
    global averagingTime
    fluxes = []
    rms = []

    log.writeComment("Using %s as flux calibrator" % (fluxcal._name))
    # prefer a planet for bootflux
    # it is done one window at a time and then averaged to get flux
    if(fluxcal._type == sources.PLANET) :
        # only rely on super windwos if we have them
        if(fluxcal.haveSuper()) :
            if(gaincal._lsbGood) :
                args = []
                args.append(globals.Variable("vis",gaincal._file,".LSB"))
                args.append(globals.Variable("ADD",fluxcal._file,".LSB"))
                args.append(globals.Variable("primary",fluxcal._name))
                args.append(globals.Variable("line","channel,1,1,%i" % (fluxcal.getSuperNumChans())))
                args.append(globals.Variable("taver",str(averagingTime)))
                args.append(globals.Variable("log","boot.%s.LSB.log" % (gaincal._file)))
                sys = log.run("bootflux",args)

                if(sys == 0) :
                    if(getFlux("boot.%s.LSB.log" % (gaincal._file),fluxes, rms)) :
                        log.writeComment("Bootflux got a flux of %f for the LSB of %s\n" % (fluxes[0],gaincal._name))
                    else :
                        fluxes.append(0.0)
                        rms.append(0.0)
                else :
                    fluxes.append(0.0)
                    rms.append(0.0)
            if(gaincal._usbGood) :
                indx = 1
                if(globals.isSci2 and globals.obsFreq() > 50.0) :
                    indx = 0
                args = []
                args.append(globals.Variable("vis",gaincal._file,".USB"))
                args.append(globals.Variable("ADD",fluxcal._file,".USB"))
                args.append(globals.Variable("primary",fluxcal._name))
                args.append(globals.Variable("line","channel,1,1,%i" % (fluxcal.getSuperNumChans())))
                args.append(globals.Variable("taver",str(averagingTime)))
                args.append(globals.Variable("log","boot.%s.USB.log" % (gaincal._file)))
                sys = log.run("bootflux",args)

                if(sys == 0) :
                    if(getFlux("boot.%s.USB.log" % (gaincal._file),fluxes, rms)) :
                        log.writeComment("Bootflux got a flux of %f for the USB of %s\n" % (fluxes[indx],gaincal._name))
                    else :
                        fluxes.append(0.0)
                        rms.append(0.0)
                else :
                    fluxes.append(0.0)
                    rms.append(0.0)
            flux,uncert = calculations.weightedAverage(fluxes, rms)
            if(flux != 0.0) :
                fluxlog.writeLog(project,gaincal._name,flux,uncert,obsFreq,"WB",obsDate,fluxcal._name)
        else :
            for window in range(1, len(gaincal._bandwidths) + 1):
                args = []
                args.append(globals.Variable("vis",gaincal._file,".w%i" % (window)))
                args.append(globals.Variable("ADD",fluxcal._file,".w%i" % (window)))
                args.append(globals.Variable("primary",fluxcal._name))
                args.append(globals.Variable("line","channel,1,1,%i" % (gaincal._numChans[window - 1])))
                args.append(globals.Variable("taver",str(averagingTime)))
                args.append(globals.Variable("log","boot.%s.w%i.log" % (gaincal._file,window)))
                sys = log.run("bootflux",args)

                if(sys == 0) :
                    if(getFlux("boot.%s.w%i.log" % (gaincal._file,window),fluxes, rms)) :
                        log.writeComment("Bootflux got a flux of %f for window %i of %s\n" % (fluxes[window - 1], window,gaincal._name))
                    else :
                        fluxes.append(0.0)
                        rms.append(0.0)
                else :
                    fluxes.append(0.0)
                    rms.append(0.0)
            flux,uncert = calculations.weightedAverage(fluxes, rms)
        if(flux == 0.0) :
            log.writeComment("Obtaining flux from catalog.")
            return getCatalogFlux(gaincal._name, obsDate, obsFreq)
        else :
            pass
            #fluxlog.writeLog(project,gaincal._name,flux,uncert,obsFreq,"NB",obsDate,fluxcal._name)
        return flux
    # if we do not have a planet then using MWC349 (if present), since is 1.0 Jy at 3mm and 1mm
    elif(fluxcal._name == "MWC349") :
        mflux = 1.0
        if(obsFreq > 120.0) :
            mflux = 1.8
        log.writeComment("Assuming a flux of %f for MWC349" % (mflux))
        if(fluxcal.haveSuper()) :
            if(gaincal._lsbGood) :
                args = []
                args.append(globals.Variable("vis",gaincal._file,".LSB"))
                args.append(globals.Variable("ADD",fluxcal._file,".LSB"))
                args.append(globals.Variable("primary",fluxcal._name))
                args.append(globals.Variable("ADD",str(mflux)))
                args.append(globals.Variable("line","channel,1,1,%i" % (fluxcal.getSuperNumChans())))
                args.append(globals.Variable("taver",str(averagingTime)))
                args.append(globals.Variable("log","boot.%s.LSB.log" % (gaincal._file)))
                sys = log.run("bootflux",args)

                if(sys == 0) :
                    if(getFlux("boot.%s.LSB.log" % (gaincal._file),fluxes, rms)) :
                        log.writeComment("Bootflux got a flux of %f for the LSB of %s\n" % (fluxes[0],gaincal._name))
                    else :
                        fluxes.append(0.0)
                        rms.append(0.0)
                else :
                    fluxes.append(0.0)
                    rms.append(0.0)
            if(gaincal._usbGood) :
                args = []
                args.append(globals.Variable("vis",gaincal._file,".USB"))
                args.append(globals.Variable("ADD",fluxcal._file,".USB"))
                args.append(globals.Variable("primary",fluxcal._name))
                args.append(globals.Variable("ADD",str(mflux)))
                args.append(globals.Variable("line","channel,1,1,%i" % (fluxcal.getSuperNumChans())))
                args.append(globals.Variable("taver",str(averagingTime)))
                args.append(globals.Variable("log","boot.%s.USB.log" % (gaincal._file)))
                sys = log.run("bootflux",args)

                if(sys == 0) :
                    if(getFlux("boot.%s.USB.log" % (gaincal._file),fluxes, rms)) :
                        log.writeComment("Bootflux got a flux of %f for the USB of %s\n" % (fluxes[0],gaincal._name))
                    else :
                        fluxes.append(0.0)
                        rms.append(0.0)
                else :
                    fluxes.append(0.0)
                    rms.append(0.0)

            flux,uncert = calculations.weightedAverage(fluxes, rms)
            if(flux != 0.0) :
                fluxlog.writeLog(project,gaincal._name,flux,uncert,obsFreq,"WB",obsDate,fluxcal._name + "(" + str(mflux) + ")")
        else :
            for window in range(1, len(gaincal._bandwidths) + 1):
                args = []
                args.append(globals.Variable("vis",gaincal._file,".w%i" % (window)))
                args.append(globals.Variable("ADD",fluxcal._file,"w%i" % (window)))
                args.append(globals.Variable("primary",fluxcal._name))
                args.append(globals.Variable("ADD",str(mflux)))
                args.append(globals.Variable("line","channel,1,1,%i" % (gaincal._numChans[window - 1])))
                args.append(globals.Variable("taver",str(averagingTime)))
                args.append(globals.Variable("log","boot.%s.w%i.log" % (gaincal._file,window)))
                sys = log.run("bootflux",args)

                if(sys == 0) :
                    if(getFlux("boot.%s.w%i.log" % (gaincal._file,window),fluxes, rms)) :
                        log.writeComment("Bootflux got a flux of %f for window %i of %s\n" % (fluxes[1],gaincal._name))
                    else :
                        fluxes.append(0.0)
                        rms.append(0.0)
                else :
                    fluxes.append(0.0)
                    rms.append(0.0)
            flux,uncert = calculations.weightedAverage(fluxes, rms)
            if(flux == 0.0) :
                log.writeComment("Obtaining flux from catalog.")
                return getCatalogFlux(gaincal._name, obsDate, obsFreq)
            else :
                pass
                #fluxlog.writeLog(project,gaincal._name,flux,uncert,obsFreq,"NB",obsDate,fluxcal._name + "(" + str(mflux) + ")")
        return flux
    else :
        # use a gain calibrator as last resort
        log.writeComment("Obtaining flux for %s from catalog." % (fluxcal._name))
        mflux = getCatalogFlux(fluxcal._name, obsDate, obsFreq)
        if(gaincal._name == fluxcal._name) :
            return mflux
        if(fluxcal.haveSuper()) :
            if(gaincal._lsbGood) :
                args = []
                args.append(globals.Variable("vis",gaincal._file,".LSB"))
                args.append(globals.Variable("ADD",fluxcal._file,".LSB"))
                args.append(globals.Variable("primary",fluxcal._name))
                args.append(globals.Variable("ADD",str(mflux)))
                args.append(globals.Variable("line","channel,1,1,%i" % (fluxcal.getSuperNumChans())))
                args.append(globals.Variable("taver",str(averagingTime)))
                args.append(globals.Variable("log","boot.%s.LSB.log" % (gaincal._file)))
                sys = log.run("bootflux",args)

                if(sys == 0) :
                    if(getFlux("boot.%s.LSB.log" % (gaincal._file),fluxes, rms)) :
                        log.writeComment("Bootflux got a flux of %f for the LSB of %s\n" % (fluxes[0],gaincal._name))
                    else :
                        fluxes.append(0.0)
                        rms.append(0.0)
                else :
                    fluxes.append(0.0)
                    rms.append(0.0)
            if(gaincal._usbGood) :
                args = []
                args.append(globals.Variable("vis",gaincal._file,".USB"))
                args.append(globals.Variable("ADD",fluxcal._file,".USB"))
                args.append(globals.Variable("primary",fluxcal._name))
                args.append(globals.Variable("ADD",str(mflux)))
                args.append(globals.Variable("line","channel,1,1,%i" % (fluxcal.getSuperNumChans())))
                args.append(globals.Variable("taver",str(averagingTime)))
                args.append(globals.Variable("log","boot.%s.USB.log" % (gaincal._file)))
                sys = log.run("bootflux",args)

                if(sys == 0) :
                    if(getFlux("boot.%s.USB.log" % (gaincal._file),fluxes, rms)) :
                        log.writeComment("Bootflux got a flux of %f for the USB of %s\n" % (fluxes[1] ,gaincal._name))
                    else :
                        fluxes.append(0.0)
                        rms.append(0.0)
                else :
                    fluxes.append(0.0)
                    rms.append(0.0)

            flux,uncert = calculations.weightedAverage(fluxes, rms)
            if(flux == 0.0) :
                log.writeComment("Obtaining flux from catalog.")
                return getCatalogFlux(gaincal._name, obsDate, obsFreq)
            else :
                fluxlog.writeLog(project,gaincal._name,flux,uncert,obsFreq,"WB",obsDate,fluxcal._name + "(" + str(mflux) + ")")
        else :
            for window in range(1, len(gaincal._bandwidths) + 1):
                mflux = getCatalogFlux(fluxcal._name, obsDate, obsFreq)
                if(gaincal._name == fluxcal._name) :
                    return mflux
                args = []
                args.append(globals.Variable("vis",gaincal._file,".w%i" % (window)))
                args.append(globals.Variable("ADD",fluxcal._file,".w%i" % (window)))
                args.append(globals.Variable("primary",fluxcal._name))
                args.append(globals.Variable("ADD",str(mflux)))
                args.append(globals.Variable("line","channel,1,1,%i" % (gaincal._numChans[window - 1])))
                args.append(globals.Variable("taver",str(averagingTime)))
                args.append(globals.Variable("log","boot.%s.w%i.log" % (gaincal._file,window)))
                sys = log.run("bootflux",args)

                if(sys == 0) :
                    if(getFlux("boot.%s.w%i.log" % (gaincal._file,window),fluxes, rms)) :
                        log.writeComment("Bootflux got a flux of %f for window %i of %s\n" % (fluxes[window - 1], window,gaincal._name))
                    else :
                        fluxes.append(0.0)
                        rms.append(0.0)
                else :
                    fluxes.append(0.0)
                    rms.append(0.0)
            flux,uncert = calculations.weightedAverage(fluxes, rms)
            if(flux == 0.0) :
                log.writeComment("Obtaining flux from catalog.")
                return getCatalogFlux(gaincal._name, obsDate, obsFreq)
            else :
                fluxlog.writeLog(project,gaincal._name,flux,uncert,obsFreq,"WB",obsDate,fluxcal._name + "(" + str(mflux) + ")")
        return flux

def getFlux(bootFile,fluxes, rms) :
    """ Method to read in the bootflux log and get the average flux
        input :
            fluxes - a list to hold the fluxes
            rms - a list to hold their rms uncertainties
        returns :
            True/False if the flux was found
    """
    # open the input file
    input = open(bootFile,'r')
    fileList = input.readlines()
    input.close()
    # get the average flux value and rms and add them to the lists
    while(len(fileList) > 0) :
        line = fileList.pop()
        if("Average Flux" in line) :
            splitLine = line.split()
            if(splitLine[2][0].isdigit() and splitLine[3][0].isdigit()) :
                fluxes.append(float(splitLine[2]))
                rms.append(float(splitLine[3]))
                return True
    return False

def getCatalogFlux(source, date, freq) :
    """ Method to get the catalog flux for the given source based on
        observation frequency (3mm/1mm) and observation date
        if an exact date match is not found then the flux is linearly interpolated
        from the closest dates
        input :
            source - the source name
            date - observation date
            freq - frequency in GHz
        returns :
            the flux of the source
    """
    scaleFactor = 1.0 #only needed if we are scaling from 3mm to 1mm
    fluxes = None
    fluxes = calflux.calflux(source,freq)
    if(len(fluxes) == 0 and freq > 130) :
        log.writeComment("No 1mm fluxes found, scaling the 3mm flux")
        fluxes = calflux.calflux(source,95.0)
        scaleFactor = 0.7
    # read in the flux file

    if(len(fluxes) == 0) :
        log.writeComment("No flux found for %s in FluxSource.cat" % (source))
        return 0.0
    low = dict()
    hi = dict()
    # sort the dates and find the 2 date bracketing the observation date
    # but only accept fluxes within 1 year of observation date

    for f in fluxes :
        fluxDate = f
        delta = date - fluxDate
        if(delta <= 1.0 and delta > 0) :
            low[fluxDate] = fluxes[f]
        elif(delta >= -1.0 and delta < 0) :
            hi[fluxDate] = fluxes[f]
        else :
            return fluxes[f]
    if(len(low) == 0) :
        if(len(hi) == 0) :
            log.writeComment("No flux found for %s in FluxSource.cat" % (source))
            return 0.0
        log.writeComment("Found flux of %d for %s in FluxSource.cat" % (hi.get(min(hi)), source))
        return hi.get(min(hi))
    elif(len(hi) == 0) :
        log.writeComment("Found flux of %d for %s in FluxSource.cat" % (low.get(max(low)), source))
        return low.get(max(low))
    minKey = max(low)
    maxKey = min(hi)
    flux = calculations.linearInterpolation(minKey,low.get(minKey),maxKey,hi.get(maxKey),date)*scaleFactor
    log.writeComment("Found flux of %d for %s in FluxSource.cat" % (flux, source))
    return flux
