import os
import logger as log
import sources
import calculations
import math
import globals
import startupTeardown
from threading import Thread
import random
import peakLocator
import continuumSubtraction
import miriadClasses
from copy import deepcopy

try :
    import preferences as p
except ImportError:
    import defaultPreferences as p

cellsize = None
imsize = None
continDone = False
regions = dict()
rmsList = dict()

"""
Module to invert uv data and produce images
Part of the CARMA data reduction pipeline
Current Version (this file): 1.1
Released with pipeline version: 1.1
Author: D. N. Friedel
"""
class invertWindowThread(Thread) :
    def __init__(self,source,window,iterate,AvgBaseline = 0,csub = False) :
        Thread.__init__(self)
        self.source = source
        self.window = window
        self.iterate = iterate
        self.avgBaseline = AvgBaseline
        self.csub = csub
    def run(self) :
        invert(self.source,self.window,self.iterate,self.csub,self.avgBaseline)

class cleanWindowThread(Thread) :
    def __init__(self,source,window,csub = False) :
        Thread.__init__(self)
        self.source = source
        self.window = window
        self.csub = csub
    def run(self) :
        cleanWindow(self.source,self.window,self.csub)

class invertContinuumThread(Thread) :
    def __init__(self,objects, obsFreq, avgBaseline ,window) :
        Thread.__init__(self)
        self.objects = objects
        self.window = window
        self.obsFreq = obsFreq
        self.avgBaseline = avgBaseline
    def run(self) :
        invertContinuum(self.objects, self.obsFreq, self.avgBaseline, [self.window], individual = True)


def invert(source,window,iterate,csub,avgBaseline = 0) :
    """ Method to invert a spectral window
        input :
            source - the source data
            window - which window is being inverted
            iterate - should we iterate to minimize the number of bad data
            csub - has continuum subtraction been done
        returns :
            none
    """
    global cellsize
    global imsize
    global continDone
    addon = ""
    noise = 0.5
    if(csub) :
        addon = ".cs"
    fileEnd = str(random.randint(1,100000))
    startChan = 1
    numChan = source._numChans[window - 1]
    if(source._bandwidths[window - 1] == 62) :
        startChan = 4
        numChan = source._numChans[window - 1] - 6
    # remove any flagged data we can
    args = []
    args.append(globals.Variable("vis",source._file, addon +".w%i" % (window)))
    args.append(globals.Variable("options","unflagged"))
    args.append(globals.Variable("out","temp.w%i.%s; sleep 3; rm -rf " % (window,fileEnd)))
    args.append(globals.Variable(None,source._file,addon + ".w%i/*; rm -rf " % (window)))
    args.append(globals.Variable(None,source._file,addon + ".w%i; sleep 3; mv " % (window)))
    args.append(globals.Variable(None,"temp.w%i.%s" % (window,fileEnd)))
    args.append(globals.Variable(None,source._file,addon + ".w%i" % (window)))
    log.run("uvcat",args,fatal=True)
    ok = False
    lastbad = 0
    # invert until the number of rejected visibilities is minimized
    while(not ok) :
        log.run("invert vis=%s.w%i map=%s.w%i.map beam=%s.w%i.beam imsize=%i cell=%f sup=0 options=systemp,double,mosaic line=channel,%i,%i,1,1 >& junk" % (source._file+addon, window, source._file+addon, window, source._file+addon, window,imsize,cellsize,numChan, startChan),[],fatal=True,logit=False)
        input = open("junk",'r')
        fileList = input.readlines()
        input.close()
        found = False
        while(len(fileList) > 0 and iterate) :
            line = fileList.pop()
            if("Theoretical rms" in line) :
                noise = float(line.split(":")[1])
        if(not found) :
            ok = True
            log.writeComment("Using an image size of %f for %s" % (imsize, source._name))
            args = []
            args.append(globals.Variable("vis",source._file,addon + ".w%i" % (window)))
            args.append(globals.Variable("map",source._file,addon + ".w%i.map" % (window)))
            args.append(globals.Variable("beam",source._file,addon + ".w%i.beam" % (window)))
            args.append(globals.Variable("imsize",str(imsize)))
            args.append(globals.Variable("cell",str(cellsize)))
            args.append(globals.Variable("sup","0"))
            args.append(globals.Variable("options","systemp,double,mosaic"))
            args.append(globals.Variable("line","channel,%i,%i,1,1" % (numChan,startChan)))
            log.run("invert",args,execute=False)
    # use mossdi to clean
    lastIter = 10000.0
    log.run("mossdi map=%s.w%i.map beam=%s.w%i.beam out=%s.w%i.clean niters=50000 cutoff=%f region=quarter" % (source._file+addon, window, source._file+addon, window, source._file+addon, window,noise*2.5),[],fatal=True,logit=False)
    # determine the synthsized beam size
    log.run("mospsf beam=%s.w%i.beam out=%s.w%i.bm""" % (source._file+addon, window, source._file+addon, window),[],fatal=True,logit=False)
    log.run("imfit in=%s.w%i.bm object=beam region=arcsec,box'('-5,-5,5,5')' > fit.log.%s" % (source._file+addon, window,fileEnd),[],fatal=True,logit=False)
    input = open("fit.log.%s" % (fileEnd),'r')
    fileList = input.readlines()
    input.close()
    majorAxis = 0.0
    minorAxis = 0.0
    positionAngle = 0.0
    while(len(fileList) > 0) :
        line = fileList.pop()
        if("Position angle" in line) :
            splitLine = line.split()
            positionAngle = float(splitLine[3])
        elif("Minor axis" in line) :
            splitLine = line.split()
            minorAxis = float(splitLine[3])
        elif("Major axis" in line) :
            splitLine = line.split()
            majorAxis = float(splitLine[3])
    log.run("restor model=%s.w%i.clean map=%s.w%i.map beam=%s.w%i.beam fwhm=%f,%f pa=%f out=%s.w%i.finalmap" % (source._file+addon, window, source._file+addon, window, source._file+addon, window, majorAxis,minorAxis,positionAngle, source._file+addon, window),[],logit=False)
    region = dict()
    region[0] = "%i,%i,%i,%i" % (int(math.floor(imsize/4)),int(math.floor((imsize/4) + (imsize/6))),int(math.floor((imsize/4) + (imsize/6))),int(math.floor((imsize/4) + (imsize/3))))
    region[1] = "%i,%i,%i,%i" % (int(math.floor((imsize/4) + (imsize/6))),int(math.floor((imsize/4) + (imsize/3))),int(math.floor((imsize/4) + (imsize/3))),int(math.floor(3*imsize/4)))
    region[2] = "%i,%i,%i,%i" % (int(math.floor((imsize/4) + (imsize/3))),int(math.floor((imsize/4) + (imsize/6))),int(math.floor(3*imsize/4)),int(math.floor((imsize/4) + (imsize/3))))
    region[3] = "%i,%i,%i,%i" % (int(math.floor((imsize/4) + (imsize/6))),int(math.floor(imsize/4)),int(math.floor((imsize/4) + (imsize/3))),int(math.floor((imsize/4) + (imsize/6))))
    rms = 1000000.0
    best = 0
    # determine the initial rms noise
    for j in range(4) :
        #log.writeLog("Locating best region for noise calculation.")
        log.run("imstat in=%s.w%i.finalmap region=box'('%s')' device=/NULL > rms.log.%s""" % (source._file+addon, window, region.get(j),fileEnd),[],logit=False)
        input = open("rms.log.%s" % (fileEnd),'r')
        fileList = input.readlines()
        input.close()
        fileList.reverse()
        while(len(fileList) > 0) :
            line = fileList.pop()
            if("Frequency" in line) :
                line = fileList.pop()
                temprms = float(line[42:51])
                if(temprms < rms) :
                    rms = temprms
                    best = j

    converged = False
    # get the clean region
    if(source.getCleanRegion() == None) :
        if(not source._mosaic) :
            if(p.preferences.get("doAutoCleanRegion")) :
                regions[window] = peakLocator.findCleanRegion("%s.w%i.finalmap" % (source._file+addon,window),rms,"%s.w%i.bm" % (source._file+addon,window),self._pointingCenters)
            elif(p.preferences.get("cleanRegion") == "quarter") :
                source.setCleanRegion("quarter")
            else :
                boxsize = int(p.preferences.get("cleanRegion"))
                source.setCleanRegion("arcsec,box'(-%f,-%f,%f,%f)'" % (boxsize,boxsize,boxsize,boxsize))
    # treat mosaics differently
        else :
            offset = [int(miriadClasses.imhead("%s.map" % source._file+addon,"crpix1")),int(miriadClasses.imhead("%s.map" % source._file+addon,"crpix2"))]
            source.setCleanRegion(calculations.calculateMosaicCleanRegion(source._pointingCenters,cellsize,imsize,offset))
    log.run("rm -rf %s.w%i.clean %s.w%i.finalmap %s.w%i.bm; sleep 3" % (source._file+addon, window, source._file+addon, window, source._file+addon, window),[],fatal=True,logit=False)

def cleanWindow(source,window,csub) :
    """ Method to clean a spectral window
        input :
            source - the source data
            window - which window to clean
            csub - has contuinuum subtraction been done
        returns :
            none
    """
    global rmsList
    lastIter = 10000.0
    addon = ""
    fileEnd = str(random.randint(1,100000))
    if(csub) :
        addon = ".cs"
    # do the initial clean
    log.run("mossdi map=%s.w%i.map beam=%s.w%i.beam out=%s.w%i.clean niters=100 region=%s" % (source._file+addon, window, source._file+addon, window, source._file+addon, window,source.getCleanRegion()),[],fatal=True,logit=False)
    # get the beam size
    log.run("mospsf beam=%s.w%i.beam out=%s.w%i.bm""" % (source._file+addon, window, source._file+addon, window),[],fatal=True,logit=False)
    log.run("imfit in=%s.w%i.bm object=beam region=arcsec,box'('-5,-5,5,5')' > fit.log.%s" % (source._file+addon, window,fileEnd),[],fatal=True,logit=False)
    input = open("fit.log.%s" % (fileEnd),'r')
    fileList = input.readlines()
    input.close()
    majorAxis = 0.0
    minorAxis = 0.0
    positionAngle = 0.0
    while(len(fileList) > 0) :
        line = fileList.pop()
        if("Position angle" in line) :
            splitLine = line.split()
            positionAngle = float(splitLine[3])
        elif("Minor axis" in line) :
            splitLine = line.split()
            minorAxis = float(splitLine[3])
        elif("Major axis" in line) :
            splitLine = line.split()
            majorAxis = float(splitLine[3])
    log.run("restor model=%s.w%i.clean map=%s.w%i.map beam=%s.w%i.beam fwhm=%f,%f pa=%f out=%s.w%i.finalmap" % (source._file+addon, window, source._file+addon, window, source._file+addon, window, majorAxis,minorAxis,positionAngle, source._file+addon, window),[],logit=False)
    region = dict()
    region[0] = "%i,%i,%i,%i" % (int(math.floor(imsize/4)),int(math.floor((imsize/4) + (imsize/6))),int(math.floor((imsize/4) + (imsize/6))),int(math.floor((imsize/4) + (imsize/3))))
    region[1] = "%i,%i,%i,%i" % (int(math.floor((imsize/4) + (imsize/6))),int(math.floor((imsize/4) + (imsize/3))),int(math.floor((imsize/4) + (imsize/3))),int(math.floor(3*imsize/4)))
    region[2] = "%i,%i,%i,%i" % (int(math.floor((imsize/4) + (imsize/3))),int(math.floor((imsize/4) + (imsize/6))),int(math.floor(3*imsize/4)),int(math.floor((imsize/4) + (imsize/3))))
    region[3] = "%i,%i,%i,%i" % (int(math.floor((imsize/4) + (imsize/6))),int(math.floor(imsize/4)),int(math.floor((imsize/4) + (imsize/3))),int(math.floor((imsize/4) + (imsize/6))))
    rms = 1000000.0
    best = 0
    # determine the rms noise
    for j in range(4) :
        log.run("imstat in=%s.w%i.finalmap region=box'('%s')' device=/NULL > rms.log.%s""" % (source._file+addon, window, region.get(j),fileEnd),[],logit=False)
        input = open("rms.log.%s" % (fileEnd),'r')
        fileList = input.readlines()
        input.close()
        fileList.reverse()
        while(len(fileList) > 0) :
            line = fileList.pop()
            if("Frequency" in line) :
                line = fileList.pop()
                temprms = float(line[42:51])
                if(temprms < rms) :
                    rms = temprms
                    best = j
    converged = False

    run = 1
    cutoff = 1000.0
    # loop over clean and restor until the noise level settles down or until 10 iterations are complete
    while((not converged) and (run <= 5)) :
        cutoff = rms * p.preferences.get("cleanThreshold")
        log.run("rm -rf %s.w%i.clean %s.w%i.finalmap; sleep 3" % (source._file+addon, window, source._file+addon, window),[],fatal=True,logit=False)
        log.run("mossdi map=%s.w%i.map beam=%s.w%i.beam out=%s.w%i.clean niters=50000 cutoff=%f region=%s" % (source._file+addon, window, source._file+addon, window, source._file+addon, window, cutoff,source.getCleanRegion()),[],fatal=True,logit=False)
        log.run("restor model=%s.w%i.clean map=%s.w%i.map beam=%s.w%i.beam fwhm=%f,%f pa=%f out=%s.w%i.finalmap" % (source._file+addon, window, source._file+addon, window, source._file+addon, window, majorAxis,minorAxis,positionAngle, source._file+addon, window),[],fatal=True,logit=False)
        log.run("imstat in=%s.w%i.finalmap region=box'('%s')' device=/NULL > rms.log.%s" % (source._file+addon, window, region.get(best),fileEnd),[],logit=False)
        input = open("rms.log.%s" % (fileEnd),'r')
        fileList = input.readlines()
        input.close()
        fileList.reverse()
        while(len(fileList) > 0) :
            line = fileList.pop()
            if("cube" in line) :
                line = fileList.pop()
                temprms = float(line[42:51])
                if(abs(temprms-rms) < 0.05 * rms) : # there must be less than a 5 % change
                    converged = True
                else :
                    rms = temprms
        run += 1
    log.writeComment("Cleaning window %i" % (window))
    if(run > 5) :
        log.writeComment("Did not reach noise level cutoff, performed 5 iterations.")

    args = []
    args.append(globals.Variable("map",source._file,addon + ".w%i.map" % (window)))
    args.append(globals.Variable("beam",source._file,addon + ".w%i.beam" % (window)))
    args.append(globals.Variable("out",source._file,addon + ".w%i.clean" % (window)))
    args.append(globals.Variable("niters","50000"))
    args.append(globals.Variable("cutoff",str(cutoff)))
    args.append(globals.Variable("region",source.getCleanRegion()))
    log.run("mossdi",args,execute=False)

    args = []
    args.append(globals.Variable("model",source._file,addon + ".w%i.clean" % (window)))
    args.append(globals.Variable("map",source._file,addon + ".w%i.map" % (window)))
    args.append(globals.Variable("beam",source._file,addon + ".w%i.beam" % (window)))
    args.append(globals.Variable("fwhm","%f,%f" % (majorAxis,minorAxis)))
    args.append(globals.Variable("pa",str(positionAngle)))
    args.append(globals.Variable("out",source._file,addon + ".w%i.finalmap" % (window)))
    log.run("restor",args,execute=False)

    args = []
    args.append(globals.Variable("in",source._file,addon + ".w%i.finalmap" % (window)))
    args.append(globals.Variable("op","xyout"))
    args.append(globals.Variable("out",source._file,addon + ".w%i.fits" % (window)))
    log.run("fits",args,logit=False)

    log.writeScript("\n")
    log.writeLog("Data reduction of %s window %i complete. Final map(%s.w%i.finalmap) has a noise level of %f Jy/beam\n" % (source._name, window, source._file+addon, window, rms))
    rmsList[window] = rms
    startupTeardown.endFile("%s.w%i.finalmap" % (source._file+addon,window))

def invertContinuum(objects, obsFreq, avgBaseline, windows, individual = False) :
    """ Method to create continuum maps
        input :
            objects - the objects
            obsFreq - observing frequency in GHz
            avgBaseline - average baseline in lambda
            windows - which windows to invert
            individual - whther we are only inverting individual windows
        returns :
            none
    """
    global cellsize
    global imsize
    global continDone
    fileEnd = str(random.randint(1,100000))
    # calculate the optimal cell size, based on the median baseline in klambda, want at least 5 pixels across the beam
    # but only if the user did not specify it in the preferences file
    if(p.preferences.get("cellSize") < 0 and cellsize == None) :
        cellsize = calculations.calcCellSize(avgBaseline,0.0)
        log.writeComment("Based on average baseline of %f,  using a cell size of %f arcseconds." % (avgBaseline, cellsize))
    elif(cellsize == None) :
        cellsize = p.preferences.get("cellSize")
        log.writeComment("Using cellsize %f" % (cellsize))
    iterate = True
    # calculate the optimal image size, based on the cell size and observing frequency
    # but only if the user did not specify an image size in the preferences file
    if(p.preferences.get("imageSize") < 0 and imsize == None) :
        imsize = calculations.calcImsize(obsFreq, cellsize)
    elif(imsize == None) :
        iterate = False
        imsize = p.preferences.get("imageSize")
    # invert each source in turn
    sourceList = deepcopy(objects._sources)
    if(not individual) :
        sourceList.append(deepcopy(objects._gaincals[0]))
    for source in sourceList :
        log.writeComment("Inverting continuum of %s" % (source._name))
        if(isinstance(source,sources.Source)) :
            if(source._mosaic) :
                if(imsize > 1000) :
                    imsize = int(imsize * 1.5)
                else :
                    imsize = int(imsize * 2.0)
                continDone = True
        invertString = ""
        invertList = []
        startChan = 1
        numChan = source._numChans[0]
        # make sure we start at the correct channel number
        if(source.haveSuper() and not individual) :
            for ending in ["LSB","USB"] :
                if(not source._lsbGood and ending == "LSB") :
                    continue
                if(not source._usbGood and ending == "USB") :
                    continue
                startChan = 1
                numChan = source.getSuperNumChans()
                invertString = invertString + ",%s.%s" % (source._file, ending)
                invertList.append([source._file,".%s" % (ending)])
                args = []
                args.append(globals.Variable("vis",source._file,".%s" % (ending)))
                args.append(globals.Variable("options","unflagged"))
                args.append(globals.Variable("out","temp." + fileEnd + "; sleep 3; rm -rf "))
                args.append(globals.Variable(None,source._file,".%s/*; rm -rf " % (ending)))
                args.append(globals.Variable(None,source._file,".%s; sleep 3; mv " % (ending)))
                args.append(globals.Variable(None,"temp." + fileEnd))
                args.append(globals.Variable(None,source._file,".%s" % (ending)))
                log.run("uvcat",args,fatal=True)
        elif(len(windows) != 0) :
            for window in windows :
                if(source._bandwidths[window - 1] == 62) :
                    startChan = 4
                    numChan = source._numChans[window - 1] - 6

                invertString = invertString + ",%s.w%i" % (source._file, window)
                invertList.append([source._file,".w%i" % (window)])
                args = []
                args.append(globals.Variable("vis",source._file,".w%i" % (window)))
                args.append(globals.Variable("options","unflagged"))
                args.append(globals.Variable("out","temp." + fileEnd + "; sleep 3; rm -rf "))
                args.append(globals.Variable(None,source._file,".w%i/*; rm -rf " % (window)))
                args.append(globals.Variable(None,source._file,".w%i; sleep 3; mv " % (window)))
                args.append(globals.Variable(None,"temp." + fileEnd))
                args.append(globals.Variable(None,source._file,".w%i" % (window)))
                log.run("uvcat",args,fatal=True)
        else :
            for window in range(globals.STARTWINDOW, globals.ENDWINDOW + 1) :
                if(source._bandwidths[window - 1] == 62) :
                    startChan = 4
                    numChan = source._numChans[window - 1] - 6
                invertString = invertString + ",%s.w%i" % (source._file, window)
                invertList.append([source._file,".w%i" % (window)])
                args = []
                args.append(globals.Variable("vis",source._file,".w%i" % (window)))
                args.append(globals.Variable("options","unflagged"))
                args.append(globals.Variable("out","temp." + fileEnd + "; sleep 3; rm -rf "))
                args.append(globals.Variable(None,source._file,".w%i/*; rm -rf " % (window)))
                args.append(globals.Variable(None,source._file,".w%i; sleep 3; mv " % (window)))
                args.append(globals.Variable(None,"temp." + fileEnd))
                args.append(globals.Variable(None,source._file,".w%i" % (window)))
                log.run("uvcat",args,fatal=True)

        invertString = invertString[1:]
        lastbad = 0
        noise = 0.5
        endings = [""]
        selects=[""]
        tempCell = cellsize
        tempImg = imsize
        if(globals.isSci2 and not individual) :
            endings.append(".short")
            if(globals.obsFreq() < 50.0) :
                selects.append(" select=uvrange'('0,2')'")
                selects.append(" select=uvrange'('2,2000000')'")
            else :
                selects.append(" select=uvrange'('0,6')'")
                selects.append(" select=uvrange'('6,6000000')'")
            endings.append(".long")
        # keep inverting until we minimize the number of rejected visibilities
        fileName = source._file
        for i in range(0,len(endings)) :
            if(endings[i].find(".short") >= 0) :
                cellsize = calculations.calcCellSize(globals.avgShortBaseline(),0.0)
                imsize = calculations.calcImsize(obsFreq, cellsize)
            elif(endings[i].find(".long") >= 0) :
                cellsize = calculations.calcCellSize(globals.avgLongBaseline(),0.0)
                imsize = calculations.calcImsize(obsFreq, cellsize)
            else :
                cellsize = tempCell
                imsize = tempImg
            fName = ""
            fEnd = ""
            ok = False
            end = endings[i]
            select = selects[i]
            if(individual) :
                fName = fileName
                fEnd = ".contin.w%i" % (windows[0])
                fileName += ".contin.w%i" % (windows[0])
            while(not ok) :
                log.run("invert vis=%s map=%s.map beam=%s.beam imsize=%i cell=%f sup=0%s options=systemp,double,mfs,mosaic line=channel,1,%i,%i >& junk" % (invertString, fileName + end, fileName + end,imsize,cellsize,select, startChan, numChan),[],fatal=True,logit=False)
                input = open("junk",'r')
                fileList = input.readlines()
                input.close()
                found = False
                while(len(fileList) > 0 and iterate) :
                    line = fileList.pop()
                    if("Visibilities rejected" in line) :
                        splitLine = line.split()
                        if(lastbad != int(splitLine[5])) :
                            cellsize = calculations.calcCellSize(avgBaseline,cellsize)
                            log.run("rm -rf %s %s; sleep 3" % (fileName + end + ".map", fileName + end + ".beam"),[],logit=False)
                            found = True
                        else :
                            cellsize *= 1.1
                            log.run("invert vis=%s map=%s.map beam=%s.beam imsize=%i cell=%f sup=0%s options=systemp,double,mfs,mosaic line=channel,1,%i,%i >& junk" % (invertString, fileName + end, fileName + end, imsize,cellsize, select, startChan, numChan),[],fatal=True,logit=False)
                    if("Theoretical rms" in line) :
                        noise = float(line.split(":")[1])
                if(not found) :
                    ok = True
                    log.writeLog("Using an image size of %f for %s" % (imsize, source._name))
                    args = []
                    args.append(globals.Variable("vis",invertList[0][0],invertList[0][1]))
                    temp = invertList.pop(0)
                    for vis in invertList :
                        args.append(globals.Variable("ADD",vis[0],vis[1]))
                    invertList.append(temp)

                    if(individual) :
                        args.append(globals.Variable("map",fName,fEnd + ".map"))
                        args.append(globals.Variable("beam",fName,fEnd + ".beam"))
                    else :
                        args.append(globals.Variable("map",source._file,end + ".map"))
                        args.append(globals.Variable("beam",source._file,end + ".beam"))
                    args.append(globals.Variable("imsize",str(imsize)))
                    args.append(globals.Variable("cell",str(cellsize)))
                    args.append(globals.Variable("sup","0"))
                    if(select != "") :
                        args.append(globals.Variable(None,select))
                    args.append(globals.Variable("options","systemp,double,mfs,mosaic"))
                    args.append(globals.Variable("line","channel,1,%i,%i" % (startChan, numChan)))
                    log.run("invert",args,execute=False)
            lastIter = 10000.0
            # treat multipoint mosics differently
            log.run("mossdi map=%s.map beam=%s.beam out=%s.clean niters=50000 cutoff=%f region=quarter" % (fileName + end, fileName + end, fileName + end,noise*2.5),[],fatal=True,logit=False)
            # calculate the correst synthesized beam
            log.run("mospsf beam=%s.beam out=%s.bm" % (fileName + end, fileName + end),[],fatal=True,logit=False)
            log.run("imfit in=%s.bm object=beam region=arcsec,box'('-5,-5,5,5')' > fit.log.%s" % (fileName + end,fileEnd),[],fatal=True,logit=False)
            input = open("fit.log.%s" % (fileEnd),'r')
            fileList = input.readlines()
            input.close()
            majorAxis = 0.0
            minorAxis = 0.0
            positionAngle = 0.0
            # get the beam parameters
            while(len(fileList) > 0) :
                line = fileList.pop()
                if("Position angle" in line) :
                    splitLine = line.split()
                    positionAngle = float(splitLine[3])
                elif("Minor axis" in line) :
                    splitLine = line.split()
                    minorAxis = float(splitLine[3])
                elif("Major axis" in line) :
                    splitLine = line.split()
                    majorAxis = float(splitLine[3])
            log.run("restor model=%s.clean map=%s.map beam=%s.beam fwhm=%f,%f pa=%f out=%s.finalmap" % (fileName + end, fileName + end, fileName + end, majorAxis,minorAxis,positionAngle, fileName + end),[],logit=False)
            # get the initial noise level
            region = dict()
            region[0] = "%i,%i,%i,%i" % (int(math.floor(imsize/4)),int(math.floor((imsize/4) + (imsize/6))),int(math.floor((imsize/4) + (imsize/6))),int(math.floor((imsize/4) + (imsize/3))))
            region[1] = "%i,%i,%i,%i" % (int(math.floor((imsize/4) + (imsize/6))),int(math.floor((imsize/4) + (imsize/3))),int(math.floor((imsize/4) + (imsize/3))),int(math.floor(3*imsize/4)))
            region[2] = "%i,%i,%i,%i" % (int(math.floor((imsize/4) + (imsize/3))),int(math.floor((imsize/4) + (imsize/6))),int(math.floor(3*imsize/4)),int(math.floor((imsize/4) + (imsize/3))))
            region[3] = "%i,%i,%i,%i" % (int(math.floor((imsize/4) + (imsize/6))),int(math.floor(imsize/4)),int(math.floor((imsize/4) + (imsize/3))),int(math.floor((imsize/4) + (imsize/6))))
            rms = 1000000.0
            best = 0
            for j in range(4) :
                #log.writeComment("Locating best region for noise calculation.")
                log.run("imstat in=%s.finalmap region=box'('%s')' device=/NULL > rms.log.%s" % (fileName + end, region.get(j),fileEnd),[],logit=False)
                input = open("rms.log.%s" % (fileEnd),'r')
                fileList = input.readlines()
                input.close()
                fileList.reverse()
                while(len(fileList) > 0) :
                    line = fileList.pop()
                    if("Frequency" in line) :
                        line = fileList.pop()
                        temprms = float(line[42:51])
                        if(temprms < rms) :
                            rms = temprms
                            best = j
        #log.writeComment("Found region %i is best: %s" % (best,  region.get(best)))
        # determine the clean region either from user input, number of pointings, or automatically from the data
            cleanReg = ""
            if(isinstance(source,sources.Source)) :
                if(source.getContinuumCleanRegion() == None) :
                    if(not source._mosaic) :
                        if(p.preferences.get("doAutoCleanRegion")) :
                            source.setContinuumCleanRegion(peakLocator.findCleanRegion("%s.finalmap" % (fileName + end),rms,"%s.bm" % (fileName),source._pointingCenters))
                            if(source.getContinuumCleanRegion() == "quarter") :
                                source.setCleanRegion("quarter")
                        elif(p.preferences.get("cleanRegion") == "quarter") :
                            source.setContinuumCleanRegion("quarter")
                            source.setCleanRegion("quarter")
                        else :
                            boxsize = int(p.preferences.get("cleanRegion"))
                            source.setContinuumCleanRegion("arcsec,box'(-%f,-%f,%f,%f)'" % (boxsize,boxsize,boxsize,boxsize))
                            source.setCleanRegion("arcsec,box'(-%f,-%f,%f,%f)'" % (boxsize,boxsize,boxsize,boxsize))
                    else :
                        if(p.preferences.get("doAutoCleanRegion")) :
                            print "AUTOCLEAN+++++++++++++++++++++++++++++"#################POINTING CENTERS
                            print source._pointingCenters
                            source.setContinuumCleanRegion(peakLocator.findCleanRegion("%s.finalmap" % (fileName + end),rms,"%s.bm" % (fileName + end),source._pointingCenters))
                            source.setCleanRegion(source.getContinuumCleanRegion())
                        else :
                            offset = [int(miriadClasses.imhead("%s.map" % fileName+end+addon,"crpix1")),int(miriadClasses.imhead("%s.map" % fileName+end+addon,"crpix2"))]
                            source.setContinuumCleanRegion(calculations.calculateMosaicCleanRegion(source._pointingCenters,cellsize,obsFreq,offset))
                            source.setCleanRegion(source.getContinuumCleanRegion())
                        log.writeLog("mosaic  " + source.getCleanRegion())
                cleanReg = source.getCleanRegion()
            else :
                # we are mapping a calibrator
                cleanReg = "quarter"
            converged = False
            log.writeComment("Cleaning image to %f times the noise level, or 100000 iterations, whichever comes first" % (p.preferences.get("cleanThreshold")))
            run = 1
            cutoff = 1000.0
            # loop over clean and restor until the noise level settles down or until 5 iterations are done
            while((not converged) and (run <= 5)) :
                cutoff = rms * p.preferences.get("cleanThreshold")
                log.run("rm -rf %s.clean %s.finalmap; sleep 3" % (fileName + end, fileName + end),[],fatal=True,logit=False)
                log.run("mossdi map=%s.map beam=%s.beam out=%s.clean niters=50000 cutoff=%f region=%s" % (fileName + end, fileName + end, fileName + end, cutoff,cleanReg),[],fatal=True,logit=False)
                log.run("restor model=%s.clean map=%s.map beam=%s.beam fwhm=%f,%f pa=%f out=%s.finalmap" % (fileName + end,fileName + end, fileName + end, majorAxis,minorAxis,positionAngle,fileName + end),[],fatal=True,logit=False)
                log.run("imstat in=%s.finalmap region=box'('%s')' device=/NULL > rms.log.%s" % (fileName + end, region.get(best),fileEnd),[],logit=False)
                input = open("rms.log.%s" % (fileEnd),'r')
                fileList = input.readlines()
                input.close()
                fileList.reverse()
                while(len(fileList) > 0) :
                    line = fileList.pop()
                    if("Frequency" in line) :
                        line = fileList.pop()
                        temprms = float(line[42:51])
                        if(abs(temprms-rms) < 0.05 * rms) : # there must be less than a 5 % change
                            converged = True
                        else :
                            rms = temprms
                run += 1
            if(run > 5) :
                log.writeComment("Did not reach noise level cutoff, performed 5 iterations.")
            args = []
            if(individual) :
                args.append(globals.Variable("map",fName,fEnd + ".map"))
                args.append(globals.Variable("beam",fName,fEnd + ".beam"))
                args.append(globals.Variable("out",fName,fEnd + ".clean"))
            else :
                args.append(globals.Variable("map",source._file,end + ".map"))
                args.append(globals.Variable("beam",source._file,end + ".beam"))
                args.append(globals.Variable("out",source._file,end + ".clean"))
            args.append(globals.Variable("niters","50000"))
            args.append(globals.Variable("cutoff",str(cutoff)))
            args.append(globals.Variable("region",cleanReg))
            log.run("mossdi",args,execute=False)

            args = []
            if(individual) :
                args.append(globals.Variable("model",fName,fEnd + ".clean"))
                args.append(globals.Variable("map",fName,fEnd + ".map"))
                args.append(globals.Variable("beam",fName,fEnd + ".beam"))
                args.append(globals.Variable("out",fName,fEnd + ".finalmap"))
            else :
                args.append(globals.Variable("model",source._file,end + ".clean"))
                args.append(globals.Variable("map",source._file,end + ".map"))
                args.append(globals.Variable("beam",source._file,end + ".beam"))
                args.append(globals.Variable("out",source._file,end + ".finalmap"))
            args.append(globals.Variable("fwhm","%f,%f" % (majorAxis,minorAxis)))
            args.append(globals.Variable("pa",str(positionAngle)))
            log.run("restor",args,execute=False)

            args = []
            if(individual) :
                args.append(globals.Variable("in",fName,fEnd + ".finalmap"))
                args.append(globals.Variable("out",fName,fEnd + ".fits"))
            else :
                args.append(globals.Variable("in",source._file,end + ".finalmap"))
                args.append(globals.Variable("out",source._file,end + ".fits"))
            args.append(globals.Variable("op","xyout"))
            log.run("fits",args,logit=False)


            log.writeComment("Data reduction of %s complete. Final map(%s.finalmap) has a noise level of %f Jy/beam" % (source._name, fileName + end, rms))
            startupTeardown.endFile("%s.finalmap" % (fileName + end))

def invertSpectra(objects, obsFreq, avgBaseline, windows) :
    """ Method to create spectral line channel maps
        input:
            objects - the objetcs
            obsFreq - observing frequency in GHz
            avgBaseline - average baseline in lambda
            windows - which windows to map
        returns :
            none
    """
    global cellsize
    global imsize
    global continDone
    global regions
    global rmsList
    # calculate the optimal cell size, based on the median baseline in klambda, want at least 5 pixels across the beam
    # but only if the user did not specify it in the preferences file
    if(p.preferences.get("cellSize") < 0 and cellsize == None) :
        cellsize = calculations.calcCellSize(avgBaseline,0.0)
        log.writeComment("Based on average baseline of %f,  using a cell size of %f arcseconds." % (avgBaseline, cellsize))
    elif(cellsize == None) :
        cellsize = p.preferences.get("cellSize")
    iterate = True
    # calculate the optimal image size from the cell size and observing frequency
    # but only if it was not specified by the user in the preferences file
    if(p.preferences.get("imageSize") < 0 and imsize == None) :
        imsize = calculations.calcImsize(obsFreq, cellsize)
    elif(imsize == None) :
        iterate = False
        imsize = p.preferences.get("imageSize")
    for source in objects._sources :
        if(source._mosaic and not continDone) :
            if(imsize > 1000) :
                imsize = int(imsize * 1.5)
            else :
                imsize = int(imsize * 2.0)
        threadList = []
        for window in windows :
            current = invertWindowThread(source,window,iterate,AvgBaseline=avgBaseline)
            threadList.append(current)
            current.start()
        for thread in threadList :
            thread.join()
        if(len(regions) > 0) :
            # get clean region
            regions[0] = source.getContinuumCleanRegion()
            source.setCleanRegion(peakLocator.compactCleanRegions(regions,[imsize,imsize]))
        del threadList[:]
        log.writeComment("Cleaning images to %f times the noise level, or 100000 iterations, whichever comes first" % (p.preferences.get("cleanThreshold")))
        log.writeComment("Note that these are threaded and may appear out of window order")
        for window in windows :
            current = cleanWindowThread(source,window)
            threadList.append(current)
            current.start()
        for thread in threadList :
            thread.join()
        regions.clear()
        # now do continuum subtraction
        if(p.preferences.get("doContinuumSubtraction")) :
            del threadList[:]
            for window in windows :
                current = continuumSubtraction.continuumSubtractionThread(source, window, 3.0*rmsList[window])
                threadList.append(current)
                current.start()
            for thread in threadList :
                thread.join()
            del threadList[:]
            for window in windows :
                current = cleanWindowThread(source,window,True)
                threadList.append(current)
                current.start()
            for thread in threadList :
                thread.join()
            del threadList[:]
            for window in windows :
                current = invertWindowThread(source,window,True)
                threadList.append(current)
                current.start()
            for thread in threadList :
                thread.join()
