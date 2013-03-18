import peakLocator
import logger as log
import calculations
import sources
import math
import globals
from copy import deepcopy
from threading import Thread

class ContinuumSubtractionThread(Thread) :
    def __init__(self,source, window, cutoff) :
        Thread.__init__(self)
        self.source = source
        self.window = window
        self.cutoff = cutoff
        self.line = line
    def run(self) :
        chans = getContinuumChannels(self.souce._file + ".w%i.final" %(window),self.source._numChans[window],self.souce._file + ".w%i.bm" %(window),self.cutoff,self._pointingCenters)
        if(chans != None) :
            removeContinuum(self.souce._file, chans,window)

def removeContinuum(file,chans,window) :
    startChan = 1
    numChan = 1
    if(source._bandwidths[window - 1] == 62) :
        startChan = 4
        numChan = source._numChans[window - 1] - 6
    if(startChan != 1) :
        args = []
        args.append(globals.Variable("vis",file,".w%i" % (window)))
        args.append(globals.Variable("line","chan,%i,%i,1,1" % (numChan,startChan)))
        args.append(globals.Variable("order","0"))
        args.append(globals.Variable("chans",chans))
        args.append(globals.Variable("out",file,".cs.w%i" % (window)))
        log.run("uvlin",args)
#        log.run("uvlin vis=%s.w%i line=chan,%i,%i,1,1 order=0 chans=%s out=%s.cs.w%i" % (file,window,numChan,startChan,chans,file,window))
    else :
        args = []
        args.append(globals.Variable("vis",file,".w%i" % (window)))
        args.append(globals.Variable("order","0"))
        args.append(globals.Variable("chans",chans))
        args.append(globals.Variable("out",file,".cs.w%i" % (window)))
        log.run("uvlin",args)
#        log.run("uvlin vis=%s.w%i order=0 chans=%s out=%s.cs.w%i" % (file,window,chans,file,window))

def avg(list) :
    """Method to average lists of values together
        input :
            list - a list of lists
        returns :
            a list of the averages of the lists (element for element)
    """
    average = dict()
    num = len(list)
    for i in range(0,len(list[0])) :
        value = 0.0
        for j in range(0,num) :
            value += list[j][i]
        average[i + 1] = value/num
    return average

def getSpec(image,point) :
    """Method to get the spectra from a cube, takes the spectra of 9 pixels and averages them
        inputs :
            image - an image cube (rotated to x,y,z)
            point - the point to take spectra from
        returns :
            dictionary of the spectra
    """
    spec = []
    for i in range(-1,2) :
        for j in range(-1,2) :
            x = point[0] + i
            y = point[1] + j
            if(x < 0 or y < 0 or x >= image.shape[0] or y >= image.shape[1]) :
                continue
            spec.append(image[x][y])
    return avg(spec)

def simplifyRegions(regions) :
    """ Method which takes a list of points and combines regions that overlap
        input :
            regions - a list of central points for each region
        returns :
            a new list of central points with overlaps removed
    """
    newRegions = []
    for r in regions :
        if(not r in newRegions) :
            newRegions.append(r)
    regions = newRegions
    while(True) :
        overlap = dict()
        newRegions = []
        remove = []
        for r in range(0,len(regions)) :
            overlap[r] = []
            for r2 in range(r+1,len(regions)) :
                if(calculations.distance(regions[r],regions[r2]) <= 1.5 * math.sqrt(2.0)) :
                    if(calculations.distance(regions[r],regions[r2]) <= 1.1 * math.sqrt(2.0)) :
                        overlap[r].append([-regions[r2][0],-regions[r2][1]])
                    else :
                        overlap[r].append(regions[r2])
            if(len(overlap[r]) == 0 and not regions[r] in remove and not regions[r] in newRegions) :
                newRegions.append(regions[r])
            else :
                for r2 in overlap[r] :
                    if(r2[0] > 0 and r2[1] > 0) :
                        mdPoint = calculations.intMidpoint(regions[r],r2)
                        if(not mdPoint in newRegions) :
                            newRegions.append(mdPoint)
                            remove.append(r2)
                    else :
                        if(not regions[r] in newRegions) :
                            newRegions.append(regions[r])
                            remove.append([abs(r2[0]),abs(r2[1])])

        # if no points were eliminated the escape
        if(len(regions) == len(newRegions)) :
            break
        regions = newRegions
    return newRegions

def getContinuumChannels(file,nchan,beamFile,cutoff,pointings) :
    """Method to determine which channels are continuum free
        inputs :
            file - the image data cube to look at
            nchan - the number of channels in the image
            beamFile - the synthesized beam
            cutoff - the rms cutoff level
    """
    continChans = ""
    chanList = [True] * nchan
    # initialize an array of continuum channels
    regions,image = peakLocator.findSpectralPeaks(file,cutoff,beamFile,pointings)
    print image[62][57][3]
    #print regions
    regions = simplifyRegions(regions)
    image = image.swapaxes(0,2) # rotate the cube
    for point in regions : ###DO THIS IN PYTHON
        chanDict = getSpec(image,point)
        foundPeak = True
        goodFit = False
        # determine which channels are continuum only
        while(foundPeak and not goodFit) :
            average = calculations.listAverage(chanDict)
            foundPeak = False
            above = 0.00001
            below = 0.00001
            for i in chanDict :
                if(chanDict[i] == -10000.0) :
                    continue
                if(abs(chanDict[i]) - cutoff >= average) : ### this may fail for very strong, narrow lines
                    chanDict[i] = -10000.0
                    chanList[i-1] = False
                    foundPeak = True
                elif(chanDict[i] > average) :
                    above += 1.0
                else :
                    below += 1.0
            if((abs(1.0 - (above/below)) < 0.05) or (abs(above - below) < 3)) :
                goodFit = True # we have a decent fit, so break out of the loop
    # if there was no peak found then break out anyway since that is the best we can do
        count = 0.0
        firstChan = False
        continChans = []
    for i in range(0,nchan) :
        if(chanList[i]) :
            count += 1.0
            # if we have not found a starting True channel
            if(not firstChan) :
                firstChan = True
                continChans.append(i + 1)
            if(i == nchan - 1) :
                continChans.append(i + 1)
            #print continChans
        else :
            if(firstChan) :
                firstChan = False
                continChans.append(i)
        #print continChans
    temp = deepcopy(continChans)
    for i in range(len(continChans) - 2,-1,-1) :
        if(i > len(continChans) - 2) :
            print i,len(continChans),"continue"
            continue
        print i,len(continChans)
        if(continChans[i] == continChans[i+1]) :
            value = continChans[i]
            continChans.remove(value)
            continChans.remove(value)
    if(len(continChans) == 0) :
        continChans = temp
    #print continChans,1.0 - (count/float(nchan)),count,nchan
    if(count/float(nchan) > 0.08) :
        return calculations.listToString(continChans)
    else :
        # not enough channels to reliably subtract the continuum
        return None

def continuumSubtraction(source,windows,cutoff) :
    threadList = []
    for window in windows :
        current = ContinuumSubtractionThread(source,window,cutoff[window-1])
        threadList.append(current)
        current.start()
    for thread in threadList :
        thread.join()

