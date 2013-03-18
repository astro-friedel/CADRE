import sources
import logger as log
import calculations
import globals
from pipeline_miriadwrap import *
import math
import flagging
try :
    import preferences as p
except ImportError:
    import defaultPreferences as p

objectPurposes = dict()

"""
Module for listobs operations
Part of the CARMA data reduction pipeline
Current Version (this file): 1.1
Released with pipeline version: 1.1
Author: D. N. Friedel
"""

def addOutput(name,purpose) :
    """ Method to add the purpose to the global dictionary
        input :
            name - the name of the object
            purpose - the purpose to add
        returns :
            none
    """
    global objectPurposes
    if(name in objectPurposes) :
        if(not (purpose in objectPurposes[name])) :
            objectPurposes[name] += purpose
    else :
        objectPurposes[name] = purpose


def orderDistance(antpos) :
    """ Method to sort the antennas based on distance from each other
        input :
            antpos - a dictionary with the antenna number as the key and the x,y,z coordinates as the value
        returns :
            a list of the antennas ordered with the most central antenna first
    """
    distances = dict()
    for k,v in antpos.iteritems() :
        avgDist = 0.0
        numAnts = 0.0
        for k1,v1 in antpos.iteritems() :
            if(k1 != k) :
                avgDist += math.sqrt(math.pow(v[0]-v1[0],2) + math.pow(v[1]-v1[1],2) + math.pow(v[2]-v1[2],2))
                numAnts += 1.0
        avgDist /= numAnts
        distances[k] = avgDist
    orderedAnts = []
    while(len(distances) > 0) :
        for k,v in distances.iteritems() :
            minDist = v
            minAnt = k
            for k1,v1 in distances.iteritems() :
                if(k1 != k) :
                    if(v1 <= v) :
                        minDist = v1
                        minAnt = k1
            orderedAnts.append(minAnt)
            del distances[minAnt]
            break
    return orderedAnts

def getBandwidths(bw) :
    """ Method to convert the bandwidth string to a list
        input :
            bw - the string containing bandwidth info
        returns :
             a list of the bandwidths (in window order)
    """
    bwlist = bw.split("-")
    for i in range(len(bwlist) - 1,-1,-1) :
        if("X" in bwlist[i] or "x" in bwlist[i]) :
            del bwlist[i]
        else :
            if("b" in bwlist[i]) :
                bwlist[i] = bwlist[i][:-2]
            bwlist[i] = int(bwlist[i])
    for i in range(0,len(bwlist)) :
        bwlist.append(bwlist[i])
    return bwlist

def addObject(obj,name,channels,bandwidths,numChans,ra,dec) :
    """ Method to add an object to the objects data
        input :
            obj - the objects data handle
            name - the name of the object to be added
            channels - the channels of the new object
            bandwidths - the bandwidths of the new object
            numChans - the number of channels of the new object
            ra,dec - the coordinates of the new object
        returns :
            none
    """
    obj.setName(name)
    obj.setCoordinates(ra,dec)
    obj.setFile(name)
    obj.setBandwidths(getBandwidths(bandwidths))
    obj.setChannels(channels)
    obj.setNumChan(numChans)

def addSource(objects,name,channels,bandwidths,numChans,ra,dec,mosaic,pointings) :
    """Method to add a source to the objects data
        input :
            obj - the objects data handle
            name - the name of the object to be added
            channels - the channels of the new object
            bandwidths - the bandwidths of the new object
            numChans - the number of channels of the new object
            ra,dec - the coordinates of the new object
            mosiac - is the data a posaic pattern
            pointings - list of mosaic pointings
        returns :
            none
    """
    addOutput(name,"S")
    source = objects.getSource(name)
    # if the source doesn't exist
    if(source == None) :
        source = sources.Source()
        globals.setScriptVar(name,"SOURCE")
        # add a new one
        addObject(source,name,channels,bandwidths,numChans,ra,dec)
        if(mosaic) :
            source.setMosaic(mosaic,pointings)
        objects.addSource(source)
    elif(mosaic) :
        source.setMosaic(mosaic,pointings)
        objects.updateSource(source)

def addGainCal(objects,name,channels,bandwidths,numChans,ra,dec,start,end) :
    """ Method to add a gain calibrator to the objects data
        input :
            obj - the objects data handle
            name - the name of the object to be added
            channels - the channels of the new object
            bandwidths - the bandwidths of the new object
            numChans - the number of channels of the new object
            ra,dec - the coordinates of the new object
            start,end - lists of the start/end times for each cycle
        returns :
            none
    """
    if(sources.isPlanet(name)) :
        log.writeComment("Ignoring %s as a gain calibrator" % (name))
        return
    addOutput(name,"G")
    gcal = objects.getGaincal(name)
    # if the source doesn't exist
    if(gcal == None) :
        gcal = sources.Gaincal()
        globals.setScriptVar(name,"GAINCAL")
        # add a new one
        addObject(gcal,name,channels,bandwidths,numChans,ra,dec)
        gcal.setTimes(start,end)
        objects.addGaincal(gcal)
    else :
        # append to the observed times
        gcal.setTimes(start,end)
        objects.updateGaincal(gcal)

def addPassCal(objects,name,channels,bandwidths,numChans,ra,dec,freqs) :
    """ Method to add a passband calibrator to the objects data
        input :
            obj - the objects data handle
            name - the name of the object to be added
            channels - the channels of the new object
            bandwidths - the bandwidths of the new object
            numChans - the number of channels of the new object
            ra,dec - the coordinates of the new object
            freqs - list of frequencies for each window
        returns :
            none
    """
    addOutput(name,"B")
    pcal = objects.getPasscal(name)
    # if the passband cal doesn't exist
    if(pcal == None) :
        pcal = sources.Passcal()
        globals.setScriptVar(name,"PASSCAL")
        # add a new one
        addObject(pcal,name,channels,bandwidths,numChans,ra,dec)
        bws = getBandwidths(bandwidths)
        for i in range(0,len(bws)):
            pcal.setHybridConf(bws[i])
            pcal.setHybridChans(bws[i],numChans[i])
        if("NOISE" in name) :
            pcal.setType(sources.NOISE)
        elif(sources.isPlanet(name)) :
            pcal.setType(sources.PLANET)
        else :
            pcal.setType(sources.QUASAR)
        pcal.setFreq(freqs)
        objects.addPasscal(pcal)
    else :
        bws = getBandwidths(bandwidths)
        if(bws != pcal.getBandwidths()) :
            pcal.setHybrid(True)
            for i in range(0,len(bws)):
                pcal.setHybridConf(bws[i])
                pcal.setHybridChans(bws[i],numChans[i])
            objects.updatePasscal(pcal)

def addFluxCal(objects,name,channels,bandwidths,numChans,ra,dec) :
    """ Method to add a flux calibrator to the objects data
        input :
            obj - the objects data handle
            name - the name of the object to be added
            channels - the channels of the new object
            bandwidths - the bandwidths of the new object
            numChans - the number of channels of the new object
            ra,dec - the coordinates of the new object
        returns :
            none
    """
    addOutput(name,"F")
    fcal = objects.getFluxcal(name)
    # if the flux cal doesn't exist
    if(fcal == None) :
        fcal = sources.Fluxcal()
        globals.setScriptVar(name,"FLUXCAL")
        # add a new one
        addObject(fcal,name,channels,bandwidths,numChans,ra,dec)
        if(sources.isPlanet(name)) :
            fcal.setType(sources.PLANET)
            sources.fluxcal[name] = True
        elif("MWC349" in name) :
            fcal.setType(sources.QUASAR)
            sources.fluxcal[name] = True
        else :
            fcal.setType(sources.QUASAR)
        objects.addFluxcal(fcal)

def addPolCal(objetcs,name,channels,bandwidths,numChans,ra,dec) :
    """ Method to add a flux calibrator to the objects data
        currently this is a placeholder until polarization is supported
        input :
            obj - the objects data handle
            name - the name of the object to be added
            channels - the channels of the new object
            bandwidths - the bandwidths of the new object
            numChans - the number of channels of the new object
            ra,dec - the coordinates of the new object
        returns :
            none
    """
    addOutput(name,"P")
    globals.setScriptVar(name,"POLCAL")

def runListobs(visFile, objects) :
    """ Method to run the miriad listobs command and then gather useful information from the output
        input :
            visFile - the main visibility file
            objects - the objects object
        returns :
            none
    """
    global objectPurposes
    handle = 0
    antpos = dict()
    oldSource = ""
    purpose = ""
    currentSource = ""
    channels = []
    bandwidths = ""
    numChans = []
    ra = 0.0
    dec = 0.0
    start = 0.0
    end = 0.0
    mosaic = False
    pointings = []
    freqs = []
    nspect = 0
    baselines = []
    # open the file
    try :
        handle = uvopen(visFile,"old")
    except :
        raise Exception,  "Could not open file: %s" % (visFile)
    # read in the preamble and intial variables
    preamble,data,flags = uvread(handle)
    nants = uvgetvri(handle,"nants",1)
    nspect = uvgetvri(handle,"nspect",1)
    goodAnts = [False] * (nants + 1)
    goodSys = [True] * (nants + 1)
    linefreq = uvgetvrd(handle,"freq",1)[0]
    positions = uvgetvrd(handle,"antpos",nants*3)
    # get the antenna positions
    for i in range(0,nants) :
        antpos[i+1] = [positions[i],positions[i + nants],positions[i + (nants * 2)]]
    for i in range(1,nants-1) :
        for j in range(i+1,nants) :
            baselines.append(math.sqrt(math.pow(linefreq*(antpos[i][1] - antpos[j][1]),2.) + math.pow(linefreq*(antpos[i][2] - antpos[j][2]),2.)))

    obsFreq = uvgetvrd(handle,'lo1',1)[0]
    obsDate = julday(preamble[2],"D")[:7]
    # start tracking variables
    uvrewind(handle)
    uvtrack(handle,"source",'u')
    uvtrack(handle,"purpose",'u')
    uvtrack(handle,"modedesc",'u')
    # scan the files
    while(uvscan(handle,'coord') == 0) :
        blnum = int(uvgetvrr(handle,"baseline",1)[0])
        ant1,ant2 = basant(blnum)
        goodAnts[ant1] = True
        goodAnts[ant2] = True
        systemp = uvgetvrr(handle,"systemp",nants*nspect)
        for i in range(0,nspect) :
            goodSys[ant1] = (systemp[nants*i + ant1 - 1] < p.preferences.get("tsysThreshold")) and goodSys[ant1]
            goodSys[ant2] = (systemp[nants*i + ant2 - 1] < p.preferences.get("tsysThreshold")) and goodSys[ant2]
        if(uvupdate(handle)) :
            # add new objects, or append if they exist
            if("S" in purpose) :
                addSource(objects,currentSource,channels,bandwidths,numChans,ra,dec,mosaic,pointings)
            if("G" in purpose) :
                addGainCal(objects,currentSource,channels,bandwidths,numChans,ra,dec,start,end)
            if("B" in purpose) :
                addPassCal(objects,currentSource,channels,bandwidths,numChans,ra,dec,freqs)
            if("F" in purpose) :
                addFluxCal(objects,currentSource,channels,bandwidths,numChans,ra,dec)
            channels = []
            bandwidths = ""
            numChans = []
            ra = 0.0
            dec = 0.0
            mosaic = False
            pointings = []
            start = []
            end = []
            freqs = []
            currentSource = uvgetvra(handle,"source")
            purpose = uvgetvra(handle,"purpose")
            bandwidths = uvgetvra(handle,"modedesc")
            ra = uvgetvrd(handle,"ra",1)[0]
            dec = uvgetvrd(handle,"dec",1)[0]
            nspect = uvgetvri(handle,"nspect",1)
            channels = uvgetvri(handle,"ischan",nspect)
            numChans = uvgetvri(handle,"nschan",nspect)
            # handle the individual object purposes (an object may have more than one purpose)
            if("S" in purpose) :
                dra = uvrdvrr(handle,"dra",0.0)
                ddec = uvrdvrr(handle,"ddec",0.0)
                if(dra != 0.0 and ddec != 0.0) :
                    if(not([dra,ddec] in pointings)) :
                        mosiac = True
                        pointings.append([dra,ddec])
            if("G" in purpose) :
                start = uvgetvrd(handle,"ut",1)[0]/(15.0*sources.degToRad)
            if("B" in purpose) :
                freqList = uvinfo(handle,"frequency")
                for i in range(0,nspect) :
                    if(freqList[channels[i] - 1] < freqList[channels[i] - 1 + numChans[i]]) :
                        freqs.append(freqList[channels[i] - 1])
                    else :
                        freqs.append(freqList[channels[i] - 1 + numChans[i]])
        if("G" in purpose):
            end = uvgetvrd(handle,"ut",1)[0]/(15.0*sources.degToRad) + uvgetvrr(handle,"inttime",1)[0]/3600.0
    if("S" in purpose) :
        addSource(objects,currentSource,channels,bandwidths,numChans,ra,dec,mosaic,pointings)
    if("G" in purpose) :
        addGainCal(objects,currentSource,channels,bandwidths,numChans,ra,dec,start,end)
    if("B" in purpose) :
        addPassCal(objects,currentSource,channels,bandwidths,numChans,ra,dec,freqs)
    if("F" in purpose) :
        addFluxCal(objects,currentSource,channels,bandwidths,numChans,ra,dec)

    avgBaseline = 0.0
    shortBaseline = 0.0
    sCounter = 0
    longBaseline = 0.0
    lCounter = 0
    cutoff = 2000.0
    if(obsFreq > 50.0) :
        cutoff = 6000.0
    for i in baselines :
        avgBaseline += i
        if(i < cutoff) :
            shortBaseline += i
            sCounter += 1
        else :
            longBaseline += i
            lCounter += 1
    avgBaseline /= len(baselines)
    if(sCounter > 0) :
        shortBaseline /= sCounter
    if(lCounter > 0) :
        longBaseline /= lCounter

    orderedAnts = orderDistance(antpos)
    # locate the best reference antenna
    refant = 9
    if(globals.isSci2) :
            refant = 20
    for i in orderedAnts :
        if(goodAnts[i] and goodSys[i]) :
            refant = i
            break
    globals.setScriptVar(str(refant),"REFANT")
    globals.setObsInfo(obsFreq, obsDate, avgBaseline, antpos, refant, shortBaseline, longBaseline)
    # print out what was found
    log.writeComment("Found the following source information")
    for key in objectPurposes :
        if("S" in objectPurposes[key]):
            log.writeComment("\t" + key + "\t\tSource")
        if("G" in objectPurposes[key]):
            log.writeComment("\t" + key + "\t\tGain Calibrator")
        if("B" in objectPurposes[key]):
            log.writeComment("\t" + key + "\t\tBandpass Calibrator")
        if("F" in objectPurposes[key]):
            log.writeComment("\t" + key + "\t\tFlux Calibrator")
        if("P" in objectPurposes[key]):
            log.writeComment("\t" + key + "\t\tPOlarization Calibrator")
