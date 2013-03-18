import logger as log
import calculations
import hybrid
import os
import decorrelation
import random
import math
import globals

"""
Module for source classes and methods
Part of the CARMA data reduction pipeline
Current Version (this file): 1.1
Released with pipeline version: 1.1
Author: D. N. Friedel
"""

# define some common constants
PLANET = "planet"
NOISE = "noise"
QUASAR = "quasar"
PRIMARY = 1
SECONDARY = 2

doDecor = False

# current available correlator bandwidths
bandwidths = [2, 8, 31, 62, 125, 250, 500]

# list of primary flux calibration sources
fluxcal = {"MERCURY" : False, "VENUS" : False, "MARS" : False, "JUPITER" : False, "SATURN" : False, "URANUS" : False, "NEPTUNE" : False, "PLUTO" : False, "MWC349" : False}

degToRad = 0.017453293;

# class for common variables for all objects
class ObjectCommons :
    def __init__(self) :
        self._name = ""            # name of the object as it is in the MIRIAD file
        self._channels = []        # starting channel number for each window
        self._bandwidths = []      # bandwidth of each window in MHz
        self._numChans = []        # number of channels in each window
        self._superwindow = []     # windows used for "super windows" for calibration and continuum
        self._lsbGood = False      # does the lsb have good data (for single sideband observations)
        self._usbGood = False      # does the usb have good data (for single sideband observations)
        self._superNumChans = 0    # number of channels in the superwide windows, this will be the same for LSB and USB
        self._file = ""            # file name of the MIRIAD file associated with this object
        self._RA = 0.0             # RA of the object
        self._DEC = 0.0            # DEC of the object

    def setName(self, id) :
        """ Method to set the gain calibrator name
            input:
                id - the name of the gain calibrator
            returns :
                none
        """
        self._name = id

    def setBandwidths(self, bw) :
        """ Method to set the bandwidth for a window
            the input value is converted from the actual value to a more readable one
            (i.e. 30.76 becomes 31)
            input :
                bw - the actual bandwidth in MHz
            return :
                none
        """
        self._bandwidths = bw
        for i in range(0,len(bw)) :
            if(bw[i] >= globals.MINCONTINUUMBW) :
                self._superwindow.append(i+1)

    def getBandwidths(self) :
        """ Method to return the bandwidths of all windows of the object
            input :
                none
            returns :
                list of window bandwidths
        """
        return self._bandwidths

    def setChannels(self, chan) :
        """ Method to set the starting channel number for a window
            input :
                chan - the starting channel
            returns :
                none
        """
        self._channels = chan

    def setNumChan(self, num) :
        """ Method to set the number of channels in a window
            input :
                num - the number of channels
            returns :
                none
        """
        self._numChans = num

    def getNumChan(self, window) :
        """ Method to get the number of channels for a given window
            input :
                window - the window to query
            returns :
                the number of channels on the window
        """
        return self._numChans[window - 1]

    def setFile(self, file) :
        """ Method to set the file name
            input :
                file - the name of the miriad file
            returns :
                none
        """
        self._file = file

    def setCoordinates(self,ra,dec) :
        """ Method to set the corrdiantes of the object
            input :
                ra - the RA of the object in radians
                dec - the DEC of the object in radians
            returns :
                none
        """
        self._RA = ra
        self._DEC = dec

    def getMaxBw(self) :
        """ Method to get the maximum bandwidth in the data
            input :
                none
            returns :
                the maximum bandwidth in MHz
        """
        return max(self._bandwidths)

    def getMinBw(self) :
        """ Method to get the minimum bandwidth in the data
            input :
                none
            returns :
                the minimum bandwidth in MHz
        """
        return min(self._bandwidths)

    def getChannelWidth(self,window) :
        """ Method to return the channel width in MHz of a given window
            input :
                window - the window number to query
            returns :
                the channel width of the requested window in MHz
        """
        return float(self._bandwidths[window-1])/float(self._numChans[window-1])

    def getBandwidth(self,window) :
        """ Method to get the bandwidth in MHz of a given window
            input :
                window - the window number to query
            return :
                the bandwidth of the requested window in MHz
        """
        return self._bandwidths[window-1]

    def isSuper(self,window) :
        """ Method to determine if a given window is a member of a
            superwindow
            input :
                window - the window number to query
            return :
                True/False if the window is/not a member of a superwindow
        """
        return (window in self._superwindow)

    def haveSuper(self) :
        """ Method to return whether the current object has superwindows
            input :
                none
            return :
                True/False if the current object has superwindows
        """
        return (len(self._superwindow) > 0)

    def getSuperNumChans(self) :
        """ Method to return the number of channels in the superwindow
            input :
                none
            return :
                the number of channels in each superwindow
        """
        return self._superNumChans

    def calcSuperChans(self) :
        """ Method to calculate the number of channels in the superwindows
            input :
                none
            return :
                none
        """
        self._superNumChans = 0
        for w in self._superwindow :
            self._superNumChans += self._numChans[w - 1]
        self._superNumChans /= 2

    def setIndividualNumChans(self,window,chans) :
        self._numChans[window-1] = chans

# Class for gain calibrators
class Gaincal (ObjectCommons):
    def __init__(self) :
        ObjectCommons.__init__(self)   # import the common variables
        self._flux = 0.0               # gain calibrator flux in Jy (calculated during reduction)
        self._start = []               # starting time for each calibration interval in hours
        self._end = []                 # ending time for each calibration interval in hours
        self._type = 1                 # type of gain calibrator (primary/secondary)

    def setFlux(self, flux) :
        """ Method to set the flux of the gain calibrator
            input:
                flux - the flux in Jy
            return :
                none
        """
        self._flux = flux

    def setTimes(self, start, end) :
        """ Method to set the starting and ending times for a calibration cycle
            input:
                start - the start time in hours
                end - the ending time in hours
            return :
                none
        """
        self._start.append(start)
        self._end.append(end)

    def setType(self, type) :
        """ Method to set the type (PRIMARY/SECONDARY)
            input :
                type - the type
            return :
                none
        """
        self._type = type

    def getFlux(self) :
        """ Method to return then flux of the gain calibrator
            input :
                none
            returns :
                the flux in Jy
        """
        return self._flux

# Class for passband calibrators
class Passcal (ObjectCommons):
    def __init__(self) :
        ObjectCommons.__init__(self)            # import the common variables
        self._type = ""                         # the type of passband (PLANET/NOISE/QUASAR)
        self._freq=[]                           # the starting frequencies of each window in GHz
        self._hybrid = False                    # was the data taken in hybrid mode
        self._hybridConf = {2:False, 8:False, 31:False, 62:False, 500:False}  # which hybrid modes were used
        self._hybridChans = {2:0, 8:0, 31:0, 62:0, 125:0, 250:0, 500:0}       # number of channels for each mode

    def setType(self, type) :
        """ Method to set the type of passband calibrator (PLANET/NOISE/QUASAR)
            input :
                type - the type of passband calibrator
            returns :
                none
        """
        self._type = type

    def setFreq(self, freq) :
        """ Method to set the starting frequency of a window
            input :
                freq - the starting frequency in GHz
            returns :
                none
        """
        self._freq = freq

    def isHybrid(self) :
        """ Method to return whether the current observation was in hybrid mode
            input :
                none
            returns :
                True/False if the observations were/not taken in hybrid mode
        """
        return self._hybrid

    def setHybrid(self, value) :
        """ Method to set whether the passband cal was taken in hybrid mode or not
            input :
                value - True/False is/isn't in hybrid mode
            returns :
                none
        """
        self._hybrid = value

    def setHybridConf(self, bw) :
        """ Method to set which bandwidths are in a hybrid observation
            input :
                bw - the bandwidth found
            returns :
                none
        """
        self._hybridConf[bw] = True

    def setHybridChans(self, bw, chans) :
        """ Method to set the channels from a hybrid observation
            input :
                bw - the bandwidth
                chans - the number of channels
            returns :
                none
        """
        self._hybridChans[bw] = chans

# Class for flux calibrators
class Fluxcal (ObjectCommons):
    def __init__(self) :
        ObjectCommons.__init__(self)      # import the common variables
        self._type = ""                   # the type of flux (PLANET/QUASAR)

    def setType(self, type) :
        """ Method to set the type of flux calibrator (PLANET/QUASAR)
            input :
                type - the type of flux calibrator
            returns :
                none
        """
        self._type = type


# Class for sources
class Source (ObjectCommons):
    def __init__(self) :
        ObjectCommons.__init__(self)      # import the common variables
        self._mosaic = False              # is this source a multipoint mosaic
        self._pointingCenters = []        # list of the pointing centers in offset arcsec from the phase center
        self._cleanRegion = None          # region to be cleaned
        self._continCleanRegion = None    # region to be cleaned in continuum maps

    def setCleanRegion(self,region) :
        """ Method to set the MIRIAD clean region
            input :
                region - the region in MIRIAD format
            returns :
                none
        """
        self._cleanRegion = region

    def getCleanRegion(self) :
        """ Method to return the MIRIAD clean region
            input :
                none
            returns :
                the miriad clean region as a string
        """
        return self._cleanRegion

    def setContinuumCleanRegion(self,region) :
        """ Method to set the MIRIAD clean region for continuum maps
            input :
                the continuum clean region as a string
            returns :
                none
        """
        self._continCleanRegion = region

    def getContinuumCleanRegion(self) :
        """ Method to return the MIRIAD clean region for continuum maps
            input :
                none
            returns :
                the continuum clean region as a string
        """
        return self._continCleanRegion

    def setMosaic(self,value,varList) :
        """ Method to set whether the source is a multipoint mosaic
            input :
                value - True/False is/isn't a multipoint mosaic
                varList - a listing of the pointing centers from varplt
            returns :
                none
        """
        self._mosaic = value
        while(len(varList) > 0) :
            line = varList.pop()
            splitLine = line.split()
            ra = float(splitLine[0])
            dec = float(splitLine[1])
            found = False
            # go through the listing and get the pointing centers, but only add to the list if they are not previously there
            for key in self._pointingCenters :
                if(ra == key[0] and dec == key[1]) :
                    found = True
                    break
            if(not found) :
                self._pointingCenters.append([ra,dec])

# class for polarization calibrators - Note: not implemented yet, just a placeholder
class Polcal (ObjectCommons) :
    def __init__(self) :
        ObjectCommons.__init__(self)

# class to hold all objects for an observation
class Objects :
    def __init__(self) :
        # set up the objects
        self._sources = []
        self._gaincals = []
        self._passcals = []
        self._fluxcals = []

    def addSource(self,src) :
        """ Method to add a source
            input :
                src - a source object
            returns :
                none
        """
        self._sources.append(src)

    def getSource(self,name) :
        """ Method to return a specific source by name
            input :
                name - the name of the source object to return
            returns :
                the source object if found, None otherwise
        """
        for src in self._sources :
            if(src._name == name) :
                return src
        return None

    def updateSource(self,source) :
        """ Method to update/replace a source object
            input :
                source - the source object to insert/update into the list
            returns :
                none
        """
        for src in self._sources :
            if(src._name == source._name) :
                src = source
                break

    def addGaincal(self,cal) :
        """ Method to add a gaincal
            input :
                cal - a gain calibrator object
            returns :
                none
        """
        self._gaincals.append(cal)

    def getGaincal(self,name) :
        """ Method to return a gaincal object by name
            input :
                name - the name of the gaincal object to return
            returns :
                the matching gaincal object, None otherwise
        """
        for gcal in self._gaincals :
            if(gcal._name == name) :
                return gcal
        return None

    def updateGaincal(self,gaincal) :
        """ Method to update/replace a gaincal object
            input :
                gaincal - the gaincal object to insert/update into the list
            returns :
                none
        """
        for gcal in self._gaincals :
            if(gcal._name == gaincal._name) :
                gcal = gaincal
                break

    def addPasscal(self,cal) :
        """ Method to add a passband cal
            input :
                cal - a passband calibrator object
            returns :
                none
        """
        self._passcals.append(cal)

    def getPasscal(self,name) :
        """ Method to return a passcal object by name
            input :
                name - the name of the passcal object to return
            returns :
                the matching passcal object, None otherwise
        """
        for pcal in self._passcals :
            if(pcal._name == name) :
                return pcal
        return None

    def updatePasscal(self,passcal) :
        """ Method to update a Passcal object
            input :
                passcal - the passcal object to update
            returns :
                none
        """
        for pcal in self._passcals :
            if(pcal._name == passcal._name) :
                pcal = passcal
                break

    def addFluxcal(self,cal) :
        """ Method to add a flux cal
            input :
                cal - a flux calibrator object
            returns :
                none
        """
        self._fluxcals.append(cal)

    def getFluxcal(self,name) :
        """ Method to get a fluxcal object by name
            inputs :
                name - the name of the fluxcal object to return
            returns :
                the matching fluxcal object, None otherwise
        """
        for fcal in self._fluxcals :
            if(fcal._name == name) :
                return fcal
        return None

    def updateFluxcal(self,fluxcal) :
        """ Method to update a specific fluxcal object
            input :
                fluxcal - the fluxcal object to update
            returns :
                none
        """
        for fcal in self._fluxcals :
            if(fcal._name == fluxcal._name) :
                fcal = fluxcal
                break

    def convertGainToSource(self,gaincal) :
        """ Method to convert a gaincal object to a source object
            input :
                gaincal - the gaincal object to convert
            returns :
                none
        """
        src = Source()
        src._name = gaincal._name
        src._channels = gaincal._channels
        src._bandwidths = gaincal._bandwidths
        src._numChans = gaincal._numChans
        src._file = gaincal._file
        src._superwindow = gaincal._superwindow
        src._lsbGood = gaincal._lsbGood
        src._usbGood = gaincal._usbGood
        src._superNumChans =  gaincal._superNumChans
        self._sources.append(src)
        del self._gaincals[self._gaincals.index(gaincal)]

    def updateIndividualNumChans(self,window,chans) :
        """ Method to update the number of channels in a given
            window for all sources
            input :
                window - the window number to update
                chans - the number of channels in the window
            returns :
                none
        """
        for o in self._sources + self._gaincals + self._passcals + self._gaincals :
            o.setIndividualNumChans(window,chans)


def isPlanet(obj) :
    """ Method to determine if the given object is a planet
        input :
            obj - the name of the object
        returns :
            True/False - is/isn't a planet
    """
    if("URANUS" in obj) :
        return True
    if("NEPTUNE" in obj) :
        return True
    if("SATURN" in obj) :
        return True
    if("JUPITER" in obj) :
        return True
    if("PLUTO" in obj) :
        return True
    if("MARS" in obj) :
        return True
    if("VENUS" in obj) :
        return True
    if("MERCURY" in obj) :
        return True
    return False

def calSort(objects) :
    """ Function to sort calibrators and flag any overlapping data
        input :
            objetcs - a list of the objects
        returns :
            none
    """
    # check for full overlap first
    done = []
    for gcal in objects._gaincals :
        for gcal2 in objects._gaincals :
            if(gcal._name != gcal2._name and not(gcal._name in done)) :
                if(gcal._start[0] >= gcal2._start[0] and gcal._end[-1] <= gcal2._end[-1]) :
                    done.append(gcal._name)
                    objects.convertGainToSource(gcal)
    # check for partial overlap
    done = []
    for gcal in objects._gaincals :
        for gcal2 in objects._gaincals :
            if(gcal._name != gcal2._name and not(gcal2._name in done)) :
                if(gcal._end[-1] > gcal2._start[0] and gcal._start[0] < gcal2._start[0]) : # we have an overlap that needs to be taken care of
                    # flag the overlaping data in the later data set
                    log.writeComment("Flag the overlaping data in the later data set")
                    args = []
                    args.append(globals.Variable("vis",gcal2._file))
                    args.append(globals.Variable("select","time'('%s,%s')'" % (calculations.unconvertTime(gcal2._start[0] - 0.0083),calculations.unconvertTime(gcal._end[-1] + 0.00833))))
                    args.append(globals.Variable("flagval","flag"))
                    sys = log.run("uvflag",args)
                    if(sys == 0) :
                        log.writeLog("Flagged %s from %s to %s to avoid primary calibrator overlap" % (gcal2._file,calculations.unconvertTime(gcal2._start[0] - 0.0083),calculations.unconvertTime(gcal._end[-1] + 0.00833)))
        done.append(gcal)

def sort(visFile, objects, refant) :
    """ Method to sort though the main data file and break it up for each source
        input :
            visFile - the main visibility file
            objects - the objects
            refant - the reference antenna
        returns :
            none
    """
    global doDecor
    # separate out the noise source before applying linelength corrections
    for pcal in objects._passcals :
        if(pcal._name == "NOISE") :
            getWindowInfo(visFile,pcal,isPasscal=True)

    # apply linelength calibration (except for the NOISE source)
    #log.writeComment("Applying line length corrections")
    #args = []
    #args.append(globals.Variable("vis",visFile))
    #log.run("linecal",args)
    # separate the source(s) from the main data
    if(len(objects._sources) > 0) :
        log.writeComment("Separating sources into individual files")
        for source in objects._sources :
            getWindowInfo(visFile,source,isSource=True)
    log.writeAll("\n")
    # separate out the gaincal file(s)
    if(len(objects._gaincals) > 0) :
        log.writeComment("Separating gain calibrators into individual files.")
        for gcal in objects._gaincals :
            if(gcal._name == "NOISE") :
                continue
            getWindowInfo(visFile,gcal,isGaincal=True)
    if(len(objects._gaincals) == 0) :
        raise Exception, "No gain calibrators located, halting data reduction"
    doDecor = decorrelation.correctDecorrelation(objects,refant)
    log.writeAll("\n")
    # separate passband calibrator(s) into file(s)
    if(len(objects._passcals) > 0) :
        log.writeComment("Separating passband calibrators into individual files.")
        for pcal in objects._passcals :
            if(pcal._name == "NOISE") :
                continue
            getWindowInfo(visFile,pcal,isPasscal=True)
    log.writeAll("\n")
    # separate the flux calibrator(s) into file(s)
    if(len(objects._fluxcals) > 0) :
        log.writeComment("Separating flux calibrators into individual files")
        for fcal in objects._fluxcals :
            if(fcal._name == "NOISE") :
                continue
            getWindowInfo(visFile,fcal,isFluxcal=True)
    log.writeAll("\n")

def findBandwidth(bw) :
    """ Method to convert the bandwidth of a window from its real value to a more readable one
        input :
            bw - the bandwidth in MHz
        returns :
            the bandwidth in MHz
    """
    for i in bandwidths :
        if(abs(bw) < i) :
            return i
    return 500

def getPrimaryFluxCal(objects) :
    """ Method to determine the primary flux calibrator
        searches the objects looking for the best one in descending order of preference
        input :
            objects - the objects object
        returns:
            the name of the primary flux calibrator
    """
    primaryFluxcal = Fluxcal()
    if(fluxcal["NEPTUNE"]) :
        for fcal in objects._fluxcals :
            if(fcal._name == "NEPTUNE") :
                primaryFluxcal = fcal
    elif(fluxcal["URANUS"]) :
        for fcal in objects._fluxcals :
            if(fcal._name == "URANUS") :
                primaryFluxcal = fcal
    elif(fluxcal["MARS"]) :
        for fcal in objects._fluxcals :
            if(fcal._name == "MARS") :
                primaryFluxcal = fcal
    elif(fluxcal["MWC349"]) :
        for fcal in objects._fluxcals :
            if(fcal._name == "MWC349") :
                primaryFluxcal = fcal
    elif(fluxcal["MERCURY"]) :
        for fcal in objects._fluxcals :
            if(fcal._name == "MERCURY") :
                primaryFluxcal = fcal
    elif(fluxcal["JUPITER"]) :
        for fcal in objects._fluxcals :
            if(fcal._name == "JUPITER") :
                primaryFluxcal = fcal
    elif(fluxcal["SATURN"]) :
        for fcal in objects._fluxcals :
            if(fcal._name == "SATURN") :
                primaryFluxcal = fcal
    elif(fluxcal["VENUS"]):
        for fcal in objects._fluxcals :
            if(fcal._name == "VENUS") :
                primaryFluxcal = fcal
    elif(fluxcal["PLUTO"]) :
        for fcal in objects._fluxcals :
            if(fcal._name == "PLUTO") :
                primaryFluxcal = fcal
    else:
        primaryFluxcal = objects._fluxcals[0]
    log.writeLog("Using %s as a flux calibrator for bootflux" % (primaryFluxcal._name))
    return primaryFluxcal

def getWindowInfo(visFile,object,isSource=False,isGaincal=False,isPasscal=False,isFluxcal=False) :
    """ Method to separate the main miriad file into smaller ones for each object
        input :
            visFile - the main MIRIAD visibility file
            object - the name of the object to separate
            source - whether the object is a source
            gaincal - whether the object is a gain calibrator
            passcal - whether the object is a passband calibrator
            fluxcal - whether the object is a flux calibrator
    """
    global doDecor
    numberOfWins = len(object._bandwidths)
    # if the output file exists then we must create a copy for the current purpose
    log.writeComment("Separating data into individual, window based, files")
    if(os.path.exists(object._file) == 1) :
        id = ""
        if(isGaincal) :
            id = globals.GAINCALEND
        elif(isPasscal) :
            id = globals.PASSCALEND
        elif(isFluxcal) :
            id = globals.FLUXCALEND
        else :
            log.writeLog("Error: File %s exists exiting" % (object._file))
            raise Exception, "Error: File %s exists exiting" % (object._file)
        print "SOURCE copy",object._file
        args = []
        args.append(globals.Variable(None,object._file))
        args.append(globals.Variable(None,object._file,"." + id))
        log.run("cp -r",args,fatal=True)
        object._file = object._file +"." + id
    else :
        print "SOURCE new",object._file
        args = []
        args.append(globals.Variable("vis",visFile))
        args.append(globals.Variable("select",object._name,"')'","-auto,source'('"))
        args.append(globals.Variable("out",object._file))
        log.run("uvcat",args,fatal=True)
        if(object._name == "NOISE") :
            args = []
            args.append(globals.Variable("vis","NOISE"))
            args.append(globals.Variable("options","noisecal"))
            args.append(globals.Variable("out","NOISE",".conj; sleep 3; rm -r"))
            args.append(globals.Variable(None,"NOISE","; sleep 3; mv"))
            args.append(globals.Variable(None,"NOISE",".conj"))
            args.append(globals.Variable(None,"NOISE"))
            log.run("uvcal",args)
    # correct for decorrelation if needed
    if(doDecor and object._name != "NOISE") :
        log.writeComment("Correcting for decorrelation on long baselines")
        fileEnd = str(random.randint(1,100000))
        args = []
        args.append(globals.Variable("vis",object._file))
        args.append(globals.Variable("options","nocal,nopass,nopol"))
        args.append(globals.Variable("delaymax","8500"))
        args.append(globals.Variable("out","temp." + fileEnd,"; sleep 3; rm -rf"))
        args.append(globals.Variable(None,object._file,"/*; rm -rf"))
        args.append(globals.Variable(None,object._file,"; sleep 3; mv temp." + fileEnd))
        args.append(globals.Variable(None,object._file))
        log.run("uvdecor",args,fatal=True)
    if(isPasscal and object.isHybrid()) :
        hybrid.splitHybrid(object, numberOfWins)
    else :
        # separate each window
        log.writeAll("\n")
        log.writeComment("Separating individual windows")
        superwin = {"LSB" : [], "USB" : []}
        for window in range(globals.STARTWINDOW, globals.ENDWINDOW + 1) :
            if(not object.isSuper(window)) :
                args = []
                args.append(globals.Variable("vis",object._file))
                args.append(globals.Variable("select","window'('%i')'" % (window)))
                args.append(globals.Variable("out",object._file,".w%i" % (window)))
                log.run("uvcat",args)
            else :
                if(window <= numberOfWins/2) :
                    superwin["LSB"].append(window)
                else :
                    superwin["USB"].append(window)
        # separate out the super windows
        lsbString = ""
        usbString = ""
        for window in superwin["LSB"] :
            lsbString += "%i," % (window)
        for window in superwin["USB"] :
            usbString += "%i," % (window)
        if(len(superwin["LSB"]) > 0) :
            object._lsbGood = True
            args = []
            args.append(globals.Variable("vis",object._file))
            args.append(globals.Variable("select","window'('%s')'" % (lsbString[:-1])))
            args.append(globals.Variable("out",object._file,".LSB"))
            log.run("uvcat",args)

        if(len(superwin["USB"]) > 0) :
            object._usbGood = True
            args = []
            args.append(globals.Variable("vis",object._file))
            args.append(globals.Variable("select","window'('%s')'" % (usbString[:-1])))
            args.append(globals.Variable("out",object._file,".USB"))
            log.run("uvcat",args)
        object.calcSuperChans()

        for window in range(globals.STARTWINDOW, globals.ENDWINDOW + 1) :
            # flag edge channels
            if(object._bandwidths[window - 1] == 62) :
                log.writeAll("\n")
                log.writeComment("Flagging edge channels in 62 MHz windows")
                args = []
                args.append(globals.Variable("vis",object._file,".w%i" % (window)))
                args.append(globals.Variable("edge","3,3"))
                args.append(globals.Variable("flagval","flag"))
                log.run("uvflag",args)

def getDistance(srcRA, srcDEC, RA, DEC) :
    """ Method to calculate the distace between the two points on the sky
        input :
            srcRA - the RA of point 1 in radians
            srcDEC -  the DEC of point 1 in radians
            RA - the RA of point 2 in radians
            DEC - the DEC of point 2 in radians
        returns :
            the distance between the points in radians
    """
    distance1 = math.sqrt((srcRA - RA)**2 + (srcDEC - DEC)**2)
    distance2 = 0
    if(srcRA > RA) :
        distance2 = math.sqrt(((srcRA - 2*math.pi)-RA)**2 + (srcDEC - DEC)**2)
    else:
        distance2 = math.sqrt(((RA - 2*math.pi)-srcRA)**2 + (srcDEC - DEC)**2)
    return min(distance1,distance2)

def sortDistances(secondaries, srcs) :
    """ Method to source objcts by distance (nearest to farthest
        input :
            secondaries - the list of objects to sort
            srcs - the reference point object from which to sort distence from
        returns :
            none
    """
    if(len(srcs) == 0) :
        return
    distances = {}
    for object in secondaries :
        distances[object] = getDistance(srcs[0]._RA,srcs[0]._DEC,object._RA,object._DEC)
    secondaries = sorted(distances, key=distances.__getitem__, reverse=True)
