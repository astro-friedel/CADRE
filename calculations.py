import time
from math import *
import datetime
import logger as log
from copy import deepcopy
from operator import itemgetter
import numpy
import globals

"""
Module of statistical tools and other calculations
Part of CARMA data reduction pipeline
Author: D. N. Friedel
"""

C = 299792458.0       # speed of light

def gregorianToNumeric(dateString) :
    """ Method to convert a date string in Gregorian format
        to a numeric format
        input :
            dateString - in the form of yymmmdd (07dec17)
        returns :
            in the form of yyyy.yyyy (2007.9616438356165)
    """
    dateTuple = time.strptime(dateString,'%y%b%d')
    year = int(dateString[0:2])
    date = 2000.0
    if(year > 50) :
        date = 1900.0
    if(year % 4 == 0) :
        date = date + float(year) + float(dateTuple.tm_yday)/366.0
    else :
        date = date + float(year) + float(dateTuple.tm_yday)/365.0
    return date

def numericalGregorianToNumeric(dateString) :
    """ Method to convert a date string in numerical Gregorian format
        to a numeric format
        input :
            dateString - in the form of yymmdd (071217)
        returns :
            in the form of yyyy.yyyy (2007.9616438356165)
    """
    dateTuple = time.strptime(dateString,'%y%m%d')
    year = int(dateString[0:2])
    date = 2000.0
    if(year % 4 == 0) :
        date = date + float(year) + float(dateTuple.tm_yday)/366.0
    else :
        date = date + float(year) + float(dateTuple.tm_yday)/365.0
    return date

def numericToConventional(date) :
    """ Method to convert numeric date to more conventional format
    """
    numdays = 365.0
    year = int(floor(date))
    if(year % 4 == 0) :
        numdays = 366.0
    day  = int(floor((date - year)*numdays))
    dateTuple = time.strptime(str(year) + " " + str(day),"%Y %j")
    month = "Dec"
    if(dateTuple.tm_mon == 1):
        month = "Jan"
    elif(dateTuple.tm_mon == 2):
        month = "Feb"
    elif(dateTuple.tm_mon == 3):
        month = "May"
    elif(dateTuple.tm_mon == 4):
        month = "Apr"
    elif(dateTuple.tm_mon == 5):
        month = "May"
    elif(dateTuple.tm_mon == 6):
        month = "Jun"
    elif(dateTuple.tm_mon == 7):
        month = "Jul"
    elif(dateTuple.tm_mon == 8):
        month = "Aug"
    elif(dateTuple.tm_mon == 9):
        month = "Sep"
    elif(dateTuple.tm_mon == 10):
        month = "Oct"
    elif(dateTuple.tm_mon == 11):
        month = "Nov"
    return "%4i-%s-%2i" % (year,month,dateTuple.tm_mday)

def convertTime(time) :
    """ Method to convert from a numeric time entry (as a string)
        to an hour based entry
        input :
            time - has the form of hhmmss.s (121500.0)
        returns :
            in the form of hh.hhhh  (12.25)
    """
    hours = float(time[0:2])
    minutes = float(time[2:4])/60.0
    seconds = float(time[4:])/3600.0
    return (hours+minutes+seconds)

def unconvertTime(time) :
    """ Method to convert from numeric hours time entry to standard format
        input :
            time - time in hours (12.25)
        returns :
            time in the form of hh:mm:ss.s (12:15:00.0)
    """
    hours = floor(time)
    mins = (time - hours) * 60.0
    minutes = floor(mins)
    seconds = (mins - minutes) * 60.0
    if(seconds > 59.99) :
        seconds -= 59.99
        minutes += 1.0
    if(minutes > 59.99) :
        minutes -= 59.99
        hours += 1.0
    if(hours >= 24) :
        hours = hours - 24
    return "%02i:%02i:%04.1f" % (hours,minutes,seconds)

def getDate() :
    """ Method to get the current date and return it as a string
        inputs :
            none
        returns :
            the current date in YYYY-MM-DDTHH:MM:SS format
    """
    dateString = ""
    date = datetime.datetime.now()
    dateString = "%4s-%02d-%02dT%02d:%02d:%02d" %(date.year,date.month,date.day,date.hour,date.minute,date.second)
    return dateString

def calcCellSize(avgBaseline,lastSize) :
    """ Method to calculate the cell size for the invert routine.
        Cell size is based on the estimated synthesized beam
        and allows for at least 5 pixels across the beam.
        input :
            avgBaseline - the average baseline in lambda
            lastSize - the last cell size calculated (used for adaptive cell sizes to minimize rejected data in invert)
        returns :
            the cell size in arcseconds
    """
    if(lastSize != 0.0) :
        return lastSize/1.1
    beamsize = 91000.0 / avgBaseline
    if(beamsize > 30.0) :
        return 8.0
    if(beamsize > 20.0) :
        return 4.0
    if(beamsize > 8.0) :
        return 2.0
    if(beamsize > 4.0) :
        return 1.0
    if(beamsize > 2.0) :
        return 0.5
    if(beamsize > 1.0) :
        return 0.25
    if(beamsize > 0.4) :
        return 0.1
    if(beamsize > 0.2) :
        return 0.05
    if(beamsize > 0.04) :
        return 0.01
    return 0.005

def calcImsize(freq, cellsize,lastSize = 0.0) :
    """ Method to calculate the image size for the invert routine.
        Calculations try to get the 10m beam to fit in the inner
        quarter of the map.
        input :
            freq - observing frequency in GHz
            cellsize - the cellsize in arcsec
            lastSize - the size of the last calculation
        returns :
            the image size in arcseconds
    """
    if(lastSize != 0.0) :
        return lastSize * 1.1
    wavelength = C / (freq * 10**9)
    dishSize = 6.1
    if(globals.isSci2) :
        dishSize = 3.5
    quarter = ((1.22 * wavelength / dishSize) * 206265 / cellsize) + 2
    # want ovro to be in inner quarter
    return int(floor(quarter))

def weightedAverage(values, uncert) :
    """ Method to calculate the wieghted average of a list of points
        input :
            values - the values to be averaged
            uncert - their uncertainties
        returns :
            the weighted average (without its uncertainty)
    """
    wts = []
    upper = 0.0
    lower = 0.0
    for i in range(0, len(values)) :
        if(values[i] != 0.0) :
            wt = 1.0/uncert[i]**2
            upper += values[i] * wt
            lower += wt
    if(upper != 0.0 and lower != 0.0) :
        return upper / lower,sqrt(1.0/lower)
    return 0.0,0.0

def createPoly(center,cellsize,freq,offset,small = False) :
    """ Method to create a polygon for a clean region for multipoint mosaics
        the polygon size and location is based on the given offset coordinates, cellsize and 10 m primary beam size
        input :
            center - the center point in offset arcseconds from the pointing center
            cellsize - the cellsize in arcseconds
            freq - the observing frequency
            offset - the offset from 0,0 of the pointing center in the map (in pixels)
        returns :
            the points of the polygon as a list
    """
    angle = 22.5
    poly = []
    if(small) :
        mainBeam = cellsize    # for small polygons - special use
    else :
        wavelength = C / (freq * 10**9)
        mainBeam = 0.51*(1.22 * wavelength /4.0) * 206265.0/cellsize  #main beam radius in pixels
    offset1 = mainBeam*sin(radians(angle))
    offset2 = sqrt(mainBeam**2 - offset1**2)
    poly.append([-int(center[0]/cellsize + offset1) + offset[0], int(center[1]/cellsize + offset2) + offset[1]])
    poly.append([-int(center[0]/cellsize + offset2) + offset[0], int(center[1] /cellsize+ offset1) + offset[1]])
    poly.append([-int(center[0]/cellsize + offset2) + offset[0], int(center[1]/cellsize - offset1) + offset[1]])
    poly.append([-int(center[0]/cellsize + offset1) + offset[0], int(center[1]/cellsize - offset2) + offset[1]])
    poly.append([-int(center[0]/cellsize - offset1) + offset[0], int(center[1]/cellsize - offset2) + offset[1]])
    poly.append([-int(center[0]/cellsize - offset2) + offset[0], int(center[1]/cellsize - offset1) + offset[1]])
    poly.append([-int(center[0]/cellsize - offset2) + offset[0], int(center[1]/cellsize + offset1) + offset[1]])
    poly.append([-int(center[0]/cellsize - offset1) + offset[0], int(center[1]/cellsize + offset2) + offset[1]])
    poly.append([-int(center[0]/cellsize + offset1) + offset[0], int(center[1]/cellsize + offset2) + offset[1]])
    return poly

def pointInPoly(poly,x,y) :
    """ Method to determine if a point is inside a ploygon
        input :
            poly - the polygon
            x - the x value
            y - the y value
        returns :
            True/False (x,y) is/isn't insode poly
    """
    inside = False
    j = 0
    for i in range(0,len(poly)) :
        j += 1
        if(j == len(poly)) :
            j = 0
        if((poly[i][1] < y and poly[j][1] >= y) or (poly[j][1] < y and poly[i][1] >= y)) :
            if((poly[i][0] + float(y - poly[i][1])/float(poly[j][1] - poly[i][1])*(poly[j][0] - poly[i][0])) < x) :
                inside = not inside
    return inside

def equals(a,b,limit) :
    """ Method to determine if 2 quantities are equal within an uncertainty
        input :
            a - the first point
            b - the second point
            limit - the acceptable tolerance
        returns :
            True/False a and b are within limit of each other
    """
    return abs(a - b) < limit

def findIntersection(start1,end1,start2,end2) :
    """ Method to find the intersection of 2 line segments
        input :
            start1 - the starting coordiantes of the first segment
            end1 - the end coordinates of first segment
            start2 - the starting coordiantes of the second segment
            end2 - the end coordiantes of the second segment
        returns :
            None = paralell
            [2,x,y] = paralell and overlaping
            [1,x,y] = intersect at x,y
            [0,x,y] intersect at x,y, but x,y outside of segments
    """
    intersection = []
    LIMIT = 0.00001
    INFINITY = 1*10**10
    a1 = 0.0
    a2 = 0.0
    y = 0.0
    x = 0.0
    if(equals(float(start1[0]),float(end1[0]),LIMIT)) :
        a1 = INFINITY
    else :
        a1 = float(start1[1] - end1[1]) / float(start1[0] - end1[0])
    if(equals(float(start2[0]),float(end2[0]),LIMIT)) :
        a2 = INFINITY
    else :
        a2 = float(start2[1] - end2[1]) / float(start2[0] - end2[0])
    b1 = float(start1[1]) - a1*float(start1[0])
    b2 = float(start2[1]) - a2*float(start2[0])

    if(equals(a1,a2,LIMIT)) :
        if(not equals(b1,b2,LIMIT)) :
            return None # we are parallel
        else :
            if(isBetween(start2[0],end2[0],start1[0]) or isBetween(start2[0],end2[0],end1[0])) :
                return [2,0,0]
        return None
    if(equals(a1,INFINITY,LIMIT)) :
        x = float(start1[0])
        y = a2 * x + b2
    elif(equals(a2,INFINITY,LIMIT)) :
        x = float(start2[0])
        y = a1 * x + b1
    else :
        x = -(b1 - b2)/(a1 - a2)
        y = a1 * x + b1
    if(isBetween(start1[0],end1[0],x) and isBetween(start2[0],end2[0],x) and isBetween(start1[1],end1[1],y) and isBetween(start2[1],end2[1],y)) :
        return [1,int(x),int(y)]
    return [0,int(x),int(y)]

def isBetween(start,end,point) :
    """ Method to determine if a point is between 2 others
        input :
            start - first point
            end - second point
            point - thrird point
        returns :
            True/False if point is between start and end
    """
    return (point <= max(start,end) and point >= min(start,end))

def mergePoly(polys,main,merge,inside) :
    """ Method to merge polygons together into a single polygon
        Note that the method used here is not generally aplicable. It assumes that
        all polygons to be merged have the same shape (i.e. all octagons)
        input :
            polys - a list of all polygons
            main - the index of the main polygon
            merge - the index of the polygon to merge
            inside - a list of points of merge that are inside of main
        returns :
            none
    """
    tempPoly = deepcopy(polys[main])
    tempMerge = deepcopy(polys[merge])
    first = -1
    all = True
    foundIn = False
    for i in range(0,len(inside)) : # find the first point that is inside of the polygon
        if(inside[i]) :
            foundIn = True
        all = all and inside[i]
        if(not all and first == -1 and foundIn) :
            first = i
    if(not foundIn) :
        first = 0
    if(all) :
        del polys[merge]
        return
    list = []
    for i in range(first,8) :
        list.append(i)
    for i in range(0,first) :
        list.append(i)
    start = None
    end = None
    remove = []
    add = []
    inpoint = -1
    outpoint = -1
    paralell = False
    for coord in list :
        if(inside[coord] and not(inside[coord + 1])) : # passing out of the polygon
            for point in range(0,len(tempPoly) - 1) :
                intersect = findIntersection(tempMerge[coord],tempMerge[coord + 1],tempPoly[point],tempPoly[point + 1])
                if(intersect != None) :
                    if(intersect[0] == 2) :
                        paralell = True
                    elif(intersect[0] == 1) :
                        start = [int(intersect[1]),int(intersect[2])]
                        add.append(coord + 1)
                        outpoint = point + 1
                        break    # there should only be one outward crossing
        elif(not inside[coord] and (inside[coord + 1])) : # passing into the polygon
            for point in range(0,len(tempPoly) - 1) :
                intersect = findIntersection(tempMerge[coord],tempMerge[coord + 1],tempPoly[point],tempPoly[point + 1])
                if(intersect != None) :
                    if(intersect[0] == 2) :
                        a = 1
                    elif(intersect[0] == 1) :
                        end = [int(intersect[1]),int(intersect[2])]
                        add.append(coord)
                        inpoint = point + 1
                        break    # there should only be one inward crossing
        elif(not inside[coord] and not(coord in add)) :# and not(coord == 8 and (0 in add)) and not(coord == 0 and (8 in add))) :
            add.append(coord)
    # now determine which points need to be removed from the main polygon and exactly where the new points need to be added
    reverse = False
    if(not paralell) :
        if(outpoint > inpoint) :
            del tempPoly[outpoint:]
            del tempPoly[0:inpoint]
            reverse = True
        else :
            if(add[0] == add[-1]) :
                add = add[:-1]
            del tempPoly[outpoint:inpoint]
    tempList = [start]
    for i in add :
        tempList.append(tempMerge[i])
    tempList.append(end)
    tempPoly[outpoint:outpoint] = tempList
    if(reverse) :
        tempPoly.append(tempPoly[0])
    polys[main] = stripNone(tempPoly)
    del polys[merge]
    return

def stripNone(poly) :
    """ Method to remove None's from a list
        input :
            poly - the list to check
        returns :
            the input list with None's removed
    """
    done = False
    while(not done) :
        try :
            indx = poly.index(None)
            del poly[indx]
        except :
            done = True
    return poly

def distance(a,b) :
    """ Method to determine the distance between 2 points
        input :
            a - one point
            b - the other point
        returns :
            the distance
    """
    return sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2)

def sortPointing(oldPointing,offset) :
    """ Method to sort pointings in increasing order of distance from the initial point
        input :
            list of pointing offsets
        returns :
            sorted list of pointing offsets
    """
    # we want the initial pointing to be the eastern most point
    pointing = []
    m1 = 0
    dist = 10000000.0
    for i in oldPointing :
        length = distance([2*offset[0],0.0],i)
        if(length < dist) :
            dist = length
            m1 = oldPointing.index(i)
    pointing.append(oldPointing[m1])
    del oldPointing[m1]
    # now sort in increasing distance
    while(len(oldPointing) > 0) :
        dist = 100000000.0
        marker = -1
        for i in oldPointing :
            length = distance(pointing[0],i)
            if(length < dist) :
                dist = length
                marker = oldPointing.index(i)
        if(marker < 0) :
            # this should never happen
            return
        pointing.append(oldPointing[marker])
        del oldPointing[marker]
    return pointing

def calculateMosaicCleanRegion(pointing,cellsize,freq,offset) :
    """ Method to calculate the clean region for a multipoint mosaic
        input :
            pointing - a list of the pointing positions in offset arcseconds from the phase center
            cellsize - the cellsize in arcseconds
            freq - observing frequency in GHz
            offset - offset from 0,0 of the phase center in the map
        returns :
            the clean region
    """
    polys = []
    temp = ""
    for i in pointing :
        temp = temp + "%i  %i\n" % (i[0],i[1])
    pointing = sortPointing(pointing,offset)
    # create the hexagonal regions to clean for each pointing
    for center in range(0,len(pointing)) :
        polys.append(createPoly(pointing[center],cellsize,freq,offset))
    # see which ones overlap and merge them one at a time
    # since all regions are the same size then they can only overlap if one has a point inside the other
    while(len(polys) > 1) :
        inside = [False,False,False,False,False,False,False,False,False]
        overlap = False
        for poly in range(1,len(polys)) :
            for coord in range(0,8) :
                if(pointInPoly(polys[0],polys[poly][coord][0],polys[poly][coord][1])) :
                    inside[coord] = True
                    overlap = True
            if(overlap) :
                inside[8] = inside[0]
                mergePoly(polys,0,poly,inside)
                break
        if(not overlap) :
            log.writeLog("regions do not overlap")
            # something is not right
            return
    doubles = True
    while(doubles) :
        doubles = False
        for i in range(0,len(polys[0]) - 1) :
            if(polys[0][i] == polys[0][i + 1]) :
                del polys[0][i]
                doubles = True
                break
    region = "polygon'('"
    for i in polys[0] :
        region = region + "%i,%i," % (max(1,min(i[0],(offset[0]*2.0))),max(1,min(i[1],(offset[1]*2.0))))
    region = region[:-1] + "')'"
    return region

def listAverage(data) :
    """ Method to calculate the average of a list of numbers
        input :
            data - the list to average
        returns :
            the average value of the list
    """
    value = 0.0
    length = 0.0
    for i in data :
        if(data[i] < -9999.0) :
            continue
        value += data[i]
        length += 1.0
    value /= length
    return value

def splitObsblockId( obsblockID ) :
    """ Method to split the input obsblockID into tokens based with dot (.) as
        the delimiter.  Valid input will be 'project.obsblock' or
        'project.obsblock.subobslock'.  Other inputs will raise
        exception.  A vector with size 3 is always returned:
        [project, obsblock, subobsblock], with subobsblock being an empty
        string if not present in the input.

        input :
            obsblockID - the fully-qualified obsblock ID
            (project.obsblock[.subobsblock]), without the trial number.
        returns :
            none
    """
    data = dict()
    pos = obsblockID.find(".mir")
    obsblockID = obsblockID[:pos]
    obpieces_ = obsblockID.split('.')
    oblen_ = len( obpieces_ )
    if ( oblen_ < 3 ) :
        raise Exception, "Obsblock ID must have at least a project.obsblock.trial"
    data["project"] = obpieces_[0]
    data["obsblock"] = obpieces_[1]
    if("_" in obpieces_[1]) :
        temp = obpieces_[1][:obpieces_[1].find("_")]
        if("A" in temp) :
            data["configuration"] = "A"
        elif("B" in temp) :
            data["configuration"] = "B"
        elif("C" in temp) :
            data["configuration"] = "C"
        elif("D" in temp) :
            data["configuration"] = "D"
        elif("E" in temp) :
            data["configuration"] = "E"
        elif("SH" in temp) :
            data["configuration"] = "SH"
        elif("SL" in temp) :
            data["configuration"] = "SL"
        else :
            data["configuration"] = "X"
    else :
        data["configuration"] = "X"
    if ( oblen_ == 3 ) :
        if ( data["project"] == "" or  data["obsblock"] == "" ) :
            raise Exception, "Empty project or obsblock not allowed."
        try :
            if(int(obpieces_[2]) > 0) :
                data["subObsblock"] = ""
                data["trial"] = int(obpieces_[2])
                return data
        except ValueError :
            raise Exception,"Obsblock ID must have a trial number"
    elif ( oblen_ == 4 ) :
        data["subObsblock"] = obpieces_[2]
        data["trial"] = int(obpieces_[3])
        return data
    raise Exception, "Invalid obsblock format"

def linearInterpolation(x1,y1,x2,y2,x3) :
    """ Method to linearly interpolate between 2 points
        input :
            x1,y1 - the first point
            x2,y2 - the second point
            x3 - the point to interpolate to (x value)
        returns :
            the y value corresponding to x3
    """
    a = (y2 - y1)/(x2 - x1)
    b = y1 - (a * x1)
    return (a * x3) + b

def freqToVel(restFreq,freq) :
    """ Method to covert from a frequency difference to velocity
        input :
            restFreq - the reference frequency
            freq - the offset frequency
        returns :
            the velocity difference (based on restFreq)
    """
    vel = 0.0
    df = (restFreq-freq) * pow(10.0,9)
    return df*(C/1000)/(restFreq * pow(10.0,9))

def velToFreq(vel,restFreq):
    """ Method to convert from a velocity offset to a frequency
        input :
            vel - the velocity offset
            restFreq - the rest frequency
        returns :
            the offset frequency
    """
    df = ((vel * restFreq * pow(10.0,9))/(C/1000))
    return restFreq - (df/pow(10.0,9))

def listToString(list) :
    """ Method to return a string containing the items of a list
        input :
            list - the list to convert
        returns :
            a string containing the items of the list
    """
    string = ""
    for i in list :
        string += str(i) + ","
    return string[:-1]

def intMidpoint(a,b) :
    """ Method for calculating the midpoint of two sets of integers
        input :
            a,b - two 2 item lists of x,y coordinates
        returns :
            the midpoint coordinate (integer values)
    """
    return [int(round((a[0] + b[0])/2)), int(round((a[1] + b[1])/2))]

def cabs(value) :
    """ Method to take the absolute value of a complex number
        input :
            value - a complex number
        returns :
            the absolute value (as fortran would)
    """
    if(not isinstance(value,complex)) :
        raise Exception, "Non complex value given to complex function"
    return sqrt(pow(value.real,2.0) + pow(value.imag,2.0))

def sortidx(values) :
    """ Method to sort a list
        input :
            values - the input list
        returns :
            a sorted list
    """
    temp = deepcopy(values)
    idx = [0] * len(values)
    temp.sort()
    temp.reverse()
    for i in range(0,len(values)) :
        idx[i] = values.index(temp[i])
    return idx

def fitPoly(x,y,order) :
    """ Method to fit a ploynomial to a set of data
        input :
            x - the x points
            y - the y points (no checking is done to make sure x and y are the same length)
            order - the order of the polynomial to fit
        returns :
            a tuple of the fit parameters (length is order + 1)
    """
    fx = numpy.array(x)
    fy = numpy.array(y)
    z = numpy.polyfit(fx,fy,order)
    return z

def getRms(xt,yt,z) :
    """ Method to calculate the rms of a fit to data points
        input :
            xt - the x points (data)
            yt - the y points (data)
            z - the fit parameters
        returns :
            the rms of the fit
    """
    p = numpy.poly1d(z)
    length = len(xt)
    value = 0.0
    rms = 0.0
    for i in range(0,length) :
        tempy = p(xt[i]) - yt[i]
        value += tempy
        rms += tempy**2
    value /= length
    return sqrt((rms - length*value**2)/(length-1))

def iround(x):
    """ Method to round a number to the nearest integer
        input :
            x - the number to round
        returns :
            the number rounded to the nearest integer
    """
    return int(round(x) - .5) + (x > 0.5)

def unwrap(x) :
    """ Method to unwrap phases, based on MIRIAD subroutines
        input :
            x - a list of the phases to unwrap
        returns :
            the phases unwrapped
    """
    theta = x[0]
    l = [0.0] * len(x)
    l[0] = theta
    for i in range(1,len(x)) :
        l[i] = x[i] - 360*(iround((x[i]-theta)/360.0))
        theta = 0.5*(l[i] + theta)
    return l
