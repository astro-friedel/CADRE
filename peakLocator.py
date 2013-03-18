import sources
import numpy
import os
import logger as log
from copy import deepcopy
import calculations
import math
from pipeline_miriadwrap import *
import miriad_functions as mfunc
import time
import threadpool as tp

AXIS_FREQ = "FREQ"
AXIS_VEL = "VEL"

beam = [0,0,0.0]    #beam size (in pixels) and pa
radToSec = 206264.806247

def imageToArray(file, region,pointings=[0.0,0.0],doMask = True,isBeam = False,cutoff=-1000.0,withHeader = False) :
    """ Method to convert a miriad image into a python array
        input :
            file - the file to convert
            region - the selection region
        returns :
            the image as a 2 or 3D array
    """
    global beam
    handle = -1
    imageArray = None
    t0 = time.time()
    header = None
    try :
        cdelt = [0.0,0.0,0.0]
        blc = [0,0,0,0,0,0,0]
        trc = [0,0,0,0,0,0,0]
        print "opening file"
        handle,num = xyopen(file,"old")
        header = mfunc.getHeader(handle)
        cdelt[0] = header.getValue("cdelt1") * radToSec
        cdelt[1] = header.getValue("cdelt2") * radToSec
        cdelt[2] = header.getValue("cdelt3")
        if(not isBeam) :
            beam[0] = math.fabs(header.getValue("bmaj") * radToSec/cdelt[0])
            beam[1] = math.fabs(header.getValue("bmin") * radToSec/cdelt[1])
            beam[2] = header.getValue("bpa")
        restFreq = header.getValue("restfreq")
        naxis = header.getValue("naxis")
        axisType = None
        freq = 0.0
        interval = 0.0
        if("VELO" in header.getValue("ctype3") or "FELO" in header.getValue("ctype3")) :
            axisType = AXIS_VEL
            freq = calculations.velToFreq(header.getValue("crval3"),restFreq)
            interval = calculations.velToFreq(cdelt[2],restFreq)-restFreq
        elif("FREQ" in header.getValue("ctype3")) :
            axisType = AXIS_FREQ
            freq = header.getValue("crval3")
            interval = cdelt[2]
        image = []
        mask = []
        xyclose(handle)
        handle = -1
        offPoints = []
        if(pointings == [0.0,0.0]) :
            offPoints.append([int(header.getValue("crpix1")),int(header.getValue("crpix2"))])
        else :
            center = [int(header.getValue("crpix1")),int(header.getValue("crpix2"))]
            print "CENTER: ",center
            for p in pointings :
                offPoints.append([int(p[0]/cdelt[0]) + center[0],int(p[1]/cdelt[1]) + center[1]])
        print "PONTINGS: ",pointings,offPoints
        # if there is only 1 plane then do not use threads, it is a waste of time
        if(header.getValue("naxis3") == 1) :
            t,img = mfunc.getPlane(file,1,header,region,int(calculations.calcImsize(freq,cdelt[1],0.0)/(2.0*1.639)),offPoints,doMask,cutoff)
            image.append(img)
        else :
            worker = tp.ThreadPool(20)
            result = dict()
            #set up all arguments for all threads
            arguments = [((file,plane,header,region,int(calculations.calcImsize(freq+((plane - 1) * interval),cdelt[1],0.0)/(2.0*1.639)),offPoints,doMask,cutoff), {}) for plane in range(1,header.getValue("naxis3") + 1)]
            # generate the requested work
            requests = tp.makeRequests(mfunc.getPlane,arguments,mfunc.getResult,result=result)
            # fill the thread pool and start working
            for req in requests:
                worker.putRequest(req)
            # wait for all threads to finish before continuing
            worker.wait()
            for plane in range(1,header.getValue("naxis3") + 1) :
                image.append(result[plane])
        imageArray = numpy.array(image)
        imageArray = imageArray.swapaxes(1,2)  # rotate the cube so that the order is z,x,y  (0,0,0) is top left of first plane
        del image
    except :
        if(handle != -1) :
            xyclose(handle)
        raise
    t1 = time.time()
    print int((t1 - t0) * 1000)
    if(withHeader) :
        return imageArray,header
    return imageArray

def getPeaks(request,result) :
    # extending lists is thread safe so no mutexes are requred
    request.result.extend(result)

def findPeaks(imageCube,beamPoints,cutoff) :
    """ Method to find all peaks in the image
        input :
            image - the image to anylize
            beamPoints - the peaks of the synthesized beam (used to remove sidelobes from the peaks)
        returns :
            a list of peak points
    """
    peaks = []
    global beam
    numPlanes = imageCube.shape[0]
    beamSize = int(math.ceil((beam[0] + beam[1])/2.0))
    width = imageCube.shape[1]
    height = imageCube.shape[2]
    if(numPlanes == 1) :
        t0 = time.time()
        image = None
        # select the proper image plane
        image = imageCube[0]
        while(True) :
            peak = findPeak(image)
            if(image[peak[0]][peak[1]] < cutoff) :
                break
            peaks.append(peak)
            for point in beamPoints :
                for x in range(-2*beamSize,2*beamSize - 1) :
                    xpoint = int(peak[0] + point[0] + x)
                    for y in range(-2*beamSize,2*beamSize - 1) :
                        ypoint = int(peak[1] + point[1] + y)
                        if(xpoint >= 0 and xpoint < width and ypoint >= 0 and ypoint < height) :
                            image[xpoint,ypoint] = 0.0
        t1 = time.time()
        print "TIME:  ",(t1-t0)*1000
    else :
        worker = tp.ThreadPool(5)
        arguments = [((imageCube[plane],beamSize,width,height,beamPoints,cutoff), {}) for plane in range(0,numPlanes)]
       # generate the requested work
        requests = tp.makeRequests(threadWork,arguments,getPeaks,result=peaks)
        # fill the thread pool and start working
        for req in requests:
            worker.putRequest(req)
        # wait for all threads to finish before continuing
        worker.wait()
    return peaks

def threadWork(image,beamSize,width,height,beamPoints,cutoff) :
    peaks = []
    while(True) :
        peak = findPeak(image)
        if(image[peak[0]][peak[1]] < cutoff) :
            break
        peaks.append(peak)
        for point in beamPoints :
            for x in range(-2*beamSize,2*beamSize - 1) :
                xpoint = int(peak[0] + point[0] + x)
                for y in range(-2*beamSize,2*beamSize - 1) :
                    ypoint = int(peak[1] + point[1] + y)
                    if(xpoint >= 0 and xpoint < width and ypoint >= 0 and ypoint < height) :
                        image[xpoint,ypoint] = 0.0
    t1 = time.time()
    return peaks

def findPeak(image) :
    """Method to find the peak value in a single plane of an image
        returns the index to the peak
        input :
            the image to be anylized
        returns :
            the index to the highest value pixel
    """
    maxValue = -100000.0
    maxIndex = [0,0]
    if(len(image.shape) == 2) :
        maxAxis = image.argmax(0)
        for pos in maxAxis :
            if(pos != 0) :
                for i in range(0,image.shape[1]) :
                    if(image[pos][i] > maxValue) :
                        maxValue = image[pos][i]
                        maxIndex = [pos,i]
        print maxValue,maxIndex
    else :
        maxAxis = image.argmax(1)
        for pos in maxAxis[0] :
            if(pos != 0) :
                for i in range(0,image.shape[1]) :
                    if(image[0][pos][i] > maxValue) :
                        maxValue = image[0][pos][i]
                        maxIndex = [pos,i]
    return maxIndex

def constructBeam(beamFile) :
    """ Method to find the peaks in the synthesized beam
        input :
            beamFile - the name of the miriad file containing the beam (single plane only)
        returns :
        a list of the peak positions of the beam (including sidelobes)
    """
    global beam
    width = int(max(beam[0],beam[1]))
    beamImage = imageToArray(beamFile,"quarter",doMask=False,isBeam=True,cutoff=0.08)  #cutoff at 8%
    center = math.ceil(beamImage.shape[1]/2)
    peakList = []
    while(True) :
        peak = findPeak(beamImage)   #peak is an index array
        if(beamImage[0][peak[0]][peak[1]] < 0.08) :
            break
        peakList.append([peak[0]-center,peak[1]-center])
        for x in range(-2*width,2*width + 1) :
            for y in range(-2*width,2*width + 1) :
                xpoint = max(0,min(peak[0] + x,beamImage.shape[1]-1))
                ypoint = max(0,min(peak[1] + y,beamImage.shape[2]-1))
                beamImage[0][xpoint, ypoint] = 0.0
    del beamImage
    return peakList

def constructCleanRegion(polys,shape) :
    """ Method to construct a clean region suitable for input into miriad, surrounding all peaks in the image
        inputs :
            polys -  list of polygons surrounding the peak points
            shape - the size of the image
        returns :
            a string containing the clean region
    """
    region = "quarter"  # the default if a suitable region cannot be found
    if(len(polys) < 1) :
        return region
    startLen = len(polys)
    done = False
    while(len(polys) > 1 and not done) :
        thisLen = len(polys)
        inside = [False,False,False,False,False,False,False,False,False]
        overlap = False
        for poly in range(0,len(polys)) :
            if(poly < len(polys) - 1) :    #python does not adapt to changing array lengths, so must do this check
                doneNow = False
                while(len(polys) > 1 and not doneNow) :
                    thisLen = len(polys)
                    inside = [False,False,False,False,False,False,False,False,False]
                    overlap = False
                    for poly2 in range(poly + 1,len(polys)) :
                        for coord in range(0,8) :
                            if(calculations.pointInPoly(polys[poly],polys[poly2][coord][0],polys[poly2][coord][1])) :
                                inside[coord] = True
                                overlap = True
                        if(overlap) :
                            inside[8] = inside[0]
                            calculations.mergePoly(polys,poly,poly2,inside)
                            thisLen = startLen
                            startLen = len(polys)
                            break
                    if(startLen == thisLen) :   #then we have at least 2 disjoint regions
                        doneNow = True
        if(startLen == thisLen) :   #then we have at least 2 disjoint regions
            done = True
    region = ""
    for poly in polys :
        region =  region +"polygon'('"
        for i in poly :
            region = region + "%i,%i," % (max(1,min(i[0],shape[1])),max(1,min(i[1],shape[2])))
        region = region[:-1] + "')',"
    return region[:-1]

def findCleanRegion(file, cutoff, beamFile,pointings) :
    """ Method to find the clean region of a map
        inputs :
            file - the miriad file name for the map
            cutoff - the minimum value to use (generally a few times the noise level)
            beamFile - the miriad file name for the beam file
        returns :
            the clean region as a string
    """
    global beam
    # convert the image to an array for easier manipulation
    image = imageToArray(file,"",pointings,False)
    numPlanes = image.shape[0]
    # convert the beam to an array for easy manipulation
    beamPoints = constructBeam(beamFile)
    peaks = findPeaks(image,beamPoints,cutoff)
    beamSize = math.ceil((beam[0] + beam[1])/2.0)
    polys = []
    # create the polygons from the peak positions
    region = "quarter"
    for peak in peaks :
        poly = calculations.createPoly([0.0,0.0],2.5*beamSize, 1.0,peak,True)
        polys.append(poly)
    if(numPlanes == 1) :
        region = constructCleanRegion(polys,image.shape)
    else :
        region = constructCleanRegion(polys,[image.shape[1],image.shape[2]])
    del image
    return region

def findSpectralPeaks(file,cutoff,beamFile,pointings) :
    """ Method to find the spectral peaks of a map
        inputs :
            file - the miriad file name for the map
            cutoff - the minimum value to use (generally a few times the noise level)
            beamFile - the miriad file name for the beam file
        returns :
            the regions to use in imspect as a list of strings
    """
    global beam
    # convert the image to an array for easier manipulation
    image = imageToArray(file,"",pointings,doMask=False)
    print "IMG1",image[62][57][3]
    numPlanes = image.shape[0]
    # convert the beam to an array for easy manipulation
    beamPoints = constructBeam(beamFile)
    peaks = findPeaks(deepcopy(image),beamPoints,cutoff)
    print "IMG2",image[62][57][3]
    regions = []
    # create the polygons from the peak positions
    width = image.shape[1]
    height = image.shape[2]
    for peak in peaks :
        regions.append(peak)
    return regions,image

def compactCleanRegions(regions,shape) :
    """ Method to compact the clean regions to make the list smaller
        input :
            regions - the list of regions
            shape - the shape of the image
        returns :
            the constructed clean regions
    """
    polys = []
    breakpt = "polygon'('"
    for window in regions :
        region = regions[window]
        while(breakpt in region) :
            poly = []
            startloc = region.find(breakpt)
            endloc = region.find("')'",startloc)
            temp = region[startloc + len(breakpt):endloc].split(",")
            for i in range(0,len(temp),2) :
                poly.append([int(temp[i]),int(temp[i+1])])
            polys.append(poly)
            region = region[endloc+2:]
    return constructCleanRegion(polys,shape)
