from pipeline_miriadwrap import *
import miriadClasses
import time
import math
import calculations
from threading import Thread
import struct
import os
import shutil
import time
import threadpool as tp
import sources
import numpy

"""
Module containing basic miriad I/O routines written in python
Part of the CARMA data reduction pipeline
Current Version (this file): 1.1
Released with pipeline version: 1.1
Author: D. N. Friedel
"""

radToSec = 206264.806247
AXIS_FREQ = "FREQ"
AXIS_VEL = "VEL"

beam = [0,0,0.0]    #beam size (in pixels) and pa

keywords = ["btype","bpa","bmin","bmaj","niters","bunit","vobs","epoch","cdelt4","cdelt3","cdelt2","cdelt1","crval4","crval3","crval2","crval1","ctype4","ctype3","ctype2","ctype1","crpix4","crpix3","crpix2","crpix1","lstep","lwidth","lstart","ltype","restfreq","telescop","object","naxis7","naxis6","naxis5","naxis4","naxis3","naxis2","naxis1","naxis","history"]

def getHeader(handle):
    """Method to get all header information from a miriad image
        inputs :
            handle - the handle to the miriad image
        returns :
            a Header object containg all header information
    """
    global keywords
    header = miriadClasses.Header()
    for keyword in keywords :
        desc,type,n = hdprobe(handle,keyword)
        if(n != 0) :
            value = ""
            type = type.lower()
            if("integer" in type) :
                if("*8" in type) :
                    value = rdhdl(handle,keyword)
                else :
                    value = rdhdi(handle,keyword)
            elif("real" in type) :
                value = rdhdr(handle,keyword)
            elif("double" in type) :
                value = rdhdd(handle,keyword)
            elif("complex" in type) :
                value = rdhdc(handle,keyword)
            elif("char" in type or "text" in type) :
                value = rdhda(handle,keyword)
            temp = miriadClasses.ValueDescPair(value,desc,type)
            header.add(keyword,temp)
    return header

def minDistance(point,pointings) :
    """ Method to determine which pointing is closest to a given point
        input :
            point - the reference point
            pointings - a list of pointings from which the one closest to point is calculated
        returns :
            the distance to the closest pointing from point
    """
    minD = 100000000.0
    for p in pointings :
        minD = min(minD,calculations.distance(p,point))
    return minD

def int2bin(n, count=24):
  """returns the binary of integer n, using count number of digits"""
  return "".join([str((n >> y) & 1) for y in range(count-1, -1, -1)])

def getPlane(file,plane,header,region,radius,pointings,doMask,cutoff=-1000.0) :
    """ Method to retrieve a plane from an image
        inputs :
            handle - handle to the miriad image
            plane - the plane number to retrieve
            xrange - the width of the image
        returns :
            the plane as a 2D array
    """
    # don't need to apply a mask if there isn't one
    if(not os.path.exists(file + "/mask") and doMask) :
        doMask = False
    t0 = time.time()
    handle = open(file + "/image")
    img = handle.read(4)
    size = struct.unpack(">l",img)[0]
    x = header.getValue("naxis1")
    y = header.getValue("naxis2")
    type = ">"
    if(size == 4) :
        for i in range(0,x) :
            type += "f"

    elif(size == 8) :
        for i in range(0,x) :
            type += "d"
    handle.seek((plane-1)*size*x*y,os.SEEK_CUR)
    #read the mask
    mask = []
    if(doMask) :
        mhandle = open(file + "/mask")
        m = mhandle.read(4)
        tmask = ""
        msize = int(x*y/31) + 2
        loc = int(math.floor(((plane-1)*x*y)/31))*4
        mhandle.seek(loc,os.SEEK_CUR)
        place = ((plane-1)*length*height) % 31
        for i in range(0,msize) :
            try :
                m = mhandle.read(4)
                num = struct.unpack(">i",m)[0]
                temp= int2bin(num,32)[::-1]
                tmask += temp[place:-1]
                place = 0
            except :
                pass
        for i in range(0,y) :
            trow = []
            for j in range(0,x) :
                if(tmask[i*x + j] == "0") :
                    trow.append(False)
                else :
                    trow.append(True)
            mask.append(trow)
        mhandle.close()
    image = []
    halfx = x/2
    halfy = y/2
    startx = 0
    starty = 0
    stopx = halfx*2
    stopy = halfy*2
    if(region == "quarter") :
        startx = int(halfx/2)
        starty = int(halfy/2)
        stopx = startx + halfx
        stopy = starty + halfy
        handle.seek(starty*x*size)
    # start reading in the image data
    for row in range(starty,stopy) :
        rawData = handle.read(size*x)
        data = list(struct.unpack(type,rawData))
        if(doMask) :
            maskRow = mask[row]
            # clip the data if masked or if there is a cutoff or if the point is outside the selected region
            if(cutoff > -999.9) :
                for point in range(startx,stopx) :
                    if(not maskRow[point] or data[point] < cutoff or minDistance([row,point],pointings) > radius) :
                        data[point] = 0.0
            else :
                for point in range(startx,stopx) :
                    if(not maskRow[point] or minDistance([row,point],pointings) > radius) :
                        data[point] = 0.0
        else :
            # clip the data if there is a cutoff or if the point is outside the selected region
            if(cutoff > -999.9) :
                for point in range(startx,stopx) :
                    if(data[point] < cutoff or minDistance([row,point],pointings) > radius) :
                        data[point] = 0.0
            else :
                for point in range(startx,stopx) :
                    if(minDistance([row,point],pointings) > radius) :
                        data[point] = 0.0
        image.append(data[startx:stopx + 1])
    handle.close()
    t1 = time.time()
    return plane,image


def getResult(request,result) :
    """ Helper method for threading
    """
    request.result[result[0]] = result[1]

def putImage(file,header,image,orig,doMask=True) :
    """ Method to write a miriad image to disk
        input :
            file - name of the output file
            header - a header object
            image - an array containing the image data to be written
            orig - name of the origianting file
            doMask - whether to apply a mask
        returns :
            none
    """

    # don't need to apply a mask if there isn't one
    if(not os.path.exists(file + "/mask") and doMask) :
        doMask = False
    if(os.path.exists(file)) :
        print "FILE EXISTS"
        return
    # do the intialization
    os.makedirs(file)
    shutil.copyfile(orig + "/header", file + "/header")
    shutil.copyfile(orig + "/history", file + "/history")
    handle = open(file + "/image",'w')
    size = struct.pack(">l",4)
    handle.write(size)
    x = header.getValue("naxis1")
    y = header.getValue("naxis2")
    type = ">f"

    # write the image
    min = 100000.0
    max = -100000.0
    for plane in image :
        for row in range(0,y) :
            data = ""
            for j in range(0,x) :
                if(plane[j][row] > max) :
                    max = plane[j][row]
                if(plane[j][row] < min) :
                    min = plane[j][row]
                data += struct.pack(type,plane[j][row])
            handle.write(data)
    handle.close()
    return

def putHeader(handle,header,h2):
    """Method to copy all header info from one miriad image
       to another
        inputs :
            handle - the handle to the miriad image
            header - the header object from the handle image
            h2 - the handle for the image to copy to
        returns :
            none
    """
    global keywords

    for item in keywords :
        type = header.getType(item)
        if(type != "none") :
            hdcopy(h2,handle,item)

def imageToArray(file, region,pointings=[0.0,0.0],doMask = True,isBeam = False,cutoff=-1000.0,withHeader = False) :
    """ Method to convert a miriad image into a python array
        input :
            file - the file to convert
            region - the selection region
        returns :
            the image as a 2 or 3D array
    """
    global beam
    global radToSec
    handle = -1
    imageArray = None
    t0 = time.time()
    header = None
    try :
        cdelt = [0.0,0.0,0.0]
        blc = [0,0,0,0,0,0,0]
        trc = [0,0,0,0,0,0,0]
        handle,num = xyopen(file,"old")
        header = getHeader(handle)
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
            for p in pointings :
                offPoints.append([int(p[0]/cdelt[0]) + center[0],int(p[1]/cdelt[1]) + center[1]])

        # if there is only 1 plane then do not use threads, it is a waste of time
        if(header.getValue("naxis3") == 1) :
            t,img = getPlane(file,1,header,region,int(calculations.calcImsize(freq,cdelt[1],0.0)/(2.0*1.639)),offPoints,doMask,cutoff)
            image.append(img)
        else :
            worker = tp.ThreadPool(20)
            result = dict()
            #set up all arguments for all threads
            arguments = [((file,plane,header,region,int(calculations.calcImsize(freq+((plane - 1) * interval),cdelt[1],0.0)/(2.0*1.639)),offPoints,doMask,cutoff), {}) for plane in range(1,header.getValue("naxis3") + 1)]
            # generate the requested work
            requests = tp.makeRequests(getPlane,arguments,getResult,result=result)
            # fill the thread pool and start working
            for req in requests:
                worker.putRequest(req)
            # wait for all threads to finish before continuing
            worker.wait()
            for plane in range(1,header.getValue("naxis3") + 1) :
                image.append(result[plane])
            del result
        imageArray = numpy.array(image)
        imageArray = imageArray.swapaxes(1,2)  # rotate the cube so that the order is z,x,y  (0,0,0) is top left of first plane
        del image
    except :
        if(handle != -1) :
            xyclose(handle)
        raise
    t1 = time.time()
    #print int((t1 - t0) * 1000)
    if(withHeader) :
        return imageArray,header
    return imageArray

def readImage(file) :
    header = None
