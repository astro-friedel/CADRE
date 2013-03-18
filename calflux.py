from pipeline_miriadwrap import *
import os
import calculations

def calflux(source,freq) :
    """ Method to get a calibrator flux from the miriad catalog
        inputs :
            source - the name of the source to search for
            freq - the observing frequency
        returns :
            the flux of the object
    """
    delfreq = 30.0
    data = dict()
    dir = os.getenv("MIR")
    input = open(dir + "/cat/FluxSource.cat",'r')
    lines = input.readlines()
    input.close()
    lines.reverse()
    found = False
    while(len(lines) > 0) :
        line = lines.pop()
        if("##" in line and source in line) :
            found = True
        elif(found and not "#" in line) :
            splitline = line.split()
            date = splitline[1]
            obsfreq = float(splitline[2])
            flux = float(splitline[3])
            if(obsfreq >= freq - delfreq and obsfreq <= freq + delfreq) :
                date = date.replace("-","")
                data[calculations.gregorianToNumeric(date[2:9])] = flux
        elif(found and "##" in line) :
            found = False
    return data
