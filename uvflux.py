from pipeline_miriadwrap import *
import calculations
import math

"""
Python version of the MIRIAD task uvflux.
Current Version (this file): 1.1
Released with pipeline version: 1.1
Author: D. N. Friedel
"""

def uvflux(file,start,step) :
    """ Python port of MIRIAD task uvflux
        input :
            file - the name of the file to operate on
            start - starting channel number
            step - step size (in channels)
    """
    doshift = False
    uvflags = "sdlcef"
    handle = uvopen(file,"old")
    uvset(handle,"data","channel",1,float(start),float(step),float(step))
    PolIndx = [0] * 14
    sources = ["  "] * 14
    isrc = -1
    nsrc = 0
    npol = 0
    fluxr = []
    fluxi = []
    amp = []
    amp2 = []
    rms2 = []
    ncnt = []
    PolMin = -9
    PolMax = 4
    for i in range(0,MAXPOL) :
        temp = [0.0] * MAXSRC
        fluxr.append(temp)
        temp = [0.0] * MAXSRC
        fluxi.append(temp)
        temp = [0.0] * MAXSRC
        amp.append(temp)
        temp = [0.0] * MAXSRC
        amp2.append(temp)
        temp = [0.0] * MAXSRC
        rms2.append(temp)
        temp = [0] * MAXSRC
        ncnt.append(temp)
    preamble,data,flags = uvread(handle)
    ipol = -20
    while(len(flags) > 0) :
        ipol = uvrdvri(handle,"pol",1)
        if(PolIndx[ipol] == 0) :
            npol += 1
            PolIndx[ipol] = npol
        ipol = PolIndx[ipol]
        t,l,update = uvprobvr(handle,"source")
        if(update) :
            source = uvgetvra(handle,"source")
            found = False
            if(isrc >= 0) :
                found = source == sources[isrc]
            if(not found) :
                if(source in sources) :
                    isrc = sources.index(source)
                    found = True
            if(not found) :
                nsrc += 1
                sources[nsrc - 1] = source
                for i in range(0,MAXPOL) :
                    fluxr[i][nsrc] = 0.0
                    fluxi[i][nsrc] = 0.0
                    amp[i][nsrc]  = 0.0
                    amp2[i][nsrc] = 0.0
                    rms2[i][nsrc] = 0.0
                    ncnt[i][nsrc] = 0
                isrc = nsrc-1
        sig2 = uvinfo(handle,"variance")[0]
        for i in range(0,len(flags)) :
            if(flags[i]) :
                fluxr[ipol][isrc] += data[i].real
                fluxi[ipol][isrc] += data[i].imag
                rms2[ipol][isrc] += sig2
                temp = calculations.cabs(data[i])
                amp[ipol][isrc]  += temp
                amp2[ipol][isrc] += temp*temp
                ncnt[ipol][isrc] += 1
        preamble,data,flags = uvread(handle)
    uvclose(handle)
    npol = 0
    p = [0] * MAXPOL
    pp = [0] * MAXPOL
    for j in range(PolMin,PolMax+1) :
        if(PolIndx[j] > 0) :
            p[npol] = j
            pp[npol] = PolIndx[j]
            npol = npol + 1
            for i in range(npol-1,1,-1) :
                print i
                if(abs(p[i]) < abs(p[i-1])) :
                    t = p[i]
                    p[i] = p[i-1]
                    p[i-1] = t
                    t = pp[i]
                    pp[i] = pp[i-1]
                    pp[i-1] = t
    PolCode = "--"
    retVal = []
    for isrc in range(0,nsrc) :
        source = sources[isrc]
        for i in range(0,npol) :
            ipol = pp[i]
            if(ncnt[ipol][isrc] > 0) :
                PolCode = polsc2p(p[i])
                fluxr[ipol][isrc] /= ncnt[ipol][isrc]
                fluxi[ipol][isrc] /= ncnt[ipol][isrc]
                vecaver  = complex(fluxr[ipol][isrc],fluxi[ipol][isrc])
                vecscat  = amp2[ipol][isrc] / (2*ncnt[ipol][isrc])- 0.5*(fluxr[ipol][isrc]**2+fluxi[ipol][isrc]**2)
                vecscat = math.sqrt(abs(vecscat))
                vecamp,vecpha = amphase(vecaver)
                scalamp = amp[ipol][isrc] / ncnt[ipol][isrc]
                scalscat = amp2[ipol][isrc] / ncnt[ipol][isrc]- (amp[ipol][isrc] / ncnt[ipol][isrc])**2
                scalscat = math.sqrt(abs(scalscat))
                sig2 = math.sqrt(rms2[ipol][isrc]/ncnt[ipol][isrc])
                retVal.append([source,PolCode,sig2,complex(fluxr[ipol][isrc],fluxi[ipol][isrc]),vecscat,scalamp,scalscat,ncnt[ipol][isrc]])
    return retVal
