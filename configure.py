import os
ok = False
try :
    import flagging
    ok = True
except :
    pass

if(not ok) :
    print "Configuring the pipeline for your system (this should only need to be run once)"
    try :
        from pipeline_miriadwrap import *
    except:
        print "Could not locate the python MIRIAD wrappers"
        print "Please install them, the install script is located in $MIR/src/scripts/python/subwrap"
        print "  and be sure the libraries (usually in $MIRLIB/python) are in your python path"
        os._exit(0)
    try :
        handle, iostat = hopen("testMiriadFile", "old")
        if(iostat != 0):
            raise Exception, "File %s not found" % (file)
        if(not hdprsnt(handle, "bandpass")):
            raise Exception,  "No bandpass present in %s" % (file)
        nfeeds = rdhdi(handle,"nfeeds",1)
        ngains = rdhdi(handle,"ngains",1)
        ntau = rdhdi(handle,"ntau")
        nchan = rdhdi(handle,"nchan0")
        nspect = rdhdi(handle,"nspect0")
        if(nfeeds <= 0 or ngains <= 0) :
            raise Exception, "Bad gain table size information"
        nants = ngains / (nfeeds+ntau)
        if(nants*(nfeeds+ntau) != ngains) :
            raise Exception, "Number of gains does equal nants*nfeeds"
        if(nchan <= 0) :
            raise Exception, "Bad number of frequencies"
        if(nspect <= 0 or nspect > nchan) :
            raise Exception, "Bad number of frequency spectral windows"
        item, iostat = haccess(handle,"freqs","read")
        if(iostat != 0) :
            raise Exception, "Error accessing the bandpass frequency table"
        off = 8
        nschan,iostat = hreadi(item,off,4)
        print nschan,iostat
        if(iostat != 0 or nschan[0] != 15) :
            os.system("mv flagging32.py flagging.py")
        else :
            os.system("mv flagging64.py flagging.py")
    except :
        print "Could not configure for your system"
        os._exit(0)



