import logger as log
import uvflux
import globals

"""
Module to determine if decorrelation correction is needed and apply it if so
Part of CARMA data reduction pipeline
Author: D. N. Friedel
"""

# import the preferences file
try :
    import preferences as p
except ImportError:
    import defaultPreferences as p

def correctDecorrelation(objects,refant) :
    """ Method to correct the miriad data for decorrelation from atmospheric phase fluctuations
        input :
            objects - the objects observed
            refant - the reference antenna
        returns :
            True/False - if it is successful or not
    """
    if(not p.preferences.get("doDecorrelation")) :
        return False
    log.writeComment("Determining if data are decorrelated on long baselines")
    win = objects._gaincals[0]._bandwidths.index(objects._gaincals[0].getMaxBw()) + 1
    numChans = objects._gaincals[0]._numChans[win - 1]
    startChan = 1
    if(objects._gaincals[0]._bandwidths[win - 1] == 62) :
        numChans -= 6
        startChan = 4
    log.run("rm -rf junk.*; sleep 3; uvcat vis=%s select=window'('%i')' out=junk.temp" % (objects._gaincals[0]._file,win),[],logit=False,fatal=False)

    log.run("mselfcal vis=junk.temp interval=%f options=phase line=chan,1,%i,%i refant=%i" % (p.preferences.get("selfcalInterval"),startChan,numChans,refant),[],logit=False,fatal=False)

    f = uvflux.uvflux("junk.temp",startChan,numChans)[0]

    if(len(f) == 0) :
        return False

    decor = f[3].real / f[6]
    if(decor > 1.25) :
        log.writeComment("Could not calculate decorrelation, no correction made")
        log.run("rm -rf junk.* flux.log",[],logit=False,fatal=False)
        return False
    if(decor >= 0.75) :
        log.writeComment("Data are not decorrelated, no corrections needed.")
        log.run("rm -rf junk.* flux.log",[],logit=False,fatal=False)
        return False
    log.run("uvdecor vis=junk.temp options=nocal,nopass,nopol delaymax=8500 out=junk.temp2",[],logit=False,fatal=False)

    log.run("mselfcal vis=junk.temp2 interval=%f options=phase line=chan,1,%i,%i refant=%i" % (p.preferences.get("selfcalInterval"),startChan,numChans,refant),[],logit=False,fatal=False)
    f = uvflux.uvflux("junk.temp2",startChan,numChans)[0]

    if(len(f) == 0) :
        return False
    decor2 = f[3].real / f[6]
    if(decor2 > 1.25 or decor2 <= decor) :
        log.writeComment("Correction for decorrelation failed.")
        log.run("rm -rf junk.* flux.log",[],logit=False,fatal=False)
        return False

    return True
