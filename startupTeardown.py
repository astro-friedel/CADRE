import os
import globals
from os.path import exists, join
import logger as log
import time
import version
import fileinput
import sys

"""
Module for the common startup/teardown routines for all pipeline data reduction
Part of the CARMA data reduction pipeline
Current Version (this file): 1.1
Released with pipeline version: 1.1
Author: D. N. Friedel
"""

# get the absolute path to the preferences file
filePath = "defaultPreferences.py"

def prepend(file_name) :
    """ Method that prepends the variable names to the output script
        input:
            file_name - the name of the file to prepend the data to
    """
    # open the file
    fobj = fileinput.FileInput(file_name, inplace=1)
    # read in the first line and save it
    first_line = fobj.readline()
    # write out the global variables
    line = "# Global variable definitions\n"
    for key in globals.scriptVarList :
        line += "set %s = \"%s\"\n" % (globals.scriptVarList[key],key)
    sys.stdout.write("%s\n\n%s\n\n" % (first_line, line))
    # write out the rest of the script
    for line in fobj:
        sys.stdout.write("%s" % line)
    sys.stdout.write("\n\n")
    # close it up
    fobj.close()
    # make the script executable
    os.system("chmod u+x %s" % (file_name))

def startup(visFile,args) :
    """ Method that starts the common processes for all data reduction
        input:
            visFile - the name of the main visibility file
    """
    # start logging system
    log.openScriptLog("reductionScript.csh")
    log.openMessageLog("reduction.notes")
    version.versionCheck()
    log.run("rm -rf sources passcals gaincals fluxcals fits",[])
    if(len(args) > 0) :
        log.run("tar",args,logit=True,fatal=True)
    log.writeAll("\n")
    # add notes to the history file about the data reduction
    args = []
    args.append(globals.Variable("in",visFile))
    args.append(globals.Variable("comment","'Data reduction by CARMA pipeline version %s started at %s with preferences:'" % (version.VERSION,time.ctime())))
    args.append(globals.Variable("file",filePath))
    args.append(globals.Variable("lines","6,100"))
    log.run("addhis",args,logit=False)

    log.writeAll("\n")
    log.writeLog("Data reduction by CARMA pipeline version %s started at %s\n" %(version.VERSION,time.ctime()))

def endFile(file) :
    """ Method to put a comment at the end of a file at the completion of a reduction run
        input:
            file - the visibility file name to write to
    """
    log.run("addhis in=%s comment='Data reduction by CARMA pipeline version %s completed at %s'" % (file,version.VERSION,time.ctime()),[],logit=False)

def transferFile(dir,fileName) :
    """ Method to transfer all completed files to subdirectories for easier viewing
        input :
            dir - the directory to transfer the files to
            fileName - string containg the file to transfer (can contain * as it is passed to os.sysem)
        returns :
            none
    """
    log.run("mkdir -p %s" % (dir),[])
    log.writeAll("\n")
    args = []
    args.append(globals.Variable(None,fileName,"*"))
    args.append(globals.Variable(None,dir + "/."))
    log.run("\mv",args)

def cleanup(objects,fail = False) :
    """ Method to close up the logs and clean up temporary files
        objetcs - objects listing
        fail - whether the data reduction failed or not
    """
    if(globals.isSci2) :
        log.run("rm -rf *.log junk antpos.list bw.sel rms.* flux.log.* fit.log* *.map",[],logit=True)
    else :
        log.run("rm -rf *.log junk antpos.list bw.sel rms.* flux.log.* fit.log* *.beam *.clean *.map *.bm",[],logit=True)

    # move the data files to directories based on purpose
    log.run("mkdir -p fits",[],logit=False)
    log.run("\mv *.fits fits/.",[],logit=False)
    for obj in objects._sources :
        transferFile("sources",obj._file)
    for obj in objects._fluxcals :
        transferFile("fluxcals",obj._file)
    for obj in objects._passcals :
        transferFile("passcals",obj._file)
    for obj in objects._gaincals :
        transferFile("gaincals",obj._file)

    # close everything up
    del objects
    log.writeLog("Data reduction by CARMA pipeline version %s completed at %s" % (version.VERSION,time.ctime()))
    log.closeLog()
    prepend(log.closeScript())
    if(not fail) :
        print "\nData reduction complete"
    else :
        print "\nData reduction failed"
