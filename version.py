import calculations
import logger as log

""" Version number for the CARMA data reduction pipeline
Part of the CARMA data reduction pipeline
Current Version (this file): 1.1
Released with pipeline version: 1.1
Author: D. N. Friedel
"""
VERSION = "1.1"
LAST_MOD_DATE = "Aug 1, 2011"

def version() :
    print ""
    print "CARMA data reduction pipeline version %s" % (VERSION)
    print "   %s" % (LAST_MOD_DATE)
    print ""

def help() :
    version()
    print ""
    print "Usage :"
    print "  python <path>/pipeline.py <DATA> [type=[c/s]]"
    print "    <DATA> = name of miriad data set to be reduced"
    print "             this can also be a tar or tar.gz file containing a miriad data set"
    print "             it is recommended to run this from the directory containing the miriad data set"
    print "    c/s = which type of reduction to do"
    print "          [c]ontinuum or [s]pectral line default is spectral line"
    print "Data reduction repferences can be given in a file called preferences.py in the"
    print "  current working directory. Use the defaultPreferences.py file as a template."
    print "Two files are created by this pipeline:"
    print "  reduction.notes - a description of what was done to the data set"
    print "  reduction.script - a standalone script that can be run to reproduce the data reduction steps taken by the pipeline"
    print "Note: This pipeline requires the most recent version of MIRIAD be installed and in your path."
    print ""

def versionCheck() :
    """ Method to make sure the proper version of listobs and gplist are installed
    """
    sys = log.run("ls -1 $MIR >/dev/null",[],logit=False)
    if(sys != 0) :
        raise Exception, "Failed to find MIRIAD. Please make sure the latest version is installed and in your path."
    log.run("more $MIR/VERSION > junk",[], logit=False)
    input = open("junk","r")
    lines = input.readlines()
    input.close()
    version = lines[0].split(".")
    if(version[0] < 4) :
        raise Exception, "Your MIRIAD version is %i, version 4 or newer is required by the pipeline, please updateyour MIRIAD distribution."
