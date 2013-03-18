import os
import globals

scriptLog = None
messageLog = None

"""
Module for general logging
Part of the CARMA data reduction pipeline
Current Version (this file): 1.1
Released with pipeline version: 1.1
Author: D. N. Friedel
"""

def writeHeader(output) :
    """ Method to write header information in the logs  in the format:
        ---------------------------------------------------------------------

        <output>

        ---------------------------------------------------------------------
        input :
            output - what to write
        returns :
            none
    """
    global scriptLog
    global messageLog
    scriptLog.write("\n\n#---------------------------------------------------------------------\n#\n")
    messageLog.write("\n\n----------------------------------------------------------------------\n\n")
    for line in output :
        scriptLog.write("#\t" + line + "\n")
        messageLog.write("\t" + line + "\n")
    scriptLog.write("#\n#---------------------------------------------------------------------\n\n")
    messageLog.write("\n----------------------------------------------------------------------\n\n")
    scriptLog.flush()
    messageLog.flush()

def writeAll(msg) :
    """ Method to write a string to all log files
        input :
            msg - the mesage to write
        returns :
            none
    """
    global scriptLog
    global messageLog
    scriptLog.write(msg)
    messageLog.write(msg)

def openScriptLog(fileName) :
    """ Method to open the script log
        input :
            fileName - the name of the log file
        returns :
            none
    """
    global scriptLog
    scriptLog = open(fileName,'w')
    scriptLog.write("#!/bin/csh -fe\n")
    scriptLog.write("\nif ( $?MIR == 0 ) then\n  echo \"Cannot find MIRIAD in your environment\"\n  exit 1\nendif\n")
    scriptLog.write("\nif ( ! -e $MAINFILE ) then\n    if(  ! -e $MAINFILE.tar ) then\n        if( ! -e $MAINFILE.tar.gz ) then\n            echo \"Can't find data file  $MAINFILE to process\"\n            exit 1\n        endif\n    endif\nendif\n\n")
    scriptLog.flush()

def writeScript(message,note = None) :
    """ Method to write 'message' to the script log and general log
        input :
            message - what to write to the logs
        returns :
            none
    """
    global scriptLog
    global messageLog
    scriptLog.write(message + "\n")
    scriptLog.flush()
    if(note == None) :
        messageLog.write("COMMAND: " + message + "\n")
    else :
        messageLog.write("COMMAND: " + note + "\n")
    messageLog.flush()

def writeComment(message) :
    """ Method to write a comment to the script and general log
        input :
            message - the comment to write
        returns :
            none
    """
    global scriptLog
    global messageLog
    scriptLog.write("# " + message + "\n")
    scriptLog.flush()
    messageLog.write(message + "\n")
    messageLog.flush()

def closeScript() :
    """ Method to close the script log
        input :
            none
        returns :
            none
    """
    global scriptLog
    file = scriptLog.name
    scriptLog.close()
    return file

def openMessageLog(fileName) :
    """ Method to open the general log
        input :
            fileName - the name of the log file
        returns :
            none
    """
    global messageLog
    messageLog = open(fileName,'a')

def writeLog(message) :
    """ Method to write to the general log
        input :
            message - whtat to write to the log
        returns :
            none
    """
    global messageLog
    messageLog.write(message + "\n")
    messageLog.flush()

def closeLog() :
    """ Method to close the message log
        input :
            none
        returns :
            none
    """
    global messageLog
    messageLog.close()

def run(command, args,  fatal=False, logit=True, execute=True) :
    """ Method to run a command line program
        input :
            command - the command to run
            args - a dictionary of the arguments to the command
            fatal - if an error occurs is it fatal
            logit - log the results of the command line operation
        returns :
            the return state of the command
    """
    commandLine = command
    logLine = command
    # parse the input arguments a form the command and log lines
    for arg in args :
        if("." + globals.GAINCALEND in arg.getArg()) :
            arg.prependPostfix("." + globals.GAINCALEND)
            arg.setArg(arg.getArg()[:arg.getArg().index("." + globals.GAINCALEND)])
        elif("." + globals.PASSCALEND in arg.getArg()) :
            arg.prependPostfix("." + globals.PASSCALEND)
            arg.setArg(arg.getArg()[:arg.getArg().index("." + globals.PASSCALEND)])
        elif("." + globals.FLUXCALEND in arg.getArg()) :
            arg.prependPostfix("." + globals.FLUXCALEND)
            arg.setArg(arg.getArg()[:arg.getArg().index("." + globals.FLUXCALEND)])
        commandLine += " "
        logLine += " "
        if(arg.getOption() == None) :
            commandLine += arg.getPrefix() + arg.getArg() + arg.getPostfix()
            if(arg.getArg() in globals.scriptVarList) :
                logLine += arg.getPrefix() + "$" + globals.scriptVarList[arg.getArg()] + arg.getPostfix()
            else :
                logLine += arg.getPrefix() + arg.getArg() + arg.getPostfix()
        elif(arg.getOption() == "ADD") :
            commandLine = commandLine[:-1]
            logLine = logLine[:-1]
            commandLine += "," + arg.getPrefix() + arg.getArg() + arg.getPostfix()
            if(arg.getArg() in globals.scriptVarList) :
                logLine += "," + arg.getPrefix() + "$" + globals.scriptVarList[arg.getArg()] + arg.getPostfix()
            else :
                logLine += "," + arg.getPrefix() + arg.getArg() + arg.getPostfix()
        else :
            commandLine += arg.getOption() + "=" + arg.getPrefix() + arg.getArg() + arg.getPostfix()
            if(arg.getOption() == "refant") :
                logLine += arg.getOption() + "=" + arg.getPrefix() + "$REFANT" + arg.getPostfix()
            elif(arg.getArg() in globals.scriptVarList) :
                logLine += arg.getOption() + "=" + arg.getPrefix() + "$" + globals.scriptVarList[arg.getArg()] + arg.getPostfix()
            else :
                logLine += arg.getOption() + "=" + arg.getPrefix() + arg.getArg() + arg.getPostfix()
    # run the command line call
    sys = 0
    if(execute) :
        sys = os.system("""%s""" % (commandLine))
    # if the call failed
    if(sys != 0 and (logit or fatal)) :
        if(logit):
            writeLog("FAILED: %s" % (logLine))
        # if the pipeline should stop
        if(fatal) : # throw an exception
            raise Exception, "Command failed: %s" % (logLine)
    else:
        if(logit) :
            writeScript(logLine,commandLine)
    return sys
