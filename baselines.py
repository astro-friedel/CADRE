import logger as log
import calculations
import globals

"""
Module to apply appropriate baseline solution to data
Part of the CARMA data reduction pipeline
Author: D. N. Friedel
"""
range = 60./365.
exact = False

def readAntpos(file) :
    """ Method to read in a baseline file and return a dictionary of antenna positions in ns
        input :
            file - the baseline file
        returns :
            dictionary containing the antenna positions
    """
    input = open(file,'r')
    antposlist = input.readlines()
    input.close()
    antposlist.reverse()
    line = antposlist.pop()
    Mantpos = dict()
    ant = 0
    while(len(antposlist) > 0) :
        line = antposlist.pop()
        if(not "#" in line) :
            ant += 1
            loc = []
            splitLine = line.split()
            loc.append(float(splitLine[0]))
            loc.append(float(splitLine[1]))
            loc.append(float(splitLine[2]))
            Mantpos[ant] = loc
    return Mantpos

def applyBaselines(dataFile, obsDate) :
    """ Method to search for the appropriate antenna position file and update the
        baselines in the data
        input :
            dataFile - the miriad data file
            obsDate - the observation date
        returns :
            none
    """
    global exact
    # get all antpos files
    log.writeComment("Applying best baseline solution to data")
    sys = log.run("ls -1 $MIRCAT/baselines/carma/antpos.* > antpos.list",[],logit=False)
    if(sys != 0) :
        log.writeComment("Failed to find antpos files, baslines not updated, update miriad via cvs with the -d option to get all new directories, this will add the antpos files")
        return
    input = open("antpos.list",'r')
    fileList = input.readlines()
    input.close()
    fileList.reverse()
    files = dict()
    # sort through the list to find all baseline solutions within 60 days
    while(len(fileList) > 0) :
        line = fileList.pop().rstrip("\n")
        date = calculations.numericalGregorianToNumeric(line[-6:])
        if(date <= obsDate + range and date >= obsDate - range) :
            files[date] = line
    if(len(files) == 0) :
        log.writeComment("Could not find appropriate updated baseline solution, no changes made.")
        return
    antpos=globals.antennaPositions()
    solution = False
    baselineFile = ""
    # iterate through all found solutions from the most recent back, to find the best one
    while(not solution and len(files) > 0) :
        baselineFile = files.get(max(files))
        del files[max(files)]
        mantpos = readAntpos(baselineFile)
        match = True
        for i in antpos :
            result,exact = comparePositions(antpos[i],mantpos[i])
            match = match and result
        if(match) :
            solution = True
    if(not solution) :
        log.writeComment("Could not find appropriate updated baseline solution, no changes made.")
        return
    # update the baselines (no harm is done if the same baseline solution is being applied to the data as is already in the data)
    if(not exact) :
        baselineFile = "$MIR" + baselineFile[baselineFile.find("/cat/baselines"):]
        log.writeComment("Updating baseline solution in %s" % (dataFile))
        args = []
        args.append(globals.Variable("vis",dataFile))
        args.append(globals.Variable("apfile",baselineFile))
        args.append(globals.Variable("out",dataFile,".temp"))
        sys=log.run("uvedit",args)
#        sys = log.run("uvedit vis=%s apfile=%s out=%s" % (dataFile,baselineFile,dataFile + ".temp"))
        if(sys == 0) :
            args = []
            args.append(globals.Variable(None,dataFile,"; sleep 3; mv "))
            args.append(globals.Variable(None,dataFile,".temp"))
            args.append(globals.Variable(None,dataFile))
            log.run("rm -rf",args)
#            log.run("rm -rf %s; sleep 3; mv %s %s" % (dataFile,dataFile + ".temp",dataFile))
    else :
        log.writeComment("Baseline solution is already the best one, no modifications done.")

def comparePositions(p1,p2) :
    """ Method to determine if the selected antpos file is in the same array as the data
        The comparison is to make sure that no antenna position has differed by more than 0.5 ns
        input :
            p1 - antenna position 1
            p2 - antenna position 2
        returns :
            True/False if the antenna is/is not in the same array position
    """
    global exact
    exact = ((p1[0] - p2[0]) + (p1[1] - p2[1]) + (p1[2] - p2[2]) == 0.0)
    return (abs(p1[0] - p2[0]) < 0.5 and abs(p1[1] - p2[1]) < 0.5 and abs(p1[2] - p2[2]) < 0.5),exact
