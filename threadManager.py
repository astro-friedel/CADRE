import threading
import globals
import time

"""
Module to handle threading for the pipeline.
This module will only allow as many threads as their are processors so that
the system does not get bogged down
It also will limit the threads if large amounts of memory are needed
(this is specifically for invert and clean)
Author: D. N. Friedel
"""

def thread(function,args,enumerator,position,watchMemory = False,imgSize = [0,0,0]) :
    """ Method to smartly handle threads
        input :
            function - the function being threaded
            args - a list of the arguments to the function
            enumerator - a list of the objects to enumerate over
            position - the position where the enumerated argument should be placed in the args list
            watchMemory - should the thread manager be cognisant of memory management
    """
    maxConcurrentThreads = globals.SystemInfo.getCPUs()
    threadList = []
    if(watchMemory) :
        #get available memory and convert to bytes (it is reported in Mb)
        availableMemory = (globals.SystemInfo.getAvailableMemory())*1024*1024
        # calculate how much memory is needed for each thread
        neededMemory = imgSize[0]*imgSize[1]*imgSize[2] * 4
        # calculate the maximum number of concurrent threads without hitting swap
        maxConcurrentThreads = max(1,min(int(availableMemory/neededMemory),maxConcurrentThreads))
    # populate the thread list and start threads up to the max limit, then wait for available threads before continuing
    for item in enumerator :
        args.insert(position - 1,item)
        current = function(args)
        threadList.append(current)
        while((threading.activeCount() - 1) >= maxConcurrentThreads) :
            time.sleep(5)
        current.start()
    # wait for the all of the threads to complete
    for thread in threadList :
        try :
            thread.join()
        except :
            pass
    del threadList
