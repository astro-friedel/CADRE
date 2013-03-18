"""
Default preferences file for CARMA data reduction pipeline.
Author: D. N. Friedel
Do not modify the format of this file, additional parameters can be added, but the preferences = line must start on line 6
"""
preferences = {"tsysThreshold" : 5000.0,           # any data with a system temp above this is flagged
               "BIMAShadowFraction" : 0.9,         # shadowing fraction for 6.1 m dishes (i.e. if this fraction of the radius is shadowed then flag it)
               "OVROShadowFraction" : 1.0,         # shadowing fraction for 10.4 m dishes (i.e. if this fraction of the radius is shadowed then flag it)
               "SZAShadowFraction" : 1.0,          # shadowing fraction for 3 m dishes (i.e. if this fraction of the radius is shadowed then flag it)
               "doBaselines" : True,               # apply the proper baseline solution
               "doDecorrelation" : False,
               "selfcalInterval" : 5.0,            # interval for self calibration averaging
               "bootfluxInterval" : 5.0,           # interval for bootflux averaging
               "amplitudeGainRange" : [0.2, 5.0],  # if the amplitude gains of an antenna are outside of this range then the corresponding data are flagged
               "maxAmplitudeGainFactor" : 3.0,     # if any antenna mean gain is more than this factor above the general mean then the antenna is flagged
               "maxGainRms" : 1.0,                 # if gain rms from an antenna is above this value then the antenna is flagged
               "cellSize" : -1.0,                  # the cell size used for invert in arcsec, a negative value lets the system decide based on the average baseline length
               "imageSize" : -1.0,                 # the image size for invert in pixels, a negative value lets the system decide based on the cell size and primary beam size
               "cleanThreshold" : 5,               # how deep clean should iterate when compared to the rms of the map the number is a scale factor (1.5 = clean to 1.5 * rms)
               "cleanRegion" : "quarter",          # what size clean region should be used, "quarter" = inner quarter, or give size (radius) in arcseconds, for mosaics the region is automatically generated
               "doContinuumSubtraction" : False,   # do automated continuum subtraction
               "doAutoCleanRegion" : False         # automatically determine the clean region
}
