import struct
import miriad_functions as mfunc

"""
Module for classes which interact with miriad files
Current Version (this file): 1.1
Released with pipeline version: 1.1
Author: D. N. Friedel
"""

# class for data/description object
class ValueDescPair :
    def __init__(self,data,description,type) :
        self.data_ = data
        self.description_ = description
        self.type_ = type

# class for holding a miriad image header
class Header :
    def __init__(self) :
        self.data_ = dict()

    def add(self,item,pair) :
        """ Method to add a item/value pair
            input :
                item - the item to add
                pair - the value to add
            returns :
                none
        """
        if(isinstance(pair,ValueDescPair)) :
            self.data_[item] = pair

    def getValue(self,item) :
        """ Method to return the value of a set
            input :
                item - the name of the item to return
            returns :
                the data value of the requested item
        """
        if(item in self.data_) :
            return self.data_[item].data_
        return -1

    def getDesc(self,item) :
        """ Method to return the description of the item
            input :
                item - the item whose description is to be returned
            returns :
                the description of the item
        """
        if(item in self.data_) :
            return self.data_[item].description_
        return ""

    def getType(self,item) :
        """ Method to return the type of the item
            input :
                item - the name of the item whose type is to be returned
            returns :
                the type of the item (real, int, etc)
        """
        if(item in self.data_) :
            return self.data_[item].type_
        return "none"

    def printHeader(self) :
        """ Method to print the full header
            input :
                none
            returns :
                none
        """
        for item in self.data_ :
            print item,": ", self.data_[item].data_,self.data_[item].description_,self.data_[item].type_

# class for holding a miriad image
class Image :
    def __init__(self,header = None,image = None, mask = None,pointings = [0.0,0.0]) :
        if(isinstance(header,Header)) :
            self.header = header
        self.image = image
        self.mask = mask
        self.pointings = pointings

    def setHeader(self,header) :
        """ Method to set the header
            input :
                header - the image header
            returns :
                none
        """
        self.header = header

    def setImage(self,image,mask = None) :
        """ Method to set the image
            inputs :
                image - the image array
                mask - the mask binary array
            returns :
                none
        """
        self.image = image
        self.mask = mask

    def getHeader(self) :
        """ Method to return the image header
            input :
                none
            returns :
                the header
        """
        return self.header

    def getImage(self) :
        """ Method to return the image and mask
            input :
                none
            returns :
                the image and mask as a tuple
        """
        return self.image,self.mask

    def getHeaderValue(self,item) :
        """ Method to return a specific header value
            input :
                item -  the item to return
            returns :
                the value of the item
        """
        return self.header.getValue(item)

    def getHeaderDesc(self,item) :
        """ Method to get the description of a specific header item
            input :
                item - the name of the item to return
            returns :
                the decription of the item
        """
        return self.header.getDesc(item)

# class to read in a miriad image to a python array
class ImageReader :
    def __init__(self,file,cutoff) :
        """ Initializer
            file - the name of the file to read in
            cutoff - only read in values above this one
        """
        handle = mfunc.xyopen(file,"old")[0]
        self.header = mfunc.getHeader(handle)
        mfunc.xyclose(handle)
        self.handle = open(file + "/image")
        img = self.handle.read(4)
        self.size = struct.unpack(">l",img)[0]
        self.x = self.header.getValue("naxis1")
        self.y = self.header.getValue("naxis2")
        self.cutoff = cutoff
        self.type = ">"
        if(self.size == 4) :
            for i in range(0,self.x) :
                self.type += "f"

        elif(size == 8) :
            for i in range(0,self.x) :
                self.type += "d"

    def setHeader(self,header) :
        """ Method to set the header
            input :
                header - the iamge header
            returns :
                none
        """
        self.header = header

    def getHeader(self) :
        """ Method to return the image header
            input :
                none
            returns :
                the image header
        """
        return self.header

    def getHeaderValue(self,item) :
        """ Method to return a specific header value
            input :
                item - the name of the item to query
            returns :
                the value of item
        """
        return self.header.getValue(item)

    def getHeaderDesc(self,item) :
        """ Method to return the description of a specific header item
            input :
                item -  the item to query
            returns :
                the description of the item
        """
        return self.header.getDesc(item)

    def getRow(self,plane,row) :
        """ Method to read in a miriad image row
            input :
                plane - the plane to read
                row - the row to read
            returns :
                a list containing the row data
        """
        self.handle.seek(plane*self.size*self.x*self.y+self.size + self.size*row*self.x)
        rawData = self.handle.read(self.size*self.x)
        data = list(struct.unpack(self.type,rawData))
        for d in range(0,len(data)) :
            if(data[d] < self.cutoff) :
                data[d] = 0.0
        return data

def imhead(file,key) :
    """ Method to read a specific header value from an image
        input :
            file - the name of the miriad image to read
            key - the name of the header item
        returns :
            the value of the item
    """
    handle,num = mfunc.xyopen(file,"old")
    desc,typeX,nr = mfunc.hdprobe(handle,key)
    mfunc.xyclose(handle)
    return desc
