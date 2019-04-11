""" Example python script for generating a gearfile with arbitrary pixel sizes
    author: Helge C. Beck helge-christoph.beck@phys.uni-goettingen.de

    The script uses the xml.ElementTree library to generate a gear xml file. Purpose is to describe the telescope geometry and the detectors to then be used in the tbsw framework.
    The back end for creating the detector planes is in the module DetLayoutGen.py. 

    Usage: 
    In the __main__ part a example is given for the creation of a gear file for a telescope with 6 mimosa26 detector planes and one Atlas FEI4 reference plane created with the detector type SiPlanesParameters. The dut plane is created with the gear type PolyPlanesParameters. 
    The base gear structure is created with the function DetLayoutGen.createGearStructure().
    Because of memory consumption all parts are written (appended) to the outfile as soon as possible. The example file for the mimosa26 telescope has a size of 100Mb on disc. Trying to keep this in memory is not adviced. Opening and closing the file for every pixel slows down the programm considerably, so a large enough buffer can be filled and then written. Care has to be taken that the telescope layers are put in between the correct opening and closing tags as done in the provided example. Because of the incremental writing a propper indentation is not done. There are some indentation and linebreaks so it is not a one line file and somewhat readable.
    
    The detector planes are handled with the class Detector from DetLayoutGen.py. 
      The Parameters of the contructor are list of the parameters used for each layer. It is assumed that no list given or to few list entries that the default parameters should be used. Missing list entries are appended at the end of the list. All lists assume that the index for a plane is the same in every parameter list. 
      The class can be used to generate detectors of two gear type representations: SiPlanesParameters and PolyPlanesParameters. The first one assumes a rectangular pixel grid whereas the second one is totaly without limitations of pixel sizes and overall geometry. The parameters to be set for each are: 
        SiPlanesParameters: uCellGroupparamsList, vCellGroupparamsList ----- List of cellGroup parameters for each layer/plane
                            detectorparams ----------- most importent the geartype and number (number of planes)
        PolyPlanesParameters: pixelPrototypeparamsList --------- List of list of prototype pixel for each layer
                              createPixelfunctionList --------- List of functions to create the pixels in each layer, has to be provided for each layer otherwise a ValueError is thrown
                              detectorparams
        The other parameter lists are common to both types and some of the parameters need to be set for a sensible description of the detectors (IDs, positions)
        Important: All parameters in the dictionaries have to be strings otherwise the Elementtree library does not parse them. 
      Parameters to be set can be found in the declared default dictionaries in the class for the different parts of the plane. The list of the respective parameters is set on construction of the detector object. Care has to be taken that the parameters for the pixel prototypes (or u/vCellGroups) are handled in a list of dictionaries ([{},{},]) because there could be multiple pixel types in one plane. 
      For PolyPlanesparameters detector types:
        Most important parameter to set is the list of function for each layer that generates the layout of the pixels the pixel matrix for each layer. It has to be set during construction. This function has to be implemented for the specific layout. An example of this function can be found in the implementation for the mimosa26 pixel -> createPixelM26(). Also a function with the basic parts is provided: createPixelBase(). Characteristics of this function:
          return: No return, all created elements are written to file directly
          Elements: name is 'pixel'
                    attributes = {"type": "0",  # global type of the pixel
                                "u": , "v",   # column, row coordinates of the pixel, pixel address in raw data, or pseudo coordinates but then a conversion processor needs to be used when reading in the raw data
                                "centreu": ,  # centre coordinate u of pixel in mm, collection point of the charges, assigned hit position in calibration and reconstruction 
                                "centrev": ,  # centre coordinate v of pixel in mm
                               }
          For conveniance there is a method called in Detector class before generating the pixel layout that creates a list of list of tuples. The list of tuples consists of the coordinates of the pixel prototype corners. Every prototype pixel is an entry in the outer list. This list is layer specific and can be accessed with self.basePointListList. Because there can be multiple protopixel in one layer the first entry of each list for the pixel is its type and not a corner point!
          The function addPixelToTH2Poly(attributes) uses the self.basePointListList to add the pixel with the attributes to a layer specific TH2Poly object. All these objects are in the last step written to a root file. 
      
      Other parameters for Detector class:
        The filename/path where the resulting gear file is written has to be provided to xmloutfile contructer parameter. In the class the write action is always append, so the file should be prepared for the detector elements. This name is also used to write the layout to a root file (.xml -> .root) if it is a PolyPlanesParameters type detector. 
      With the default parameters 6 mimosa26 planes are created with ID=0 and at positionX/Y/Z=0.0. These two parameters are the ones that have to be set for every plane. For every plane the default parameters are assigned first and then with the provided parameters overwritten. So not provided parameters stay the default ones. 
      The Detector class does not return anything. All elements are written to the outputfile with the writeFile function.
    
    Parameters for the telescope setup have to be provided in some way. 

"""
# update text TODO
import xml.etree.ElementTree as ET
#from copy import deepcopy 
import os
import sys

import tbsw

"""
# Example (incomplete) of the basics that should be in a createPixel function for creating a PolyPlaneParameters detector layout with the Detector class.

def createPixelBase(self):
  pixelList = []
  attribute = {}
  attribute["type"] = "0" # or specific for each pixel 
  layoutStartX = 61 
  layoutEndX = 78
  layoutStartY = 0
  layoutEndY = 92
  pitchx = 0.25
  pitchy = 0.05
  count = 0
  # Go through the range of pixel indices, or similar
  for x in range(layoutStartX, layoutEndX+1): # +1 for inclusive
    for y in range(layoutStartY, layoutEndY+1):
      attribute["u"] = str(x)
      attribute["v"] = str(y)
      #calculate centre of pixel in layout, in tbsw the centre of the layout will be shifted to 0,0
      centeru = 
      centerv = 
      attribute["centeru"] = str(centeru)
      attribute["centerv"] = str(centerv)
      pixelList.append(ET.Element('pixel', attribute))
      
      # test TH2Poly layout, Detector class writes a root file with the PolyPlaneParameters type layers as TH2Poly objects. Use this function to construct the bins/pixel in the TH2Poly
      self.addPixelToTH2Poly(attribute)

      # write to xml file multiple pixel so that io actions do not slow down the script 
      count += 1
      if count == 500 or (x == layoutEndX and y == layoutEndY):
        writePixelListToFile(self.xmloutfile, pixelList)
        pixelList[:] = []
        count = 0
"""

def writePixelListToFile(xmlfilename, pixelList):
  with open(xmlfilename, "a") as f:
    for element in pixelList:
      tbsw.DetLayoutGen.indent(element)
      f.write(ET.tostring(element)+"\n")
    f.flush()
    os.fsync(f.fileno())
    pixelList[:] = []


def createPixelM26(self):
  pixelList = []
  attributes = {}
  attributes["type"] = "0"
  npixelsU = 1151
  npixelsV = 575
  pitch = 0.0184
  shiftu = -0.5*npixelsU*pitch+0.5*pitch # 0.5 pitch because prototype has centre at 0|0 so need to correct for the half pitch to centre the whole layout at 0|0
  shiftv = -0.5*npixelsV*pitch+0.5*pitch
  count = 0
  for i in range(npixelsU):
    for j in range(npixelsV):
      attributes["u"] = str(i)
      attributes["v"] = str(j)
     
      centreu = pitch*i+shiftu 
      centrev = pitch*j+shiftv 
      attributes["centreu"] = str(centreu)
      attributes["centrev"] = str(centrev)
      pixelList.append(ET.Element('pixel', attributes))

      self.addPixelToTH2Poly(attributes)

      count += 1
      if count == 1000 or (i == npixelsU-1 and j == npixelsV-1):
        writePixelListToFile(self.xmloutfile, pixelList)
        pixelList[:] = []
        count = 0

# creates the pixel layout for a mixture of rectangular and hexagonal pixel used by H. C. Beck
def createPixel2Ddiamond(self):
   
  createFEI42DDiamond(self)
  createRect2DDiamond(self)
  createHex2DDiamond(self)


def createFEI42DDiamond(self):
  pixelList = []
  attribute = {}
  attribute["type"] = "0"
  layoutStartX = 61 
  layoutEndX = 78
  layoutStartY = 0
  layoutEndY = 92
  metalStartX = 64
  metalEndX = 76
  metalStartY = 12
  metalEndY = 76
  pitchx = 0.25
  pitchy = 0.05
  shiftx = 2.25 #global shift of FEI4 area into the layout: half x size
  shifty = 2.325 # half y size
  count = 0
  for x in range(layoutStartX, layoutEndX+1): # +1 for inclusive
    for y in range(layoutStartY, layoutEndY+1):
      if x >= metalStartX and x < metalEndX and y >= metalStartY and y < metalEndY:
        continue;
      attribute["u"] = str(x)
      attribute["v"] = str(y)
      centeru = 0.5*pitchx + (x-layoutStartX)*pitchx - shiftx
      centerv = 0.5*pitchy + (y-layoutStartY)*pitchy - shifty
      attribute["centeru"] = str(centeru)
      attribute["centerv"] = str(centerv)
      pixelList.append(ET.Element('pixel', attribute))
      
      # test TH2Poly layout
      self.addPixelToTH2Poly(attribute)

      count += 1
      if count == 500 or (x == layoutEndX and y == layoutEndY):
        writePixelListToFile(self.xmloutfile, pixelList)
        pixelList[:] = []
        count = 0

def createRect2DDiamond(self):
  pixelList = []
  attribute = {}
  
  count = 0
  metalStartX = 64
  metalEndX = 76
  metalStartY = 12
  hexStartY = 43
  rectperiod = 4 # 4*0.125 is one 0.5 block
  pitchx = 0.125
  pitchy = 0.1
  shiftx = 1.500 # shift rect area middle to 0,0, is also x middle in total layout
  shifty = 0.775 + 0.950 # shift rect area middle to 0,0, then shift to total layout middle
  for x in range(metalStartX, metalEndX):
    for y in range(metalStartY, hexStartY):
      attribute["type"] = "1"
      attribute["u"] = str(x)
      attribute["v"] = str(y)
      centeru = 0.0
      centerv = 0.0
      # v coordinate first because small rect will shift this a bit but everything else is ok
      if (y-metalStartY)%2 == 0: # even row, first row later
        centerv = (y-metalStartY)/2*pitchy - shifty
      else: #odd row
        centerv = (y+1-metalStartY)/2*pitchy - shifty
      
      # first row is special
      if y == metalStartY: # first row only metal dot on solder bump no trace -> smallest rectangle
        attribute["type"] = "8"
        if (x-metalStartX)%2 == 0: # right column
          centeru = 0.25*pitchx + rectperiod*pitchx*((x-metalStartX)/2) - shiftx
        else: # left column
          centeru = -0.25*pitchx + rectperiod*pitchx*((x+1-metalStartX)/2) - shiftx
        centerv = 0.25*pitchy - shifty
      # not first row, normal pixel
      elif (x-metalStartX)%2 == 0: #right column
        if (y-metalStartY)%2 == 1: # right column, odd row
          centeru = pitchx + rectperiod*pitchx*((x-metalStartX)/2) - shiftx
        elif ((x-metalStartX)/2)%2 == 0: #even right dcolumn, even row
          centeru = 2.0*pitchx + rectperiod*pitchx*((x-metalStartX)/2) - shiftx
        else: # odd right dcolumn, even row
          attribute["type"] = "2" # small rect
          centeru = 0.25*pitchx + rectperiod*pitchx*((x-metalStartX)/2) - shiftx
          #centerv += 0.25*pitchy # no shift because center is not in upper or lower half of rectangle
      else : #left column
        if ((x+1-metalStartX)/2)%2 == 0: #even left dcolumn
          if (y-metalStartY)%2 == 1: # even left dcolumn, odd row
            centeru = -pitchx + rectperiod*pitchx*((x+1-metalStartX)/2) - shiftx
          else: # second row
            centeru = -2.0*pitchx + rectperiod*pitchx*((x+1-metalStartX)/2) - shiftx
        else: # odd left dcolumn
          if (y-metalStartY)%2 == 1: # odd row, small rect
            attribute["type"] = "2"
            centeru = -0.25*pitchx + rectperiod*pitchx*((x+1-metalStartX)/2) - shiftx
            #centerv -= 0.25*pitchy
          else: # even row
            centeru = -pitchx + rectperiod*pitchx*((x+1-metalStartX)/2) - shiftx
             
      attribute["centeru"] = str(centeru)
      attribute["centerv"] = str(centerv)
      pixelList.append(ET.Element('pixel', attribute))

      self.addPixelToTH2Poly(attribute)

      count += 1
      if count == 200 or (x == metalEndX-1 and y == hexStartY-1):
        writePixelListToFile(self.xmloutfile, pixelList)
        pixelList[:] = []
        count = 0

def createHex2DDiamond(self):
  pixelList = []
  attribute = {}
  
  count = 0
  metalStartX = 64
  metalEndX = 76
  hexStartY = 43
  metalEndY = 76
  hexperiod = 4
  periodlengthx = 0.5
  pitchx = 0.1138
  pitchy = 0.1333
  hexRectX = 0.0793
  lastRowRectY = 0.09166
  shiftx = 1.500 # same local and total layout middle
  shifty = 0.825 - 0.675 # already gap of 0.25 between rect and hex calculated
  for x in range(metalStartX, metalEndX):
    for y in range(hexStartY, metalEndY):
      attribute["type"] = "3" # hex, 4 recthex, 5 small rect
      attribute["u"] = str(x)
      attribute["v"] = str(y)
      centeru = 0.0
      centerv = 0.0

      if (y-hexStartY)%2 == 0: # first row
        centerv = 0.5*pitchy + 0.75*(y-hexStartY)/2*pitchy - shifty
      else : # second row
        centerv = 0.5*pitchy + 0.75*(y-1-hexStartY)/2*pitchy - shifty

      if (x-metalStartX)%2 == 0: # right column
        if (y-hexStartY)%4 == 0: #first row, hex
          centeru = hexRectX + pitchx + periodlengthx*(x-metalStartX)/2 - shiftx
        elif (y-hexStartY)%4 == 1: # second row, recthex
          attribute["type"] = "4"
          centeru = hexRectX + periodlengthx*(x-metalStartX)/2 - shiftx
        elif (y-hexStartY)%4 == 2: # third row, hex
          centeru = hexRectX + 0.5*pitchx + periodlengthx*(x-metalStartX)/2 - shiftx
        elif (y-hexStartY)%4 == 3: # forth row
          if (x-metalStartX)/2%2 == 0: # even right dcolumn, small rect
            attribute["type"] = "5"
            centeru = 0.5*hexRectX + periodlengthx*(x-metalStartX)/2 - shiftx
          else: # odd right dcolumn, hex
            centeru = hexRectX + 1.5*pitchx + periodlengthx*(x-metalStartX)/2 - shiftx
      else: # left column
        if (y-hexStartY)%4 == 0: #first row, hex
          centeru = -hexRectX - pitchx + periodlengthx*(x+1-metalStartX)/2 - shiftx
        elif (y-hexStartY)%4 == 1: # second row, recthex
          attribute["type"] = "9"
          centeru = -hexRectX + periodlengthx*(x+1-metalStartX)/2 - shiftx
        elif (y-hexStartY)%4 == 2: # third row, hex
          centeru = -hexRectX - 0.5*pitchx + periodlengthx*(x+1-metalStartX)/2 - shiftx
        elif (y-hexStartY)%4 == 3: # forth row
          if (x+1-metalStartX)/2%2 == 0: # even left dcolumn, small rect
            attribute["type"] = "5"
            centeru = -0.5*hexRectX + periodlengthx*(x+1-metalStartX)/2 - shiftx
          else: # odd left dcolum, hex
            centeru = -hexRectX - 1.5*pitchx + periodlengthx*(x+1-metalStartX)/2 - shiftx
      if y == metalEndY-3: # second pixel in last recthex no trace but assign to cut off hex in last row, center ok from general sorting
        attribute["type"] = "6"
      elif y == metalEndY-2: # first row in smaller rectangles in last rows
        attribute["type"] = "7"
        if (x-metalStartX)%2 == 0: #right column
          centeru = 0.5*hexRectX + periodlengthx*(x-metalStartX)/2 - shiftx
        else: # left column
          centeru = -0.5*hexRectX + periodlengthx*(x+1-metalStartX)/2 - shiftx
        centerv = centerv - 0.25*pitchy + 0.25*lastRowRectY
      elif y == metalEndY-1: # last row no metal
        attribute["type"] = "7"
        if (x-metalStartX)%2 == 0: #right column
          centeru = 0.5*hexRectX + periodlengthx*(x-metalStartX)/2 - shiftx
        else: # left column
          centeru = -0.5*hexRectX + periodlengthx*(x+1-metalStartX)/2 - shiftx
        centerv = centerv -(0.75+0.25)*pitchy + 0.75*lastRowRectY # normaly is assigned to next row so substract that and correct for small rect center

      attribute["centeru"] = str(centeru)
      attribute["centerv"] = str(centerv)
      pixelList.append(ET.Element('pixel', attribute))

      self.addPixelToTH2Poly(attribute) 

      count += 1
      if count == 200 or (x == metalEndX-1 and y == metalEndY-1):
        writePixelListToFile(self.xmloutfile, pixelList)
        pixelList[:] = []
        count = 0


if __name__ == '__main__':
  
  # Stuff that has to be individually provided for the specific geometry
  outfile = 'gear_dia_planar_batch2_alltype.xml'

  #batch2 2d dia
  telPositionZ = [0, 150, 300, 653, 743, 898]
  dutPositionZ = [403, 494] # dia, ref
  dutID = [21, 22] #dia, ref

  # pixel prototypes
  protoFEI4 = dict(type="0", distu="0.3", distv="0.01", points="-0.125 -0.025, -0.125 0.025, 0.125 0.025, 0.125 -0.025")

  protoRect = dict(type="1", distu="0.14", distv="0.12", points="-0.0625 -0.05, -0.0625 0.05, 0.0625 0.05, 0.0625 -0.05")

  protoRectHalf = dict(type="2", distu="0.1", distv="0.12", points="-0.03125 -0.05, -0.03125 0.05, 0.03125 0.05, 0.03125 -0.05")

  protoRectHalfHalf = dict(type="8", distu="0.1", distv="0.12", points="-0.03125 -0.025, -0.03125 0.025, 0.03125 0.025, 0.03125 -0.025")

  protoHex = dict(type="3", distu="0.14", distv="0.12", points="0.0 -0.0666, -0.0569 -0.0333, -0.0569 0.0333, 0.0 0.0666, 0.0569 0.0333, 0.0569 -0.0333")

  protoRectHex = dict(type="4", distu="0.162", distv="0.12", points="0.0 -0.0666, -0.0793 -0.0666, -0.0793 0.0666, 0.0 0.0666, 0.0569 0.0333, 0.0569 -0.0333")

  protoRectHexFlipped = dict(type="9", distu="0.162", distv="0.12", points="0.0 -0.0666, -0.0569 -0.0333, -0.0569 0.0333, 0.0 0.0666, 0.0793 0.0666, 0.0793 -0.0666")

  protoHexRectSmall = dict(type="5", distu="0.1", distv="0.12", points="-0.03965 -0.0333, -0.03965 0.0333, 0.03965 0.0333, 0.03965 -0.0333")

  protoHexCut = dict(type="6", distu="0.14", distv="0.12", points="-0.0569 0.05833, -0.0569 -0.0333, 0.0 -0.0666, 0.0569 -0.0333, 0.0569 0.05833")

  protoHexRectLastRows = dict(type="7", distu="0.1", distv="0.1", points="-0.03965 -0.022915, -0.03965 0.022915, 0.03965 0.022915, 0.03965 -0.022915")

  # end of input parameters


  # base mimosa 26 telescope with Atlas FEI4 reference plane in standard gear/xml format
  f = open(outfile, "w") # reset the file to empty because use append mode later
  f.close()
  # create the basic xml file structure, for paramters see DetLayoutGen.py
  gearBaseTags = tbsw.DetLayoutGen.createGearStructure()
  # write the opening tags of the basic xml structure. This is needed because the xml nodes are itterativly created and written to file to reduce memory consumption
  with open(outfile, "a") as f:
    f.write(gearBaseTags[0])
    f.flush()
    os.fsync(f.fileno())

  spList = []
  uCGpList = []
  vCGpList = []
  # Detector in SiPlaneParameters description
  # Mimosa26 planes, use default parameters
  for i in range(len(telPositionZ)):
    spList.append({"ID":str(i), "positionZ":str(telPositionZ[i])})
    uCGpList.append([{}])
    vCGpList.append([{}])

  # Atlas FEi4 plane
  spList.append({"ID":str(dutID[1]), "positionZ":str(dutPositionZ[1]), "thickness":"1.0"})
  uCGpList.append([{"minCell":"0", "maxCell":"80", "pitch":"0.25"}])
  vCGpList.append([{"minCell":"0", "maxCell":"336", "pitch":"0.05"}])

  dparams = {"number":"7"}
  # create and initialise the Detector object
  detector = tbsw.DetLayoutGen.Detector(xmloutfile=outfile, sensparamsList=spList, detectorparams=dparams, uCellGroupparamsList=uCGpList, vCellGroupparamsList=vCGpList)
  # create the detector with the planes and write it to outfile
  detector.createDetectorElement()

  # 2D dia plane in full PolyGear description
  cpfList = [createPixel2Ddiamond]
  dparams = {"geartype":"PolyPlanesParameters", "number":"1"}
  lpList = [{"radLength":"121.3068", "atomicNumber":"6", "atomicMass":"12"}]
  spList = [{"ID":str(dutID[0]), "positionZ":str(dutPositionZ[0]), "thickness":"0.824", "radLength":"121.3068", "atomicNumber":"6", "atomicMass":"12"}]
  pppList = []
  ppp = []
  ppp.append(protoFEI4)
  ppp.append(protoRect)
  ppp.append(protoRectHalf)
  ppp.append(protoRectHalfHalf)
  ppp.append(protoHex)
  ppp.append(protoRectHex)
  ppp.append(protoRectHexFlipped)
  ppp.append(protoHexRectSmall)
  ppp.append(protoHexCut)
  ppp.append(protoHexRectLastRows)
  pppList.append(ppp)
  # There are two errors raised on __init__ of detector both if something with the createPixelfunctionList parameter is not correct, no list or to few list items. It is not needed to use this try-except because the exception in Detector class cancels the program then, but it is neater.
  try: 
    detector = tbsw.DetLayoutGen.Detector(xmloutfile=outfile, ladderparamsList=lpList, sensparamsList=spList, createPixelfunctionList=cpfList, detectorparams=dparams, pixelPrototypeparamsList=pppList)
  except ValueError as ex:
    sys.exit(str(ex)+ " Aborting!")
  detector.createDetectorElement()

  # write the closing tags of the basic xml structure
  with open(outfile, "a") as f:
    f.write(gearBaseTags[1])
    f.flush()
    os.fsync(f.fileno())
