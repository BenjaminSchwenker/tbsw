""" Example python script for generating a gearfile with arbitrary pixel sizes
    author: Helge C. Beck helge-christoph.beck@phys.uni-goettingen.de

    The script uses the xml.ElementTree library to generate a gear xml file. Purpose is to describe the telescope geometry and the detectors to then be used in the tbsw framework.

    Usage: 
    In the __main__ part a example is given for the creation of a gear file for a telescope with 6 mimosa26 detector planes and one Atlas FEI4 reference plane. The last telescope plane is created with the gear type PolyPlanesParameters. 
    The base gear structure is (at the moment hard coded) in the function createGearStructure(). Parameter to change there is the B-field.
    All parameter values have to be strings otherwise the parsing does not work.
    Because of memory consumption all parts are written (appended) to the outfile as soon as possible. The example file for the mimosa26 telescope has a size of 100Mb on disc. Trying to keep this in memory is not adviced. Opening and closing the file for every pixel slows down the programm considerably, so a large enough buffer can be filled and then written. Care has to be taken that the telescope layers are put in between the correct opening and closing tags as done in the provided example. This is not a very intuitive way of doing it put for the mimosa26 this reduces the memory consumption from 50% (capped) to less then 1% on 8Gb RAM machine.
    The detector planes are handled with the class Detector. 
      The Parameters of the contructor are list of the parameters used for each layer. It is assumed that no list given or to few list entries that the default parameters should be used. Missing list entries are appended at the end of the list. All lists assume that the index for a plane is the same in every parameter list. 
      The class can be used to generate detectors of two gear type representations: SiPlanesParameters and PolyPlanesParameters. The first one assumes a rectangular pixel grid whereas the second one is totaly without limitations of pixel sizes and overall geometry. The parameters to be set for each are: 
        SiPlanesParameters: uCellGroupparamsList, vCellGroupparamsList ----- List of cellGroup parameters for each layer/plane
                            detectorparams ----------- most importent the geartype and number (number of planes)
        PolyPlanesParameters: pixelPrototypeparamsList --------- List of list of prototype pixel for each layer
                              detectorparams
        The other parameter lists are common to both types and some of the parameters need to be set for a sensible description of the detectors (IDs, positions)
      Parameters to be set can be found in the declared default dictionaries in the class for the different parts of the plane. The list of the respective parameters can be set on construction of the planes object or with the setXParams(dict) functions. Care has to be taken that the parameters for the pixel prototypes are handled in a list of dictionaries ([{},{},]) because there could be multiple pixel types in one plane. 
      Pixel prototypes centre is at 0.0. The pixel with this type will shift their centre to a new position shifting the full pixel with it updating the pixel corners. 
      Most important parameter to set is the function that generates the layout of the pixels the pixel matrix. It is the first parameter of the constructer (after self) and has to be set during construction (at the moment no set function provided). This function has to be implemented for the specific layout. An example of this function can be found in the implementation for the mimosa26 pixel -> createPixelM26(). Characteristics of this function:
        return: No return, all created elements are written to file directly
        Elements: name is 'pixel'
                  attributes = {"type": "0",  # global type of the pixel
                                "u": , "v",   # column, row coordinates of the pixel, pixel address in raw data, or pseudo coordinates but then a conversion processor needs to be used when reading in the raw data
                                "centreu": ,  # centre coordinate u of pixel in mm, collection point of the charges, assigned hit position in calibration and reconstruction 
                                "centrev": ,  # centre coordinate v of pixel in mm
                               }
        For conveniance there is a method called in Detector class before generating the pixel layout that creates a list of list of tuples. The list of tuples consists of the coordinates of the pixel prototype corners. Every prototype pixel is an entry in the outer list. This list is (should be) plane object specific and can be accessed with self.basePointListList. 
      Other parameters for Detector class:
        The filename/path where the resulting gear file is written has to be provided to xmloutfile contructer parameter. In the class the write action is always append, so the file should be prepared for the plane elements.
        With the default parameters (including now the createPixelM26 function that has to be set manualy at construction) a mimosa26 plane is created with ID=0 and at positionX/Y/Z=0.0. These two parameters are the ones that have to be set for every plane. Either on construction or with the setSensitiveParams() function. Not every parameter has to be set. For every plane the default parameters are assigned first and then with the provided parameters overwritten. So not provided parameters stay the default ones. 
        The Detector class does not return anything. All elements are written to the outputfile with the writeFile function.
    The path/name of the output gear xml file has to be provided to the Plane class constructor. 
    Parameters for the telescope setup have to be provided in some way. 

"""
import xml.etree.ElementTree as ET
from copy import deepcopy 
import os
import array

from ROOT import TFile, TGraph, TH2Poly

# pixel prototypes
protoFEI4 = dict(type="0", distu="0.3", distv="0.01", points="-0.125 -0.025, -0.125 0.025, 0.125 0.025, 0.125 -0.025")

protoRect = dict(type="1", distu="0.14", distv="0.12", points="-0.0625 -0.05, -0.0625 0.05, 0.0625 0.05, 0.0625 -0.05")

protoRectHalf = dict(type="2", distu="0.1", distv="0.12", points="-0.03125 -0.05, -0.03125 0.05, 0.03125 0.05, 0.03125 -0.05")

protoRectHalfHalf = dict(type="8", distu="0.1", distv="0.12", points="-0.03125 -0.025, -0.03125 0.025, 0.03125 0.025, 0.03125 -0.025")

protoHex = dict(type="3", distu="0.14", distv="0.12", points="0.0 -0.0666, -0.0569 -0.0333, -0.0569 0.0333, 0.0 0.0666, 0.0569 0.0333, 0.0569 -0.0333")

protoRectHex = dict(type="4", distu="0.162", distv="0.12", points="0.0 -0.0666, -0.0793 -0.0666, -0.0793 0.0666, 0.0 0.0666, 0.0569 0.0333, 0.0569 -0.0333")

protoRectHexFlipped = dict(type="9", distu="0.162", distv="0.12", points="0.0 -0.0666, -0.0569 -0.0333, -0.0569 0.0333, 0.0 0.0666, 0.0793 0.0666, 0.0793 -0.0666")

protoHexRectSmall = dict(type="5", distu="0.1", distv="0.12", points="-0.03965 -0.0333, -0.03965 0.0333, 0.03965 0.0333, 0.03965 -0.0333")

protoHexCut = dict(type="6", distu="0.14", distv="0.12", points="-0.0569 0.05833, -0.0569 -0.0333, 0.0 0.0666, 0.0569 -0.0333, 0.0569 0.05833")

protoHexRectLastRows = dict(type="7", distu="0.1", distv="0.1", points="-0.03965 -0.022915, -0.03965 0.022915, 0.03965 0.022915, 0.03965 -0.022915")

protoFEI4points = [[-0.125, -0.025],[-0.125, 0.025], [0.125, 0.025], [0.125, -0.025]]
protoRectpoints = [[-0.0625, -0.05], [-0.0625, 0.05], [0.0625, 0.05], [0.0625, -0.05]]
protoRectHalfpoints = [[-0.03125, -0.05], [-0.03125, 0.05], [0.03125, 0.05], [0.03125, -0.05]]
protoRectHalfHalfpoints = [[-0.03125, -0.025], [-0.03125, 0.025], [0.03125, 0.025], [0.03125, -0.025]]
protoHexpoints = [[0.0, -0.0666], [-0.0569, -0.0333], [-0.0569, 0.0333], [0.0, 0.0666], [0.0569, 0.0333], [0.0569, -0.0333]]
protoRectHexpoints = [[0.0, -0.0666], [-0.0793, -0.0666], [-0.0793, 0.0666], [0.0, 0.0666], [0.0569, 0.0333], [0.0569, -0.0333]]
protoRectHexFlippedpoints = [[0.0, -0.0666], [-0.0569, -0.0333], [-0.0569, 0.0333], [0.0, 0.0666], [0.0793, 0.0666], [0.0793, -0.0666]]
protoHexRectSmallpoints = [[-0.03965, -0.0333], [-0.03965, 0.0333], [0.03965, 0.0333], [0.03965, -0.0333]]
protoHexCutpoints = [[-0.0569, 0.05833], [-0.0569, -0.0333], [0.0, -0.0666], [0.0569, -0.0333], [0.0569, 0.05833]]
protoHexRectLastRowpoints = [[-0.03965, -0.022915], [-0.03965, 0.022915], [0.03965, 0.022915], [0.03965, -0.022915]]

def createPixelM26(self):
  pixelList = []
  attributes = {}
  attributes["type"] = "0"
  npixelsU = 1151
  npixelsV = 575
  basepoints = self.basePointListList[0]
  npoints = len(basepoints)
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
      count += 1
      if count == 1000 or (i == npixelsU-1 and j == npixelsV-1):
        with open(self.xmloutfile, "a") as f:
          for element in pixelList:
            f.write(ET.tostring(element))
          f.flush()
          os.fsync(f.fileno())
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
  metalStartY = 12#13
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
      gpixel = TGraph(4)
      gpixel.SetName("{},{}".format(x,y))
      i = 0
      for point in protoFEI4points:
        gpixel.SetPoint(i, point[0]+centeru, point[1]+centerv)
        i +=1
      layout.AddBin(gpixel)

      count += 1
      if count == 500 or (x == layoutEndX and y == layoutEndY):
        with open(self.xmloutfile, "a") as f:
          for element in pixelList:
            f.write(ET.tostring(element))
          f.flush()
          os.fsync(f.fileno())
        pixelList[:] = []
        count = 0

def createRect2DDiamond(self):
  pixelList = []
  attribute = {}
  
  count = 0
  metalStartX = 64
  metalEndX = 76
  metalStartY = 12#13
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

      gpixel = TGraph(4)
      gpixel.SetName("{},{}".format(x, y))
      i = 0
      if attribute["type"] == "1":
        for point in protoRectpoints:
          gpixel.SetPoint(i, point[0]+centeru, point[1]+centerv)
          i += 1
      elif attribute["type"] == "2":
        for point in protoRectHalfpoints:
          gpixel.SetPoint(i, point[0]+centeru, point[1]+centerv)
          i += 1
      elif attribute["type"] == "8":
        for point in protoRectHalfHalfpoints:
          gpixel.SetPoint(i, point[0]+centeru, point[1]+centerv)
          i += 1
      layout.AddBin(gpixel)

      count += 1
      if count == 200 or (x == metalEndX-1 and y == hexStartY-1):
        with open(self.xmloutfile, "a") as f:
          for element in pixelList:
            f.write(ET.tostring(element))
          f.flush()
          os.fsync(f.fileno())
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

      
      i = 0
      if attribute["type"] == "3":
        gpixel = TGraph(6)
        gpixel.SetName("{},{}".format(x, y))
        for point in protoHexpoints:
          gpixel.SetPoint(i, point[0]+centeru, point[1]+centerv)
          i += 1
        layout.AddBin(gpixel)
      elif attribute["type"] == "4":
        gpixel = TGraph(6)
        gpixel.SetName("{},{}".format(x, y))
        for point in protoRectHexpoints:
          gpixel.SetPoint(i, point[0]+centeru, point[1]+centerv)
          i += 1
        layout.AddBin(gpixel)
      elif attribute["type"] == "9":
        gpixel = TGraph(6)
        gpixel.SetName("{},{}".format(x, y))
        for point in protoRectHexFlippedpoints:
          gpixel.SetPoint(i, point[0]+centeru, point[1]+centerv)
          i += 1
        layout.AddBin(gpixel)
      elif attribute["type"] == "5":
        gpixel = TGraph(4)
        gpixel.SetName("{},{}".format(x, y))
        for point in protoHexRectSmallpoints:
          gpixel.SetPoint(i, point[0]+centeru, point[1]+centerv)
          i += 1
        layout.AddBin(gpixel)
      elif attribute["type"] == "6":
        gpixel = TGraph(5)
        gpixel.SetName("{},{}".format(x, y))
        for point in protoHexCutpoints:
          gpixel.SetPoint(i, point[0]+centeru, point[1]+centerv)
          i += 1
        layout.AddBin(gpixel)
      elif attribute["type"] == "7":
        gpixel = TGraph(4)
        gpixel.SetName("{},{}".format(x, y))
        for point in protoHexRectLastRowpoints:
          gpixel.SetPoint(i, point[0]+centeru, point[1]+centerv)
          i += 1
        layout.AddBin(gpixel)
      

      count += 1
      if count == 200 or (x == metalEndX-1 and y == metalEndY-1):
        with open(self.xmloutfile, "a") as f:
          for element in pixelList:
            f.write(ET.tostring(element))
          f.flush()
          os.fsync(f.fileno())
        pixelList[:] = []
        count = 0


#def createPixel3DRectangle(self):


#def createPixel3DHexagonal(self):


class Detector:
  """Base class which holds all parameters and functions for creating a detector with multiple planes which have the same gear type in a gear file """

  #detector parameters
  ddefault = dict(name="SiPlanes", geartype="SiPlanesParameters", ID="250", number="8")

  # ladder paramaters
  ldefault = dict(sizeU="200", sizeV="100", thickness="0.07", 
                  radLength="93.660734", atomicNumber="14", atomicMass="28")

  #sensitive parameters
  sdefault = dict(ID="0", PixType="0", 
              positionX="0.00", positionY="0.00", positionZ="0.00", 
              thickness="0.07", 
              radLength="93.660734", atomicNumber="14", atomicMass="28", 
              alpha="0", beta="0", gamma="0", 
              rotation1="1.0", rotation2="0.0", rotation3="0.0", rotation4="1.0",)
  
  # uCellGroup default
  uCellGroupdefault = dict(minCell="0", maxCell="1151", pitch="0.0184")

  #vCellGroup default 
  vCellGroupdefault = dict(minCell="0", maxCell="575", pitch="0.0184")

  #pixel prototyp, within dist pixels are counted as neigbouring, default 0.02 a bit larger then 0.0184 pixel pitch
  protopixdefault = dict(type="0", distu="0.02", distv="0.02", points="-0.0092 -0.0092, -0.0092 0.0092, 0.0092 0.0092, 0.0092 -0.0092")

  #pixel matrix
  #createPixel

  def __init__(self, xmloutfile, createPixelfunctionList=None, ladderparamsList=None, sensparamsList=None, pixelPrototypeparamsList=None, detectorparams=None, uCellGroupparamsList=None, vCellGroupparamsList=None):
    self.xmloutfile = xmloutfile
    self.setDetectorParams(detectorparams)
    
    if ladderparamsList is not None:
      if len(ladderparamsList) < self.nlayers:
        for i in range(self.nlayers-len(ladderparamsList)):
          ladderparamsList.append(deepcopy(ldefault))
    else:
      ladderparamsList = [deepcopy(self.ldefault) for i in range(self.nlayers)]
    self.ladderparamslist = ladderparamsList

    if sensparamsList is not None:
      if len(sensparamsList) < self.nlayers:
        for i in range(self.nlayers-len(sensparamsList)):
          sensparamsList.append(deepcopy(sdefault))
    else:
      sensparamsList = [deepcopy(self.sdefault) for i in range(self.nlayers)]
    self.sensparamslist = sensparamsList

    if self.gearType == "PolyPlanesParameters":
      self.createPixelFunctionsList = createPixelfunctionList
      if pixelPrototypeparamsList is not None:
        if len(pixelPrototypeparamsList) < self.nlayers:
          for i in range(self.nlayers-len(pixelPrototypeparamsList)):
            pixelPrototypeparamsList.append(deepcopy(protopixdefault))
      else:
        pixelPrototypeparamsList = [[deepcopy(self.protopixdefault)] for i in range(self.nlayers)]
      self.pixelPrototypeParamsListList = pixelPrototypeparamsList

    elif self.gearType == "SiPlanesParameters":
      if uCellGroupparamsList is not None:
        if len(uCellGroupparamsList) < self.nlayers:
          for i in range(self.nlayers-len(uCellGroupparamsList)):
            uCellGrouparamsList.append(deepcopy(uCellGroupdefault))
      else:
        uCellGroupparamsList = [deepcopy(self.uCellGroupdefault) for i in range(self.nlayers)]
      self.uCellGroupParamsList = uCellGroupparamsList

      if vCellGroupparamsList is not None:
        if len(vCellGroupparamsList) < self.nlayers:
          for i in range(self.nlayers-len(vCellGroupparamsList)):
            vCellGrouparamsList.append(deepcopy(vCellGroupdefault))
      else:
        vCellGroupparamsList = [deepcopy(self.vCellGroupdefault) for i in range(self.nlayers)]
      self.vCellGroupParamsList = vCellGroupparamsList
      
  def createLadder(self):
    ladder = ET.Element('ladder', self.ladderparams)
    return ladder

  def createSensitive(self):
    attributes = deepcopy(self.sensparams)
    sensitive = ET.Element('sensitive', attributes)
    return sensitive

  def createPixelPrototype(self):
    for i in range(self.nprotopixel):
      self.writeFile(ET.Element('pixelProtoype', self.protopixlistdic[i]))

  def createPixelMatrix(self):
    matrixstart = '<pixelMatrix>'
    matrixend = '</pixelMatrix>'
    self.writeFile(matrixstart)
    self.createPixel(self)
    self.writeFile(matrixend)

  def setDetectorParams(self, detectorparamsdic):
    self.detectorparams = deepcopy(self.ddefault)
    if detectorparamsdic is not None:
      self.detectorparams.update(detectorparamsdic)
    self.gearType = self.detectorparams["geartype"]
    self.nlayers = int(self.detectorparams["number"])
  
  def setLadderParams(self, ladderparamsdic):
    self.ladderparams = deepcopy(self.ldefault)
    if ladderparamsdic is not None:
      self.ladderparams.update(ladderparamsdic)

  def setSensitiveParams(self, sensparamsdic):
    self.sensparams = deepcopy(self.sdefault)
    if sensparamsdic is not None:
      self.sensparams.update(sensparamsdic)
    
  def setPixelPrototypeParams(self, protopixdic):
    self.protopixlistdic = []
    if protopixdic is not None:
      self.nprotopixel = len(protopixdic)
      self.protopixlistdic = [deepcopy(self.protopixdefault) for i in range(self.nprotopixel)]
      for i in range(self.nprotopixel):
        self.protopixlistdic[i].update(protopixdic[i])
    else:
      self.nprotopixel = 1
      self.protopixlistdic.append(self.protopixdefault)
    self.setPrototypeBasePoints()

  def setCellGroupParams(self, uCellGroupdict, vCellGroupdict):
    self.uCellGroupparams = deepcopy(self.uCellGroupdefault)
    if uCellGroupdict is not None:
      self.uCellGroupparams.update(uCellGroupdict)
    self.vCellGroupparams = deepcopy(self.vCellGroupdefault)
    if vCellGroupdict is not None:
      self.vCellGroupparams.update(vCellGroupdict)

  def setPrototypeBasePoints(self):
    self.basePointListList = []
    for i in range(self.nprotopixel):
      pointList = []
      for point in self.protopixlistdic[i]["points"].split(','):
        pointList.append(tuple(x for x in point.split()))

      self.basePointListList.append(pointList)
  
  def createDetectorElement(self):
    detector = ET.Element('detector', name=self.detectorparams["name"], geartype=self.detectorparams["geartype"])
    siplanesID = ET.SubElement(detector, 'siplanesID', ID=self.detectorparams["ID"])
    siplanesNumber = ET.SubElement(detector, 'siplanesNumber', number=self.detectorparams["number"])
    stringdetector = ET.tostring(detector)
    layerInsertPos = stringdetector.find('</detector>')
    layerstart, layerend = stringdetector[:layerInsertPos], stringdetector[layerInsertPos:]
    layerstart += '<layers>'
    layerend = '</layers>' + layerend
    self.writeFile(layerstart)

    for i in range(self.nlayers):
      self.createLayerElement(i) 

    self.writeFile(layerend)

  def createLayerElement(self, layernumber):
    layer = ET.Element('layer')
    self.setLadderParams(self.ladderparamslist[layernumber])
    self.setSensitiveParams(self.sensparamslist[layernumber])
    ladder = self.createLadder()
    sensitive = self.createSensitive()
    layer.append(ladder)
    layer.append(sensitive)
    stringlayer = ET.tostring(layer)
    pixelInsertPos = stringlayer.find('</layer>')
    layerstart, layerend = stringlayer[:pixelInsertPos], stringlayer[pixelInsertPos:]
    self.writeFile(layerstart)

    if self.gearType == "SiPlanesParameters":
      self.setCellGroupParams(self.uCellGroupParamsList[layernumber], self.vCellGroupParamsList[layernumber])
      self.createPixel = self.createRectPixel
      self.createPixel()
    elif self.gearType == "PolyPlanesParameters":
      self.createPixel = self.createPixelFunctionsList[layernumber]
      self.setPixelPrototypeParams(self.pixelPrototypeParamsListList[layernumber])
      self.createPixelPrototype()
      self.createPixelMatrix()

    self.writeFile(layerend)
  
  def createRectPixel(self):
    ucells = ET.Element('uCellGroup', self.uCellGroupparams);
    vcells = ET.Element('vCellGroup', self.vCellGroupparams);
    self.writeFile(ucells);
    self.writeFile(vcells);

  def writeFile(self, writeobject):  # write object is either an ET.Element or a string
    writestring = ""
    if ET.iselement(writeobject):
      writestring = ET.tostring(writeobject)
    else:
      writestring = writeobject
    with open(self.xmloutfile, "a") as f:
      f.write(writestring)
      f.flush()
      os.fsync(f.fileno())


outfile = 'gear_dia_planar_batch2_alltype.xml'

#batch2 2d dia
telPositionZ = [0, 150, 300, 653, 743, 898]
dutPositionZ = [403, 494] # dia, ref
dutID = [21, 22] #dia, ref

rootoutfile = 'gear_dia_planar_batch2_alltype.root'
layout = TH2Poly()
layout.SetName("FullPlanarTest")

def createGearStructure():
  gearTag = ET.Element('gear')
  globalTag = ET.SubElement(gearTag, 'global', detectorName="EUTelescope")
  bfieldTag = ET.SubElement(gearTag, 'BField', type="ConstantBField", x="0", y="0", z="0")

  stringgearTag = ET.tostring(gearTag)
  layerInsertPos = stringgearTag.find('</gear>')
  openingTags, closingTags = stringgearTag[:layerInsertPos], stringgearTag[layerInsertPos:]
  openingTags += '<detectors>'
  closingTags = '</detectors>' + closingTags
  gearTag.clear()
  return [openingTags, closingTags]
  

if __name__ == '__main__':
  # base mimosa 26 telescope with Atlas FEI4 reference plane in standard gear format
  f = open(outfile, "w") # reset the file to empty because use append mode later
  f.close()
  gearBaseTags = createGearStructure()

  with open(outfile, "a") as f:
    f.write(gearBaseTags[0])
    f.flush()
    os.fsync(f.fileno())

  spList = []
  uCGpList = []
  vCGpList = []
  for i in range(len(telPositionZ)):
    spList.append({"ID":str(i), "positionZ":str(telPositionZ[i])})
    uCGpList.append({})
    vCGpList.append({})

  spList.append({"ID":str(dutID[1]), "positionZ":str(dutPositionZ[1]), "thickness":"1.0"})
  uCGpList.append({"minCell":"0", "maxCell":"80", "pitch":"0.25"})
  vCGpList.append({"minCell":"0", "maxCell":"336", "pitch":"0.05"})

  dparams = {"number":"7"}
  detector = Detector(xmloutfile=outfile, sensparamsList=spList, detectorparams=dparams, uCellGroupparamsList=uCGpList, vCellGroupparamsList=vCGpList)
  detector.createDetectorElement()

  # 2D dia plane in full PolyGear description
  cpfList = [createPixel2Ddiamond]
  dparams = {"geartype":"PolyPlanesParameters", "number":"1"}
  spList = [{"ID":str(dutID[0]), "positionZ":str(dutPositionZ[0]), "thickness":"0.824", "radLength":"121.3068", "atomicNumber":"6", "atomicMass":"12"}]
  pppList = []
  ppp = []
  ppp.append(protoFEI4)
  ppp.append(protoRect)
  ppp.append(protoRectHalf)
  ppp.append(protoHex)
  ppp.append(protoRectHex)
  ppp.append(protoHexRectSmall)
  ppp.append(protoHexCut)
  ppp.append(protoHexRectLastRows)
  pppList.append(ppp)
  detector = Detector(xmloutfile=outfile, createPixelfunctionList=cpfList, sensparamsList=spList, detectorparams=dparams, pixelPrototypeparamsList=pppList)
  detector.createDetectorElement()

  rootfile = TFile(rootoutfile, "RECREATE")

  for x in range(61,79):
    for y in range(93):
      layout.Fill("{},{}".format(x,y),1)
  layout.Draw()
  layout.Write()
  rootfile.Write()
  rootfile.Close()

  with open(outfile, "a") as f:
    f.write(gearBaseTags[1])
    f.flush()
    os.fsync(f.fileno())
