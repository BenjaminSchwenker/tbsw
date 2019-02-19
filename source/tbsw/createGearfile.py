""" Example python script for generating a gearfile with arbitrary pixel sizes
    author: Helge C. Beck helge-christoph.beck@phys.uni-goettingen.de

    The script uses the xml.ElementTree library to generate a gear xml file. Purpose is to describe the telescope geometry and the detectors to then be used in the tbsw framework.

    Usage: 
    In the __main__ part a example is given for the creation of a gear file for just a telescope with 6 mimosa26 detector planes. 
    The base gear structure is (at the moment hard coded) in the function createGearStructure(). Parameters to change there are the number of planes and potentialy the B-field.
    All parameter values have to be strings otherwise the parsing does not work.
    Because of memory consumption all parts are written (appended) to the outfile as soon as possible. The example file for the mimosa26 telescope has a size of 1.1Gb on disc. Trying to keep this in memory is not adviced. Opening and closing the file for every pixel slows down the programm considerably, so a large enough buffer can be filled and then written. Care has to be taken that the telescope layers are put in between the correct opening and closing tags as done in the provided example. This is not a very intuitive way of doing it put for the mimosa26 this reduces the memory consumption from 50% (capped) to less then 1% on 8Gb RAM machine.  
    The detector planes are handled with the class Plane. 
      Parameters to be set can be found in the declared default dictionaries in the class for the different parts of the plane. The parameters can be set on construction of the planes object or with the setXParams(dict) functions. Care has to be taken that the parameters for the pixel prototypes are handled in a list of dictionaries ([{},{},]) because there could be multiple pixel types in one plane. 
      Most important parameter to set is the function that generates the layout of the pixels the pixel matrix. It is the first parameter of the constructer (after self) and has to be set during construction (at the moment no set function provided). This function has to be implemented for the specific layout. An example of this function can be found in the implementation for the mimosa26 pixel -> createPixelM26(). Characteristics of this function:
        return: No return, all created elements are written to file directly
        Elements: name is 'pixel'
                  attributes = {"type": "0",  # global type of the pixel
                                "u": , "v",   # column, row coordinates of the pixel, pixel address in raw data
                                "points": ,   # string of corner points of the pixel in mm, actual position in layout, centre of sensitiv area at (0,0). Comma seperated x y values that are seperated by whitespace: "x0 y0, x1, y1, ..."
                                "adjacent": , # list of u, v addresses that are considered adjacent to the pixel, same style as points
                                "centreu": ,  # centre coordinate u of pixel in mm, collection point of the charges, assigned hit position in calibration and reconstruction 
                                "centrev": ,  # centre coordinate v of pixel in mm
                               }
        For conveniance there is a method called in Plane class before generating the pixel layout that creates a list of list of tuples. The list of tuples consists of the coordinates of the pixel prototype corners. Every prototype pixel is an entry in the outer list. This list is (should be) plane object specific and can be accessed with self.basePointListList. 
      Other parameters for Plane class:
        The filename/path where the resulting gear file is written has to be provided to xmloutfile contructer parameter. In the class the write action is always append, so the file should be prepared for the plane elements.
        With the default parameters (including now the createPixelM26 function that has to be set manualy at construction) a mimosa26 plane is created with ID=0 and at positionX/Y/Z=0.0. These two parameters are the ones that have to be set for every plane. Either on construction or with the setSensitiveParams() function. Not every parameter has to be set. For every plane the default parameters are assigned first and then with the provided parameters overwritten. So not provided parameters stay the default ones. 
        The Plane class does not return anything. All elements are written to the outputfile with the writeFile function.
    The path/name of the output gear xml file has to be provided to the Plane class constructor. 
    Parameters for the telescope setup have to be provided in some way. 

"""
import xml.etree.ElementTree as ET
from copy import deepcopy 
import os

def createPixelM26(self):
  pixelList = []
  attributes = {}
  attributes["type"] = "0"
  npixelsU = 1151
  npixelsV = 575
  basepoints = self.basePointListList[0]
  count = 0
  for i in range(npixelsU):
    for j in range(npixelsV):
      attributes["u"] = str(i)
      attributes["v"] = str(j)
     
      npoints = len(basepoints)
      pitch = 0.018402778
      shiftu = -0.5*npixelsU*pitch
      pointsu = [float(point[0])+pitch*i+shiftu for point in basepoints]
      shiftv = -0.5*npixelsV*pitch
      pointsv = [float(point[1])+pitch*j+shiftv for point in basepoints]
      pointstring = ','.join(["{} {}".format(float(point[0])+pitch*i+shiftu, float(point[1])+pitch*j+shiftv) for point in basepoints])
      centreu = sum(pointsu)/ npoints
      centrev = sum(pointsv)/ npoints
      adjacentstring = ""
      for adju in range(i-1, i+2):
        if adju < 0 or adju > npixelsU:
          continue
        for adjv in range(j-1, j+2):
          if adjv < 0 or adjv > npixelsV:
            continue
          if adju == i and adjv == j:
            continue
          adjacentstring += "{} {}, ".format(adju, adjv)

      attributes["points"] = pointstring
      attributes["adjacent"] = adjacentstring
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


#def createPixel3DRectangle(self):


#def createPixel3DHexagonal(self):


class Plane:
  """Base class which holds all parameters and functions for creating a plane in a gear file """

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
  
  #pixel prototyp
  protopixdefault = dict(type="0", npoints="4", points="0.0 0.0, 0.0 0.018402778, 0.018402778 0.018402778, 0.018402778 0.0")

  #pixel matrix
  #createPixel

  def __init__(self, createPixelfunction, xmloutfile, ladderparams=None, sensparams=None, pixelProtoypeparams=None):
    self.createPixel = createPixelfunction 
    self.xmloutfile = xmloutfile
    self.setLadderParams(ladderparams)
    self.setSensitiveParams(sensparams)
    self.setPixelPrototypeParams(pixelProtoypeparams)

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

  def setPrototypeBasePoints(self):
    self.basePointListList = []
    for i in range(self.nprotopixel):
      pointList = []
      for point in self.protopixlistdic[i]["points"].split(','):
        pointList.append(tuple(x for x in point.split()))

      self.basePointListList.append(pointList)
    
  def createLayerElement(self):
    layer = ET.Element('layer')
    ladder = self.createLadder()
    sensitive = self.createSensitive()
    layer.append(ladder)
    layer.append(sensitive)
    stringlayer = ET.tostring(layer)
    pixelInsertPos = stringlayer.find('</layer>')
    layerstart, layerend = stringlayer[:pixelInsertPos], stringlayer[pixelInsertPos:]
    self.writeFile(layerstart)

    self.createPixelPrototype()
    self.createPixelMatrix()

    self.writeFile(layerend)

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


outfile = 'gearTest.xml'

telPositionZ = [0, 100, 200, 500, 600, 700]

def createGearStructure():
  gearTag = ET.Element('gear')
  globalTag = ET.SubElement(gearTag, 'global', detectorName="EUTelescope")
  bfieldTag = ET.SubElement(gearTag, 'BField', type="ConstantBField", x="0", y="0", z="0")
  detectorsTag = ET.SubElement(gearTag, 'detectors')
  detectorTag = ET.SubElement(detectorsTag, 'detector', name="SiPlanes", gearType="SiPlanesParameters")
  siPlanesIDTag = ET.SubElement(detectorTag, 'siplanesID', ID="250")
  siPlanesNumber = ET.SubElement(detectorTag, 'siplanesNumber', number="6")

  stringgearTag = ET.tostring(gearTag)
  layerInsertPos = stringgearTag.find('</detector>')
  openingTags, closingTags = stringgearTag[:layerInsertPos], stringgearTag[layerInsertPos:]
  openingTags += '<layers>'
  closingTags = '</layers>' + closingTags
  gearTag.clear()
  return [openingTags, closingTags]
  

if __name__ == '__main__':
  # base mimosa 26 telescope
  f = open(outfile, "w") # reset the file to empty because use append mode later
  f.close()
  gearBaseTags = createGearStructure()

  with open(outfile, "a") as f:
    f.write(gearBaseTags[0])
    f.flush()
    os.fsync(f.fileno())

  for i in range(len(telPositionZ)):
    layer = Plane(createPixelM26, xmloutfile=outfile, sensparams={"ID":str(i), "positionZ":str(telPositionZ[i])})
    layer.createLayerElement()

  with open(outfile, "a") as f:
    f.write(gearBaseTags[1])
    f.flush()
    os.fsync(f.fileno())
