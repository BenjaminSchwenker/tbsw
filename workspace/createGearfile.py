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
                                "u": , "v",   # column, row coordinates of the pixel, pixel address in raw data
                                "adjacent": , # list of u, v addresses that are considered adjacent to the pixel, comma seperated u v value pares that are separated by whitespaces: "u0 v0, u1 v1, ..."
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
      pitch = 0.0184
      shiftu = -0.5*npixelsU*pitch+0.5*pitch # 0.5 pitch because prototype has centre at 0|0 so need to correct for the half pitch to centre the whole layout at 0|0
      shiftv = -0.5*npixelsV*pitch+0.5*pitch
      centreu = pitch*i+shiftu 
      centrev = pitch*j+shiftv 
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

  #pixel prototyp
  protopixdefault = dict(type="0", npoints="4", points="-0.0092 -0.0092, -0.0092 0.0092, 0.0092 0.0092, 0.0092 -0.0092")

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


outfile = 'gearTest.xml'

telPositionZ = [0, 100, 200, 350, 500, 600, 700, 800]

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
  planesNo = 0
  for i in range(len(telPositionZ)-1):
    if i != 3:
      spList.append({"ID":str(planesNo), "positionZ":str(telPositionZ[i])})
      uCGpList.append({})
      vCGpList.append({})
    if i == 3:
      spList.append({"ID":str(21), "positionZ":str(telPositionZ[i])})
      uCGpList.append({"minCell":"0", "maxCell":"80", "pitch":"0.25"})
      vCGpList.append({"minCell":"0", "maxCell":"336", "pitch":"0.05"})
      planesNo -= 1
    planesNo += 1

  dparams = {"number":"7"}
  detector = Detector(xmloutfile=outfile, sensparamsList=spList, detectorparams=dparams, uCellGroupparamsList=uCGpList, vCellGroupparamsList=vCGpList)
  detector.createDetectorElement()

  # 8th plane in full PolyGear description
  cpfList = [createPixelM26]
  dparams = {"geartype":"PolyPlanesParameters", "number":"1"}
  spList = [{"ID":"7", "positionZ":str(telPositionZ[7])}]
  #pppList = [[{}]] use default 
  detector = Detector(xmloutfile=outfile, createPixelfunctionList=cpfList, sensparamsList=spList, detectorparams=dparams)
  detector.createDetectorElement()

  with open(outfile, "a") as f:
    f.write(gearBaseTags[1])
    f.flush()
    os.fsync(f.fileno())
