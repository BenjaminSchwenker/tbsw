""" Backend class and functions for creation of geometry xml files
    author: Helge C. Beck helge-christoph.beck@phys.uni-goettingen.de

    The script uses the xml.ElementTree library to generate a xml nodes. Purpose is to describe the telescope geometry and the detectors to then be used in the tbsw framework.

    Usage:
    The base xml/gear structure is in the function createGearStructure().
    The Detector class generates a detector xml tag with potentially multiple layers of detectors which are described by the same gear/detector parameters.
      The Parameters of the contructor are list of the parameters used for each layer. It is assumed that no list given or to few list entries that the default parameters should be used. Missing list entries are appended at the end of the list. All lists assume that the index for a plane is the same in every parameter list. 
      The class can be used to generate detectors of two gear type representations: SiPlanesParameters and PolyPlanesParameters. The first one assumes a rectangular pixel grid whereas the second one is totaly without limitations of pixel sizes and overall geometry. The parameters to be set for each are: 
        SiPlanesParameters: uCellGroupparamsList, vCellGroupparamsList ----- List of cellGroup parameters for each layer/plane
                            detectorparams ----------- most importent the geartype and number (number of planes)
        PolyPlanesParameters: pixelPrototypeparamsList --------- List of list of prototype pixel for each layer
                              createPixelfunctionList --------- List of functions for each layer to create the pixel, if not provided or not for every layer will throw a ValueError
                              detectorparams
        All parameter values have to be strings otherwise the parsing does not work.
        The other parameter lists are common to both types and some of the parameters need to be set for a sensible description of the detectors (IDs, positions)
      Parameters to be set can be found in the declared default dictionaries in the class. The list of the respective parameters is set on construction of the detector object. Care has to be taken that the parameters for the pixel prototypes (or u/vCellGroups) are handled in a list of dictionaries ([{},{},]) because there could be multiple pixel types in one plane. 
      For PolyPlanesParameters detector types:
        Most important parameter to set is the list of functions that generate the layout of the pixels the pixel matrix for each layer. It has to be set during construction. This function has to be implemented for the specific layout (in a different file, see createGearfile.py). Characteristics of this function:
        return: No return, all created elements have to be written to file directly.
        First parameter is self of the Detector class so that parameters can be accessed like the xmloutfile name.
        Elements: name is 'pixel'
                  attributes = {"type": "0",  # global type of the pixel
                                "u": , "v",   # column, row coordinates of the pixel, pixel address in raw data, or pseudo coordinates but then a conversion processor needs to be used when reading in the raw data
                                "centeru": ,  # centre coordinate u of pixel in mm, collection point of the charges, assigned hit position in calibration and reconstruction 
                                "centerv": ,  # centre coordinate v of pixel in mm
                               }
        For conveniance there is a method called in Detector class before generating the pixel layout that creates a list of list of tuples. The list of tuples consists of the coordinates of the pixel prototype corners. Every prototype pixel is an entry in the outer list. This list is layer specific and can be accessed with self.basePointListList. Be aware that the provided string values from the prototype dictionaries are converted to floats, precision loss could happen. The first entry for each pixel is its type and not a point!
        Use the function addPixelToTH2Poly(attributes) in the createPixel() function to add the pixel to a TH2Poly object which is written to a root file at the end of the detector creation. Name of the root file from the xmloutfile. 
      Other parameters for Detector class:
        The filename/path where the resulting gear file is written has to be provided to xmloutfile contructer parameter. In the class the write action is always append, so the file should be prepared for the detector elements.
        With the default parameters 6 mimosa26 planes are created with ID=0 and at positionX/Y/Z=0.0. These two parameters are the ones that have to be set for every plane. Not every parameter has to be set. For every plane the default parameters are assigned first and then with the provided parameters overwritten. So not provided parameters stay the default ones. 
        The Detector class does not return anything. All elements are written to the outputfile with the writeFile function.

        basic usage of the module:
        import DetLayoutGen
        #create new .xml file to append everything to
        gearTags = DetLayoutGen.createGearStructure(parameters)
        #write (append) gearTags[0] to file (opening tags)
        detector = DetLayoutGen.Detector(all the parameters you need)
        detector.createDetectorElement()
        #write (append) gearTags[1] to file (closing tags)

"""
import xml.etree.ElementTree as ET
from copy import deepcopy 
import os

from ROOT import TFile, TH2Poly, TGraph

# better styling of xml outfile, like actual linebreaks! -> currently linebreaks are in but not a staggered indentation because there is no information about the indentation where pieces are added to. 
# list of local class variables for reference
class Detector:
  """Base class which holds all parameters and functions for creating a detector with multiple planes which have the same gear type in a gear file 
  Default parameters are for a 6 layer Mimosa26 telescope, all layers at the same position!"""
  
  # gear/detector types available.
  dettypes = dict(SquareDet="SiPlanesParameters", PolyDet="PolyPlanesParameters")

  #detector parameters. Do not be confused by the dettypes and geartype names. Historically gear was used until the swap to only xml files for geometry description. For backwards compatibility the parameter is still called geartype even if it is a detector type now.
  ddefault = dict(name="SiPlanes", geartype=dettypes["SquareDet"], ID="250", number="6")

  # ladder paramaters
  ldefault = dict(sizeU="200", sizeV="100", thickness="0.07", 
                  radLength="93.660734", atomicNumber="14", atomicMass="28")

  # sensitive parameters
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
  # prototype for Mimosa26 pixel with pitch=0.0184mm
  protopixdefault = dict(type="0", distu="0.02", distv="0.02", points="-0.0092 -0.0092, -0.0092 0.0092, 0.0092 0.0092, 0.0092 -0.0092")
  
  # example of the list created by setPrototypeBasePoints() for the default protopixel, first entry is type
  #self.basePointListList = [[0, (-0.0092, -0.0092), (-0.0092, 0.0092), (0.0092, 0.0092), (0.0092, -0.0092)]]

  # variables which are local/individual for every detector object created, local to remove overwriting between the objects. Just for the overview. Do not put back in!

  #self.xmloutfile --- file where the output is written to. Can be accessed by the createPixel(self) function!

  #self.nlayers --- number of layers determined from detectorparams
  #self.detType --- type of detector layout, from detectorparams, possible options in dettypes dictionary 
  #self.ladderparamslist --- ladder parameters for each layer, appended with default values if list to short
  #self.sensparamslist --- sensitive parameters for each layer, "-"
  #self.createPixelFunctionsList --- functions that create pixels for each layer, from createPixelfunctionList, if not or not complete provided raises an exception ValueError
  #self.pixelPrototypeParamsListList --- list of pixel prototypes for each layer, append default if not complete
  #self.uCellGroupParamsListList --- List of uCellGroups for each layer
  #self.vCellGroupParamsListList --- List of vCellGroups for each layer
  #self.uCellGroupParamsList --- uCell parameters for the current layer (for SquareDet type), from uCellGroupparamsList, appended by default if not complete
  #self.vCellGroupParamsList --- vCell parameters for the current layer "-"
  #self.nucellGroups --- number of ucellGroups in the current layer
  #self.nvcellGroups --- number of vcellGroups in the current layer
  #self.protopixlistdic --- List of parameters for each protopixel in the current layer
  #self.nprotopixel
  #self.basePointListList --- list of points for every protopixel in the current layer
  #self.createPixel --- object that gets passed the next used function from the createPixelFunctionsList
  #self.layout --- List of TH2Poly objects for each layer, can be filled with bins with addPixelToTH2Poly(), is written to root file at the end of detector creation (for PolyDet type)
  #self.layer --- current layer number 
  #self.detectorparams --- detector parameters for the class object

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

    if self.detType == self.dettypes["PolyDet"]:
      if createPixelfunctionList is not None:
        if len(createPixelfunctionList) < self.nlayers:
          raise ValueError("ValueError: createPixelfunctionList parameter has less list items ({}) then number of layers ({})".format(len(createPixelfunctionList), self.nlayers))
        self.createPixelFunctionsList = createPixelfunctionList
      else:
        raise ValueError("ValueError: Using {} detector type and not providing a list of functions to create the pixels is not forseen!".format(self.detType))
      if pixelPrototypeparamsList is not None:
        if len(pixelPrototypeparamsList) < self.nlayers:
          for i in range(self.nlayers-len(pixelPrototypeparamsList)):
            pixelPrototypeparamsList.append([deepcopy(self.protopixdefault)])
      else:
        pixelPrototypeparamsList = [[deepcopy(self.protopixdefault)] for i in range(self.nlayers)]
      self.pixelPrototypeParamsListList = pixelPrototypeparamsList

    elif self.detType == self.dettypes["SquareDet"]:
      if uCellGroupparamsList is not None:
        if len(uCellGroupparamsList) < self.nlayers:
          for i in range(self.nlayers-len(uCellGroupparamsList)):
            uCellGroupparamsList.append([deepcopy(uCellGroupdefault)])
      else:
        uCellGroupparamsList = [[deepcopy(self.uCellGroupdefault)] for i in range(self.nlayers)]
      self.uCellGroupParamsListList = uCellGroupparamsList

      if vCellGroupparamsList is not None:
        if len(vCellGroupparamsList) < self.nlayers:
          for i in range(self.nlayers-len(vCellGroupparamsList)):
            vCellGroupparamsList.append([deepcopy(vCellGroupdefault)])
      else:
        vCellGroupparamsList = [[deepcopy(self.vCellGroupdefault)] for i in range(self.nlayers)]
      self.vCellGroupParamsListList = vCellGroupparamsList
      
  def createLadder(self):
    ladder = ET.Element('ladder', self.ladderparams)
    return ladder

  def createSensitive(self):
    sensitive = ET.Element('sensitive', self.sensparams)
    return sensitive

  def createPixelPrototype(self):
    for i in range(self.nprotopixel): 
      self.writeFile(ET.Element('pixelPrototype', self.protopixlistdic[i]))

  def createPixelMatrix(self):
    matrixstart = '<pixelMatrix>'#indentation
    matrixend = '</pixelMatrix>'
    self.writeFile(matrixstart)
    self.createPixel(self)
    self.writeFile(matrixend)

  def setDetectorParams(self, detectorparamsdic):
    self.detectorparams = deepcopy(self.ddefault)
    if detectorparamsdic is not None:
      self.detectorparams.update(detectorparamsdic)
    self.detType = self.detectorparams["geartype"]
    if self.detType not in self.dettypes.values():
      raise ValueError("The provided det/gear type ({}) is not supported (or not in the list of dettypes).".format(self.detType))
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
    self.uCellGroupParamsList = []
    if uCellGroupdict is not None:
      self.nucellGroups = len(uCellGroupdict)
      self.uCellGroupParamsList = [deepcopy(self.uCellGroupdefault) for i in range(self.nucellGroups)]
      for i in range(self.nucellGroups):
        self.uCellGroupParamsList[i].update(uCellGroupdict[i])
    else:
      self.nucellGroups = 1
      self.uCellGroupParamsList.append(self.uCellGroupdefault)

    self.vCellGroupParamsList = []
    if vCellGroupdict is not None:
      self.nvcellGroups = len(vCellGroupdict)
      self.vCellGroupParamsList = [deepcopy(self.vCellGroupdefault) for i in range(self.nvcellGroups)]
      for i in range(self.nvcellGroups):
        self.vCellGroupParamsList[i].update(vCellGroupdict[i])
    else:
      self.nvcellGroups = 1
      self.vCellGroupParamsList.append(self.vCellGroupdefault)
    
  def setPrototypeBasePoints(self):
    self.basePointListList = []
    for i in range(self.nprotopixel):
      pointList = []
      pointList.append(int(float(self.protopixlistdic[i]["type"])))
      for point in self.protopixlistdic[i]["points"].split(','):
        pointList.append(tuple(float(x) for x in point.split()))
      self.basePointListList.append(pointList)
  
  def createDetectorElement(self):
    detector = ET.Element('detector', name=self.detectorparams["name"], geartype=self.detectorparams["geartype"])
    siplanesID = ET.SubElement(detector, 'siplanesID', ID=self.detectorparams["ID"])
    siplanesNumber = ET.SubElement(detector, 'siplanesNumber', number=self.detectorparams["number"])
    indent(detector)
    stringdetector = ET.tostring(detector)
    layerInsertPos = stringdetector.find('</detector>')
    layersstart, layersend = stringdetector[:layerInsertPos], stringdetector[layerInsertPos:]
    layersstart += '<layers>'# indentation?
    layersend = '</layers>' + layersend
    self.writeFile(layersstart)

    self.layout = [] # list of TH2Poly objects for each layer, filled with addPixelToTH2Poly() if put in to createPixel() function

    for i in range(self.nlayers):
      self.layer = i
      self.createLayerElement(i) 

    self.writeFile(layersend)
    self.layer = -1
    self.writeLayout()

  def createLayerElement(self, layernumber):
    layer = ET.Element('layer')
    self.setLadderParams(self.ladderparamslist[layernumber])
    self.setSensitiveParams(self.sensparamslist[layernumber])
    ladder = self.createLadder()
    sensitive = self.createSensitive()
    layer.append(ladder)
    layer.append(sensitive)
    indent(layer)
    stringlayer = ET.tostring(layer)
    pixelInsertPos = stringlayer.find('</layer>')
    layerstart, layerend = stringlayer[:pixelInsertPos], stringlayer[pixelInsertPos:]
    self.writeFile(layerstart)

    self.layout.append(TH2Poly())
    self.layout[layernumber].SetName("layer"+str(layernumber))

    if self.detType == self.dettypes["SquareDet"]:
      self.setCellGroupParams(self.uCellGroupParamsListList[layernumber], self.vCellGroupParamsListList[layernumber])
      self.createPixel = self.createRectPixel
      self.createPixel()
    elif self.detType == self.dettypes["PolyDet"]:
      self.createPixel = self.createPixelFunctionsList[layernumber]
      self.setPixelPrototypeParams(self.pixelPrototypeParamsListList[layernumber])
      self.createPixelPrototype()
      self.createPixelMatrix()

    self.writeFile(layerend)
  
  def createRectPixel(self):
    if self.nucellGroups == self.nvcellGroups:
      for i in range(self.nucellGroups):
        ucells = ET.Element('uCellGroup', self.uCellGroupParamsList[i]);
        vcells = ET.Element('vCellGroup', self.vCellGroupParamsList[i]);
        self.writeFile(ucells);
        self.writeFile(vcells);
    else:
      raise ValueError("ValueError: not equal number of uCellGroups ({}) and vCellGroups ({}). Have to be equal!".format(self.nucellGroups, self.nvcellGroups))

  def addPixelToTH2Poly(self, attributes):
    pointsList = [plist for plist in self.basePointListList if plist[0] == int(float(attributes["type"]))][0]
    gpixel = TGraph(len(pointsList)-1)
    gpixel.SetName(attributes["u"] + "," + attributes["v"])
    i = 0
    for point in pointsList[1:]:
      gpixel.SetPoint(i, point[0]+float(attributes["centeru"]), point[1]+float(attributes["centerv"]))
      i += 1
    self.layout[self.layer].AddBin(gpixel)

  def writeLayout(self):
    if self.detType != self.dettypes["PolyDet"]:
      return
    index = self.xmloutfile.find(".xml")
    if index == -1:
      print("Can not parse xmloutfile name to a rootfile (.xml -> .root) because no .xml found. Not writing Th2Poly layout!")
      return
    rootfilename = self.xmloutfile[0:self.xmloutfile.find(".xml")] + ".root"
    rootfile = TFile(rootfilename, "RECREATE")
    detTypedir = rootfile.mkdir(self.detType)
    detTypedir.cd()
    layoutdir = detTypedir.mkdir("layouts")
    layoutdir.cd()
    for polyhist in self.layout:
      if polyhist.GetNumberOfBins > 0:
        polyhist.Draw()
        polyhist.Write()

    detTypedir.cd()
    testdir = detTypedir.mkdir("test")
    testdir.cd()
    for polyhist in self.layout:
      # fill every bin with 1
      nobins = polyhist.GetNumberOfBins()
      if nobins > 0:
        for i in range(nobins):
          polyhist.Fill(polyhist.GetBinName(i), 1)
        polyhist.Draw("COLZ")
        polyhist.Write()
        
    rootfile.Write()
    rootfile.Close()

  def writeFile(self, writeobject):  # write object is either an ET.Element or a string
    writestring = ""
    if ET.iselement(writeobject):
      indent(writeobject)
      writestring = ET.tostring(writeobject)
    else:
      writestring = writeobject
    with open(self.xmloutfile, "a") as f:
      f.write(writestring)
      f.flush()
      os.fsync(f.fileno())

# creates the basic gear/xml layout file structure. Returns this structur in two parts so that the detectors generated with the Detector class can be inserted in between when writing the xml file.
def createGearStructure(detName="EUTelescope", bFieldtype="ConstantBField", bField=[0,0,0]):
  gearTag = ET.Element('gear')
  globalTag = ET.SubElement(gearTag, 'global', detectorName=detName)
  bfieldTag = ET.SubElement(gearTag, 'BField', type=bFieldtype, x=str(bField[0]), y=str(bField[1]), z=str(bField[2]))
  indent(gearTag)
  stringgearTag = ET.tostring(gearTag)
  layerInsertPos = stringgearTag.find('</gear>')
  openingTags, closingTags = stringgearTag[:layerInsertPos], stringgearTag[layerInsertPos:]
  openingTags += '<detectors>'
  closingTags = '</detectors>' + closingTags
  gearTag.clear()
  return [openingTags, closingTags]


# function for indenting a Elementtree from: http://effbot.org/zone/element-lib.htm#prettyprint
def indent(elem, level=0):
  i = "\n" + level*"  "
  if len(elem):
    if not elem.text or not elem.text.strip():
      elem.text = i + "  "
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
    for elem in elem:
      indent(elem, level+1) 
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
  else:
    if level and (not elem.tail or not elem.tail.strip()):
      elem.tail = i
