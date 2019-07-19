""" Backend classes and functions for creation of geometry xml files
    
    
    The xml.ElementTree library is used to generate all xml nodes. Purpose is to describe
    the telescope geometry and the detectors to then be used in the tbsw framework.

    Usage:
    
    A list of detectors must be created. All detectors are objects created from either the 
    class SquareDetector or the class PolyDetector. Both classes inherit from the class 
    BaseDetector.   

    The function WriteGearfile(outfile, sensors) writes the output gearfile.  
    
        
    author: Helge C. Beck helge-christoph.beck@phys.uni-goettingen.de
    author: Benjamin Schwenker benjamin.schwenker@phys.uni-goettingen.de
"""

import xml.etree.ElementTree as ET
from copy import deepcopy 
import os

from ROOT import TFile, TH2Poly, TGraph


class BaseDetector(object):
  """
  BaseDetector class implements an interface define a detector in a geomtry xml file. 
   
  sensitiveParams -----  Dictionary of parameters describing the sensitive part of the detector 
  supportParams   -----  Dictionary of parameters describing the support parts of the sensor (outside acceptance) 
  gearType        -----  Gear Type (str) 

  The unit for lenght is mm. The unit of angles are degrees. Initial values for the sensitiveParams and supportParams
  dictionary are set in the constructor. All values can be overriden by the member functions sensParams() and suppParams(). 
  
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  def __init__(self):
    """  
    Initialize with default values for the sensitive and parameter dictionaries. 
    """    
    self.gearType = None
    
    self.sensitiveParams = {
                             "ID": "0",
                             "positionX": "0.00",  
                             "positionY": "0.00", 
                             "positionZ": "0.00",  
                             "alpha": "0",          
                             "beta": "0",              
                             "gamma": "0",         
                             "rotation1": "1.0", 
                             "rotation2": "0.0",  
                             "rotation3": "0.0", 
                             "rotation4": "1.0",
                             "PixType": "0",
                             "thickness": "0.07",  
                             "radLength": "93.660734", 
                             "atomicNumber": "14",
                             "atomicMass": "28", 
                           }
      
    self.supportParams = {
                           "sizeU": "200", 
                           "sizeV": "100",  
                           "thickness": "0.07",  
                           "radLength": "93.660734", 
                           "atomicNumber": "14", 
                           "atomicMass": "28",   
                         }
    
  def sensParams(self, params):
    """  
    Updates dictionary of sensitive parameters.  
    """    
    self.sensitiveParams.update(params)

  def suppParams(self, params):
    """  
    Updates dictionary of support parameters.  
    """    
    return self.supportParams.update(params)
    

class SquareDetector(BaseDetector):
  """
  SquareDetector class implements an interface to define a detector with a pixel matrix having
  only rectangular pixel cells.  

        
  uCellGroups ---- List of cellGroup dictionaries describing the segmentation of sensor along u axis
  vCellGroups ---- List of cellGroup dictionaries describing the segmentation of sensor along u axis 

  The cellGroup dictionary is expected to look like: 

  cellGroup = {"minCell": 0, "maxCell": 335, "pitch": 0.05 }
  
  It tells us that pixels with numbers [0, 335] have a pitch of 0.05mm. There can be multiple cellGroups with non
  overlapping numbers to reflect a change in the pixel pitch along one sensor direction. A 2d grid is defined from 
  the list of cellGroups for the u and v directions.   
  
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  def __init__(self, sensParams={}, suppParams={}):
    """  
    """ 
    # Set default values for sensitive and support params
    BaseDetector.__init__(self)

    # Set gearType variable
    self.gearType = "SiPlanesParameters" 

    # Update params
    self.suppParams(suppParams)
    self.sensParams(sensParams)
   
    # By default, there are no pixels defined 
    self.uCellGroups = []
    self.vCellGroups = []  
    
       
  def addUCellGroup(self, cellGroup):
    """
    Add cell group for u axis 
    """
    self.uCellGroups.append(cellGroup) 

  def addVCellGroup(self, cellGroup):
    """
    Add cell group for v axis 
    """
    self.vCellGroups.append(cellGroup) 

class PolyDetector(BaseDetector):
  """
  PolyDetector class implements an interface to define a detector with a pixel matrix having
  polygonal shaped pixel cells.  

  pixelShapes      ----- List of pixelShape dictionaries describing all different pixel types 
  generatePixels   ----- Generator function to create pixel cells placed in the sensitive area
  
  The pixelShape dictionary is expected to look like
 
  pixelShape = {
                 "type": 0,         # int, unique for the detector   
                 "distu": 0.3,      # float, extension along u direction in mm 
                 "distv": 0.01,     # float, extension along v direction in mm 
                 "points": [-0.125, -0.025, -0.125, 0.025, 0.125, 0.025, 0.125, -0.025]
               } 
   
  The type (int) of the pixelShape is a unique integer value. The distu and distv parameters meaure the extension of the pixel cell 
  in the u and v direction and are used to define the neighborhood of the pixel. The points are a list of edge points of the polygonal
  outer shape of the pixel. The center of the pixel at 0,0 can be the geometric center or the charge collection center. 

  The generatePixels generator function yields a sequence of all pixels placed on the sensitive area. A placed pixel is defined by 
  a dictionary of the following form. 

  pixel = {
            "type": "0",  # global type of the pixel
            "u": , "v",   # column, row coordinates of the pixel, pixel address in raw data
            "centeru": ,  # centre coordinate u of pixel in mm
            "centerv": ,  # centre coordinate v of pixel in mm
          } 
  
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  
  def __init__(self, pixelShapes=[], generatePixels=None):
    """  
    :@pixelShapes: list of pixel shape dictionaries 
    """ 
    
    BaseDetector.__init__(self)
    self.pixelShapes = []
    self.pixelShapes.extend(pixelShapes)
    self.gearType = "PolyPlanesParameters"     
     
    def generateNoPixels():
      pass
    
    if not generatePixels == None: 
      self.generatePixels = generatePixels
    else: 
      self.generatePixels = generateNoPixels
   
  def addPixelShapes(self, pixelShapes):
    """
    Extends the list of pixel shapes with new shapes 
    """
    self.pixelShapes.extend(pixelShapes)


def WriteGearfile(xmloutfile, sensors=[]):
  """Create a geometry xml file with path xmloutfile for the specified list of sensors.
  """

  # reset the file to empty because use append mode later
  f = open(xmloutfile, "w") 
  f.close()
  # create the basic xml file structure, for parameters see DetLayoutGen.py
  gearBaseTags = createGearStructure()
  # write the opening tags of the basic xml structure. This is needed because the 
  # xml nodes are itterativly created and written to file to reduce memory consumption
  with open(xmloutfile, "a") as f:
    f.write(gearBaseTags[0])
    f.flush()
    os.fsync(f.fileno())

  # Filter the sensors 
  square_dets = []
  poly_dets = []
  for sensor in sensors: 
    if sensor.gearType == "SiPlanesParameters":
      square_dets.append(sensor)
    elif sensor.gearType == "PolyPlanesParameters": 
      poly_dets.append(sensor)
    else: 
      print("Sensor with unkown gearType={}. Skipping it!".format(sensor.gearType)) 
   
  # Write all SquareDets into xmloutfile    
  WriteDetsIntoGearfile(xmloutfile, name="SiPlanes", sensors=square_dets)
  
  #  Write all PolyDets into xmloutfile 
  WriteDetsIntoGearfile(xmloutfile, name="PolyPlanes", sensors=poly_dets)
  
  # write the closing tags of the basic xml structure
  with open(xmloutfile, "a") as f:
    f.write(gearBaseTags[1])
    f.flush()
    os.fsync(f.fileno())
 
def WriteDetsIntoGearfile(xmloutfile, name, sensors=[]):
  
  if len(sensors) == 0: 
    return 
  
  nlayers = len(sensors)
  gearType = sensors[0].gearType
  dparams = {"name": name, "ID": 250, "geartype": gearType, "number":str(nlayers)}
  
  detector = ET.Element('detector', name=dparams["name"], geartype=dparams["geartype"])
  siplanesID = ET.SubElement(detector, 'siplanesID', ID=str(dparams["ID"]))
  siplanesNumber = ET.SubElement(detector, 'siplanesNumber', number=str(dparams["number"]))
  indent(detector)
  stringdetector = ET.tostring(detector)
  layerInsertPos = stringdetector.find('</detector>')
  layersstart, layersend = stringdetector[:layerInsertPos], stringdetector[layerInsertPos:]
  layersstart += '<layers>'# indentation?
  layersend = '</layers>' + layersend
  writeFile(layersstart, xmloutfile)
  
  # Dictionary of TH2Poly objects for each layer, filled with addPixelToTH2Poly()
  layout = {} 
  
  for sensor in sensors:
    createLayerElement(sensor, layout, xmloutfile) 
   
  writeFile(layersend, xmloutfile)
  #TODO fix this function
  #writeLayout(xmloutfile, gearType, layout)
  

def createPixelMatrix(generatePixels, xmloutfile):
  matrixstart = '<pixelMatrix>'#indentation
  matrixend = '</pixelMatrix>'
  writeFile(matrixstart, xmloutfile)
  generatePixels()
  writeFile(matrixend, xmloutfile)
  
  
def setPrototypeBasePoints():
  basePointListList = []
  for i in range(self.nprotopixel):
    pointList = []
    pointList.append(int(float(self.protopixlistdic[i]["type"])))
    for point in self.protopixlistdic[i]["points"].split(','):
      pointList.append(tuple(float(x) for x in point.split()))
    basePointListList.append(pointList)
  return basePointList
  
  
def addPixelToTH2Poly(attributes):
  pointsList = [plist for plist in self.basePointListList if plist[0] == int(float(attributes["type"]))][0]
  gpixel = TGraph(len(pointsList)-1)
  gpixel.SetName(attributes["u"] + "," + attributes["v"])
  i = 0
  for point in pointsList[1:]:
    gpixel.SetPoint(i, point[0]+float(attributes["centeru"]), point[1]+float(attributes["centerv"]))
    i += 1
  self.layout[self.layer].AddBin(gpixel)

def writeLayout(xmloutfile, gearType, layout):
  if gearType != "PolyPlanesParameters":
    return
  index = xmloutfile.find(".xml")
  if index == -1:
    print("Can not parse xmloutfile name to a rootfile (.xml -> .root) because no .xml found. Not writing Th2Poly layout!")
    return
  rootfilename = xmloutfile[0:xmloutfile.find(".xml")] + ".root"
  rootfile = TFile(rootfilename, "RECREATE")
  detTypedir = rootfile.mkdir(gearType)
  detTypedir.cd()
  layoutdir = detTypedir.mkdir("layouts")
  layoutdir.cd()
  for polyhist in layout:
    if polyhist.GetNumberOfBins() > 0:
      polyhist.Draw()
      polyhist.Write()

  detTypedir.cd()
  testdir = detTypedir.mkdir("test")
  testdir.cd()
  for sensorID, polyhist in layout.iteritems():
    # fill every bin with 1
    nobins = polyhist.GetNumberOfBins
    if nobins > 0:
      for i in range(1, nobins+1): # bins start at 1 and end at including nobins
        polyhist.Fill(polyhist.GetBinName(i), 1)
      polyhist.Draw("COLZ")
      polyhist.Write()
        
  rootfile.Write()
  rootfile.Close()

  
def createRectPixel(uCellGroups, vCellGroups, xmloutfile):
  for uCellGroup in uCellGroups:
    ucells = ET.Element('uCellGroup', {str(k) : str(v) for k, v in uCellGroup.iteritems()})
    writeFile(ucells, xmloutfile)
  for vCellGroup in vCellGroups:
    vcells = ET.Element('vCellGroup', {str(k) : str(v) for k, v in vCellGroup.iteritems()})
    writeFile(vcells, xmloutfile)

def createLayerElement(sensor, layout, xmloutfile):
  layer = ET.Element('layer')
  ladder = ET.Element('ladder', {str(k) : str(v) for k, v in sensor.supportParams.iteritems()})  
  sensitive = ET.Element('sensitive', {str(k) : str(v) for k, v in sensor.sensitiveParams.iteritems()})
  layer.append(ladder)
  layer.append(sensitive)
  indent(layer)
  stringlayer = ET.tostring(layer)
  pixelInsertPos = stringlayer.find('</layer>')
  layerstart, layerend = stringlayer[:pixelInsertPos], stringlayer[pixelInsertPos:]
  writeFile(layerstart, xmloutfile)
               
  sensorID = int(sensor.sensitiveParams["ID"]) 
  layout[sensorID] = TH2Poly()
  layout[sensorID].SetName("layer"+str(sensorID))
  
  if sensor.gearType == "SiPlanesParameters":
     createRectPixel(sensor.uCellGroups, sensor.vCellGroups, xmloutfile)
  elif sensor.gearType == "PolyPlanesParameters":
    for pixelShape in sensor.pixelShapes: 
      writeFile(ET.Element('pixelPrototype', pixelShape), xmloutfile)
    createPixelMatrix(sensor.generatePixels, xmloutfile)
  writeFile(layerend, xmloutfile)

def writeFile(writeobject, xmloutfile):  
  writestring = ""
  if ET.iselement(writeobject):
    indent(writeobject)
    writestring = ET.tostring(writeobject)
  else:
    writestring = writeobject
  with open(xmloutfile, "a") as f:
    f.write(writestring)
    f.flush()
    os.fsync(f.fileno())


def createGearStructure(detName="EUTelescope", bFieldtype="ConstantBField", bField=[0,0,0]):
  """Creates the basic xml structure of the gearfile. 
     Returns this structur in two parts so that the detectors 
     generated with the Detector class can be inserted in
     between when writing the xml file.
  """
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
