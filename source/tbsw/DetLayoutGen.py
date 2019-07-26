""" Backend classes and functions for creation of geometry xml files
    
    
    The xml.ElementTree library is used to generate all xml nodes. Purpose is to describe
    the telescope geometry and the detectors to then be used in the tbsw framework.

    Usage:
    
    A list of detectors must be created. All detectors are objects created from either the 
    class SquareDetector or the class PolyDetector. Both classes inherit from the class 
    BaseDetector.   

    The function WriteGearfile(outfile, sensors) writes the output gearfile.  
    
    The function WriteLayoutRootfile(outfile, sensors) writes a root file with histograms to plot the pixel matrices of sensors.   
        
    author: Helge C. Beck helge-christoph.beck@phys.uni-goettingen.de
    author: Benjamin Schwenker benjamin.schwenker@phys.uni-goettingen.de
"""

import xml.etree.ElementTree as ET
from copy import deepcopy 
import os
import array

from ROOT import TFile, TCanvas, TH2Poly, TH2F, TGraph, TLine


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

  The groups are sorted internaly accending in minCell. 

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
    Add cell group for u axis and sort the groups accending in minCell
    """
    self.uCellGroups.append(cellGroup) 
    self.uCellGroups = sorted(self.uCellGroups, key=lambda k: k['minCell'])

  def addVCellGroup(self, cellGroup):
    """
    Add cell group for v axis and sort the groups accending in minCell
    """
    self.vCellGroups.append(cellGroup) 
    self.vCellGroups = sorted(self.vCellGroups, key=lambda k: k['minCell'])

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
                 "points": [(-0.125, -0.025), (-0.125, 0.025), (0.125, 0.025), (0.125, -0.025)]
               } 
   
  The type (int) of the pixelShape is a unique integer value. type == -1 is used for none described areas of the layout. 
  The distu and distv parameters measure the extension of the pixel cell in the u and v direction and are used to define the neighborhood of the pixel. 
  -> The distu and distv parameters give the distance from the pixel center in which gets checked for neighboring pixel aka. other pixel centers. TODO
  The points are a list of edge points of the polygonal outer shape of the pixel. They polygon is built with the points in order of the list. 
  The center of the pixel at 0,0 can be the geometric center or the charge collection center. In tbsw it is used as the hit position of the pixel.

  The generatePixels generator function yields a sequence of all pixels placed on the sensitive area. A placed pixel is defined by 
  a dictionary of the following form. 

  pixel = {
            "type": 0,     # global type of the pixel
            "u": 0,        # uCell (column) address of pixel in raw data 
            "v": 0,        # vCell (row) address of the pixel in raw data
            "centeru": 0,  # center coordinate u of pixel in mm
            "centerv": 0,  # center coordinate v of pixel in mm
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
  """
  Create a geometry xml file with path xmloutfile for the specified list of sensors.
  """

  # reset the file to empty because use append mode later
  f = open(xmloutfile, "w") 
  f.close()
  # create the basic xml file structure
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
  
  for sensor in sensors:
    createLayerElement(sensor, xmloutfile) 
    print("Sensor with ID=" + str(sensor.sensitiveParams["ID"]) + " written to xml file")
  
  writeFile(layersend, xmloutfile)
  
  
def WriteLayoutRootfile(outfile, sensors=[], maxPixel=30000, fillTest=True): 
  """
  Create a root file containing TH2Poly histograms representing the matrix layout 
  of the PolyDet sensors and TH2F histograms representing the SquareDet sensors. 
  Limit the number of pixels to maxPixel per sensor to reduce processing time.
  Filling TH2Poly for testing takes also very long, so disable it with fillTest=False.
  """
  # Dictionary of TH2F objects for each square det layer, filed with 
  layoutSquareDet = {} 
  # Dictionary of TH2Poly objects for each poly det layer, filled with addPixelToTH2Poly()
  layoutPolyDet = {}
  for sensor in sensors:
    if sensor.gearType == "SiPlanesParameters":
      sensorID, squarehist = createTH2F(sensor)
      layoutSquareDet[sensorID] = squarehist
    elif sensor.gearType == "PolyPlanesParameters":
      sensorID, polyhist = createTH2Poly(sensor, maxPixel)
      layoutPolyDet[sensorID] = polyhist 
    
  writeLayout(outfile, layoutSquareDet, layoutPolyDet, fillTest)


def createTH2F(sensor):
  sensorID = int(sensor.sensitiveParams["ID"])

  lowEdgesU, lowEdgesV = generateSquarePixelBins(sensor.uCellGroups, sensor.vCellGroups)
  
  squarehist = TH2F("layer"+str(sensorID), "", len(lowEdgesU)-1, array.array('d', lowEdgesU), len(lowEdgesV)-1, array.array('d', lowEdgesV))
  return sensorID, squarehist


def createTH2Poly(sensor, maxPixel):
  sensorID = int(sensor.sensitiveParams["ID"])              
  polyhist = TH2Poly()
  polyhist.SetName("layer"+str(sensorID))
  
  pixelCount = 0
  for pixel in sensor.generatePixels():
    edges = [pixelShape["points"] for pixelShape in sensor.pixelShapes if pixelShape["type"] == int(pixel["type"])][0]   
    addPixelToTH2Poly(pixel, edges, polyhist)
    pixelCount += 1
    if pixelCount >= maxPixel: 
      print("TH2Poly histogram for Sensor with ID={} truncated after maxPixel={}".format(sensorID, maxPixel))
      break
      
  return sensorID, polyhist


def generateSquarePixelBins(uCellGroups, vCellGroups):
  """
  Generate list of lower bin edges aka pixel edges for display in TH2F. 
  Layout starts at 0,0 (in tbsw then centred but not here).
  """
  lowEdgesU = [0,]
  lowEdgesV = [0,]
  
  i = 1
  for vCell in vCellGroups:
    pitchV = vCell["pitch"]
    for v in range(vCell["minCell"], vCell["maxCell"] + 1):
      lowEdgesV.append(lowEdgesV[i-1] + pitchV)
      i += 1
  
  i = 1
  for uCell in uCellGroups:
    pitchU = uCell["pitch"]
    for u in range(uCell["minCell"], uCell["maxCell"] + 1):
      lowEdgesU.append(lowEdgesU[i-1] + pitchU)
      i += 1
  
  return lowEdgesU, lowEdgesV


def createPixelMatrix(sensor, xmloutfile):
  matrixstart = '<pixelMatrix>'#indentation
  matrixend = '</pixelMatrix>'
  writeFile(matrixstart, xmloutfile)
  
  pixelList = []
  count = 0
  for pixel in sensor.generatePixels():
    pixelList.append(ET.Element('pixel', {str(k) : str(v) for k, v in pixel.iteritems()} )) 
    count += 1
    if count == 500:
      writePixelListToFile(xmloutfile, pixelList)
      pixelList[:] = []
      count = 0  
  
  writePixelListToFile(xmloutfile, pixelList)
  writeFile(matrixend, xmloutfile)

def writePixelListToFile(xmlfilename, pixelList):
  with open(xmlfilename, "a") as f:
    for element in pixelList:
      indent(element)
      f.write(ET.tostring(element)+"\n")
    f.flush()
    os.fsync(f.fileno())
    pixelList[:] = []


def addPixelToTH2Poly(pixel, edges, polyhist): 
  gpixel = TGraph(len(edges))
  gpixel.SetName(str(pixel["u"]) + "," + str(pixel["v"]))
  
  for i, point in enumerate(edges):
    gpixel.SetPoint(i, float(point[0])+float(pixel["centeru"]), float(point[1])+float(pixel["centerv"]))
  polyhist.AddBin(gpixel)

def writeLayout(rootfilename, layoutSquareDet, layoutPolyDet, fillTest):
  rootfile = TFile(rootfilename, "RECREATE")
  layoutdir = rootfile.mkdir("layouts")
  layoutdir.cd()
  for sensorID, squarehist in layoutSquareDet.iteritems():
    nbinsU = squarehist.GetNbinsX()
    nbinsV = squarehist.GetNbinsY()
    print("Drawing sensor with ID=" + str(sensorID) + " with " + str(nbinsU) + " times " + str(nbinsV) + " bins.")
    # Generating a grid of lines representing the bins/pixels
    canvas = TCanvas("c"+str(sensorID))
    canvas.cd()
    squarehist.Draw()
    line = TLine()
    canvas.Update()
    umin = canvas.GetUxmin()
    umax = canvas.GetUxmax()
    vmin = canvas.GetUymin()
    vmax = canvas.GetUymax()
    lowEdgesU = array.array('d', [0.0]*(nbinsU+1)) 
    lowEdgesV = array.array('d', [0.0]*(nbinsV+1))
    squarehist.GetXaxis().GetLowEdge(lowEdgesU)
    squarehist.GetYaxis().GetLowEdge(lowEdgesV)
    for lowEdgeU in lowEdgesU:
      line.DrawLine(lowEdgeU, vmin, lowEdgeU, vmax)
    for lowEdgeV in lowEdgesV:
      line.DrawLine(umin, lowEdgeV, umax, lowEdgeV)
    canvas.Draw()
    canvas.Write()
    squarehist.Write()

  for sensorID, polyhist in layoutPolyDet.iteritems():
    if polyhist.GetNumberOfBins() > 0:
      print("Drawing sensor " + str(sensorID) + " with " + str(polyhist.GetNumberOfBins()) + " bins.")
      polyhist.Draw()
      polyhist.Write()
  
  if fillTest:
    print("Warning: No test filling for SquareDet sensors!")
    rootfile.cd()
    testdir = rootfile.mkdir("test")
    testdir.cd()
    for sensorID, polyhist in layoutPolyDet.iteritems():
      # fill every bin with 1
      nobins = polyhist.GetNumberOfBins()
      if nobins > 0:
        for i in range(1, nobins+1): # bins start at 1 and end at including nobins
          polyhist.Fill(polyhist.GetBinName(i), 1)
        polyhist.Draw("COLZ")
        polyhist.Write()
        print("Test filling of sensor " + str(sensorID) + " done.")
        
  rootfile.Write()
  rootfile.Close()

  
def createRectPixel(uCellGroups, vCellGroups, xmloutfile):
  for uCellGroup in uCellGroups:
    ucells = ET.Element('uCellGroup', {str(k) : str(v) for k, v in uCellGroup.iteritems()})
    writeFile(ucells, xmloutfile)
  for vCellGroup in vCellGroups:
    vcells = ET.Element('vCellGroup', {str(k) : str(v) for k, v in vCellGroup.iteritems()})
    writeFile(vcells, xmloutfile)

def createLayerElement(sensor, xmloutfile):
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
               
  if sensor.gearType == "SiPlanesParameters":
     createRectPixel(sensor.uCellGroups, sensor.vCellGroups, xmloutfile)
  elif sensor.gearType == "PolyPlanesParameters":
    for pixelShape in sensor.pixelShapes: 
      # We need to serialize the list of edge points into a string 
      tmpPixelShape = deepcopy(pixelShape)
      tmpPixelShape["points"] = [ "{} {}".format(str(point[0]), str(point[1])) for point in tmpPixelShape["points"] ]
      tmpPixelShape["points"] = " ".join(tmpPixelShape["points"])  
      writeFile(ET.Element('pixelPrototype', {str(k) : str(v) for k, v in tmpPixelShape.iteritems()} ), xmloutfile)
    createPixelMatrix(sensor, xmloutfile)
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
  """
  Creates the basic xml structure of the gearfile. 
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



def indent(elem, level=0):
  """
  Function for indenting a Elementtree from: http://effbot.org/zone/element-lib.htm#prettyprint.
  """
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
