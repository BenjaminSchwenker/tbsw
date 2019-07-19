""" Example python script for generating a xml geometry description (aka gearfile) of a tracking telescope for data 
    analysis with tbsw. The script uses the xml.ElementTree library to generate the xml file.  

    Usage: 
    python createGearfile.py  
    
    In the __main__ part, an example is given for the creation of a gear file for a telescope with 6 mimosa26 detector
    planes and one Atlas FEI4 reference plane created with the detector type SiPlanesParameters. The dut plane is 
    created with the gear type PolyPlanesParameters.
      
    author: Helge C. Beck helge-christoph.beck@phys.uni-goettingen.de
    author: Benjamin Schwenker benjamin.schwenker@phys.uni-goettingen.de
"""

import xml.etree.ElementTree as ET 
import os
import sys
import tbsw
import dummy

class BaseDetector(object):
  """
  BaseDetector class implements an interface define a detector in a geomtry xml file. 
   
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
  
  def getSensParams(self): 
    """  
    Return dictionary of sensitive parameters.
    All length units are in millimeters and angles are in degree.  
    """    
    return self.sensitiveParams

  def getSuppParams(self): 
    """  
    Return dictionary of support parameters.
    All length units are in millimeters and angles are in degree.  
    """    
    return self.supportParams
    

class SquareDetector(BaseDetector):
  """
  SquareDetector class implements an interface to define a detector with a pixel matrix having
  rectangular pixel cells.  

  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  def __init__(self, sensParams={}, suppParams={}):
    """  
    """ 
    # Set default values for sensitive and support params
    BaseDetector.__init__()

    # Set gearType variable
    self.gearType = "SiPlanesParameters" 

    # Update params
    self.suppParams(suppParams)
    self.sensParams(sensParams)

    # By default, there are no pixels defined 
    self.uCGpList = list()  
    self.vCGpList = list()  
    
       
  def getUCellGroups(self):
    """
    Returns list of cell groups for u axis 
    """
    return self.uCGpList 

  def getVCellGroups(self):
    """
    Returns list of cell groups for v axis 
    """
    return self.vCGpList 

class PolyDetector(BaseDetector):
  """
  PolyDetector class implements an interface to define a detector with a pixel matrix having
  polygonal shaped pixel cells.  
   
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  
  def __init__(self, pixelParams=[], createPixels=None):
    """  
    :@pixelPrototypeParams: pixel prototype parameter dictionary
    """ 
    
    BaseDetector.__init__()
    self.pixelParams = list().extend(pixelParams)
    self.gearType = "PolyPlanesParameters"     
     
    def createEmptyPixelList():
      return list()
    
    if not createPixels == None: 
      self.createPixelListFunc = createPixels
    else: 
      self.createPixelFunc = createEmptyPixelList
   
  def pixelParams(self):
    """
    Returns list of pixel params. 
    """
    return self.pixelParams
 
if __name__ == '__main__':
  
  # Name of output gearfile 
  outfile = 'gear_dia_planar_batch2_alltype.xml'
  
  # In this example, we construct a tracking telescope consisting of six mimosa26 planes, a FEI4 
  # timing plane and a 2d diamond sensor with matrix layout that is a mixture of polygonal and 
  # rectangular pixels.   
  telescope = list()
   
  # Define z positions (in mm) and sensor ID for the mimosa26 reference sensors. 
  telPositionZ = [0, 150, 300, 653, 743, 898]
  telID = [0, 1, 2, 3, 4, 5]
   
  # Create six mimosa 26 reference sensors 
  for i, z in zip(telID, telPositionZ):
    # Create detector and override default sensor ID and z position
    m26 = SquareDetector( sensitiveParams={"ID": i, "positionZ": z} )
    # Define pixel numbering and segmentation in u direction     
    m26.getUCellGroups().append( dict("minCell": 0, "maxCell": 1151, "pitch": 0.0184) )
    # Define pixel numbering and segmentation in v direction   
    m26.getVCellGroups().append( dict("minCell": 0, "maxCell": 575, "pitch": 0.0184)  )
    # Add sensor to telescope 
    telescope.append(m26)
  
  # Create Atlas FEi4 reference sensor 
  fei4 = SquareDetector( sensitiveParams={"ID": 22, "positionZ": 494, "thickness": 1.0} )    
  fei4.getUCellGroups() = [{"minCell": 0, "maxCell": 79, "pitch": 0.25 }] 
  fei4.getVCellGroups() = [{"minCell": 0, "maxCell": 335, "pitch": 0.05 }]
  telescope.append(fei4)
  
  # Create a costum DUT sensor with hexagonal pixel cells 
  dut = PolyDetector( createPixel2Ddiamond  )
  dut.sensParams( {"ID":21, "positionZ": 403, "thickness": 0.824, "radLength": 121.3068, "atomicNumber": 6, "atomicMass": 12} )
  dut.suppParams( {"radLength": 121.3068, "atomicNumber": 6, "atomicMass": 12} ) 
  
  # Define the outer shape of all pixels to be placed on the sensitive area of the DUT 
  pixelDUT = dict("type": 0, "distu": 0.3, "distv": 0.01, "points": [-0.125, -0.025, -0.125, 0.025, 0.125, 0.025, 0.125, -0.025]) 
  # Add a =list of the pixel outer shapes to the DUT 
  dut.pixelParams() = [ pixelDUT ]   
  
  # Generator to create the placed pixels on the sensitive area one by one 
  # of placed pixels 
  def generateDUTPixels():
    attributes = {}
    attributes["type"] = 0
    npixelsU = 1151
    npixelsV = 575
    pitchU = 0.25
    pitchV = 0.05
    for i in range(npixelsU):
      for j in range(npixelsV):
        attributes["u"] = i
        attributes["v"] = j
        attributes["centreu"] = pitch*i 
        attributes["centrev"] = pitch*j 
        yield attributes

  dut.generatePixels = createDUTPixels 
  telescope.append(dut)
  
  # There are two errors raised on __init__ of detector both if something with the 
  # createPixelfunctionList parameter is not correct, no list or to few list items. 
  # It is not needed to use this try-except because the exception in Detector class 
  # cancels the program then, but it is neater.
  try: 
    detector = tbsw.DetLayoutGen.CreateDetector(xmloutfile=outfile, sensors=telescope)
  except ValueError as ex:
    sys.exit(str(ex)+ " Aborting!")
  
  #########################################################################################
  # The following code should be moved into the redesigned tbsw.DetLayoutGen.CreateDetector
  # method and fully hidden from the user. 
  #########################################################################################
  
  # reset the file to empty because use append mode later
  f = open(outfile, "w") 
  f.close()
  # create the basic xml file structure, for parameters see DetLayoutGen.py
  gearBaseTags = tbsw.DetLayoutGen.createGearStructure()
  # write the opening tags of the basic xml structure. This is needed because the 
  # xml nodes are itterativly created and written to file to reduce memory consumption
  with open(outfile, "a") as f:
    f.write(gearBaseTags[0])
    f.flush()
    os.fsync(f.fileno())
  
  # create and initialise the Detector object
  detector = tbsw.DetLayoutGen.Detector(xmloutfile=outfile, sensparamsList=spList, detectorparams=dparams, uCellGroupparamsList=uCGpList, vCellGroupparamsList=vCGpList)
  # create the detector with the planes and write it to outfile
  detector.createDetectorElement()
  
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
