""" Example python script for generating a xml geometry description (aka gearfile) of a tracking telescope for data 
    analysis with tbsw. 

    Usage: 
    python createGearfile.py  
    
    In the __main__ part, an example is given for the creation of a gear file for a telescope with 6 mimosa26 detector
    planes and one Atlas FEI4 reference plane created with the detector type SiPlanesParameters. The dut plane has 
    hexagonal pixel cells and is created with the gear type PolyPlanesParameters.
      
    author: Helge C. Beck helge-christoph.beck@phys.uni-goettingen.de
    author: Benjamin Schwenker benjamin.schwenker@phys.uni-goettingen.de
"""

import os
import sys
import tbsw


 
if __name__ == '__main__':
  
  # Name of output gearfile 
  outfile = 'gear_diamond.xml'
  
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
    m26 = tbsw.DetLayoutGen.SquareDetector( sensParams={"ID": i, "positionZ": z} )
    # Define pixel numbering and segmentation in u direction     
    m26.addUCellGroup( {"minCell": 0, "maxCell": 1151, "pitch": 0.0184} )
    # Define pixel numbering and segmentation in v direction   
    m26.addVCellGroup( {"minCell": 0, "maxCell": 575, "pitch": 0.0184}  )
    # Add sensor to telescope 
    telescope.append(m26)
    
  
  # Create Atlas FEi4 reference sensor 
  fei4 = tbsw.DetLayoutGen.SquareDetector( sensParams={"ID": 22, "positionZ": 494, "thickness": 1.0} )    
  fei4.addUCellGroup( {"minCell": 0, "maxCell": 79, "pitch": 0.25 } )
  fei4.addVCellGroup( {"minCell": 0, "maxCell": 335, "pitch": 0.05 } )
  telescope.append(fei4)
  
  # Create a costum DUT sensor with hexagonal pixel cells 
  dut = tbsw.DetLayoutGen.PolyDetector( )
  dut.sensParams( {"ID":21, "positionZ": 403, "thickness": 0.824, "radLength": 121.3068, "atomicNumber": 6, "atomicMass": 12} )
  dut.suppParams( {"radLength": 121.3068, "atomicNumber": 6, "atomicMass": 12} ) 
  
  # Define the outer shape of all pixels to be placed on the sensitive area of the DUT 
  pixelDUT = {"type": 0, "distu": 0.3, "distv": 0.01, "points": [-0.125, -0.025, -0.125, 0.025, 0.125, 0.025, 0.125, -0.025]} 
  # Add a list of the pixel outer shapes to the DUT 
  dut.addPixelShapes( [] )   
  
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

  dut.generatePixels = generateDUTPixels
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
