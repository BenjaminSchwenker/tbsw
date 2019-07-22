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
  outfile = 'gear_telescope.xml'
  
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
    
    # Uncomment the following lines to see the documentation and 
    # a full listing of default values for the SquareDetector class. 
    # print(m26.__doc__) 
    # print(m26.sensitiveParams)
    # print(m26.supportParams)
    
    # Define the segmentation and numbering of pixel cells pixel in u direction     
    # This creates cells numbered 0, 1151 along u axis with pitch 18.4 micron
    m26.addUCellGroup( {"minCell": 0, "maxCell": 1151, "pitch": 0.0184} )
    # Define pixel numbering and segmentation in v direction   
    m26.addVCellGroup( {"minCell": 0, "maxCell": 575, "pitch": 0.0184}  )
    # Add sensor to telescope 
    telescope.append(m26)

  # Create Atlas FEi4 reference sensor 
  fei4 = tbsw.DetLayoutGen.SquareDetector( sensParams={"ID": 22, "positionZ": 494, "thickness": 1.0} )    
  # In case the local u axis is parallel to the global y axis,
  # we can define a disrcete rotation matrix D 
  #  
  # D(1,1) =  rotation1  ,  D(1,2) =  rotation2  , D(1,3) =  0
  # D(2,1) =  rotation3  ,  D(2,2) =  rotation4  , D(2,3) =  0
  # D(3,1) =  0          ,  D(3,2) =  0          , D(3,3) =  sign
  #   
  # The rotation parameters are either +/-1 or 0. The component D(3,3) is determined by 
  # det(D)=1 condition. The matrix D reflects the different possibilities to 
  # install a detector in the beam.
  fei4.sensParams({"rotation1": 0, "rotation2": 1, "rotation3": 1, "rotation4": 0})

  # Of course it is also possible to add additional contiouus rotaions around the local 
  # sensor axes. Angles are specified in degree. 
  fei4.sensParams({"alpha": 10.1, "beta": 0, "gamma": 0})
   
  # Define the rectangular pixel matrix 
  fei4.addUCellGroup( {"minCell": 0, "maxCell": 79, "pitch": 0.25 } )
  fei4.addVCellGroup( {"minCell": 0, "maxCell": 335, "pitch": 0.05 } )
  # If the pixel pitch increases along the v axis, we simply can add more 
  # cellGroups: 
  fei4.addVCellGroup( {"minCell": 336, "maxCell": 435, "pitch": 0.15 } ) 
  telescope.append(fei4)
  
  # Create a costum DUT sensor with hexagonal pixel cells 
  dut = tbsw.DetLayoutGen.PolyDetector( )
  dut.sensParams( {"ID":21, "positionZ": 403, "thickness": 0.824, "radLength": 121.3068, "atomicNumber": 6, "atomicMass": 12} )
  dut.suppParams( {"radLength": 121.3068, "atomicNumber": 6, "atomicMass": 12} ) 
  
  # Define the outer shape of all pixels to be placed on the sensitive area of the DUT 
  pixelDUT = {"type": 0, "distu": 0.3, "distv": 0.01, "points": [-0.125, -0.025, -0.125, 0.025, 0.125, 0.025, 0.125, -0.025]} 
  # Add a list of the pixel outer shapes to the DUT 
  dut.addPixelShapes( [pixelDUT] )   
  
  # Generator to create the placed pixels on the sensitive area one by one 
  # of placed pixels 
  # TODO put a slightly more complicated examplep here with some hexagonal pixels. 
  # BUT the example should not be get much longer to keep this example self contained. 
  def generateDUTPixels():
    attributes = {}
    attributes["type"] = 0
    npixelsU = 50
    npixelsV = 50
    pitchU = 0.25
    pitchV = 0.05
    for i in range(npixelsU):
      for j in range(npixelsV):
        attributes["u"] = i
        attributes["v"] = j
        attributes["centreu"] = pitchU*i 
        attributes["centrev"] = pitchV*j 
        yield attributes
   
  dut.generatePixels = generateDUTPixels
  telescope.append(dut)
  
  
  try: 
    tbsw.DetLayoutGen.WriteGearfile(xmloutfile=outfile, sensors=telescope)
  except IOError as (errno, strerror):
    print "I/O error({0}): {1}".format(errno, strerror)
  except ValueError:
    print "Could not convert data to an integer."
  except:
    print "Unexpected error:", sys.exc_info()[0]
    raise

  print("Successfully created gear file \"{}\".".format(outfile))
  

  
