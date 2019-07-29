""" Example python script for generating a xml geometry description (aka gearfile) of a tracking
    telescope for data  analysis with tbsw. 

    Usage: 
    python createGearfile.py  
    
    In the __main__ part, an example is given for the creation of a gear file for a telescope with 6
    mimosa26 detector planes and one Atlas FEI4 trigger reference plane. As device under test (DUT), 
    a diamond sensor with hexagonal pixel cells is placed in the center of the telescope. 
      
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
  
  # List of all sensors placed in the beam
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
  
  # Create Atlas FEI4 reference sensor 
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
  # fei4.sensParams({"rotation1": 0, "rotation2": 1, "rotation3": 1, "rotation4": 0})

  # Of course it is also possible to add additional continuous rotations around the local 
  # sensor axes. Angles are specified in degree. 
  # fei4.sensParams({"alpha": 10.1, "beta": 0, "gamma": 0})
   
  # Define the rectangular pixel matrix 
  # In this example, the u pitch changes at uCell 160 from 0.25 to 0.5mm 
  fei4.addUCellGroup( {"minCell": 0, "maxCell": 159, "pitch": 0.25 } )
  fei4.addUCellGroup( {"minCell": 160, "maxCell": 319, "pitch": 0.5 } )
  # and the v pitch changes at vCell 336 from 0.05 to 0.1mm 
  fei4.addVCellGroup( {"minCell": 0, "maxCell": 335, "pitch": 0.05 } )
  fei4.addVCellGroup( {"minCell": 336, "maxCell": 671, "pitch": 0.1 } )
  # There is no limitation for the number of cell groups along u and v axis. 
  telescope.append(fei4)
  
  # Create a custom DUT sensor with hexagonal pixel cells 
  dut = tbsw.DetLayoutGen.PolyDetector( )
  # In order to look at the documentation, outcomment the next line.
  # print(dut.__doc__) 

  # Override parameters for sensitive and support. In particular, we want 
  # to set the atomic number and mass to be diamond. 
  dut.sensParams( {"ID":21, "positionZ": 403, "thickness": 0.824, "radLength": 121.3068, "atomicNumber": 6, "atomicMass": 12} )
  dut.suppParams( {"radLength": 121.3068, "atomicNumber": 6, "atomicMass": 12} ) 
  
  # Define the outer shape of all pixels to be placed on the sensitive area of the DUT
  pixelDUTList = []
  # Rectangular pixel
  pixelDUTRect = {"type": 0, "distu": 0.3, "distv": 0.07, "points": [(-0.125, -0.025), (-0.125, 0.025), (0.125, 0.025), (0.125, -0.025)]} 
  # Hexagonal pixel
  pixelDUTHex = {"type": 1, "distu": 0.14, "distv": 0.12, "points": [(0.0, -0.1332), (-0.1138, -0.0666), (-0.1138, 0.0666), (0.0, 0.1332), (0.1138, 0.0666), (0.1138, -0.0666)]}
  pixelDUTList.append(pixelDUTRect)
  pixelDUTList.append(pixelDUTHex)

  # Add a list of the pixel outer shapes to the DUT 
  dut.addPixelShapes( pixelDUTList )   
  
  # Generator to create the placed pixels on the sensitive area one by one 
  # This example DUT consists of rectangular pixel (type == 0) and 
  # staggered hexagonal pixel (type == 1).
  # More complicated layouts with more pixel types can be created. 
  def generateDUTPixels():
    attributes = {}
    # Placing rectangular pixels type == 0
    # Note that the pixel type will be searched in the list of pixelShapes
    # added to a PolyDet instance. 
    attributes["type"] = 0
    npixelsU = 50
    npixelsVrect = 50
    pitchUrect = 0.25
    pitchVrect = 0.05
    # Placing pixels means to define the mapping of hit pixel 
    # addresses in raw data to a geometric position on the sensor. 
    for i in range(npixelsU):
      for j in range(npixelsVrect):
        # The pixels u/v adress in raw data
        attributes["u"] = i 
        attributes["v"] = j 
        # The pixels u/v center position in mm 
        attributes["centeru"] = pitchUrect*i 
        attributes["centerv"] = pitchVrect*j 
        yield attributes
   
    # Placing hexagonal pixels type == 1
    attributes["type"] = 1
    npixelsVhex = 25
    pitchUhex = 2.*0.1138
    pitchVhex = 0.4 # pitch between 2 rows of pixel
    gapRectHexV = 0.2 # distance to last row rectangular pixel
    offsetHexV = pitchVrect*npixelsVrect + gapRectHexV
    for i in range(npixelsU):
      for j in range(npixelsVhex):
        attributes["u"] = i
        attributes["v"] = j + npixelsVrect
        if j%2 == 0:
          attributes["centeru"] = pitchUhex*i
          attributes["centerv"] = pitchVhex*j/2 + offsetHexV
        else: # stagger the hexagons every second row
          attributes["centeru"] = pitchUhex*i - pitchUhex/2. 
          attributes["centerv"] = pitchVhex*(j-1)/2 + pitchVhex/2. + offsetHexV
        yield attributes

  dut.generatePixels = generateDUTPixels
  telescope.append(dut)
  
  # Try to create the gearfile needed for tbsw
  try: 
    print("Start writting gear file \"{}\" ...".format(outfile))
    tbsw.DetLayoutGen.WriteGearfile(xmloutfile=outfile, sensors=telescope)
  except IOError as (errno, strerror):
    print "I/O error({0}): {1}".format(errno, strerror)
  except ValueError:
    print "Could not convert data to an integer."
  except:
    print "Unexpected error:", sys.exc_info()[0]
    raise
   
  print("Successfully created gear file \"{}\".".format(outfile))
  
  # Try to create the rootfile containing 2d histograms for visualizing the pixel
  # matrix of sensors. We will visualize one m26 sensor, the fei4 sensor and 
  # the dut sensor.  
  # PolyDetectors are visualized as TH2Poly histograms. This can take a long time if 
  # the sensor has very many pixel cells. The maximum number of pixels that will be 
  # created is controlled by parameter maxPixel. Default value is maxPixel=30000   
  # but maxPixel=-1 creates all pixel. 
  # For PolyDetectors, the layout can be tested by filling a one in each histogram 
  # bin. The test can take long for very many pixels and can be disabled with parameter
  # fillTest=False.   
  # SquareDet sensors are visualised as TH2F histograms with different bins sizes. 
  # Lines are drawn at every bin edge. 
  rootfilename = outfile[0:outfile.find(".xml")] + ".root"
  sensorVisualise = [fei4, dut, telescope[0]]
  
  try:
    print("Start writting layout root file \"{}\" ...".format(rootfilename)) 
    tbsw.DetLayoutGen.WriteLayoutRootfile(outfile=rootfilename, sensors=sensorVisualise, fillTest=True, maxPixel=30000)
  except IOError as (errno, strerror):
    print "I/O error({0}): {1}".format(errno, strerror)
  except ValueError:
    print "Could not convert data to an integer."
  except:
    print "Unexpected error:", sys.exc_info()[0]
    raise
  
  print("Successfully created layout root file \"{}\".".format(rootfilename))
  
