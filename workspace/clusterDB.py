"""
This is an example script to demonstrate how to visualize the contents of clusterDBs for resolutions studies. 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

import tbsw 
import yaml
import os


# Read experiment db
expdb = yaml.load(open("expdb.yaml", "r"),  Loader=yaml.FullLoader)  

for expName in expdb.keys(): 
  
  if not os.path.isfile("localDB/{:s}/clusterDB-PXD.root".format(expName)):
    print("Cannot find " + "localDB/{:s}/clusterDB-PXD.root".format(expName) )
    continue
  
  print("Analyzing clusterDB {:} with phi={}, theta={}:".format(expName, expdb[expName]["phi"], expdb[expName]["theta"]))
  
  pxdDB = tbsw.clusterDB.ClusterDB("localDB/{:s}/clusterDB-PXD.root".format(expName))
  
  for pixelType in range(2):
    
    # You can querry the ClusterDB to get estimates for the spatial 
    # resolution (sigma's) for a specific cluster shape. 
    # You can use regular expression for shapes to get the sigmas  
    # on different levels of granularity. 
    #
    # 0) Example: '^E{:d}P{:d}.{:d}.{:d}D0.0$'.format(etaBin, vPeriod, uPeriod, pixelType)
    #    
    # You can specify the etaBin, vPeriod and uPeriod and pixelType to distinguish clusters
    # from different regions of the sensor. You also request cluster with exactly one firing pixel. 
    #
    # 1) Example: '^E[0-9]+P[0-9]+.[0-9]+.{:d}'.format(pixelType) 
    #    
    # You can specify the pixelType to distinguish between different a region with a specific  
    # pixel pitch. All values for etaBin, vPeriod, uPeriod and all arrangements of firing
    # pixels are matched.   
    # 
    # 2) Example: '^E[0-9]+P[0-9]+.[0-9]+.1D0.0D0.1$'
    #
    # You can select the pixel type (=1) and the detailed arrangement of 
    # fired pixels in the cluster (exactly two pixels in the same column).
    
    # Select a shape pattern
    shape = '^E[0-9]+P[0-9]+.[0-9]+.{:d}'.format(pixelType) 
           
    sigU, sigUError = pxdDB.getSigmaU(shape)
    sigV, sigVError = pxdDB.getSigmaV(shape)
    frac = pxdDB.getFraction(shape) 
    
    print("  PixelType={:d}".format(pixelType))
    print("  Fraction={:.3f}%".format(frac))
    print("  SigmaU={:.4f}+/- {:.5f} mm".format(sigU, sigUError))
    print("  SigmaV={:.4f}+/- {:.5f} mm".format(sigV, sigVError))
    
    
    





