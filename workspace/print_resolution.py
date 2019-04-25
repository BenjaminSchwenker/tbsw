"""
This is an example script prints the extracted cluster resolutions (sigmas) from clusterDBs in the localDB folder.

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""


if __name__ == '__main__':
  import tbsw 
  import glob
  import os
  import argparse
  
  parser = argparse.ArgumentParser(description="Perform calibration and reconstruction of a test beam run")
  parser.add_argument('--dbpath', dest='dbpath', default='', type=str, help='Location of clusterDB file')
  parser.add_argument('--output', dest='output', default='', type=str, help='Location of outputs')
  parser.add_argument('--size', dest='size', default='', type=str, help='Cluster size string')
  args = parser.parse_args()
  
  # Make folder for output plots 
  os.mkdir(args.output)
   
  # Create a db object 
  sensorDB = tbsw.clusterDB.ClusterDB(args.dbpath)
    
  print("Analyzing clusterDB {:} with thetaU={:.1f} degree, thetaV={:.1f} degree:".format(args.dbpath, sensorDB.getThetaU(), sensorDB.getThetaV() ))
  
  # The coverage efficiency gives the fraction of training clusters that could be succefully 
  # calibrated. Only for these clusters we can make statements for sigma's etc. 
  
  coverage = sensorDB.getCoverage()
  telsigmaU = sensorDB.getTelescopeSigmaU()
  telsigmaV = sensorDB.getTelescopeSigmaV()
  telrho = sensorDB.getTelescopeRho() 
  
  print("ClusterDB coverage={:.1f}%".format(coverage))
  print("Telescope sigmaU={:.5f}mm".format(telsigmaU))
  print("Telescope sigmaV={:.5f}mm".format(telsigmaV))
  print("Telescope rho={:.5f}".format(telrho))
  
  for clusterType in sensorDB.getClusterTypes():
    sensorDB.plotClusterType(clusterType=clusterType, imagePath="{}/{}.png".format(args.output, clusterType))

  for pixelType in sensorDB.getPixelTypes():
    for uPeriod in sensorDB.getPeriodsU(): 
      for vPeriod in sensorDB.getPeriodsV(): 
        
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
        if not args.size in ["", "D0.0D1.0$" , "D0.0$",  "D0.0D0.1$"]: 
          print("You selected a non default cluster size string. Be sure you know what you are doing!!") 
        
        # Select a shape pattern
        shape = '^E[0-9]+P[0-9]+.[0-9]+.{:d}{:s}'.format(pixelType, args.size) 
            
        sigU, sigUError = sensorDB.getSigmaU(shape)
        sigV, sigVError = sensorDB.getSigmaV(shape)
        frac = sensorDB.getFraction(shape) 
        
        print("  PixelType={:d}".format(pixelType))
        print("  Fraction={:.3f}%".format(frac))
        print("  PositionU={:.4f} mm".format(sensorDB.getPositionU(shape)))
        print("  PositionV={:.4f} mm".format(sensorDB.getPositionV(shape)))
        print("  SigmaU={:.4f}+/- {:.5f} mm".format(sigU, sigUError))
        print("  SigmaV={:.4f}+/- {:.5f} mm".format(sigV, sigVError))
        print("  Rho={:.4f}".format(sensorDB.getRho(shape)))
        
