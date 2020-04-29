"""
This is an example script prints the extracted cluster resolutions (sigmas) from clusterDBs in the localDB folder.


Usage: 

python print_resolution.py --dbpath=<path-to-localDB/clusterDB-{detname}.root --output=<path-to-plots>


Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

def str2bool(v):
  if v.lower() in ('yes', 'true', 'on', 't', 'y', '1'):
    return True
  elif v.lower() in ('no', 'false', 'off', 'f', 'n', '0'):
    return False
  else:
    raise argparse.ArgumentTypeError('Boolean value expected')

if __name__ == '__main__':
  import tbsw 
  import glob
  import os
  import argparse
  
  parser = argparse.ArgumentParser(description="Perform calibration and reconstruction of a test beam run")
  parser.add_argument('--dbpath', dest='dbpath', default='', type=str, help='Location of clusterDB file')
  parser.add_argument('--output', dest='output', default='', type=str, help='Location of outputs')
  parser.add_argument('--size', dest='size', default='', type=str, help='Cluster size string')
  parser.add_argument('--clustype', dest='cltype', default=True, type=str2bool, help='Cluster descriptor type poly True (default) or False')
  parser.add_argument('--shape', dest='shape', default=None, type=str, help='Print results for query shape')
  args = parser.parse_args()
  
  # Make folder for output plots 
  os.mkdir(args.output)
   
  # Create a db object 
  sensorDB = tbsw.clusterDB.ClusterDB(args.dbpath)
    
  print("Analyzing clusterDB {:} with thetaU={:.3f} degree, thetaV={:.3} degree.".format(args.dbpath, sensorDB.getThetaU(), sensorDB.getThetaV() ))
  
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
  
  # You can querry the ClusterDB to get estimates for the spatial 
  # resolution (sigma's) for a specific cluster shape. 
  # You can use regular expression for shapes to get the sigmas  
  # on different levels of granularity. 
  #
  # 0) Example: '^E{:d}P{:d}.{:d}.{:d}D0.0.{:d}$'.format(etaBin, vPeriod, uPeriod, pixelType, pixelType)
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
  # The clustersize strings are a sorted sequence of firing pixel cells seperated by the literal 'D'. 
  # Pixel cells are sorted with increasing pixel cell center coordinates along the v and u coordintes, 
  # i.e. starting with the lower left cell. The firing pixel cells are written in the format
  # 'D{centerV}.{centerU}.{pixeltype}'. A call to sensorDB.getClusterTypes() produces a list of all 
  # clusterType strings.  
  
  # Plot hit estimators for all available cluster types
  for typeID, clusterType in enumerate(sensorDB.getClusterTypes()):
    sensorDB.plotClusterType(clusterType=clusterType, imagePath="{}/typeID_{:d}.png".format(args.output, typeID), poly=args.cltype)
  
  # Print hit estimator for selected cluster type
  for pixelType in sensorDB.getPixelTypes(): 
    for uPeriod in sensorDB.getPeriodsU(): 
      for vPeriod in sensorDB.getPeriodsV(): 
            
        # Select a shape pattern
        shape = '^E[0-9]+P{:d}.{:d}.{:d}{:s}'.format(vPeriod, uPeriod, pixelType, args.size) 
            
        sigU, sigUError = sensorDB.getSigmaU(shape)
        sigV, sigVError = sensorDB.getSigmaV(shape)
        frac = sensorDB.getFraction(shape) 
        
        print("  ClusterType={:s}".format(shape))
        print("  PixelType={:d}".format(pixelType))
        print("  Fraction={:.3f}%".format(frac))
        print("  SigmaU={:.5f}+/- {:.5f} mm".format(sigU, sigUError))
        print("  SigmaV={:.5f}+/- {:.5f} mm".format(sigV, sigVError))
        print("  Rho={:.4f}".format(sensorDB.getRho(shape)))



  if not args.shape==None: 
    sigU, sigUError = sensorDB.getSigmaU(args.shape)
    sigV, sigVError = sensorDB.getSigmaV(args.shape)
    frac = sensorDB.getFraction(args.shape) 
        
    print("  ClusterType={:s}".format(args.shape))
    print("  Fraction={:.3f}%".format(frac))
    print("  PositionU={:.4f} mm".format(sensorDB.getPositionU(args.shape)))
    print("  PositionV={:.4f} mm".format(sensorDB.getPositionV(args.shape)))
    print("  SigmaU={:.5f}+/- {:.5f} mm".format(sigU, sigUError))
    print("  SigmaV={:.5f}+/- {:.5f} mm".format(sigV, sigVError))
    print("  Rho={:.4f}".format(sensorDB.getRho(args.shape)))


        
