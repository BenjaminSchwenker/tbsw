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
  
  sigU, sigUError = pxdDB.getSigmaU()
  sigV, sigVError = pxdDB.getSigmaV()
  
  print("  SigmaU={:.4f}+/- {:.5f}mm".format(sigU, sigUError))
  print("  SigmaU={:.4f}+/- {:.5f}mm".format(sigV, sigVError))
    
    
    





