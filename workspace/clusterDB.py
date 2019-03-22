"""
This is an example script to demonstrate how to visualize the contents of clusterDBs for resolutions studies. 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

import tbsw 
import os
import glob
import yaml
import shutil
import ROOT


fullpath = os.getcwd() 

# Read experiment db
db = yaml.load(open("expdb.yaml", "r"),  Loader=yaml.FullLoader)  

# Create a folder for clusterDB plots 
if os.path.isdir("resolution_study"):
  shutil.rmtree("resolution_study")

os.mkdir("resolution_study")

for expName in db.keys(): 

  if not os.path.isfile("{:}/localDB/{:s}/clusterDB-PXD.root".format(fullpath, expName)):
    print("Cannot find " + "{:}/localDB/{:s}/clusterDB-PXD.root".format(fullpath, expName) )
    continue
  
  # Create folder for experiment specific plots
  os.mkdir("{:s}/resolution_study/{:s}".format(fullpath, expName))

  # Go into this folder for plotting
  os.chdir("{:s}/resolution_study/{:s}".format(fullpath, expName))

  print("Plotting {:}".format(expName))
  
  # and plot 
  tbsw.DQMplots.plot_clusterDB_parameters('{:}/localDB/{:s}/clusterDB-PXD.root'.format(fullpath, expName) )   
   
  os.chdir("{:s}".format(fullpath))
    
    
    





