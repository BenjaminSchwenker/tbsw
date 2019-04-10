"""
This is an example script to plot pxd angle scan data from nov 18 test beam 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

import tbsw 
import yaml
import os
import ROOT 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
  
  import argparse

  parser = argparse.ArgumentParser(description="Plot cluster sigmas from clusterDB")
  parser.add_argument('--dut', dest='dut', default='W40_IF', type=str, help='DUT name')
  parser.add_argument('--size', dest='size', default='', type=str, help='Cluster size string')
  parser.add_argument('--pixeltype', dest='pixeltype', default=0, type=int, help='PixelType')
  args = parser.parse_args()
  
  if not args.size in ["", "D0.0D1.0$" , "D0.0$",  "D0.0D0.1$"]: 
    print("You selected a non default cluster size string. Be sure you know what you are doing!!") 
  
  # Read experiment db
  expdb = yaml.load(open("expdb.yaml", "r"),  Loader=yaml.FullLoader)  

  # Keep some numbers to gain overview
  summary_dict = {  'dut': [],
                    'exp': [],
                    'shape' : [],
                    'coverage' : [],
                    'sigma_u' : [],
                    'sigma_v' : [],
                    'rho' : [],
                    'theta'     : [],
                    'phi'   : [],
                    'pixeltype' : [],
                    'geoid' : [],
                    'energy': [],
                  }


  for k in expdb.keys(): 
  
    if not os.path.isfile("localDB/{:s}/clusterDB-PXD.root".format(k)):
      print("Cannot find " + "localDB/{:s}/clusterDB-PXD.root".format(k) )
      continue
      
    pxdDB = tbsw.clusterDB.ClusterDB("localDB/{:s}/clusterDB-PXD.root".format(k))
  
    for pixelType in pxdDB.getPixelTypes():
      
      # Build complete cluster shape string
      shape = '^E[0-9]+P[0-9]+.[0-9]+.{:d}{:s}'.format(pixelType, args.size) 
      
      #print("full shape={}".format(shape))
      
      sigU, sigUError = pxdDB.getSigmaU(shape)
      sigV, sigVError = pxdDB.getSigmaV(shape)
      rho = pxdDB.getRho(shape)
      cov = pxdDB.getCoverage() 
      
      
      # Fill summary dict
      summary_dict['dut'].append(expdb[k]["dut"])  
      summary_dict['exp'].append(k)            
      summary_dict['phi'].append(expdb[k]["phi"])       
      summary_dict['theta'].append(expdb[k]["theta"])       
      summary_dict['pixeltype'].append(pixelType)             
      summary_dict['shape'].append(shape)       
      summary_dict['coverage'].append(cov)       
      summary_dict['sigma_u'].append(sigU*1000) # in um      
      summary_dict['sigma_v'].append(sigV*1000) # in um       
      summary_dict['rho'].append(rho )        
      summary_dict['geoid'].append( expdb[k]["geoid"] )   
      summary_dict['energy'].append( expdb[k]["energy"] )                 
      

  # Create a pandas data 
  df = pd.DataFrame(summary_dict)
  
  print ("DUT={}".format( args.dut)  )
  print ("PixelType={}".format( args.pixeltype)  )
    
  df2 = df[(df['dut'] == args.dut)  & (df['pixeltype'] == args.pixeltype) ][['exp', 'energy', 'geoid','theta', 'phi', 'sigma_u', 'sigma_v', 'rho']]  

  print(df2)

  x = range(df2['exp'].shape[0])
    
  plt.plot(x, df2['sigma_u'], 'ro',  label='U')
  plt.title("DUT={} PixelType={} Size={}".format(args.dut,args.pixeltype, args.size))
  # You can specify a rotation for the tick labels in degrees or with keywords.
  plt.xticks(x, df2['exp'], rotation='vertical')
  plt.ylabel("cluster sigma /microns")
  # Pad margins so that markers don't get clipped by the axes
  plt.margins(0.2)
  # Tweak spacing to prevent clipping of tick-labels
  plt.subplots_adjust(bottom=0.15)
  
  plt.legend(loc='upper center', shadow=True, fontsize='x-large')

  #plt.show()
    
  plt.savefig("SigmaU_{}_pixeltype_{}_size_{}.png".format(args.dut,args.pixeltype, args.size))

  plt.clf()  
    
  plt.plot(x, df2['sigma_v'], 'ro',  label='V')
  plt.title("DUT={} PixelType={} Size={}".format(args.dut,args.pixeltype, args.size))
  # You can specify a rotation for the tick labels in degrees or with keywords.
  plt.xticks(x, df2['exp'], rotation='vertical')
  plt.ylabel("cluster sigma /microns")
  # Pad margins so that markers don't get clipped by the axes
  plt.margins(0.2)
  # Tweak spacing to prevent clipping of tick-labels
  plt.subplots_adjust(bottom=0.15)
  
  plt.legend(loc='upper center', shadow=True, fontsize='x-large')

  #plt.show()
    
  plt.savefig("SigmaV_{}_pixeltype_{}_size_{}.png".format(args.dut,args.pixeltype, args.size))

  plt.clf()    


