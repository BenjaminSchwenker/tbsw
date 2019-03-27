"""
This is an example script to plot pxd angle scan data from nov 18 test beam 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

import tbsw 
import yaml
import os
import ROOT 
import pandas as pd
import matplotlib.pyplot as plt


# Read experiment db
expdb = yaml.load(open("expdb.yaml", "r"),  Loader=yaml.FullLoader)  



# Keep some numbers to gain overview
summary_dict = {  'dut': [],
                  'shape' : [],
                  'coverage' : [],
                  'sigma_u' : [],
                  'sigma_v' : [],
                  'rho' : [],
                  'theta'     : [],
                  'phi'   : [],
                  'pixeltype' : [],
                }


for k in expdb.keys(): 
  
  if not os.path.isfile("localDB/{:s}/clusterDB-PXD.root".format(k)):
    print("Cannot find " + "localDB/{:s}/clusterDB-PXD.root".format(k) )
    continue
  
  pxdDB = tbsw.clusterDB.ClusterDB("localDB/{:s}/clusterDB-PXD.root".format(k))
  
  for pixelType in pxdDB.getPixelTypes():
    shape = '^E[0-9]+P[0-9]+.[0-9]+.{:d}'.format(pixelType) 
             
    sigU, sigUError = pxdDB.getSigmaU(shape)
    sigV, sigVError = pxdDB.getSigmaV(shape)
    rho = pxdDB.getRho(shape)
    cov = pxdDB.getCoverage() 

    
    # Fill summary dict
    summary_dict['dut'].append(expdb[k]["dut"])       
    summary_dict['phi'].append(expdb[k]["phi"])       
    summary_dict['theta'].append(expdb[k]["theta"])       
    summary_dict['pixeltype'].append(pixelType)             
    summary_dict['shape'].append(shape)       
    summary_dict['coverage'].append(cov)       
    summary_dict['sigma_u'].append(sigU*1000) # in um      
    summary_dict['sigma_v'].append(sigV*1000) # in um       
    summary_dict['rho'].append(rho )             


# Create a pandas data 
df = pd.DataFrame(summary_dict)


for dut in set(df["dut"]): 
  for pixeltype in set(df["pixeltype"]): 
    
    print ("DUT={}".format( dut)  )
    print ("PixelType={}".format( pixeltype)  )
    
    df2 = df[(df['dut'] == dut)  & (df['pixeltype'] == pixeltype) ][['theta', 'phi', 'sigma_u', 'sigma_v', 'rho']]  

    print(df2)
    
    #pivot_table = df.pivot(index='phi', columns='theta', values='sigma_u')
    #fig = plt.figure(figsize=(12,12))
    #ax = fig.add_subplot(111)
    #ax.set_xlabel('thetaU / degree', size = 20)
    #ax.set_ylabel('thetaV / degree', size = 20)
    #ax.set_title('Number of labels kind={:d}'.format(pixelkind), size = 20)
    #ax = sns.heatmap(pivot_table, mask=pivot_table.isnull(), annot=True, fmt="d", linewidths=.5, square = True,   cmap = 'Blues_r', cbar_kws={'label': '#labels'})
    #ax.invert_yaxis()
    #fig.savefig(resultdir + '/Label_Heatmap_kind_{:d}.png'.format(pixelkind), dpi=100)
    #fig.clf() 
    #plt.close(fig)  

