import re
import ROOT
import math

shape = "E0P0.0.0D0.0"

def get_labeltype(label=None): 
  return_type = ""
  for tok in re.split('D',label):
    addr = re.split('\.',tok)
    if len(addr) == 1:
     return_type += tok    
    else: 
      return_type += "D" + addr[0] + '.' + addr[1]   
    
  return return_type

class ClusterDB(object):
  """
  The ClusterDB class implements an interface to extract cluster resolutions and cluster positions from 
  a clusterDB root file produced by the GoeClusterCalibrator processor in tbsw.
   
  :author: benjamin.schwenker@phys.uni-goettinge.de  
  """
  def __init__(self, dbpath):
    """  
    :@dbpath: path to clusterDB root file
    """    
    self.dbpath = dbpath
    
    
  def getSigmaU(self): 
    """
    Get pair of sigma (and sigma error) along the u axis as weighted mean of 
    all calibrated cluster shapes in the db. 
    """
    dbfile = ROOT.TFile( self.dbpath, 'READ' ) 
    histo_weights  = dbfile.Get("hDB_Weight")
    histo_sigmas2  = dbfile.Get("hDB_Sigma2_U")
    
    mean2 = 0.0
    error2= 0.0
    norm = 0.0
    
    for bin in range(1,histo_sigmas2.GetNbinsX()+1):        
      weight = histo_weights.GetBinContent(bin)
      norm += weight   
      mean2 += weight*histo_sigmas2.GetBinContent(bin)
      error2 += weight*weight*histo_sigmas2.GetBinError(bin)
    
    dbfile.Close()
  
    if norm >0 and mean2>0:
      return math.sqrt(mean2/norm) , math.sqrt(error2/norm/norm)
    else: 
      return -1, -1; 

  def getSigmaV(self): 
    """
    Get pair of sigma (and sigma error) along the v axis as weighted mean of 
    all calibrated cluster shapes in the db. 
    """
    dbfile = ROOT.TFile( self.dbpath, 'READ' ) 
    histo_weights  = dbfile.Get("hDB_Weight")
    histo_sigmas2  = dbfile.Get("hDB_Sigma2_V")
    
    mean2 = 0.0
    error2= 0.0
    norm = 0.0
    
    for bin in range(1,histo_sigmas2.GetNbinsX()+1):        
      weight = histo_weights.GetBinContent(bin)
      norm += weight   
      mean2 += weight*histo_sigmas2.GetBinContent(bin)
      error2 += weight*weight*histo_sigmas2.GetBinError(bin)
    
    dbfile.Close()
  
    if norm >0 and mean2>0:
      return math.sqrt(mean2/norm) , math.sqrt(error2/norm/norm)
    else: 
      return -1, -1; 

    
