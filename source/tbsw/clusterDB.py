import re
import ROOT
import math

example = 'E0P0.0.1D0.0'

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
    self.shapes = []
    
    dbfile = ROOT.TFile( self.dbpath, 'READ' ) 
    histo = dbfile.Get("hDB_Weight")
    
    for bin in range(1,histo.GetNbinsX()+1):
      shape = histo.GetXaxis().GetBinLabel(bin)
      self.shapes.append(shape) 
     
    dbfile.Close()
    
  def getNClusters(self, shape=''): 
    """
    Get number of clusters used in calibration for a given selection of cluster shapes. 
    """
    dbfile = ROOT.TFile( self.dbpath, 'READ' ) 
    histo  = dbfile.Get("hDB_Weight")
     
    nClusters = 0.0
    
    regex = re.compile(shape)       
    selected_shapes = filter(regex.match, self.shapes)
    
    for bin in range(1,histo.GetNbinsX()+1):  
      if histo.GetXaxis().GetBinLabel(bin) in selected_shapes:         
        nClusters += histo.GetBinContent(bin)   
      
    dbfile.Close()
    
    return nClusters
    
  def getFraction(self, shape=''): 
    """
    Get fraction of clusters for a given selection of cluster shapes. 
    """
    dbfile = ROOT.TFile( self.dbpath, 'READ' ) 
    histo  = dbfile.Get("hDB_Weight")
     
    fraction = 0.0
    norm = 0.0    
     
    regex = re.compile(shape)       
    selected_shapes = filter(regex.match, self.shapes)
    
    for bin in range(1,histo.GetNbinsX()+1):  
      weight = histo.GetBinContent(bin)
      norm += weight 
      if histo.GetXaxis().GetBinLabel(bin) in selected_shapes:         
        fraction += weight    
      
    dbfile.Close()
    
    if norm==0: 
      return 0
    else:  
      return 100*fraction/norm
    

  def getSigmaU(self, shape=''): 
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
    
    regex = re.compile(shape)       
    selected_shapes = filter(regex.match, self.shapes)
    
    for bin in range(1,histo_sigmas2.GetNbinsX()+1):    
      if histo_sigmas2.GetXaxis().GetBinLabel(bin) in selected_shapes:        
        weight = histo_weights.GetBinContent(bin)
        norm += weight   
        mean2 += weight*histo_sigmas2.GetBinContent(bin)
        error2 += weight*weight*histo_sigmas2.GetBinError(bin)
    
    dbfile.Close()
  
    if norm >0 and mean2>0:
      return math.sqrt(mean2/norm) , math.sqrt(error2/norm/norm)
    else: 
      return float('nan'), float('nan')

  def getSigmaV(self, shape=''): 
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
    
    regex = re.compile(shape)       
    selected_shapes = filter(regex.match, self.shapes)
     
    for bin in range(1,histo_sigmas2.GetNbinsX()+1):  
      if histo_sigmas2.GetXaxis().GetBinLabel(bin) in selected_shapes:        
        weight = histo_weights.GetBinContent(bin)
        norm += weight   
        mean2 += weight*histo_sigmas2.GetBinContent(bin)
        error2 += weight*weight*histo_sigmas2.GetBinError(bin)
    
    dbfile.Close()
    
    if norm >0 and mean2>0:
      return math.sqrt(mean2/norm) , math.sqrt(error2/norm/norm)
    else: 
      return float('nan'), float('nan') 
    
    
