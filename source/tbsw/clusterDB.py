import re
import ROOT
import math


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
    self.pixelTypes = []
    self.uPeriods = []
    self.vPeriods = []  
    
    dbfile = ROOT.TFile( self.dbpath, 'READ' ) 
    
    self.coverage = dbfile.Get("hDB_Coverage").GetBinContent(1)
    self.thetaU = math.atan(dbfile.Get("DB_angles")[0])*180/math.pi
    self.thetaV = math.atan(dbfile.Get("DB_angles")[1])*180/math.pi 

    if dbfile.Get("DB_telcov"): 
      self.telsigmaU = math.sqrt(dbfile.Get("DB_telcov")[0])
      self.telsigmaV = math.sqrt(dbfile.Get("DB_telcov")[1]) 
      self.telrho =  dbfile.Get("DB_telcov")[2]/self.telsigmaU/self.telsigmaV  
    else: 
      self.telsigmaU = float('nan')
      self.telsigmaV = float('nan')
      self.telrho = float('nan') 

    histo = dbfile.Get("hDB_Weight")
    for bin in range(1,histo.GetNbinsX()+1):
      shape = histo.GetXaxis().GetBinLabel(bin)
      self.shapes.append(shape) 
      
      # Parse the shape to obtain pixelType and periods 
      header = shape.split('D')[0].split('P')[1].split('.')
      self.pixelTypes.append(int(header[2]))
      self.uPeriods.append(int(header[1]))
      self.vPeriods.append(int(header[0]))
    
    self.pixelTypes = list(set(self.pixelTypes))  
    self.uPeriods = list(set(self.uPeriods))  
    self.vPeriods = list(set(self.vPeriods))        
    
    dbfile.Close()
    
  def getCoverage(self):
    """
    Get coverage efficiency in percent for clusterDB.
    Coverage efficiency is defined as fraction of clusters from training that could be calibrated. 
    """
    return self.coverage  

  def getThetaU(self):
    """
    Get mean track incidence angle thetaU into sensor. 
    """
    return self.thetaU 

  def getThetaV(self):
    """
    Get mean track incidence angle thetaV into sensor. 
    """
    return self.thetaV
    
  def getPixelTypes(self): 
    """
    Get list of all pixel types found in the clusterDB
    """
    return self.pixelTypes
  
  def getPeriodsU(self): 
    """
    Get list of u periods found in the clusterDB
    """
    return self.uPeriods
  
  def getPeriodsV(self): 
    """
    Get list of u periods found in the clusterDB
    """
    return self.vPeriods
 
  def getTelescopeSigmaU(self):
    """
    Get get the unbiased mean telelscope sigmaU. 
    """  
    return self.telsigmaU

  def getTelescopeSigmaV(self):
    """
    Get get the unbiased mean telelscope sigmaU. 
    """  
    return self.telsigmaV

  def getTelescopeRho(self):
    """
    Get get the unbiased mean telelscope rho. 
    """  
    return self.telrho
  
  def getSelectedShapes(self, shape=''): 
    """
    Get get list of selected shapes 
    """
    regex = re.compile(shape)       
    return filter(regex.match, self.shapes)

  def getNClusters(self, shape=''): 
    """
    Get number of training clusters for selected shapes. 
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
    Get fraction selected shapes among all shapes. 
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
    Get cluster sigma (and sigma error) along the u axis as weighted mean 
    over all selected shapes. 
    """
    dbfile = ROOT.TFile( self.dbpath, 'READ' ) 
    histo_weights  = dbfile.Get("hDB_Weight")
    histo_sigmas2  = dbfile.Get("hDB_Sigma2_U")
    
    mean2 = 0.0
    error2= 0.0
    norm = 0.0
    
    regex = re.compile(shape)       
    selected_shapes = filter(regex.match, self.shapes)
    
    for bin in range(1,histo_weights.GetNbinsX()+1):    
      if histo_weights.GetXaxis().GetBinLabel(bin) in selected_shapes:        
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
    Get cluster sigma (and sigma error) along the v axis as weighted mean  
    over all selected shapes.
    """
    dbfile = ROOT.TFile( self.dbpath, 'READ' ) 
    histo_weights  = dbfile.Get("hDB_Weight")
    histo_sigmas2  = dbfile.Get("hDB_Sigma2_V")
    
    mean2 = 0.0
    error2= 0.0
    norm = 0.0
    
    regex = re.compile(shape)       
    selected_shapes = filter(regex.match, self.shapes)
     
    for bin in range(1,histo_weights.GetNbinsX()+1):  
      if histo_weights.GetXaxis().GetBinLabel(bin) in selected_shapes:        
        weight = histo_weights.GetBinContent(bin)
        norm += weight   
        mean2 += weight*histo_sigmas2.GetBinContent(bin)
        error2 += weight*weight*histo_sigmas2.GetBinError(bin)
    
    dbfile.Close()
    
    if norm >0 and mean2>0:
      return math.sqrt(mean2/norm) , math.sqrt(error2/norm/norm)
    else: 
      return float('nan'), float('nan') 
  
  def getRho(self, shape=''): 
    """
    Get uv correlation coefficient as weighted mean over  
    all selected shapes. 
    """
    dbfile = ROOT.TFile( self.dbpath, 'READ' ) 
    histo_weights  = dbfile.Get("hDB_Weight")
    histo_cov  = dbfile.Get("hDB_Cov_UV")
    
    mean_cov = 0.0
    norm = 0.0
    
    regex = re.compile(shape)       
    selected_shapes = filter(regex.match, self.shapes)
     
    for bin in range(1,histo_weights.GetNbinsX()+1):  
      if histo_weights.GetXaxis().GetBinLabel(bin) in selected_shapes:        
        weight = histo_weights.GetBinContent(bin)
        norm += weight   
        mean_cov += weight*histo_cov.GetBinContent(bin)
      
    dbfile.Close()
    
    sigmaU, _ = self.getSigmaU(shape)
    sigmaV, _ = self.getSigmaV(shape)
    
    if norm >0 and sigmaU>0 and sigmaV>0 :
      return mean_cov/sigmaU/sigmaV/norm
    else: 
      return float('nan') 
  
  def getPositionU(self, shape=''): 
    """
    Get pair of position (and position error) along the u axis as weighted mean of 
    all calibrated cluster shapes in the db. 
    """
    dbfile = ROOT.TFile( self.dbpath, 'READ' ) 
    histo_weights  = dbfile.Get("hDB_Weight")
    histo_position  = dbfile.Get("hDB_U")
    
    mean = 0.0
    norm = 0.0
    
    regex = re.compile(shape)       
    selected_shapes = filter(regex.match, self.shapes)
    
    for bin in range(1,histo_weights.GetNbinsX()+1):    
      if histo_weights.GetXaxis().GetBinLabel(bin) in selected_shapes:        
        weight = histo_weights.GetBinContent(bin)
        norm += weight   
        mean += weight*histo_position.GetBinContent(bin)
        
    dbfile.Close()
  
    if norm >0:
      return mean/norm 
    else: 
      return float('nan')

  def getPositionV(self, shape=''): 
    """
    Get pair of position (and position error) along the v axis as weighted mean of 
    all calibrated cluster shapes in the db. 
    """
    dbfile = ROOT.TFile( self.dbpath, 'READ' ) 
    histo_weights  = dbfile.Get("hDB_Weight")
    histo_position  = dbfile.Get("hDB_V")
    
    mean = 0.0
    norm = 0.0
    
    regex = re.compile(shape)       
    selected_shapes = filter(regex.match, self.shapes)
    
    for bin in range(1,histo_weights.GetNbinsX()+1):    
      if histo_weights.GetXaxis().GetBinLabel(bin) in selected_shapes:        
        weight = histo_weights.GetBinContent(bin)
        norm += weight   
        mean += weight*histo_position.GetBinContent(bin)
        
    dbfile.Close()
  
    if norm >0:
      return mean/norm 
    else: 
      return float('nan')


