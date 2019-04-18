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
   
  def getClusterTypes(self):
    """
    Get list of unique cluster types found in the clusterDB.
    """ 
    clusterTypes = set( ['P'+re.split('P',shape)[1] for shape in self.shapes ] ) 
    return list(clusterTypes)
    
  def getCoverage(self):
    """
    Get coverage (in percent) efficiency in percent for clusterDB.
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
    Get number of training clusters used to train selected shapes. 
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
    Get fraction (in percent) selected shapes among all clusters seen in training. 
    All clusters includes those which could not be calibrated. 
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
      return fraction*self.getCoverage()/norm  
    

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
  
  def plotClusterType(self, clusterType, imagePath, pitchU, pitchV):
    """
    Create an image of the pixel cells of the given cluster type 
    overlaid with 68% error ellipses of position estimators. 
    """
     
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    from matplotlib.patches import Rectangle
    import numpy as np
    
    def get_cells(clusterType):
      uCells = list()
      vCells = list()
      digits = re.split('D',clusterType)[1:] 
      for digit in digits:
        uCells.append( int(re.split('\.',digit)[1]) )  
        vCells.append( int(re.split('\.',digit)[0]) )      
      return uCells, vCells
    
    def create_pixels(ucells, vcells, pitchU, pitchV):
      size = len(ucells)
      pixels = []
      for i in range(size): 
        pix = Rectangle( (-pitchU/2. + ucells[i]*pitchU, -pitchV/2. + vcells[i]*pitchV), pitchU, pitchV,)		
        pix.set_alpha(0.1)
        pix.set_facecolor('blue')
        pix.set_edgecolor('blue')
        pixels.append(pix)
      return pixels
      
    def create_error_ellipse(u, v, cov):
      vals, vecs = np.linalg.eigh(cov)
      order = vals.argsort()[::-1]
      vals = vals[order]
      vecs = vecs[:,order]
      theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
      w, h = 2 * math.sqrt(2.30) * np.sqrt(vals)
      ell = Ellipse( (u, v), width=w, height=h, angle=theta,)
      ell.set_alpha(0.4)
      ell.set_facecolor('none')
      return ell
      
    # Sum the fractions for all shape belonging to the cluster type
    # Note that we do not want to match shapes with more digits
    sum_prob = self.getFraction('^E[0-9]'+clusterType+'$') 
    
    # Compute the sigmaU, sigmaV and rho for the weighted avarage 
    # of the covariance matrices of all contributing shapes. 
    av_sigU, _ = self.getSigmaU('^E[0-9]'+clusterType+'$')
    av_sigV, _ = self.getSigmaV('^E[0-9]'+clusterType+'$')
    av_rho = self.getRho('^E[0-9]'+clusterType+'$')    
    
    # Compute list of uCells and vCells 
    uCells, vCells = get_cells(clusterType)
       
    # Prepare new figure 
    fig = plt.figure(0)
    
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_xlim(-pitchU, (np.max(uCells)+1)*pitchU)
    ax.set_ylim(-pitchV, (np.max(vCells)+1)*pitchV)
    ax.set_xlabel('offset u / mm')
    ax.set_ylabel('offset v / mm')
    ax.set_title(clusterType)
    ax.text(0.5, 0.9, 'prob={:.2f}% \n $\sigma_u$={:.1f}$\mu$m, $\sigma_v$={:.1f}$\mu$m, $\\rho$={:.2f}'.format(sum_prob, 1000*av_sigU, 1000*av_sigV, av_rho),
          style='italic',
          bbox={'facecolor':'red', 'alpha':0.5, 'pad':10},
          horizontalalignment='center',
          verticalalignment='center',
          transform = ax.transAxes)
    
    # Draw the cluster outline
    pixels = create_pixels(uCells, vCells, pitchU, pitchV)      
    for pixel in pixels:
      pixel.set_clip_box(ax.bbox)
      ax.add_artist(pixel)
      
    # Loop over all shapes   
    for shape in self.getSelectedShapes('^E[0-9]'+clusterType+'$') :
      prob = self.getFraction('^'+shape+'$')
      sigU, _ = self.getSigmaU('^'+shape+'$')
      sigV, _ = self.getSigmaV('^'+shape+'$')
      rho = self.getRho('^'+shape+'$')           
      posU = self.getPositionU('^'+shape+'$')
      posV = self.getPositionV('^'+shape+'$') 
      cov = np.array( [ [sigU*sigU, rho*sigU*sigV], [rho*sigU*sigV, sigV*sigV]] )
      
      # Draw 68% error ellipse 
      ell = create_error_ellipse(posU, posV, cov)
      ell.set_clip_box(ax.bbox)
      ell.set_color('black')
      ax.add_artist(ell)
      ax.scatter(posU, posV, color='black')
        
    fig.savefig(imagePath)
    fig.clf()


