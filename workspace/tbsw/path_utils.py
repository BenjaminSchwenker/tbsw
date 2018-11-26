"""
Some helper code to define processor paths

:author: ulf.stolzenberg@phys.uni-goettinge.de  
"""

from tbsw import *
 

def add_rawinputprocessor(path,rawfile):
  """
  Adds raw input processor to the path
  """ 

  rawinput_processor = Processor(name="rawinput_processor",proctype="EudaqInputProcessor")
  rawinput_processor.param("DetectorName","EUTelescope")
  rawinput_processor.param("FileNames",rawfile)
  path.add_processor(rawinput_processor)

  return path 


def add_M26unpacker(path):
  """
  Adds M26 unpacker to the path
  """ 

  M26unpacker = Processor(name="M26Unpacker",proctype="NIUnpacker")
  M26unpacker.param("InputCollectionName","NI")
  M26unpacker.param("OutputCollectionName","zsdata_m26")
  path.add_processor(M26unpacker)

  return path 


def add_M26clusterizer(path):
  """
  Adds M26 clusterizer to the path
  """  
    
  m26clust = Processor(name="M26Clusterizer",proctype="PixelClusterizer")   
  m26clust.param("NoiseDBFileName","localDB/NoiseDB-M26.root")
  m26clust.param("SparseDataCollectionName","zsdata_m26")
  m26clust.param("ClusterCollectionName","zscluster_m26")
  m26clust.param("SparseClusterCut",0)
  m26clust.param("SparseSeedCut", 0)
  m26clust.param("SparseZSCut", 0)   
  path.add_processor(m26clust)  

  return path


def add_trackfinder_loosecut(path, beamenergy):
  """
  Adds trackfinding with loose cut criteria to the path
  """  
    
  trackfinder_loosecut = Processor(name="AlignTF_LC",proctype="FastTracker")
  trackfinder_loosecut.param("InputHitCollectionNameVec","hit_m26")
  trackfinder_loosecut.param("AlignmentDBFileName","localDB/alignmentDB.root")
  trackfinder_loosecut.param("ExcludeDetector", "3")
  trackfinder_loosecut.param("MaxTrackChi2", 10000000)
  trackfinder_loosecut.param("MaximumGap", 0)
  trackfinder_loosecut.param("MinimumHits",6)
  trackfinder_loosecut.param("OutlierChi2Cut", 100000000)
  trackfinder_loosecut.param("ParticleCharge","-1")
  trackfinder_loosecut.param("ParticleMass","0.000511")
  trackfinder_loosecut.param("ParticleMomentum", beamenergy)
  trackfinder_loosecut.param("BackwardPass_FirstPlane","-1")
  trackfinder_loosecut.param("BackwardPass_SecondPlane","-1")
  trackfinder_loosecut.param("ForwardPass_FirstPlane","-1")
  trackfinder_loosecut.param("ForwardPass_SecondPlane","-1")
  trackfinder_loosecut.param("SingleHitSeeding", "0 1")
  trackfinder_loosecut.param("MaxResidualU","0.4")
  trackfinder_loosecut.param("MaxResidualV","0.4")
  path.add_processor(trackfinder_loosecut)

  return path


def add_trackletfinder_loosecut(path, beamenergy, excludeplanes, minhits):
  """
  Adds trackfinding with loose cut criteria to the path
  """  
    
  trackfinder_loosecut = Processor(name="AlignTF_LC",proctype="FastTracker")
  trackfinder_loosecut.param("InputHitCollectionNameVec","hit_m26")
  trackfinder_loosecut.param("AlignmentDBFileName","localDB/alignmentDB.root")
  trackfinder_loosecut.param("ExcludeDetector", excludeplanes)
  trackfinder_loosecut.param("MaxTrackChi2", 10000000)
  trackfinder_loosecut.param("MaximumGap", 0)
  trackfinder_loosecut.param("MinimumHits",minhits)
  trackfinder_loosecut.param("OutlierChi2Cut", 100000000)
  trackfinder_loosecut.param("ParticleCharge","-1")
  trackfinder_loosecut.param("ParticleMass","0.000511")
  trackfinder_loosecut.param("ParticleMomentum", beamenergy)
  trackfinder_loosecut.param("BackwardPass_FirstPlane","-1")
  trackfinder_loosecut.param("BackwardPass_SecondPlane","-1")
  trackfinder_loosecut.param("ForwardPass_FirstPlane","-1")
  trackfinder_loosecut.param("ForwardPass_SecondPlane","-1")
  trackfinder_loosecut.param("SingleHitSeeding", "0 1")
  trackfinder_loosecut.param("MaxResidualU","0.4")
  trackfinder_loosecut.param("MaxResidualV","0.4")
  path.add_processor(trackfinder_loosecut)

  return path


def add_prealigner(path):
  """
  Adds prealigner (no z alignment) to the path
  """  
  prealigner = Processor(name="PreAligner",proctype="KalmanAligner")
  prealigner.param("AlignmentDBFileName","localDB/alignmentDB.root")
  prealigner.param('ErrorsShiftX' , '0 10 10 0 10 10 0')
  prealigner.param('ErrorsShiftY' , '0 10 10 0 10 10 0')
  prealigner.param('ErrorsShiftZ' , '0 0 0 0 0 0 0')
  prealigner.param('ErrorsAlpha'  , '0 0 0 0 0 0 0')
  prealigner.param('ErrorsBeta'   , '0 0 0 0 0 0 0')
  prealigner.param('ErrorsGamma'  , '0 0.01 0.01 0 0.01 0.01 0')
  path.add_processor(prealigner) 

  return path 


def add_trackfinder_tightcut(path, beamenergy):
  """
  Adds trackfinding with tight cut criteria to the path
  """  

  trackfinder_tightcut = Processor(name="AlignTF_TC",proctype="FastTracker")
  trackfinder_tightcut.param("InputHitCollectionNameVec","hit_m26")
  trackfinder_tightcut.param("AlignmentDBFileName","localDB/alignmentDB.root")
  trackfinder_tightcut.param("ExcludeDetector", "3")
  trackfinder_tightcut.param("MaxTrackChi2", 50)
  trackfinder_tightcut.param("MaximumGap", 0)
  trackfinder_tightcut.param("MinimumHits",6)
  trackfinder_tightcut.param("OutlierChi2Cut", 10)
  trackfinder_tightcut.param("ParticleCharge","-1")
  trackfinder_tightcut.param("ParticleMass","0.000511")
  trackfinder_tightcut.param("ParticleMomentum", beamenergy)
  trackfinder_tightcut.param("BackwardPass_FirstPlane","-1")
  trackfinder_tightcut.param("BackwardPass_SecondPlane","-1")
  trackfinder_tightcut.param("ForwardPass_FirstPlane","-1")
  trackfinder_tightcut.param("ForwardPass_SecondPlane","-1")
  trackfinder_tightcut.param("SingleHitSeeding", "0 1")
  trackfinder_tightcut.param("MaxResidualU","0.4")
  trackfinder_tightcut.param("MaxResidualV","0.4")
  path.add_processor(trackfinder_tightcut)

  return path

def add_trackletfinder_tightcut(path, beamenergy, excludeplanes, minhits):
  """
  Adds trackfinding with tight cut criteria to the path
  """  

  trackfinder_tightcut = Processor(name="AlignTF_TC",proctype="FastTracker")
  trackfinder_tightcut.param("InputHitCollectionNameVec","hit_m26")
  trackfinder_tightcut.param("AlignmentDBFileName","localDB/alignmentDB.root")
  trackfinder_tightcut.param("ExcludeDetector", excludeplanes)
  trackfinder_tightcut.param("MaxTrackChi2", 50)
  trackfinder_tightcut.param("MaximumGap", 0)
  trackfinder_tightcut.param("MinimumHits",minhits)
  trackfinder_tightcut.param("OutlierChi2Cut", 10)
  trackfinder_tightcut.param("ParticleCharge","-1")
  trackfinder_tightcut.param("ParticleMass","0.000511")
  trackfinder_tightcut.param("ParticleMomentum", beamenergy)
  trackfinder_tightcut.param("BackwardPass_FirstPlane","-1")
  trackfinder_tightcut.param("BackwardPass_SecondPlane","-1")
  trackfinder_tightcut.param("ForwardPass_FirstPlane","-1")
  trackfinder_tightcut.param("ForwardPass_SecondPlane","-1")
  trackfinder_tightcut.param("SingleHitSeeding", "0 1")
  trackfinder_tightcut.param("MaxResidualU","0.4")
  trackfinder_tightcut.param("MaxResidualV","0.4")
  path.add_processor(trackfinder_tightcut)

  return path


def add_aligner(path):
  """
  Adds alignment to the path
  """  

  aligner = Processor(name="Aligner",proctype="KalmanAligner")
  aligner.param("AlignmentDBFileName","localDB/alignmentDB.root")
  aligner.param('ErrorsShiftX' , '0 10 10 0 10 10 0' )
  aligner.param('ErrorsShiftY' , '0 10 10 0 10 10 0')
  aligner.param('ErrorsShiftZ' , '0 10 10 0 10 10 0')
  aligner.param('ErrorsAlpha'  , '0 0 0 0 0 0 0')
  aligner.param('ErrorsBeta'   , '0 0 0 0 0 0 0')
  aligner.param('ErrorsGamma'  , '0 0.01 0.01 0 0.01 0.01 0')
  path.add_processor(aligner) 

  return path 



def add_M26hitmaker(path, hitmakertype):
  """
  Adds M26 hitmaker to the path
  """  
  if hitmakertype=="goe":
    m26goehitmaker = Processor(name="M26GoeHitMaker",proctype="GoeHitMaker")
    m26goehitmaker.param("ClusterCollection","zscluster_m26")
    m26goehitmaker.param("HitCollectionName","hit_m26")
    m26goehitmaker.param("ClusterDBFileName","localDB/clusterDB-M26.root")
    path.add_processor(m26goehitmaker)  

  else:
    m26coghitmaker = Processor(name="M26CogHitMaker",proctype="CogHitMaker")
    m26coghitmaker.param("ClusterCollection","zscluster_m26")
    m26coghitmaker.param("HitCollectionName","hit_m26")
    m26coghitmaker.param("SigmaUCorrections", "0.698 0.31 0.315")
    m26coghitmaker.param("SigmaVCorrections", "0.698 0.31 0.315")
    path.add_processor(m26coghitmaker)  
    
  return path 



def add_clustercalibrator(path):
  """
  Adds M26 cluster calibration to the path
  """ 

  cluster_calibrator = Processor(name="M26ClusterCalibrator",proctype="GoeClusterCalibrator")
  cluster_calibrator.param("ClusterDBFileName","localDB/clusterDB-M26.root")
  cluster_calibrator.param("AlignmentDBFileName","localDB/alignmentDB.root")
  cluster_calibrator.param("MinClusters", "2000")
  cluster_calibrator.param("MinVarianceU", 2e-06)
  cluster_calibrator.param("MinVarianceV", 2e-06)
  cluster_calibrator.param("IgnoreIDs",11)
  path.add_processor(cluster_calibrator)

  return path 


def append_tracklett_aligner(Env, path, gearfile, nevents_cali, beamenergy, hitmakertype, trackletttype):
  """
  Adds M26 cluster calibration to the path
  """ 

  if trackletttype=="triplet":
    xerrors="0 10 0 0 0 0 0"
    yerrors="0 10 0 0 0 0 0"
    gammaerrors="0 0.01 0 0 0 0 0"
    minhits=3
    excludeplanes="3 4 5 6"

  elif trackletttype=="quadruplet":
    xerrors="0 10 10 0 0 0 0"
    yerrors="0 10 10 0 0 0 0"
    gammaerrors="0 0.01 0.01 0 0 0 0"
    minhits=4
    excludeplanes="3 5 6"

  elif trackletttype=="quintet":
    xerrors="0 10 10 0 10 0 0"
    yerrors="0 10 10 0 10 0 0"
    gammaerrors="0 0.01 0.01 0 0.01 0 0"
    minhits=5
    excludeplanes="3 6"

  else:
    xerrors="0 10 10 0 10 10 0"
    yerrors="0 10 10 0 10 10 0"
    gammaerrors="0 0.01 0.01 0 0.01 0.01 0"
    minhits=6
    excludeplanes="3"

  # Create path for tracklett prealigner
  tracklett_prealigner_path = Env.create_path(hitmakertype+'hits_'+trackletttype+'_prealigner_path')
  tracklett_prealigner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  nevents_cali, 'LCIOInputFiles': "tmp.slcio"  }) 

  tracklett_prealigner_path=add_M26hitmaker(tracklett_prealigner_path,hitmakertype) 
  tracklett_prealigner_path=add_trackletfinder_loosecut(tracklett_prealigner_path, beamenergy, excludeplanes, minhits)

  prealigner_tracklett = Processor(name="Tracklett_PreAligner",proctype="KalmanAligner")
  prealigner_tracklett.param("AlignmentDBFileName","localDB/alignmentDB.root")
  prealigner_tracklett.param('ErrorsShiftX' , xerrors)
  prealigner_tracklett.param('ErrorsShiftY' , yerrors)
  prealigner_tracklett.param('ErrorsShiftZ' , '0 0 0 0 0 0 0')
  prealigner_tracklett.param('ErrorsAlpha'  , '0 0 0 0 0 0 0')
  prealigner_tracklett.param('ErrorsBeta'   , '0 0 0 0 0 0 0')
  prealigner_tracklett.param('ErrorsGamma'  , gammaerrors)
  tracklett_prealigner_path.add_processor(prealigner_tracklett) 

  # Finished with tracklett prealigner
  path.append(tracklett_prealigner_path) 
    

  # Create path for tracklett aligner
  tracklett_aligner_path = Env.create_path(hitmakertype+'hits_'+trackletttype+'_aligner_path')
  tracklett_aligner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  nevents_cali, 'LCIOInputFiles': "tmp.slcio"  }) 

  tracklett_aligner_path=add_M26hitmaker(tracklett_aligner_path,hitmakertype) 
  tracklett_aligner_path=add_trackletfinder_tightcut(tracklett_aligner_path, beamenergy, excludeplanes, minhits)

  aligner_tracklett = Processor(name="Tracklett_Aligner",proctype="KalmanAligner")
  aligner_tracklett.param("AlignmentDBFileName","localDB/alignmentDB.root")
  aligner_tracklett.param('ErrorsShiftX' , xerrors )
  aligner_tracklett.param('ErrorsShiftY' , yerrors )
  aligner_tracklett.param('ErrorsShiftZ' , '0 0 0 0 0 0 0')
  aligner_tracklett.param('ErrorsAlpha'  , '0 0 0 0 0 0 0')
  aligner_tracklett.param('ErrorsBeta'   , '0 0 0 0 0 0 0')
  aligner_tracklett.param('ErrorsGamma'  , gammaerrors)
  tracklett_aligner_path.add_processor(aligner_tracklett) 

  # Finished with tracklett aligner
  path.append(tracklett_aligner_path) 


  # Create path for tracklett dqm
  tracklett_dqm_path = Env.create_path(hitmakertype+'hits_'+trackletttype+'_dqm_path')
  tracklett_dqm_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  nevents_cali, 'LCIOInputFiles': "tmp.slcio"  }) 

  tracklett_dqm_path=add_M26hitmaker(tracklett_dqm_path,hitmakertype) 
  tracklett_dqm_path=add_trackletfinder_tightcut(tracklett_dqm_path, beamenergy, excludeplanes, minhits)

  tracklett_dqm = Processor(name="TelescopeDQM", proctype="TrackFitDQM") 
  tracklett_dqm.param("AlignmentDBFileName","localDB/alignmentDB.root")
  tracklett_dqm.param("RootFileName","TelescopeDQM_"+hitmakertype+'hits_'+trackletttype+".root")
  tracklett_dqm_path.add_processor(tracklett_dqm)  

  # Finished with tracklett dqm
  path.append(tracklett_dqm_path)


  # Create path for tracklett correlator
  tracklett_correlator_path = Env.create_path(hitmakertype+'hits_'+trackletttype+'_correlator_path')
  tracklett_correlator_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  nevents_cali, 'LCIOInputFiles': "tmp.slcio"  }) 

  tracklett_correlator_path=add_M26hitmaker(tracklett_correlator_path,"cog") 
  tracklett_correlator_path=add_trackletfinder_tightcut(tracklett_correlator_path, beamenergy, excludeplanes, 3)

  tracklett_correlator = Processor(name="Tracklett_Correlator", proctype="TriplettCorrelator") 
  tracklett_correlator.param("InputHitCollectionNameVec","hit_m26")
  tracklett_correlator.param("OutputRootFileName","Correlator_"+hitmakertype+'hits_'+trackletttype+".root")
  tracklett_correlator.param("TrackCollectionName","tracks")
  tracklett_correlator.param("UpdateAlignment","true")
  tracklett_correlator.param("NewAlignment","false")
  tracklett_correlator.param("AlignmentDBFileName", "localDB/alignmentDB.root")
  tracklett_correlator_path.add_processor(tracklett_correlator)

  # Finished with tracklett dqm
  path.append(tracklett_correlator_path)

  return path


def create_anglereco_path(Env, rawfile, gearfile, numberofevents, usesinglehitseeding, useclusterdb, beamenergy, mcdata):
  """
  Returns a list of tbsw path objects to reconstruct a test beam run 
  """


  reco_path = Env.create_path('reco')

  if not mcdata:
    reco_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : numberofevents }) 
    reco_path=add_rawinputprocessor(reco_path,rawfile) 
    reco_path=add_M26unpacker(reco_path) 

  else:
    reco_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : numberofevents, 'LCIOInputFiles': rawfile })

  reco_path=add_M26clusterizer(reco_path);

  if useclusterdb:
    reco_path=add_M26hitmaker(reco_path,"goe")
  else:
    reco_path=add_M26hitmaker(reco_path,"cog")

  downstream_TF=Processor(name="downstream_TF",proctype="FastTracker")
  downstream_TF.param("InputHitCollectionNameVec","hit_m26")
  downstream_TF.param("AlignmentDBFileName","localDB/alignmentDB.root")
  downstream_TF.param("HitCollectionName","unusedhits_down")
  downstream_TF.param("OutputTrackCollectionName","down_tracks")
  downstream_TF.param("ExcludeDetector", "0 1 2 3")
  downstream_TF.param("MaxTrackChi2", 20)
  downstream_TF.param("MaximumGap", 0)
  downstream_TF.param("MinimumHits", 3)
  downstream_TF.param("OutlierChi2Cut", 10)
  downstream_TF.param("ParticleCharge","-1")
  downstream_TF.param("ParticleMass","0.000511")
  downstream_TF.param("ParticleMomentum", beamenergy)
  downstream_TF.param("MaxResidualU","0.4")
  downstream_TF.param("MaxResidualV","0.4")

  if usesinglehitseeding:
    downstream_TF.param("BackwardPass_FirstPlane","-1")
    downstream_TF.param("BackwardPass_SecondPlane","-1")
    downstream_TF.param("ForwardPass_FirstPlane","-1")
    downstream_TF.param("ForwardPass_SecondPlane","-1")
    downstream_TF.param("SingleHitSeeding", "6")

  else:
    downstream_TF.param("BackwardPass_FirstPlane","5")
    downstream_TF.param("BackwardPass_SecondPlane","6")
    downstream_TF.param("ForwardPass_FirstPlane","-1")
    downstream_TF.param("ForwardPass_SecondPlane","-1")
    downstream_TF.param("SingleHitSeeding", "")
  reco_path.add_processor(downstream_TF)


  upstream_TF=Processor(name="upstream_TF",proctype="FastTracker")
  upstream_TF.param("InputHitCollectionNameVec","hit_m26")
  upstream_TF.param("AlignmentDBFileName","localDB/alignmentDB.root")
  upstream_TF.param("HitCollectionName","unusedhits_up")
  upstream_TF.param("OutputTrackCollectionName","up_tracks")
  upstream_TF.param("ExcludeDetector", "3 4 5 6")
  upstream_TF.param("MaxTrackChi2", 20)
  upstream_TF.param("MaximumGap", 0)
  upstream_TF.param("MinimumHits", 3)
  upstream_TF.param("OutlierChi2Cut", 10)
  upstream_TF.param("ParticleCharge","-1")
  upstream_TF.param("ParticleMass","0.000511")
  upstream_TF.param("ParticleMomentum", beamenergy)
  upstream_TF.param("MaxResidualU","0.4")
  upstream_TF.param("MaxResidualV","0.4")

  if usesinglehitseeding:
    upstream_TF.param("BackwardPass_FirstPlane","-1")
    upstream_TF.param("BackwardPass_SecondPlane","-1")
    upstream_TF.param("ForwardPass_FirstPlane","-1")
    upstream_TF.param("ForwardPass_SecondPlane","-1")
    upstream_TF.param("SingleHitSeeding", "0")

  else:
    upstream_TF.param("BackwardPass_FirstPlane","-1")
    upstream_TF.param("BackwardPass_SecondPlane","-1")
    upstream_TF.param("ForwardPass_FirstPlane","0")
    upstream_TF.param("ForwardPass_SecondPlane","1")
    upstream_TF.param("SingleHitSeeding", "")
  reco_path.add_processor(upstream_TF)

  x0imageproducer=Processor(name="x0imageproducer",proctype="X0ImageProducer")
  x0imageproducer.param("DownStreamTrackCollection","down_tracks")
  x0imageproducer.param("UpStreamTrackCollection","up_tracks")
  x0imageproducer.param("AlignmentDBFileName","localDB/alignmentDB.root")
  x0imageproducer.param("DUTPlane","3")
  x0imageproducer.param("VertexFitSwitch","false")
  x0imageproducer.param("ToyBetheHeitlerSwitch","true")
  x0imageproducer.param("ToyScatteringSwitch","false")
  x0imageproducer.param("ToyRecoError","-1")
  x0imageproducer.param("MaxDist","0.5")
  x0imageproducer.param("RootFileName","X0.root")
  reco_path.add_processor(x0imageproducer)
    
  return [ reco_path ]


# Processor settings and sequence during telescope calibration
def create_x0analysis_calibration_longtelescope_path(Env, rawfile, gearfile, nevents_cali, useclusterdb, beamenergy, mcdata):
  """
  Returns a list of tbsw path objects to calibrate the tracking telescope
  """
  
  # Calibrations are organized in a sequence of calibration paths. 
  # The calibration paths are collected in a list for later execution
  calpath = []

  # Create path for noisy pixel masking
  mask_path = Env.create_path('mask_path')

  if not mcdata:
    mask_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 1000000 }) 
    mask_path=add_rawinputprocessor(mask_path, rawfile) 
    mask_path=add_M26unpacker(mask_path) 

  else:
    mask_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': rawfile }) 

  m26hotpixelkiller = Processor(name="M26HotPixelKiller",proctype="HotPixelKiller",)
  m26hotpixelkiller.param("InputCollectionName", "zsdata_m26")
  m26hotpixelkiller.param("MaxOccupancy", 0.01)
  m26hotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-M26.root")
  m26hotpixelkiller.param("OfflineZSThreshold", 0)
  mask_path.add_processor(m26hotpixelkiller)

  # Finished with path for noisy pixel masking
  calpath.append(mask_path)

  # Create path for detector level creation of clusters
  clusterizer_path = Env.create_path('clusterizer_path')


  if not mcdata:
    clusterizer_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali }) 
    clusterizer_path=add_rawinputprocessor(clusterizer_path,rawfile) 
    clusterizer_path=add_M26unpacker(clusterizer_path) 

  else:
    clusterizer_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': rawfile }) 

  clusterizer_path=add_M26clusterizer(clusterizer_path)

  lciooutput = Processor(name="LCIOOutput",proctype="LCIOOutputProcessor")
  lciooutput.param("LCIOOutputFile","tmp.slcio")
  lciooutput.param("LCIOWriteMode","WRITE_NEW")
  clusterizer_path.add_processor(lciooutput)  

  # Finished with path for clusters
  calpath.append(clusterizer_path) 


  # Create path for m26 mc cluster calibration
  if mcdata:
    mc_clustercal_path = Env.create_path('mc_clustercal_path')
    mc_clustercal_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': "tmp.slcio" })  

    mc_clustercal_path=add_M26hitmaker(mc_clustercal_path,"cog") 

    m26clustercalibrationformc = Processor(name="M26ClusterCalibrationFromMC", proctype="GoeClusterCalibratorForMC")
    m26clustercalibrationformc.param("HitCollection","hit_m26")
    m26clustercalibrationformc.param("SimTrackerHitCollection","SimTrackerHits")
    m26clustercalibrationformc.param("ClusterDBFileName","localDB/clusterDB-M26-MC.root")
    m26clustercalibrationformc.param("MaxResidualU", 0.1)
    m26clustercalibrationformc.param("MaxResidualV", 0.1)
    m26clustercalibrationformc.param("MinClusters", 500) 
    mc_clustercal_path.add_processor(m26clustercalibrationformc)  

    # Finished with mc cluster calibration
    calpath.append(mc_clustercal_path)


  # Create path for correlator
  correlator_path = Env.create_path('correlator_path')
  correlator_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  nevents_cali, 'LCIOInputFiles': "tmp.slcio"  }) 

  correlator_path=add_M26hitmaker(correlator_path,"cog") 
  
  hitdqm = Processor(name="RawDQM",proctype="RawHitDQM")
  hitdqm.param("InputHitCollectionNameVec","hit_m26")  
  hitdqm.param("RootFileName","RawDQM.root")
  correlator_path.add_processor(hitdqm)  
  
  correlator = Processor(name="TelCorrelator", proctype="Correlator")
  correlator.param("InputHitCollectionNameVec","hit_m26")
  correlator.param("AlignmentDBFileName", "localDB/alignmentDB.root")
  correlator.param("OutputRootFileName","XCorrelator.root")
  correlator.param("ReferencePlane","0")
  correlator.param("ParticleCharge","-1")
  correlator.param("ParticleMass","0.000511")
  correlator.param("ParticleMomentum", beamenergy)
  correlator.param("NewAlignment","true")
  correlator.param("UpdateAlignment","true")
  correlator_path.add_processor(correlator)  

  # Finished with correlator
  calpath.append(correlator_path) 

  # Append tracklett aligners
  calpath=append_tracklett_aligner(Env, calpath, gearfile, nevents_cali, beamenergy, "cog", "triplet")
  calpath=append_tracklett_aligner(Env, calpath, gearfile, nevents_cali, beamenergy, "cog", "quadruplet")
  calpath=append_tracklett_aligner(Env, calpath, gearfile, nevents_cali, beamenergy, "cog", "quintet")

  # Create path for pre alignment with loose cut track sample 
  prealigner_path = Env.create_path('prealigner_path')
  prealigner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  
  prealigner_path=add_M26hitmaker(prealigner_path,"cog")
  prealigner_path=add_trackfinder_loosecut(prealigner_path, beamenergy)
  prealigner_path=add_prealigner(prealigner_path)
  
  # Finished with path for prealigner
  calpath.append(prealigner_path)  
  calpath.append(prealigner_path) 
  calpath.append(prealigner_path) 
  calpath.append(prealigner_path) 

  # Create path for some track based dqm using current calibrations
  dqm_loose_path = Env.create_path('dqm_loose_path')
  dqm_loose_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  

  dqm_loose_path=add_M26hitmaker(dqm_loose_path,"cog")
  dqm_loose_path=add_trackfinder_loosecut(dqm_loose_path, beamenergy)
  
  teldqm_loose = Processor(name="TelescopeDQM_loose", proctype="TrackFitDQM") 
  teldqm_loose.param("AlignmentDBFileName","localDB/alignmentDB.root")
  teldqm_loose.param("RootFileName","TelescopeDQM_loose.root")
  dqm_loose_path.add_processor(teldqm_loose) 

  calpath.append(dqm_loose_path) 

  # Create path for alignment with tight cut track sample 
  aligner_path = Env.create_path('aligner_path')
  aligner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  
  aligner_path=add_M26hitmaker(aligner_path,"cog")
  aligner_path=add_trackfinder_tightcut(aligner_path, beamenergy)
  aligner_path=add_aligner(aligner_path)   
  
  # Finished with path for aligner
  # Repeat this 3x
  calpath.append(aligner_path)  
  calpath.append(aligner_path)
  calpath.append(aligner_path)
  calpath.append(aligner_path)
  calpath.append(aligner_path)

  # Create path for some track based dqm using current calibrations
  dqm_path = Env.create_path('dqm_path')
  dqm_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  

  dqm_path=add_M26hitmaker(dqm_path,"cog")
  dqm_path=add_trackfinder_tightcut(dqm_path, beamenergy)
  
  teldqm = Processor(name="TelescopeDQM", proctype="TrackFitDQM") 
  teldqm.param("AlignmentDBFileName","localDB/alignmentDB.root")
  teldqm.param("RootFileName","TelescopeDQM.root")
  dqm_path.add_processor(teldqm) 

  # Finished with path for teldqm
  calpath.append(dqm_path)


  # Create path for first part of cluster calibration
  cluster_calibration1_path = Env.create_path('cluster_calibration_path')
  cluster_calibration1_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': "tmp.slcio" }) 

  cluster_calibration1_path=add_M26hitmaker(cluster_calibration1_path,"cog")
  cluster_calibration1_path=add_trackfinder_tightcut(cluster_calibration1_path, beamenergy)
  cluster_calibration1_path=add_clustercalibrator(cluster_calibration1_path)

  print(calpath)

  # Finished with first part of cluster calibration
  if useclusterdb:
    calpath.append(cluster_calibration1_path)

  print(calpath)


  # Create path for 6 hit pre alignment with loose cut track sample 
  prealigner2_path = Env.create_path('prealigner2_path')
  prealigner2_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" }) 

  prealigner2_path=add_M26hitmaker(prealigner2_path,"goe")
  prealigner2_path=add_trackfinder_loosecut(prealigner2_path, beamenergy)
  prealigner2_path=add_prealigner(prealigner2_path)
  
  # Finished with path for 6 hit prealigner
  if useclusterdb:
    calpath.append(prealigner2_path)  


  # Create path for second part of cluster calibration
  cluster_calibration2_path = Env.create_path('cluster_calibration2_path')
  cluster_calibration2_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': "tmp.slcio" }) 

  cluster_calibration2_path=add_M26hitmaker(cluster_calibration2_path,"goe")
  cluster_calibration2_path=add_trackfinder_tightcut(cluster_calibration2_path, beamenergy)
  cluster_calibration2_path=add_clustercalibrator(cluster_calibration2_path)
  # Finished with second part of cluster calibration
  # Repeat this 7x
  if useclusterdb:
    calpath.append(cluster_calibration2_path)
    print(calpath)
    calpath.append(cluster_calibration2_path)
    print(calpath)
    calpath.append(cluster_calibration2_path)
    calpath.append(cluster_calibration2_path)
    calpath.append(cluster_calibration2_path)
    calpath.append(cluster_calibration2_path)
    calpath.append(cluster_calibration2_path)
    print(calpath)

  # Create path for correlator
  correlator2_path = Env.create_path('correlator2_path')
  correlator2_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  nevents_cali, 'LCIOInputFiles': "tmp.slcio"  }) 

  correlator2_path=add_M26hitmaker(correlator2_path,"goe") 
  
  hitdqm2 = Processor(name="RawDQM",proctype="RawHitDQM")
  hitdqm2.param("InputHitCollectionNameVec","hit_m26")  
  hitdqm2.param("RootFileName","RawDQM2.root")
  correlator2_path.add_processor(hitdqm2)  
  
  correlator2 = Processor(name="TelCorrelator2", proctype="Correlator")
  correlator2.param("InputHitCollectionNameVec","hit_m26")
  correlator2.param("AlignmentDBFileName", "localDB/alignmentDB.root")
  correlator2.param("OutputRootFileName","XCorrelator.root")
  correlator2.param("ReferencePlane","0")
  correlator2.param("ParticleCharge","-1")
  correlator2.param("ParticleMass","0.000511")
  correlator2.param("ParticleMomentum", beamenergy)
  correlator2.param("NewAlignment","true")
  correlator2.param("UpdateAlignment","true")
  correlator2_path.add_processor(correlator2)  

  # Finished with correlator
  if useclusterdb:
    calpath.append(correlator2_path) 


  # Append tracklett aligners
  calpath=append_tracklett_aligner(Env, calpath, gearfile, nevents_cali, beamenergy, "goe", "triplet")
  calpath=append_tracklett_aligner(Env, calpath, gearfile, nevents_cali, beamenergy, "goe", "quadruplet")
  calpath=append_tracklett_aligner(Env, calpath, gearfile, nevents_cali, beamenergy, "goe", "quintet")

  # Create path for pre alignment with loose cut track sample 
  prealigner3_path = Env.create_path('prealigner3_path')
  prealigner3_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  

  prealigner3_path=add_M26hitmaker(prealigner3_path,"goe")
  prealigner3_path=add_trackfinder_loosecut(prealigner3_path, beamenergy)
  prealigner3_path=add_prealigner(prealigner3_path)
  
  # Finished with path for prealigner
  if useclusterdb:
    calpath.append(prealigner3_path)  

  # Create second path for alignment with tight cut track sample 
  aligner2_path = Env.create_path('aligner2_path')
  aligner2_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  

  aligner2_path=add_M26hitmaker(aligner2_path,"goe")
  aligner2_path=add_trackfinder_tightcut(aligner2_path, beamenergy)
  aligner2_path=add_aligner(aligner2_path)   
  
  # Finished with path for aligner
  # Repeat this 3x
  if useclusterdb:
    calpath.append(aligner2_path)  
    calpath.append(aligner2_path)
    calpath.append(aligner2_path)


  # Creeate second path for some track based dqm using current calibrations
  dqm2_path = Env.create_path('dqm2_path')
  dqm2_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  

  dqm2_path=add_M26hitmaker(dqm2_path,"goe")
  dqm2_path=add_trackfinder_tightcut(dqm2_path, beamenergy)
  
  teldqm = Processor(name="TelescopeDQM", proctype="TrackFitDQM") 
  teldqm.param("AlignmentDBFileName","localDB/alignmentDB.root")
  teldqm.param("RootFileName","TelescopeDQM2.root")
  dqm2_path.add_processor(teldqm)  
  
  # Finished with path for teldqm
  if useclusterdb:
    calpath.append(dqm2_path)

  return calpath


def create_x0analysis_calibration_path(Env, rawfile, gearfile, nevents_cali, useclusterdb, beamenergy, mcdata):
  """
  Returns a list of tbsw path objects to calibrate the tracking telescope
  """

  # Calibrations are organized in a sequence of calibration paths. 
  # The calibration paths are collected in a list for later execution
  calpath = []

  # Create path for noisy pixel masking
  mask_path = Env.create_path('mask_path')
  mask_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': rawfile }) 

  if not mcdata:
    mask_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali })
    mask_path=add_rawinputprocessor(mask_path,rawfile) 
    mask_path=add_M26unpacker(mask_path) 

  else:
    mask_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': rawfile }) 

  m26hotpixelkiller = Processor(name="M26HotPixelKiller",proctype="HotPixelKiller",)
  m26hotpixelkiller.param("InputCollectionName", "zsdata_m26")
  m26hotpixelkiller.param("MaxOccupancy", 0.01)
  m26hotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-M26.root")
  m26hotpixelkiller.param("OfflineZSThreshold", 0)
  mask_path.add_processor(m26hotpixelkiller)

  # Finished with path for noisy pixel masking
  calpath.append(mask_path)

  # Create path for detector level creation of clusters
  clusterizer_path = Env.create_path('clusterizer_path')

  if not mcdata:
    clusterizer_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali })
    clusterizer_path=add_rawinputprocessor(mask_path,rawfile) 
    clusterizer_path=add_M26unpacker(mask_path) 

  else:
    clusterizer_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': rawfile }) 

  clusterizer_path=add_M26clusterizer(clusterizer_path)

  lciooutput = Processor(name="LCIOOutput",proctype="LCIOOutputProcessor")
  lciooutput.param("LCIOOutputFile","tmp.slcio")
  lciooutput.param("LCIOWriteMode","WRITE_NEW")
  clusterizer_path.add_processor(lciooutput)  

  # Finished with path for clusters
  calpath.append(clusterizer_path) 


  # Create path for m26 mc cluster calibration
  if mcdata:
    mc_clustercal_path = Env.create_path('mc_clustercal_path')
    mc_clustercal_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': "tmp.slcio" })  

    mc_clustercal_path=add_M26hitmaker(mc_clustercal_path,"cog") 

    m26clustercalibrationformc = Processor(name="M26ClusterCalibrationFromMC", proctype="GoeClusterCalibratorForMC")
    m26clustercalibrationformc.param("HitCollection","hit_m26")
    m26clustercalibrationformc.param("SimTrackerHitCollection","SimTrackerHits")
    m26clustercalibrationformc.param("ClusterDBFileName","localDB/clusterDB-M26-MC.root")
    m26clustercalibrationformc.param("MaxResidualU", 0.1)
    m26clustercalibrationformc.param("MaxResidualV", 0.1)
    m26clustercalibrationformc.param("MinClusters", 500) 
    mc_clustercal_path.add_processor(m26clustercalibrationformc)  

    # Finished with mc cluster calibration
    calpath.append(mc_clustercal_path)


  # Create path for correlator
  correlator_path = Env.create_path('correlator_path')
  correlator_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  nevents_cali, 'LCIOInputFiles': "tmp.slcio"  }) 

  correlator_path=add_M26hitmaker(correlator_path,"cog") 
  
  hitdqm = Processor(name="RawDQM",proctype="RawHitDQM")
  hitdqm.param("InputHitCollectionNameVec","hit_m26")  
  hitdqm.param("RootFileName","RawDQM.root")
  correlator_path.add_processor(hitdqm)  
  
  correlator = Processor(name="TelCorrelator", proctype="Correlator")
  correlator.param("InputHitCollectionNameVec","hit_m26")
  correlator.param("AlignmentDBFileName", "localDB/alignmentDB.root")
  correlator.param("OutputRootFileName","XCorrelator.root")
  correlator.param("ReferencePlane","0")
  correlator.param("ParticleCharge","-1")
  correlator.param("ParticleMass","0.000511")
  correlator.param("ParticleMomentum", beamenergy)
  correlator_path.add_processor(correlator)  

  # Finished with correlator
  calpath.append(correlator_path) 


  # Create path for pre alignment with loose cut track sample 
  prealigner_path = Env.create_path('prealigner_path')
  prealigner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  
  prealigner_path=add_M26hitmaker(prealigner_path,"cog")
  prealigner_path=add_trackfinder_loosecut(prealigner_path, beamenergy)
  prealigner_path=add_prealigner(prealigner_path)
  
  # Finished with path for prealigner
  calpath.append(prealigner_path)  

  # Create path for alignment with tight cut track sample 
  aligner_path = Env.create_path('aligner_path')
  aligner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  
  aligner_path=add_M26hitmaker(aligner_path,"cog")
  aligner_path=add_trackfinder_tightcut(aligner_path, beamenergy)
  aligner_path=add_aligner(aligner_path)   
  
  # Finished with path for aligner
  # Repeat this 3x
  calpath.append(aligner_path)  
  calpath.append(aligner_path)
  calpath.append(aligner_path)

   
  # Creeate path for some track based dqm using current calibrations
  dqm_path = Env.create_path('dqm_path')
  dqm_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  

  dqm_path=add_M26hitmaker(dqm_path,"cog")
  dqm_path=add_trackfinder_tightcut(dqm_path, beamenergy)
  
  teldqm = Processor(name="TelescopeDQM", proctype="TrackFitDQM") 
  teldqm.param("AlignmentDBFileName","localDB/alignmentDB.root")
  teldqm.param("RootFileName","TelescopeDQM.root")
  dqm_path.add_processor(teldqm)  
  
  # Finished with path for teldqm
  calpath.append(dqm_path)

  # Create second path that will only be added in case the Use_clusterDB flag is enabled
  calpath2 = []

  # Create path for first part of cluster calibration
  cluster_calibration1_path = Env.create_path('cluster_calibration_path')
  cluster_calibration1_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': "tmp.slcio" }) 

  cluster_calibration1_path=add_M26hitmaker(cluster_calibration1_path,"cog")
  cluster_calibration1_path=add_trackfinder_tightcut(cluster_calibration1_path, beamenergy)
  cluster_calibration1_path=add_clustercalibrator(cluster_calibration1_path)

  # Finished with first part of cluster calibration
  calpath2.append(cluster_calibration1_path)


  # Create path for pre alignment2 with loose cut track sample 
  prealigner2_path = Env.create_path('prealigner2_path')
  prealigner2_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" }) 

  prealigner2_path=add_M26hitmaker(prealigner2_path,"goe")
  prealigner2_path=add_trackfinder_loosecut(prealigner2_path, beamenergy)
  prealigner2_path=add_prealigner(prealigner2_path)
  
  # Finished with path for second prealigner
  calpath2.append(prealigner2_path)   


  # Create path for second part of cluster calibration
  cluster_calibration2_path = Env.create_path('cluster_calibration2_path')
  cluster_calibration2_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': "tmp.slcio" }) 

  cluster_calibration2_path=add_M26hitmaker(cluster_calibration2_path,"goe")
  cluster_calibration2_path=add_trackfinder_tightcut(cluster_calibration2_path, beamenergy)
  cluster_calibration2_path=add_clustercalibrator(cluster_calibration2_path)

  # Finished with second part of cluster calibration
  # Repeat this 7x
  calpath2.append(cluster_calibration2_path)
  calpath2.append(cluster_calibration2_path)
  calpath2.append(cluster_calibration2_path)
  calpath2.append(cluster_calibration2_path)
  calpath2.append(cluster_calibration2_path)
  calpath2.append(cluster_calibration2_path)
  calpath2.append(cluster_calibration2_path)


  # Create path for second correlator
  correlator2_path = Env.create_path('correlator2_path')
  correlator2_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  nevents_cali, 'LCIOInputFiles': "tmp.slcio"  }) 

  correlator2_path=add_M26hitmaker(correlator2_path,"goe")
  
  hitdqm = Processor(name="RawDQM",proctype="RawHitDQM")
  hitdqm.param("InputHitCollectionNameVec","hit_m26")  
  hitdqm.param("RootFileName","RawDQM2.root")
  correlator2_path.add_processor(hitdqm)  
  
  correlator = Processor(name="TelCorrelator", proctype="Correlator")
  correlator.param("InputHitCollectionNameVec","hit_m26")
  correlator.param("AlignmentDBFileName", "localDB/alignmentDB.root")
  correlator.param("OutputRootFileName","XCorrelator2.root")
  correlator.param("ReferencePlane","0")
  correlator.param("ParticleCharge","-1")
  correlator.param("ParticleMass","0.000511")
  correlator.param("ParticleMomentum", beamenergy)
  correlator2_path.add_processor(correlator)  

  # Finished with correlator
  calpath2.append(correlator2_path) 


  # Create path for pre alignment with loose cut track sample 
  prealigner3_path = Env.create_path('prealigner3_path')
  prealigner3_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  

  prealigner3_path=add_M26hitmaker(prealigner3_path,"goe")
  prealigner3_path=add_trackfinder_loosecut(prealigner3_path, beamenergy)
  prealigner3_path=add_prealigner(prealigner3_path)
  
  # Finished with path for prealigner
  calpath2.append(prealigner3_path)  

  # Create second path for alignment with tight cut track sample 
  aligner2_path = Env.create_path('aligner2_path')
  aligner2_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  

  aligner2_path=add_M26hitmaker(aligner2_path,"goe")
  aligner2_path=add_trackfinder_tightcut(aligner2_path, beamenergy)
  aligner2_path=add_aligner(aligner2_path)   
  
  # Finished with path for aligner
  # Repeat this 3x
  calpath2.append(aligner2_path)  
  calpath2.append(aligner2_path)
  calpath2.append(aligner2_path)


  # Creeate second path for some track based dqm using current calibrations
  dqm2_path = Env.create_path('dqm2_path')
  dqm2_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  

  dqm2_path=add_M26hitmaker(dqm2_path,"goe")
  dqm2_path=add_trackfinder_tightcut(dqm2_path, beamenergy)
  
  teldqm = Processor(name="TelescopeDQM", proctype="TrackFitDQM") 
  teldqm.param("AlignmentDBFileName","localDB/alignmentDB.root")
  teldqm.param("RootFileName","TelescopeDQM2.root")
  dqm2_path.add_processor(teldqm)  
  
  # Finished with path for teldqm
  calpath2.append(dqm2_path)


  if useclusterdb:
    calpath.extend(calpath2)

  return calpath



def create_x0sim_path(Env, rawfile_air, rawfile_alu, gearfile_air, gearfile,  nevents_air, nevents_alu, beamenergy):
  """
  Returns a list of tbsw path objects to simulate a test beam run without any DUT
  """
  
  sim_air = Env.create_path('sim_air')
  sim_alu = Env.create_path('sim_alu')

  sim_air.set_globals(params={'GearXMLFile': gearfile_air , 'MaxRecordNumber' : nevents_air}) 
  sim_alu.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_alu})   

  infosetter = Processor(name="InfoSetter", proctype='EventInfoSetter')
  infosetter.param("RunNumber","0")
  infosetter.param("DetectorName","EUTelescope") 
  infosetter2 = Processor(name="InfoSetter", proctype='EventInfoSetter')
  infosetter2.param("DetectorName","EUTelescope")  
  infosetter2.param("RunNumber","1")  

  sim_air.add_processor(infosetter)
  sim_alu.add_processor(infosetter2)

  gun = Processor(name="ParticleGun",proctype="ParticleGunGenerator")
  gun.param("BeamIntensity","60000")
  gun.param("BeamMomentum", str(beamenergy))
  gun.param("BeamMomentumSigma", 0.001)
  gun.param("BeamTimeWindow", 0.0001)
  gun.param("BeamVertexX","0")
  gun.param("BeamVertexY","0")
  gun.param("BeamVertexZ","-10")  
  gun.param("BeamVertexXSigma","20")
  gun.param("BeamVertexYSigma","10")  
  gun.param("BeamSlopeXSigma","0.0035") 
  gun.param("CorrelationVertexXvsMomentum","0.005")  
  gun.param("CorrelationVertexXvsSlopeX","0.1") 
  gun.param("CorrelationVertexYvsMomentum","0.005") 
  gun.param("CorrelationVertexYvsSlopeY","0.1") 
  gun.param("PDG","11")
  gun.param("ParticleCharge","-1")
  gun.param("ParticleMass","0.000511")
  sim_air.add_processor(gun)
  sim_alu.add_processor(gun)

  fastsim = Processor(name="FastSim",proctype="FastSimulation")
  fastsim.param("MCParticleCollectionName","MCParticles")
  fastsim.param("SimTrackerHitCollectionName","SimTrackerHits")
  fastsim.param("AlignmentDBFileName","localDB/alignmentDB.root")
  fastsim.param("ScatterModel","0")
  fastsim.param("DoEnergyLossStraggling","false")
  fastsim.param("DoFractionalBetheHeitlerEnergyLoss","false")
  sim_air.add_processor(fastsim)
  sim_alu.add_processor(fastsim)

  m26digi = Processor(name="M26Digitizer",proctype="SiPixDigitizer")
  m26digi.param("DigitCollectionName","zsdata_m26")  
  m26digi.param("FrontEndType","1") 
  m26digi.param("ComparatorThrehold","1100")
  m26digi.param("ElectronicNoise","300")
  m26digi.param("ElectronGroupSize","100")
  m26digi.param("FilterIDs","0 1 2 3 4 5")
  m26digi.param("IntegrationWindow","true")
  m26digi.param("StartIntegration","0")
  m26digi.param("MaxSegmentLength","0.005")
  m26digi.param("SiBulkDoping","10")
  m26digi.param("TopPlaneVoltage","-5")
  m26digi.param("ZSThreshold","0")
  m26digi.param("StopIntegration","100000")
  m26digi.param("uSideBorderLength","3")
  m26digi.param("vSideBorderLength","3")
  m26digi.param("NoiseFraction","0.000001")
  sim_air.add_processor(m26digi)
  sim_alu.add_processor(m26digi)

  lciooutput = Processor(name="LCIOOutput",proctype="LCIOOutputProcessor")
  lciooutput.param("LCIOOutputFile",rawfile_air)
  lciooutput.param("LCIOWriteMode","WRITE_NEW")  
  sim_air.add_processor(lciooutput)
  
  lciooutput2 = Processor(name="LCIOOutput",proctype="LCIOOutputProcessor")
  lciooutput2.param("LCIOOutputFile",rawfile_alu)
  lciooutput2.param("LCIOWriteMode","WRITE_NEW")  
  sim_alu.add_processor(lciooutput2)
    
  simpath = [ sim_air,
              sim_alu, 
            ]

  return simpath

