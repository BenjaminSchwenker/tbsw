"""
Some helper code to define processor paths for X/X0 studies

:author: benjamin.schwenker@phys.uni-goettinge.de  
:author: ulf.stolzenberg@phys.uni-goettinge.de  
"""

import tbsw 

def add_geometry(path, applyAlignment="true", overrideAlignment="true"): 
  """
  Adds geometry processor to the path
  """   
  
  geo = tbsw.Processor(name="Geo",proctype="Geometry")
  geo.param("AlignmentDBFilePath", "localDB/alignmentDB.root")
  geo.param("ApplyAlignment", applyAlignment)
  geo.param("OverrideAlignment", overrideAlignment)
  path.add_processor(geo)
   
  return path

def add_rawinputprocessor(path,rawfile):
  """
  Adds raw input processor to the path
  """ 

  rawinput_processor = tbsw.Processor(name="rawinput_processor",proctype="EudaqInputProcessor")
  rawinput_processor.param("DetectorName","EUTelescope")
  rawinput_processor.param("FileNames",rawfile)
  path.add_processor(rawinput_processor)

  return path 

def add_csvrawinputprocessor(path,rawfile, runno=0):
  """
  Adds raw input processor to the path
  """ 

  csvinput = tbsw.Processor(name="CSVInputProcessor",proctype="AsciiInputProcessor")
  csvinput.param("RawHitCollectionName", "zsdata")
  csvinput.param("RunNumber", runno)
  csvinput.param('FileNames', rawfile)
  csvinput.param("DetectorName","EUTelescope")
  path.add_processor(csvinput)
  
  return path 


def add_M26unpacker(path):
  """
  Adds M26 unpacker to the path
  """ 

  M26unpacker = tbsw.Processor(name="M26Unpacker",proctype="NIUnpacker")
  M26unpacker.param("InputCollectionName","NI")
  M26unpacker.param("OutputCollectionName","zsdata_m26")
  path.add_processor(M26unpacker)

  return path 

def add_CSVunpacker(path, ids=[1, 2, 3, 4, 5, 6], colname="zsdata_m26"): 
  """
  Adds csv unpacker to the path
  """ 
  csvunpacker = tbsw.Processor(name="CSVUnpacker",proctype="HitsFilterProcessor")   
  csvunpacker.param("InputCollectionName", "zsdata")
  csvunpacker.param("OutputCollectionName", colname)
  csvunpacker.param("FilterIDs"," ".join([ str(sensorID) for sensorID in ids ]))
  path.add_processor(csvunpacker)
   
  return path    

def add_M26clusterizer(path):
  """
  Adds M26 clusterizer to the path
  """  
    
  m26clust = tbsw.Processor(name="M26Clusterizer",proctype="PixelClusterizer")   
  m26clust.param("NoiseDBFileName","localDB/NoiseDB-M26.root")
  m26clust.param("SparseDataCollectionName","zsdata_m26")
  m26clust.param("ClusterCollectionName","zscluster_m26")
  m26clust.param("SparseClusterCut",0)
  m26clust.param("SparseSeedCut", 0)
  m26clust.param("SparseZSCut", 0)   
  path.add_processor(m26clust)  

  return path


def add_trackfinder(path, beamenergy, excludeplanes="3", minhits=6, maxTrkChi2=20, maxOutlierChi2=10, maxgap=0):
  """
  Adds trackfinding to the path
  """  
  #fixme: add gap to deal with additional sensors (timing planes)
  maxgap=1  
  trackfinder = tbsw.Processor(name="TF",proctype="FastTracker")
  trackfinder.param("InputHitCollectionNameVec","hit_m26")
  trackfinder.param("ExcludeDetector", excludeplanes)
  trackfinder.param("MaxTrackChi2", maxTrkChi2)  # loosecut version 10,000,000  
  trackfinder.param("MaximumGap", maxgap)
  trackfinder.param("MinimumHits",minhits)
  trackfinder.param("OutlierChi2Cut", maxOutlierChi2)  # loosecut version 10,000,000
  trackfinder.param("ParticleCharge","-1")
  trackfinder.param("ParticleMass","0.000511")
  trackfinder.param("ParticleMomentum", beamenergy)
  trackfinder.param("BackwardPass_FirstPlane","-1")
  trackfinder.param("BackwardPass_SecondPlane","-1")
  trackfinder.param("ForwardPass_FirstPlane","-1")
  trackfinder.param("ForwardPass_SecondPlane","-1")
  trackfinder.param("SingleHitSeeding", "0 1")
  trackfinder.param("MaxResidualU","0.4")
  trackfinder.param("MaxResidualV","0.4")
  path.add_processor(trackfinder)
  
  return path

def add_aligner(path, xerrors='0 10 10 0 10 10 0', yerrors='0 10 10 0 10 10 0', zerrors='0 10 10 0 10 10 0', gammaerrors='0 0.01 0.01 0 0.01 0.01 0'):
  """
  Adds aligner to the path
  """  
  aligner = tbsw.Processor(name="PreAligner",proctype="KalmanAligner")
  aligner.param('ErrorsShiftX' , xerrors)
  aligner.param('ErrorsShiftY' , yerrors)
  aligner.param('ErrorsShiftZ' , zerrors)
  aligner.param('ErrorsAlpha'  , '0 0 0 0 0 0 0')
  aligner.param('ErrorsBeta'   , '0 0 0 0 0 0 0')
  aligner.param('ErrorsGamma'  , gammaerrors)
  aligner.param('pValueCut'  , "0.0005")
  path.add_processor(aligner) 
  
  return path 

def add_M26hitmaker(path, hitmakertype):
  """
  Adds M26 hitmaker to the path
  """  
  if hitmakertype=="goe":
    m26goehitmaker = tbsw.Processor(name="M26GoeHitMaker",proctype="GoeHitMaker")
    m26goehitmaker.param("ClusterCollection","zscluster_m26")
    m26goehitmaker.param("HitCollectionName","hit_m26")
    m26goehitmaker.param("ClusterDBFileName","localDB/clusterDB-M26.root")
    m26goehitmaker.param("UseCenterOfGravityFallback","false")
    path.add_processor(m26goehitmaker)  
  else:
    m26coghitmaker = tbsw.Processor(name="M26CogHitMaker",proctype="CogHitMaker")
    m26coghitmaker.param("ClusterCollection","zscluster_m26")
    m26coghitmaker.param("HitCollectionName","hit_m26")
    m26coghitmaker.param("SigmaUCorrections", "0.62118 0.28235 0.31373")
    m26coghitmaker.param("SigmaVCorrections", "0.62118 0.28235 0.31373")
    path.add_processor(m26coghitmaker)  
    
  return path 


def add_clustercalibrator(path, use_outerplanes=False):
  """
  Adds M26 cluster calibration to the path
  """ 
  
  cluster_calibrator = tbsw.Processor(name="M26ClusterCalibrator",proctype="GoeClusterCalibrator")
  cluster_calibrator.param("ClusterDBFileName","localDB/clusterDB-M26.root")
  cluster_calibrator.param("MinClusters", "1000")
  
  if use_outerplanes:
    cluster_calibrator.param("SelectPlanes","0 1 2 4 5 6")
  else:
    cluster_calibrator.param("SelectPlanes","1 2 4 5")

  path.add_processor(cluster_calibrator)
  
  return path 


def append_tracklett_aligner(Env, paths, gearfile, nevents, beamenergy, hitmakertype, trackletttype):
  """
  Adds tracklet finder and aligner to the path
  """ 
  
  if trackletttype=="triplet":
    xerrors="0 10 0 0 0 0 0"
    yerrors="0 10 0 0 0 0 0"
    zerrors="0 10 0 0 0 0 0"
    gammaerrors="0 0.01 0 0 0 0 0"
    minhits=3
    excludeplanes="3 4 5 6"
  elif trackletttype=="quadruplet":
    xerrors="0 10 10 0 0 0 0"
    yerrors="0 10 10 0 0 0 0"
    zerrors="0 10 10 0 0 0 0"
    gammaerrors="0 0.01 0.01 0 0 0 0"
    minhits=4
    excludeplanes="3 5 6"
  elif trackletttype=="quintet":
    xerrors="0 10 10 0 10 0 0"
    yerrors="0 10 10 0 10 0 0"
    zerrors="0 10 10 0 10 0 0"
    gammaerrors="0 0.01 0.01 0 0.01 0 0"
    minhits=5
    excludeplanes="3 6"
  else:
    xerrors="0 10 10 0 10 10 0"
    yerrors="0 10 10 0 10 10 0"
    zerrors="0 10 10 0 10 10 0"
    gammaerrors="0 0.01 0.01 0 0.01 0.01 0"
    minhits=6
    excludeplanes="3"
  
  # Create path for tracklett prealigner 
  tracklett_prealigner_path = Env.create_path(hitmakertype+'hits_'+trackletttype+'_prealigner_path')
  tracklett_prealigner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  nevents, 'LCIOInputFiles': "tmp.slcio"}) 
  tracklett_prealigner_path=add_geometry(tracklett_prealigner_path, applyAlignment="true", overrideAlignment="true")
  tracklett_prealigner_path=add_M26hitmaker(tracklett_prealigner_path,hitmakertype) 
  tracklett_prealigner_path=add_trackfinder(tracklett_prealigner_path, beamenergy, excludeplanes, minhits, maxTrkChi2=10000000, maxOutlierChi2=10000000)
  tracklett_prealigner_path=add_aligner(tracklett_prealigner_path, xerrors, yerrors, '0 0 0 0 0 0 0', gammaerrors)
  
  # Perform the prealignment 2x
  for i in range(2):
    paths.append(tracklett_prealigner_path)      

  # Create path for tracklett aligner
  tracklett_aligner_path = Env.create_path(hitmakertype+'hits_'+trackletttype+'_aligner_path')
  tracklett_aligner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  nevents, 'LCIOInputFiles': "tmp.slcio"  }) 
  tracklett_aligner_path=add_geometry(tracklett_aligner_path, applyAlignment="true", overrideAlignment="true")
  tracklett_aligner_path=add_M26hitmaker(tracklett_aligner_path,hitmakertype) 
  tracklett_aligner_path=add_trackfinder(tracklett_aligner_path, beamenergy, excludeplanes, minhits, maxTrkChi2=20, maxOutlierChi2=10)
  tracklett_aligner_path=add_aligner(tracklett_aligner_path, xerrors, yerrors, zerrors, gammaerrors)
  
  # Perform the final alignment 3x
  for i in range(3):
    paths.append(tracklett_aligner_path) 
    
  # Create path for tracklett dqm
  tracklett_dqm_path = Env.create_path(hitmakertype+'hits_'+trackletttype+'_dqm_path')
  tracklett_dqm_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  nevents, 'LCIOInputFiles': "tmp.slcio"  })
  tracklett_dqm_path=add_geometry(tracklett_dqm_path, applyAlignment="true", overrideAlignment="true") 
  tracklett_dqm_path=add_M26hitmaker(tracklett_dqm_path,hitmakertype) 
  tracklett_dqm_path=add_trackfinder(tracklett_dqm_path, beamenergy, excludeplanes, minhits, maxTrkChi2=20, maxOutlierChi2=10)
  tracklett_dqm = tbsw.Processor(name="TelescopeDQM", proctype="TrackFitDQM") 
  tracklett_dqm.param("RootFileName","TelescopeDQM_"+hitmakertype+'hits_'+trackletttype+".root")
  tracklett_dqm_path.add_processor(tracklett_dqm)  
  paths.append(tracklett_dqm_path)
  
  # Create path for tracklett correlator
  tracklett_correlator_path = Env.create_path(hitmakertype+'hits_'+trackletttype+'_correlator_path')
  tracklett_correlator_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  nevents, 'LCIOInputFiles': "tmp.slcio"  }) 
  tracklett_correlator_path=add_geometry(tracklett_correlator_path, applyAlignment="true", overrideAlignment="true")
  tracklett_correlator_path=add_M26hitmaker(tracklett_correlator_path,"cog") 
  tracklett_correlator_path=add_trackfinder(tracklett_correlator_path, beamenergy, excludeplanes, minhits , maxTrkChi2=20, maxOutlierChi2=10)
  tracklett_correlator = tbsw.Processor(name="Tracklett_Correlator", proctype="TriplettCorrelator") 
  tracklett_correlator.param("InputHitCollectionNameVec","hit_m26")
  tracklett_correlator.param("OutputRootFileName","Correlator_"+hitmakertype+'hits_'+trackletttype+".root")
  tracklett_correlator.param("TrackCollectionName","tracks")
  tracklett_correlator_path.add_processor(tracklett_correlator)
  paths.append(tracklett_correlator_path)

  return paths


def create_anglereco_path(Env, rawfile, gearfile, numberofevents, usesinglehitseeding, useclusterdb, beamenergy, mcdata, csvdata=False):
  """
  Returns a list of tbsw path objects to reconstruct a test beam run 
  """
  
  reco_path = Env.create_path('reco')

  if not mcdata:
    reco_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : numberofevents }) 
    if csvdata:
      reco_path=add_csvrawinputprocessor(reco_path, rawfile)
      reco_path=add_CSVunpacker(reco_path)  
    else: 
      reco_path=add_rawinputprocessor(reco_path, rawfile)
      reco_path=add_M26unpacker(reco_path)  
  else:
    reco_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : numberofevents, 'LCIOInputFiles': rawfile })
  
  reco_path=add_geometry(reco_path, applyAlignment="true", overrideAlignment="true")
  reco_path=add_M26clusterizer(reco_path);

  if useclusterdb:
    reco_path=add_M26hitmaker(reco_path,"goe")
  else:
    reco_path=add_M26hitmaker(reco_path,"cog")

  downstream_TF=tbsw.Processor(name="downstream_TF",proctype="FastTracker")
  downstream_TF.param("InputHitCollectionNameVec","hit_m26")
  downstream_TF.param("HitCollectionName","unusedhits_down")
  downstream_TF.param("OutputTrackCollectionName","down_tracks")
  downstream_TF.param("ExcludeDetector", "0 1 2 3")
  downstream_TF.param("MaxTrackChi2", 20)
  downstream_TF.param("MaximumGap", 1)    #fixme: add gap to deal with additional sensors (timing planes)
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


  upstream_TF=tbsw.Processor(name="upstream_TF",proctype="FastTracker")
  upstream_TF.param("InputHitCollectionNameVec","hit_m26")
  upstream_TF.param("HitCollectionName","unusedhits_up")
  upstream_TF.param("OutputTrackCollectionName","up_tracks")
  upstream_TF.param("ExcludeDetector", "3 4 5 6")
  upstream_TF.param("MaxTrackChi2", 20)
  upstream_TF.param("MaximumGap", 1)   #fixme: add gap to deal with additional sensors (timing planes)
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

  x0imageproducer=tbsw.Processor(name="x0imageproducer",proctype="X0ImageProducer")
  x0imageproducer.param("DownStreamTrackCollection","down_tracks")
  x0imageproducer.param("UpStreamTrackCollection","up_tracks")
  x0imageproducer.param("DUTPlane","3")
  x0imageproducer.param("VertexFitSwitch","false")
  x0imageproducer.param("ToyBetheHeitlerSwitch","true")
  x0imageproducer.param("ToyScatteringSwitch","false")
  x0imageproducer.param("ToyRecoError","-1")
  x0imageproducer.param("MaxDist","0.5")
  x0imageproducer.param("RootFileName","X0.root")
  reco_path.add_processor(x0imageproducer)
    
  return [ reco_path ]



def create_x0analysis_calibration_paths(Env, rawfile, gearfile, nevents, useClusterDB, beamenergy, mcdata, isLongTelescope=True, use_outerplanes=False, csvdata=False):
  """
  Returns a list of tbsw path objects to calibrate a tracking telescope
  """
  
  # Calibrations are organized in a sequence of calibration paths. 
  # The calibration paths are collected in a list for later execution
  calpaths = []

  # Create path for noisy pixel masking
  mask_path = Env.create_path('mask_path')

  if not mcdata:
    mask_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 1000000, 'Verbosity': "MESSAGE3" }) 
    if csvdata:
      mask_path=add_csvrawinputprocessor(mask_path, rawfile)
      mask_path=add_CSVunpacker(mask_path)  
    else: 
      mask_path=add_rawinputprocessor(mask_path, rawfile)
      mask_path=add_M26unpacker(mask_path)  
  else:
    mask_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 1000000, 'LCIOInputFiles': rawfile }) 
    
  
  mask_path=add_geometry(mask_path, applyAlignment="true", overrideAlignment="true")

  m26hotpixelkiller = tbsw.Processor(name="M26HotPixelKiller",proctype="HotPixelKiller")
  m26hotpixelkiller.param("InputCollectionName", "zsdata_m26")
  m26hotpixelkiller.param("MaxNormedOccupancy", 5)
  m26hotpixelkiller.param("MinNormedOccupancy", -1)
  m26hotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-M26.root")
  m26hotpixelkiller.param("OfflineZSThreshold", 0)
  mask_path.add_processor(m26hotpixelkiller)
  
  # Finished with path for noisy pixel masking
  calpaths.append(mask_path)
  
  # Create path for detector level creation of clusters
  clusterizer_path = Env.create_path('clusterizer_path')
  
  if not mcdata:
    clusterizer_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents, 'Verbosity': "MESSAGE3" }) 
    if csvdata:
      clusterizer_path=add_csvrawinputprocessor(clusterizer_path, rawfile)
      clusterizer_path=add_CSVunpacker(clusterizer_path)  
    else: 
      clusterizer_path=add_rawinputprocessor(clusterizer_path, rawfile)
      clusterizer_path=add_M26unpacker(clusterizer_path)  
  else:
    clusterizer_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents, 'LCIOInputFiles': rawfile }) 
    

  clusterizer_path=add_geometry(clusterizer_path, applyAlignment="true", overrideAlignment="true")
  clusterizer_path=add_M26clusterizer(clusterizer_path)

  lciooutput = tbsw.Processor(name="LCIOOutput",proctype="LCIOOutputProcessor")
  lciooutput.param("LCIOOutputFile","tmp.slcio")
  lciooutput.param("LCIOWriteMode","WRITE_NEW")
  clusterizer_path.add_processor(lciooutput)  
  
  # Finished with path for clusters
  calpaths.append(clusterizer_path) 
  
  # Create path for m26 mc cluster calibration
  if mcdata:
    mc_clustercal_path = Env.create_path('mc_clustercal_path')
    mc_clustercal_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents, 'LCIOInputFiles': "tmp.slcio" })  
    mc_clustercal_path=add_geometry(mc_clustercal_path, applyAlignment="true", overrideAlignment="true")
    mc_clustercal_path=add_M26hitmaker(mc_clustercal_path,"cog") 

    m26clustercalibrationformc = tbsw.Processor(name="M26ClusterCalibrationFromMC", proctype="GoeClusterCalibratorForMC")
    m26clustercalibrationformc.param("HitCollection","hit_m26")
    m26clustercalibrationformc.param("SimTrackerHitCollection","SimTrackerHits")
    m26clustercalibrationformc.param("ClusterDBFileName","localDB/clusterDB-M26-MC.root")
    m26clustercalibrationformc.param("MaxResidualU", 0.1)
    m26clustercalibrationformc.param("MaxResidualV", 0.1)
    m26clustercalibrationformc.param("MinClusters", 500) 
    mc_clustercal_path.add_processor(m26clustercalibrationformc)  

    # Finished with mc cluster calibration
    calpaths.append(mc_clustercal_path)


  # Create path for correlator
  correlator_path = Env.create_path('correlator_path')
  correlator_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  nevents, 'LCIOInputFiles': "tmp.slcio", 'Verbosity': "MESSAGE3"  }) 
  correlator_path=add_geometry(correlator_path, applyAlignment="false", overrideAlignment="true")
  correlator_path=add_M26hitmaker(correlator_path,"cog") 
  
  hitdqm = tbsw.Processor(name="RawDQM",proctype="RawHitDQM")
  hitdqm.param("InputHitCollectionNameVec","hit_m26")  
  hitdqm.param("RootFileName","RawDQM.root")
  correlator_path.add_processor(hitdqm)  
  
  correlator = tbsw.Processor(name="TelCorrelator", proctype="Correlator")
  correlator.param("InputHitCollectionNameVec","hit_m26")
  correlator.param("OutputRootFileName","XCorrelator.root")
  correlator.param("ReferencePlane","0")
  correlator.param("ParticleCharge","-1")
  correlator.param("ParticleMass","0.000511")
  correlator.param("ParticleMomentum", beamenergy)
  correlator_path.add_processor(correlator)  

  # Finished with correlator
  calpaths.append(correlator_path) 

  if isLongTelescope: 
    # Append tracklett aligners
    calpaths=append_tracklett_aligner(Env, calpaths, gearfile, nevents, beamenergy, "cog", "triplet")
    calpaths=append_tracklett_aligner(Env, calpaths, gearfile, nevents, beamenergy, "cog", "quadruplet")
    calpaths=append_tracklett_aligner(Env, calpaths, gearfile, nevents, beamenergy, "cog", "quintet")

  calpaths=append_tracklett_aligner(Env, calpaths, gearfile, nevents, beamenergy, "cog", "full")
  
  if useClusterDB:   
    # The code below produces cluster calibration constants
    # (clusterDB). IF you only want to use CoG hits, this part
    # is not needed. 
    
    # Compute a first iteration of clusterDB where track fits are done using Center of Gravity 
    cluster_calibration1_path = Env.create_path('cluster_calibration_path')
    cluster_calibration1_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents, 'LCIOInputFiles': "tmp.slcio" }) 
    cluster_calibration1_path=add_geometry(cluster_calibration1_path, applyAlignment="true", overrideAlignment="true")
    cluster_calibration1_path=add_M26hitmaker(cluster_calibration1_path,"cog")
    cluster_calibration1_path=add_trackfinder(cluster_calibration1_path, beamenergy, excludeplanes="3", minhits=6, maxTrkChi2=20, maxOutlierChi2=10)
    cluster_calibration1_path=add_clustercalibrator(cluster_calibration1_path, use_outerplanes)
     
    calpaths.append(cluster_calibration1_path)
    
    # Perform a track based alignment using calibrated clusters
    calpaths=append_tracklett_aligner(Env, calpaths, gearfile, nevents, beamenergy, "goe", "full")
    
    # Compute a cluster calibration iteratively where track fit are done using calibrated clusters
    # This needs to be iterated multiple times
    cluster_calibration2_path = Env.create_path('cluster_calibration2_path')
    cluster_calibration2_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents, 'LCIOInputFiles': "tmp.slcio" }) 
    cluster_calibration2_path=add_geometry(cluster_calibration2_path, applyAlignment="true", overrideAlignment="true")
    cluster_calibration2_path=add_M26hitmaker(cluster_calibration2_path,"goe")
    cluster_calibration2_path=add_trackfinder(cluster_calibration2_path, beamenergy, excludeplanes="3", minhits=6, maxTrkChi2=20, maxOutlierChi2=10)
    cluster_calibration2_path=add_clustercalibrator(cluster_calibration2_path, use_outerplanes)
    # Finished with second part of cluster calibration
    # Repeat this 7x
    for i in range(7):
      calpaths.append(cluster_calibration2_path)
    
    # Perform a track based alignment using calibrated clusters
    calpaths=append_tracklett_aligner(Env, calpaths, gearfile, nevents, beamenergy, "goe", "full")
    
    # Create tracking DQM  using cluster calibration
    dqm2_path = Env.create_path('dqm2_path')
    dqm2_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" }) 
    dqm2_path=add_geometry(dqm2_path, applyAlignment="true", overrideAlignment="true") 
    dqm2_path=add_M26hitmaker(dqm2_path,"goe")
    dqm2_path=add_trackfinder(dqm2_path, beamenergy, excludeplanes="3", minhits=6, maxTrkChi2=20, maxOutlierChi2=10)
    
    teldqm = tbsw.Processor(name="TelescopeDQM", proctype="TrackFitDQM") 
    teldqm.param("RootFileName","TelescopeDQM2.root")
    dqm2_path.add_processor(teldqm)  
    
    calpaths.append(dqm2_path)
    
  return calpaths

def create_x0sim_path(Env, name, rawfile, gearfile, nevents, beamenergy):
  """
  Returns a list of tbsw path objects to simulate a test beam run 
  """
  
  x0sim = Env.create_path(name)
  x0sim.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents, 'Verbosity': "MESSAGE3"})   
  
  infosetter = tbsw.Processor(name="InfoSetter", proctype='EventInfoSetter')
  infosetter.param("RunNumber","0")
  infosetter.param("DetectorName","EUTelescope") 
  x0sim.add_processor(infosetter)
  
  geo_noalign = tbsw.Processor(name="Geo",proctype="Geometry")
  geo_noalign.param("AlignmentDBFilePath", "localDB/alignmentDB.root")
  geo_noalign.param("ApplyAlignment", "false")
  geo_noalign.param("OverrideAlignment", "true")
  x0sim.add_processor(geo_noalign)
   
  gun = tbsw.Processor(name="ParticleGun",proctype="ParticleGunGenerator")
  gun.param("BeamIntensity","20000")
  gun.param("BeamMomentum", str(beamenergy))
  gun.param("BeamMomentumSigma", 0.001)
  gun.param("BeamTimeWindow", 0.0001)
  gun.param("BeamVertexX","0")
  gun.param("BeamVertexY","0")
  gun.param("BeamVertexZ","-10")  
  gun.param("BeamVertexXSigma","10")
  gun.param("BeamVertexYSigma","5")  
  gun.param("BeamSlopeXSigma","0.0035") 
  gun.param("BeamSlopeYSigma","0.0035") 
  gun.param("CorrelationVertexXvsMomentum","0.005")  
  gun.param("CorrelationVertexXvsSlopeX","0.1") 
  gun.param("CorrelationVertexYvsMomentum","0.005") 
  gun.param("CorrelationVertexYvsSlopeY","0.1") 
  gun.param("PDG","11")
  gun.param("ParticleCharge","-1")
  gun.param("ParticleMass","0.000511")
  x0sim.add_processor(gun)
  
  fastsim = tbsw.Processor(name="FastSim",proctype="FastSimulation")
  fastsim.param("MCParticleCollectionName","MCParticles")
  fastsim.param("SimTrackerHitCollectionName","SimTrackerHits")
  fastsim.param("ScatterModel","0")
  fastsim.param("DoEnergyLossStraggling","false")
  fastsim.param("DoFractionalBetheHeitlerEnergyLoss","false")
  x0sim.add_processor(fastsim)
 
  m26digi = tbsw.Processor(name="M26Digitizer",proctype="SiPixDigitizer")
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
  x0sim.add_processor(m26digi)
  
  lciooutput = tbsw.Processor(name="LCIOOutput",proctype="LCIOOutputProcessor")
  lciooutput.param("LCIOOutputFile",rawfile)
  lciooutput.param("LCIOWriteMode","WRITE_NEW")  
  x0sim.add_processor(lciooutput)
    
  return x0sim

