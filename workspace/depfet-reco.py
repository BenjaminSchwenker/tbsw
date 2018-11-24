"""
Script for processing depfet testbeam November 2018.

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *

def add_unpackers(path):
  """
  Adds unpackers to the path
  """  
  
  m26unpacker = Processor(name="M26Unpacker",proctype="NIUnpacker")   
  m26unpacker.param("InputCollectionName", "NI")
  m26unpacker.param("OutputCollectionName","zsdata_m26")
  path.add_processor(m26unpacker)
  
  fei4unpacker = Processor(name="FEI4Unpacker",proctype="PyBARUnpacker")
  fei4unpacker.param("InputCollectionName","PyBAR")
  fei4unpacker.param("OutputCollectionName", "zsdata_fei4")
  fei4unpacker.param("SensorID","21")
  path.add_processor(fei4unpacker)   
  
  pxdunpacker = Processor(name="PXDUnpacker",proctype="DEPFETUnpacker")
  pxdunpacker.param('InputCollectionName', 'DEPFET')
  pxdunpacker.param('OutputCollectionName','zsdata_pxd')
  path.add_processor(pxdunpacker)
  
  h5unpacker = Processor(name="H5Unpacker",proctype="DEPFETUnpacker")
  h5unpacker.param('InputCollectionName','DEPFE5')
  h5unpacker.param('OutputCollectionName','zsdata_h5')
  h5unpacker.param("Mapping","Hybrid5")  
  path.add_processor(h5unpacker)   
  
  return path


def add_clusterizers(path):
  """
  Adds clusterizers to the path
  """  
    
  m26clust = Processor(name="M26Clusterizer",proctype="PixelClusterizer")   
  m26clust.param("NoiseDBFileName","localDB/NoiseDB-M26.root")
  m26clust.param("SparseDataCollectionName","zsdata_m26")
  m26clust.param("ClusterCollectionName","zscluster_m26")
  m26clust.param("SparseClusterCut",0)
  m26clust.param("SparseSeedCut", 0)
  m26clust.param("SparseZSCut", 0)   
  path.add_processor(m26clust)  

  fei4clust = Processor(name="FEI4Clusterizer",proctype="PixelClusterizer")   
  fei4clust.param("NoiseDBFileName","localDB/NoiseDB-FEI4.root")
  fei4clust.param("SparseDataCollectionName","zsdata_fei4")
  fei4clust.param("ClusterCollectionName","zscluster_fei4")
  fei4clust.param("SparseClusterCut",0)
  fei4clust.param("SparseSeedCut", 0)
  fei4clust.param("SparseZSCut", 0)   
  path.add_processor(fei4clust)  

  pxdclust = Processor(name="PXDClusterizer",proctype="PixelClusterizer")   
  pxdclust.param("NoiseDBFileName","localDB/NoiseDB-PXD.root")
  pxdclust.param("SparseDataCollectionName","zsdata_pxd")
  pxdclust.param("ClusterCollectionName","zscluster_pxd")
  pxdclust.param("SparseClusterCut",4)
  pxdclust.param("SparseSeedCut", 4)
  pxdclust.param("SparseZSCut", 4)   
  path.add_processor(pxdclust)  
  
  h5clust = Processor(name="H5Clusterizer",proctype="PixelClusterizer")   
  h5clust.param("NoiseDBFileName","localDB/NoiseDB-H5.root")
  h5clust.param("SparseDataCollectionName","zsdata_h5")
  h5clust.param("ClusterCollectionName","zscluster_h5")
  h5clust.param("SparseClusterCut",4)
  h5clust.param("SparseSeedCut", 4)
  h5clust.param("SparseZSCut", 4)   
  path.add_processor(h5clust)  

  return path

def add_hitmakers(path):
  """
  Adds center of gravity hitmakers to the path
  """  

  m26hitmaker = Processor(name="M26CogHitMaker",proctype="CogHitMaker")
  m26hitmaker.param("ClusterCollection","zscluster_m26")
  m26hitmaker.param("HitCollectionName","hit_m26")
  m26hitmaker.param("SigmaU1",0.0037)
  m26hitmaker.param("SigmaU2",0.0033)
  m26hitmaker.param("SigmaU3",0.0050)
  m26hitmaker.param("SigmaV1",0.0037)
  m26hitmaker.param("SigmaV2",0.0033)
  m26hitmaker.param("SigmaV3",0.0050)
  path.add_processor(m26hitmaker)

  fei4hitmaker = Processor(name="FEI4CogHitMaker",proctype="CogHitMaker")
  fei4hitmaker.param("ClusterCollection","zscluster_fei4")
  fei4hitmaker.param("HitCollectionName","hit_fei4")
  fei4hitmaker.param("SigmaU1",0.072)
  fei4hitmaker.param("SigmaU2",0.072)
  fei4hitmaker.param("SigmaU3",0.072)
  fei4hitmaker.param("SigmaV1",0.0144)
  fei4hitmaker.param("SigmaV2",0.0144)
  fei4hitmaker.param("SigmaV3",0.0144)
  path.add_processor(fei4hitmaker)

  pxdhitmaker = Processor(name="PXDCogHitMaker",proctype="CogHitMaker")
  pxdhitmaker.param("ClusterCollection","zscluster_pxd")
  pxdhitmaker.param("HitCollectionName","hit_pxd")
  pxdhitmaker.param("SigmaU1",0.0134)
  pxdhitmaker.param("SigmaU2",0.0077)
  pxdhitmaker.param("SigmaU3",0.0077)
  pxdhitmaker.param("SigmaV1",0.024)
  pxdhitmaker.param("SigmaV2",0.014)
  pxdhitmaker.param("SigmaV3",0.014)
  path.add_processor(pxdhitmaker)
  
  h5hitmaker = Processor(name="H5CogHitMaker",proctype="CogHitMaker")
  h5hitmaker.param("ClusterCollection","zscluster_h5")
  h5hitmaker.param("HitCollectionName","hit_h5")
  h5hitmaker.param("SigmaU1",0.0134)
  h5hitmaker.param("SigmaU2",0.0077)
  h5hitmaker.param("SigmaU3",0.0077)
  h5hitmaker.param("SigmaV1",0.0134)
  h5hitmaker.param("SigmaV2",0.0077)
  h5hitmaker.param("SigmaV3",0.0077)
  path.add_processor(h5hitmaker)   
  
  return path

def add_hitmakersDB(path):
  """
  Add cluster shape hitmakers to the path (requiring clusterDBs)
  """  
  
  m26goehitmaker = Processor(name="M26GoeHitMaker",proctype="GoeHitMaker")   
  m26goehitmaker.param("ClusterCollection","zscluster_m26")
  m26goehitmaker.param("HitCollectionName","hit_m26")
  m26goehitmaker.param("ClusterDBFileName","localDB/clusterDB-M26.root")
  path.add_processor(m26goehitmaker)  
    
  fei4goehitmaker = Processor(name="FEI4GoeHitMaker",proctype="GoeHitMaker")   
  fei4goehitmaker.param("ClusterCollection","zscluster_fei4")
  fei4goehitmaker.param("HitCollectionName","hit_fei4")
  fei4goehitmaker.param("ClusterDBFileName","localDB/clusterDB-FEI4.root")
  path.add_processor(fei4goehitmaker) 
  
  pxdgoehitmaker = Processor(name="PXDGoeHitMaker",proctype="GoeHitMaker")   
  pxdgoehitmaker.param("ClusterCollection","zscluster_pxd")
  pxdgoehitmaker.param("HitCollectionName","hit_pxd")
  pxdgoehitmaker.param("ClusterDBFileName","localDB/clusterDB-PXD.root")
  path.add_processor(pxdgoehitmaker)   
  
  h5goehitmaker = Processor(name="H5GoeHitMaker",proctype="GoeHitMaker")   
  h5goehitmaker.param("ClusterCollection","zscluster_h5")
  h5goehitmaker.param("HitCollectionName","hit_h5")
  h5goehitmaker.param("ClusterDBFileName","localDB/clusterDB-H5.root")
  path.add_processor(h5goehitmaker)   

  return path

def add_clustercalibrators(path):
  """
  Add cluster calibration processors to create clusterDB's
  """
  
  m26clustdb = Processor(name="M26ClusterCalibrator",proctype="GoeClusterCalibrator")   
  m26clustdb.param("AlignmentDBFileName","localDB/alignmentDB.root")
  m26clustdb.param("ClusterDBFileName","localDB/clusterDB-M26.root")  
  m26clustdb.param("MinClusters","200")
  m26clustdb.param("IgnoreIDs","6 7 21")
  path.add_processor(m26clustdb)  
    
  pxdclustdb = Processor(name="PXDClusterCalibrator",proctype="GoeClusterCalibrator")   
  pxdclustdb.param("AlignmentDBFileName","localDB/alignmentDB.root")
  pxdclustdb.param("ClusterDBFileName","localDB/clusterDB-PXD.root")  
  pxdclustdb.param("MinClusters","200")
  pxdclustdb.param("IgnoreIDs","0 1 2 3 4 5 7 21")
  path.add_processor(pxdclustdb)  
    
  h5clustdb = Processor(name="H5ClusterCalibrator",proctype="GoeClusterCalibrator")   
  h5clustdb.param("AlignmentDBFileName","localDB/alignmentDB.root")
  h5clustdb.param("ClusterDBFileName","localDB/clusterDB-H5.root")  
  h5clustdb.param("MinClusters","200")
  h5clustdb.param("IgnoreIDs","0 1 2 3 4 5 6 21")
  path.add_processor(h5clustdb)  

  fei4clustdb = Processor(name="FEI4ClusterCalibrator",proctype="GoeClusterCalibrator")   
  fei4clustdb.param("AlignmentDBFileName","localDB/alignmentDB.root")
  fei4clustdb.param("ClusterDBFileName","localDB/clusterDB-FEI4.root")  
  fei4clustdb.param("MinClusters","200")
  fei4clustdb.param("IgnoreIDs","0 1 2 3 4 5 6 7")
  path.add_processor(fei4clustdb)  
  
  return path

def create_calibration_path(Env, rawfile, gearfile, energy, useClusterDB):
  """
  Returns a list of tbsw path objects needed to calibrate the tracking telescope
  """

  # Calibrations are organized in a sequence of calibration paths. 
  # The calibration paths are collected in a list for later execution
  calpaths = []
  
  # Create path for detector level masking of hot channels 
  mask_path = Env.create_path('mask_path')
  mask_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1 })  
  
  rawinput = Processor(name="RawInputProcessor",proctype="EudaqInputProcessor")
  rawinput.param('FileNames', rawfile)
  mask_path.add_processor(rawinput)
  
  mask_path = add_unpackers(mask_path)
   
  m26hotpixelkiller = Processor(name="M26HotPixelKiller",proctype="HotPixelKiller",)
  m26hotpixelkiller.param("InputCollectionName", "zsdata_m26")
  m26hotpixelkiller.param("MaxOccupancy", 0.001)
  m26hotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-M26.root")
  m26hotpixelkiller.param("OfflineZSThreshold", 0)
  mask_path.add_processor(m26hotpixelkiller)
   
  fei4hotpixelkiller = Processor(name="FEI4HotPixelKiller", proctype="HotPixelKiller")
  fei4hotpixelkiller.param("InputCollectionName", "zsdata_fei4")
  fei4hotpixelkiller.param("MaxOccupancy", 0.001)
  fei4hotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-FEI4.root")
  fei4hotpixelkiller.param("OfflineZSThreshold", 0)
  mask_path.add_processor(fei4hotpixelkiller)
   
  pxdhotpixelkiller = Processor(name="PXDHotPixelKiller", proctype="HotPixelKiller")
  pxdhotpixelkiller.param("InputCollectionName", "zsdata_pxd")
  pxdhotpixelkiller.param("MaxOccupancy", 0.001)
  pxdhotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-PXD.root")
  pxdhotpixelkiller.param("OfflineZSThreshold", 0)
  mask_path.add_processor(pxdhotpixelkiller)  
  
  h5hotpixelkiller = Processor(name="H5HotPixelKiller", proctype="HotPixelKiller")
  h5hotpixelkiller.param("InputCollectionName", "zsdata_h5")
  h5hotpixelkiller.param("MaxOccupancy", 0.001)
  h5hotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-H5.root")
  h5hotpixelkiller.param("OfflineZSThreshold", 0)
  mask_path.add_processor(h5hotpixelkiller)  
  
  # Add path for masking
  calpaths.append(mask_path)  

  # Create path for detector level creation of clusters
  clusterizer_path = Env.create_path('clusterizer_path')
  clusterizer_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  -1}) 
  
  clusterizer_path.add_processor(rawinput)  
  clusterizer_path = add_unpackers(clusterizer_path) 
  clusterizer_path = add_clusterizers(clusterizer_path)    
   
  lciooutput = Processor(name="LCIOOutput",proctype="LCIOOutputProcessor")
  lciooutput.param("LCIOOutputFile","tmp.slcio")
  lciooutput.param("LCIOWriteMode","WRITE_NEW")
  clusterizer_path.add_processor(lciooutput)  
   
  # Finished with path for clusterizers
  calpaths.append(clusterizer_path)   
  
  # Create path for pre alignmnet and dqm based on hits
  correlator_path = Env.create_path('correlator_path')
  correlator_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 200000, 'LCIOInputFiles': "tmp.slcio" })  
  
  correlator_path = add_hitmakers(correlator_path) 
  
  hitdqm = Processor(name="RawDQM",proctype="RawHitDQM")
  hitdqm.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_pxd hit_h5")  
  hitdqm.param("RootFileName","RawDQM.root")
  correlator_path.add_processor(hitdqm)  
   
  correlator = Processor(name="TelCorrelator", proctype="Correlator")
  correlator.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_pxd hit_h5")
  correlator.param("AlignmentDBFileName", "localDB/alignmentDB.root")
  correlator.param("OutputRootFileName","XCorrelator.root")
  correlator.param("ReferencePlane","0")
  correlator.param("ParticleCharge","-1")
  correlator.param("ParticleMass","0.000511")
  correlator.param("ParticleMomentum", energy)
  correlator_path.add_processor(correlator)  
  
  # Finished with path for hit based pre alignment
  calpaths.append(correlator_path)  
  
  # Create path for pre alignment with loose cut track sample 
  prealigner_path = Env.create_path('prealigner_path')
  prealigner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 200000, 'LCIOInputFiles': "tmp.slcio" })  
  
  prealigner_path = add_hitmakers(prealigner_path)
   
  trackfinder_loosecut = Processor(name="AlignTF_LC",proctype="FastTracker")
  trackfinder_loosecut.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_pxd hit_h5")
  trackfinder_loosecut.param("AlignmentDBFileName","localDB/alignmentDB.root")
  trackfinder_loosecut.param("ExcludeDetector", "")
  trackfinder_loosecut.param("MaxTrackChi2", 10000000)
  trackfinder_loosecut.param("MaximumGap", 1)
  trackfinder_loosecut.param("MinimumHits",7)
  trackfinder_loosecut.param("OutlierChi2Cut", 100000000)
  trackfinder_loosecut.param("ParticleCharge","-1")
  trackfinder_loosecut.param("ParticleMass","0.000511")
  trackfinder_loosecut.param("ParticleMomentum", energy)
  trackfinder_loosecut.param("SingleHitSeeding", "0")
  trackfinder_loosecut.param("MaxResidualU","0.5")
  trackfinder_loosecut.param("MaxResidualV","0.5")
  prealigner_path.add_processor(trackfinder_loosecut)
 
  prealigner = Processor(name="PreAligner",proctype="KalmanAligner")
  prealigner.param("AlignmentDBFileName","localDB/alignmentDB.root")
  prealigner.param('ErrorsShiftX' , '0 10 10 10 10 10 0 10 10')
  prealigner.param('ErrorsShiftY' , '0 10 10 10 10 10 0 10 10')
  prealigner.param('ErrorsShiftZ' , '0 0 0 0 0 0 0 0 0')
  prealigner.param('ErrorsAlpha'  , '0 0 0 0 0 0 0 0 0')
  prealigner.param('ErrorsBeta'   , '0 0 0 0 0 0 0 0 0')
  prealigner.param('ErrorsGamma'  , '0 0.01 0.01 0.01 0.01 0.01 0 0.01 0.01')
  prealigner_path.add_processor(prealigner)  

  # Finished with path for prealigner
  calpaths.append(prealigner_path)  
  
  # Create path for alignment with tight cut track sample 
  aligner_path = Env.create_path('aligner_path')
  aligner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 200000, 'LCIOInputFiles': "tmp.slcio" })  
  
  aligner_path = add_hitmakers(aligner_path)
  
  trackfinder_tightcut = Processor(name="AlignTF_TC",proctype="FastTracker")
  trackfinder_tightcut.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_pxd hit_h5")
  trackfinder_tightcut.param("AlignmentDBFileName","localDB/alignmentDB.root")
  trackfinder_tightcut.param("ExcludeDetector", "")
  trackfinder_tightcut.param("MaxTrackChi2", 100)
  trackfinder_tightcut.param("MaximumGap", 1)
  trackfinder_tightcut.param("MinimumHits",7)
  trackfinder_tightcut.param("OutlierChi2Cut", 20)
  trackfinder_tightcut.param("ParticleCharge","-1")
  trackfinder_tightcut.param("ParticleMass","0.000511")
  trackfinder_tightcut.param("ParticleMomentum", energy)
  trackfinder_tightcut.param("SingleHitSeeding", "0")
  trackfinder_tightcut.param("MaxResidualU","0.4")
  trackfinder_tightcut.param("MaxResidualV","0.4")
  aligner_path.add_processor(trackfinder_tightcut)
   
  aligner = Processor(name="Aligner",proctype="KalmanAligner")
  aligner.param("AlignmentDBFileName","localDB/alignmentDB.root")
  aligner.param('ErrorsShiftX' , '0 10 10 10 10 10 0 10 10' )
  aligner.param('ErrorsShiftY' , '0 10 10 10 10 10 0 10 10')
  aligner.param('ErrorsShiftZ' , '0 10 10 10 10 10 0 10 10')
  aligner.param('ErrorsAlpha'  , '0 0 0 0.01 0 0 0 0.01 0.01')
  aligner.param('ErrorsBeta'   , '0 0 0 0.01 0 0 0 0.01 0.01')
  aligner.param('ErrorsGamma'  , '0 0.01 0.01 0.01 0.01 0.01 0 0.01 0.01')
  aligner_path.add_processor(aligner)    
  
  # Finished with path for aligner
  # Repeat this 3x
  calpaths.append(aligner_path)  
  calpaths.append(aligner_path)
  calpaths.append(aligner_path)
   
  # Creeate path for some track based dqm using current calibrations
  dqm_path = Env.create_path('dqm_path')
  dqm_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 200000, 'LCIOInputFiles': "tmp.slcio" })  

  dqm_path = add_hitmakers(dqm_path)
  dqm_path.add_processor(trackfinder_tightcut)

  teldqm = Processor(name="TelescopeDQM", proctype="TrackFitDQM") 
  teldqm.param("AlignmentDBFileName","localDB/alignmentDB.root")
  teldqm.param("RootFileName","TelescopeDQM.root")
  dqm_path.add_processor(teldqm)  
  
  # Finished with path for teldqm
  calpaths.append(dqm_path)
  
  if useClusterDB: 
    # The code below produces cluster calibration constants
    # (clusterDB). IF you only want to use CoG hits, this part
    # is not needed.
    
    # Creeate path for first iteration for computing clusterDBs for all sensors 
    preclustercal_path = Env.create_path('preclustercal_path')
    preclustercal_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1, 'LCIOInputFiles': "tmp.slcio" })  
    preclustercal_path = add_hitmakers(preclustercal_path) 
    preclustercal_path.add_processor(trackfinder_tightcut)      
    preclustercal_path = add_clustercalibrators(preclustercal_path)
    
    # Finished with path for pre cluster calibration 
    calpaths.append(preclustercal_path)
    
    # Create path for alignment with tight cut track sample and cluster DB
    aligner_db_path = Env.create_path('aligner_db_path')
    aligner_db_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 200000, 'LCIOInputFiles': "tmp.slcio" })  
    
    aligner_db_path = add_hitmakersDB(aligner_db_path) 
    aligner_db_path.add_processor(trackfinder_tightcut) 
    aligner_db_path.add_processor(aligner)   
    
    # Finished with path for alignemnt with hits from pre clusterDB 
    # Repeat this 2x
    for i in range(2):
      calpaths.append(aligner_db_path) 
    
    # Creeate path for next iterations for computing clusterDBs for all sensors 
    clustercal_path = Env.create_path('clustercal_path')
    clustercal_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1, 'LCIOInputFiles': "tmp.slcio" })   
    clustercal_path = add_hitmakersDB(clustercal_path) 
    clustercal_path.add_processor(trackfinder_tightcut) 
    clustercal_path = add_clustercalibrators(clustercal_path)
     
    # Finished with path for pre cluster calibration
    # Repeat this 6x
    for i in range(6): 
      calpaths.append(clustercal_path)
             
    # Finished with path for alignemnt with hits from final clusterDB 
    # Repeat this 2x
    for i in range(2):
      calpaths.append(aligner_db_path) 
    
    # Creeate path for dqm using cluster calibrations
    dqm_db_path = Env.create_path('dqm_db_path')
    dqm_db_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 200000, 'LCIOInputFiles': "tmp.slcio" })
    
    dqm_db_path = add_hitmakersDB(dqm_db_path)   
    dqm_db_path.add_processor(trackfinder_tightcut)  
    
    teldqm_db = Processor(name="TelescopeDQM_DB", proctype="TrackFitDQM") 
    teldqm_db.param("AlignmentDBFileName","localDB/alignmentDB.root")
    teldqm_db.param("RootFileName","TelescopeDQM_DB.root")
    dqm_db_path.add_processor(teldqm_db)  
    
    # Finished with path for dqm with cluster calibration
    calpaths.append(dqm_db_path)
    
  return calpaths


def create_reco_path(Env, rawfile, gearfile, energy, useClusterDB):
  """
  Returns a list of tbsw path objects for reconstruciton of a test beam run 
  """
  
  reco_path = Env.create_path('reco_path')
  reco_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1 }) 
  
  rawinput = Processor(name="RawInputProcessor",proctype="EudaqInputProcessor")
  rawinput.param('FileNames', rawfile) 
  reco_path.add_processor(rawinput)
  
  # Create path for all reconstruction up to hits
  reco_path = add_unpackers(reco_path)    
  reco_path = add_clusterizers(reco_path)    
  
  if useClusterDB: 
    reco_path = add_hitmakersDB(reco_path)   
  else: 
    reco_path = add_hitmakers(reco_path) 

  trackfinder = Processor(name="TrackFinder",proctype="FastTracker")
  trackfinder.param("InputHitCollectionNameVec","hit_m26 hit_fei4")
  trackfinder.param("AlignmentDBFileName","localDB/alignmentDB.root")
  trackfinder.param("ExcludeDetector", "3 8")
  trackfinder.param("MaxTrackChi2", "100")
  trackfinder.param("MaximumGap", "1")
  trackfinder.param("MinimumHits","6")
  trackfinder.param("OutlierChi2Cut", "20")
  trackfinder.param("ParticleCharge","-1")
  trackfinder.param("ParticleMass","0.000511")
  trackfinder.param("ParticleMomentum", energy)
  trackfinder.param("SingleHitSeeding", "0")
  trackfinder.param("MaxResidualU","0.4")
  trackfinder.param("MaxResidualV","0.4")
  reco_path.add_processor(trackfinder)  

  hybrid_analyzer = Processor(name="HybridAnalyzer",proctype="PixelDUTAnalyzer")
  hybrid_analyzer.param("HitCollection","hit_h5")  
  hybrid_analyzer.param("DigitCollection","zsdata_h5")
  hybrid_analyzer.param("AlignmentDBFileName","localDB/alignmentDB.root")
  hybrid_analyzer.param("DUTPlane","8")
  hybrid_analyzer.param("ReferencePlane","7")
  hybrid_analyzer.param("MaxResidualU","0.2")
  hybrid_analyzer.param("MaxResidualU","0.2")
  hybrid_analyzer.param("RootFileName","Histos-H5.root")
  reco_path.add_processor(hybrid_analyzer)  

  pxd_analyzer = Processor(name="PXDAnalyzer",proctype="PixelDUTAnalyzer")
  pxd_analyzer.param("HitCollection","hit_pxd")  
  pxd_analyzer.param("DigitCollection","zsdata_pxd")
  pxd_analyzer.param("AlignmentDBFileName","localDB/alignmentDB.root")
  pxd_analyzer.param("DUTPlane","3")
  pxd_analyzer.param("ReferencePlane","7")
  pxd_analyzer.param("MaxResidualU","0.2")
  pxd_analyzer.param("MaxResidualU","0.2")
  pxd_analyzer.param("RootFileName","Histos-PXD.root")
  reco_path.add_processor(pxd_analyzer)   
  
  return [ reco_path ]  
  
  
def calibrate(params):
  
  rawfile, steerfiles, gearfile, energy, caltag, useClusterDB = params
   
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
  CalObj = Calibration(steerfiles=steerfiles, name=caltag + '-cal') 
  
  # Create list of calibration paths
  calpaths = create_calibration_path(CalObj, rawfile, gearfile, energy, useClusterDB)
  
  # Run the calibration steps 
  CalObj.calibrate(paths=calpaths,ifile=rawfile,caltag=caltag)  
   
  
def reconstruct(params):
  
  rawfile, steerfiles, gearfile, energy, caltag, useClusterDB = params 
   
  # Reconsruct the rawfile using caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=caltag + '-reco' )
  
  # Create reconstuction path
  recopath = create_reco_path(RecObj, rawfile, gearfile, energy, useClusterDB)  
  
  # Run the reconstuction  
  RecObj.reconstruct(path=recopath,ifile=rawfile,caltag=caltag) 

if __name__ == '__main__':
  
  import argparse
  parser = argparse.ArgumentParser(description="Perform reconstruction of a test beam run")
  parser.add_argument('--rawfile', dest='rawfile', default='/home/benjamin/Desktop/run000020.raw', type=str, help='Location of rawfile to process')
  parser.add_argument('--gearfile', dest='gearfile', default='gear_desy_W11OF2.xml', type=str, help='Location of gearfile')
  parser.add_argument('--energy', dest='energy', default=5.0, type=float, help='Beam energy in GeV')
  parser.add_argument('--steerfiles', dest='steerfiles', default='steering-files/depfet-tb/', type=str, help='Path to steerfiles')
  parser.add_argument('--caltag', dest='caltag', default='', type=str, help='Name of calibration tag to use')
  parser.add_argument('--special', dest='special', default='', type=str, help='Added to name of new caltag to distinguish different ways to process the same run')
  args = parser.parse_args()
  
  # Use cluster calibration
  useClusterDB = True
  
  if args.caltag=='':
    print("Compute a new calibration tag directly from the rawfile {}".format(args.rawfile))
    args.caltag = os.path.splitext(os.path.basename(args.rawfile))[0] + args.special
    calibrate( (args.rawfile, args.steerfiles, args.gearfile, args.energy, args.caltag, useClusterDB) )   
  else: 
    print("Use existing caltag {}".format(args.caltag))
  
  reconstruct( (args.rawfile, args.steerfiles, args.gearfile, args.energy, args.caltag, useClusterDB) ) 
  



