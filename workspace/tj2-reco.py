#!/usr/bin/env python
# coding: utf8
"""
Script for processing tj2 testbeam June/July 2022 at Desy. 

This script shows the calibration and reconstruction process for testbeam data.

Usage: 

python3 tj2-reco.py --rawfile /home/benjamin/desy_tb/run388.txt --runno 388

python3 histo-plotter-tj2.py --ifile=root-files/Histos-<something>.root

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw.tbsw import Simulation, Processor, Calibration, Reconstruction
import os

maxRecordNrLong  = -1
maxRecordNrShort = 200000

energy = 4.0
mass = 0.000511
useClusterDB = False

def add_rawinput(path):
  rawinput = Processor(name="CorryInputProcessor",proctype="CorryInputProcessor")
  rawinput.param('FileNames', rawfile)
  rawinput.param('SensorIDs', "0 1 2 3 4 5 22 21")
  rawinput.param('SensorNames', "MIMOSA26_0 MIMOSA26_1 MIMOSA26_2 MIMOSA26_3 MIMOSA26_4 MIMOSA26_5 Monopix2_0 FEI4_0")
  rawinput.param('RawHitCollectionName', "rawdata")
  rawinput.param("RunNumber", runno)
  path.add_processor(rawinput)

  return path

def add_unpackers(path):
  """
  Adds unpackers to the path
  """  

  m26unpacker = Processor(name="TelUnpacker", proctype="HitsFilterProcessor") 
  m26unpacker.param("InputCollectionName", "rawdata")
  m26unpacker.param("OutputCollectionName", "zsdata_m26")
  m26unpacker.param("FilterIDs", "0 1 2 3 4 5")
  path.add_processor(m26unpacker)

  fei4unpacker = Processor(name="FEI4Unpacker",proctype="HitsFilterProcessor")
  fei4unpacker.param("InputCollectionName","rawdata")
  fei4unpacker.param("OutputCollectionName", "zsdata_fei4")
  fei4unpacker.param("FilterIDs","21")
  path.add_processor(fei4unpacker)   

  tj2unpacker = Processor(name="FEI4Unpacker",proctype="HitsFilterProcessor")
  tj2unpacker.param("InputCollectionName","rawdata")
  tj2unpacker.param("OutputCollectionName", "zsdata_tj2")
  tj2unpacker.param("FilterIDs","22")
  path.add_processor(tj2unpacker)   
  
  return path

def add_pixelmaskers(path):
  """
  Adds hot/dead pixel masking processors to the path
  """

  m26hotpixelkiller = Processor(name="M26HotPixelKiller",proctype="HotPixelKiller")
  m26hotpixelkiller.param("InputCollectionName", "zsdata_m26")
  m26hotpixelkiller.param("MaxNormedOccupancy", 5)
  m26hotpixelkiller.param("MinNormedOccupancy", -1)   # deactivated
  m26hotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-M26.root")
  m26hotpixelkiller.param("OfflineZSThreshold", 0)
  m26hotpixelkiller.param("MaxOccupancy", 6e-03)
  m26hotpixelkiller.param("MaskNormalized", "false")
  m26hotpixelkiller.param("MaskAsAscii", "true")
  path.add_processor(m26hotpixelkiller)
   
  fei4hotpixelkiller = Processor(name="FEI4HotPixelKiller", proctype="HotPixelKiller")
  fei4hotpixelkiller.param("InputCollectionName", "zsdata_fei4")
  fei4hotpixelkiller.param("MaxNormedOccupancy", 5)
  fei4hotpixelkiller.param("MinNormedOccupancy", -1) 
  fei4hotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-FEI4.root")
  fei4hotpixelkiller.param("OfflineZSThreshold", 0)
  path.add_processor(fei4hotpixelkiller)
   
  tj2hotpixelkiller = Processor(name="TJ2HotPixelKiller", proctype="HotPixelKiller")
  tj2hotpixelkiller.param("InputCollectionName", "zsdata_tj2")
  tj2hotpixelkiller.param("MaxNormedOccupancy", 5)
  tj2hotpixelkiller.param("MinNormedOccupancy", 0.0)   # mask dead pixels 
  tj2hotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-TJ2.root")
  tj2hotpixelkiller.param("OfflineZSThreshold", 0)
  tj2hotpixelkiller.param("MaskNormalized", "false")
  tj2hotpixelkiller.param("MaxOccupancy", 0.01) # 6e-03)
  path.add_processor(tj2hotpixelkiller)  
  
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

  tj2clust = Processor(name="TJ2Clusterizer",proctype="PixelClusterizer")   
  tj2clust.param("NoiseDBFileName","localDB/NoiseDB-TJ2.root")
  tj2clust.param("SparseDataCollectionName","zsdata_tj2")
  tj2clust.param("ClusterCollectionName","zscluster_tj2")
  tj2clust.param("SparseClusterCut",0)
  tj2clust.param("SparseSeedCut", 0)
  tj2clust.param("SparseZSCut", 1)   
  path.add_processor(tj2clust)  
  
  return path

def add_hitmakers(path):
  """
  Adds center of gravity hitmakers to the path
  """  

  m26hitmaker = Processor(name="M26CogHitMaker",proctype="CogHitMaker")
  m26hitmaker.param("ClusterCollection","zscluster_m26")
  m26hitmaker.param("HitCollectionName","hit_m26")
  m26hitmaker.param("SigmaUCorrections", "0.698 0.31 0.315")  
  m26hitmaker.param("SigmaVCorrections", "0.698 0.31 0.315")
  path.add_processor(m26hitmaker)

  fei4hitmaker = Processor(name="FEI4CogHitMaker",proctype="CogHitMaker")
  fei4hitmaker.param("ClusterCollection","zscluster_fei4")
  fei4hitmaker.param("HitCollectionName","hit_fei4")
  fei4hitmaker.param("SigmaUCorrections", "1.0 0.5 0.3")  
  fei4hitmaker.param("SigmaVCorrections", "1.0 0.5 0.3")
  path.add_processor(fei4hitmaker)

  tj2hitmaker = Processor(name="TJ2CogHitMaker",proctype="CogHitMaker")
  tj2hitmaker.param("ClusterCollection","zscluster_tj2")
  tj2hitmaker.param("HitCollectionName","hit_tj2")
  tj2hitmaker.param("SigmaUCorrections", "0.8 0.3 0.3")  
  tj2hitmaker.param("SigmaVCorrections", "0.8 0.3 0.3")
  path.add_processor(tj2hitmaker)
  
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
  
  tj2goehitmaker = Processor(name="TJ2GoeHitMaker",proctype="GoeHitMaker")   
  tj2goehitmaker.param("ClusterCollection","zscluster_tj2")
  tj2goehitmaker.param("HitCollectionName","hit_tj2")
  tj2goehitmaker.param("ClusterDBFileName","localDB/clusterDB-TJ2.root")
  tj2goehitmaker.param("UseCenterOfGravityFallback","true")
  tj2goehitmaker.param("SigmaUCorrections", "0.8 0.3 0.3")  
  tj2goehitmaker.param("SigmaVCorrections", "0.8 0.3 0.3")
  path.add_processor(tj2goehitmaker)   
  
  return path

def add_clustercalibrators(path):
  """
  Add cluster calibration processors to create clusterDB's
  """
  
  m26clustdb = Processor(name="M26ClusterCalibrator",proctype="GoeClusterCalibrator")   
  m26clustdb.param("ClusterDBFileName","localDB/clusterDB-M26.root")  
  m26clustdb.param("MinClusters","200")
  m26clustdb.param("SelectPlanes","1 2 4 5")
  path.add_processor(m26clustdb)  
    
  tj2clustdb = Processor(name="TJ2ClusterCalibrator",proctype="GoeClusterCalibrator")   
  tj2clustdb.param("ClusterDBFileName","localDB/clusterDB-TJ2.root")  
  tj2clustdb.param("MinClusters","200")
  tj2clustdb.param("MaxEtaBins","7")
  tj2clustdb.param("SelectPlanes","3")
  path.add_processor(tj2clustdb)  
    
  fei4clustdb = Processor(name="FEI4ClusterCalibrator",proctype="GoeClusterCalibrator")   
  fei4clustdb.param("ClusterDBFileName","localDB/clusterDB-FEI4.root")  
  fei4clustdb.param("MinClusters","200")
  fei4clustdb.param("SelectPlanes","6")
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
  mask_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : maxRecordNrLong })
  
  mask_path = add_rawinput(mask_path)


  geo = Processor(name="Geo",proctype="Geometry")
  geo.param("AlignmentDBFilePath", "localDB/alignmentDB.root")
  geo.param("ApplyAlignment", "true")
  geo.param("OverrideAlignment", "true")
  mask_path.add_processor(geo)
  
  mask_path = add_unpackers(mask_path)
   
  mask_path = add_pixelmaskers(mask_path)
  
  # Add path for masking
  calpaths.append(mask_path)  

  # Create path for detector level creation of clusters
  clusterizer_path = Env.create_path('clusterizer_path')
  clusterizer_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  maxRecordNrLong})
  
  clusterizer_path = add_rawinput(clusterizer_path)

  clusterizer_path.add_processor(geo)
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
  correlator_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : maxRecordNrShort, 'LCIOInputFiles': "tmp.slcio" })
  correlator_path.add_processor(geo)
  correlator_path = add_hitmakers(correlator_path) 
  
  hitdqm = Processor(name="RawDQM",proctype="RawHitDQM")
  hitdqm.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_tj2")  
  hitdqm.param("RootFileName","RawDQM.root")
  correlator_path.add_processor(hitdqm)  
   
  correlator = Processor(name="TelCorrelator", proctype="Correlator")
  correlator.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_tj2")
  correlator.param("OutputRootFileName","XCorrelator.root")
  correlator.param("ReferencePlane","0")
  correlator.param("ParticleCharge","-1")
  correlator.param("ParticleMass", mass)
  correlator.param("ParticleMomentum", energy)
  correlator_path.add_processor(correlator)  
  
  # Finished with path for hit based pre alignment
  calpaths.append(correlator_path)  
  
  # Create path for pre alignment with loose cut track sample 
  prealigner_path = Env.create_path('prealigner_path')
  prealigner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : maxRecordNrShort, 'LCIOInputFiles': "tmp.slcio" })
  prealigner_path.add_processor(geo)
  prealigner_path = add_hitmakers(prealigner_path)
   
  trackfinder_loosecut = Processor(name="AlignTF_LC",proctype="FastTracker")
  trackfinder_loosecut.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_tj2")
  trackfinder_loosecut.param("ExcludeDetector", "")
  # HitQualitySelection==1: use cog hits and goehits for tracking
  # HitQualitySelection==0: only use goehits for tracking
  #trackfinder_loosecut.param("HitQualitySelection", 1)
  trackfinder_loosecut.param("MaxTrackChi2", 10000000)
  trackfinder_loosecut.param("MaximumGap", 1)
  trackfinder_loosecut.param("MinimumHits",6)
  trackfinder_loosecut.param("OutlierChi2Cut", 100000000)
  trackfinder_loosecut.param("ParticleCharge","-1")
  trackfinder_loosecut.param("ParticleMass",mass)
  trackfinder_loosecut.param("ParticleMomentum", energy)
  trackfinder_loosecut.param("SingleHitSeeding", "0")
  trackfinder_loosecut.param("MaxResidualU","0.5")
  trackfinder_loosecut.param("MaxResidualV","0.5")
  prealigner_path.add_processor(trackfinder_loosecut)
 
  prealigner = Processor(name="PreAligner",proctype="KalmanAligner")
  prealigner.param('ErrorsShiftX' , '0 10 10 10 10 10 10 0')
  prealigner.param('ErrorsShiftY' , '0 10 10 10 10 10 10 0')
  prealigner.param('ErrorsShiftZ' , '0 0 0 0 0 0 0 0')
  prealigner.param('ErrorsAlpha'  , '0 0 0 0 0 0 0 0')
  prealigner.param('ErrorsBeta'   , '0 0 0 0 0 0 0 0')
  prealigner.param('ErrorsGamma'  , '0 0.01 0.01 0.01 0.01 0.01 0.01 0')
  prealigner_path.add_processor(prealigner)  

  # Finished with path for prealigner
  calpaths.append(prealigner_path)  
  

  # Create path for alignment with tight cut track sample 
  alignert_path = Env.create_path('alignert_path')
  alignert_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : maxRecordNrShort, 'LCIOInputFiles': "tmp.slcio" })
  alignert_path.add_processor(geo)
  alignert_path = add_hitmakers(alignert_path)
  
  alignert_path.add_processor(trackfinder_loosecut)
  tcorrelator = Processor(name="MyTriplettCorrelator",proctype="TriplettCorrelator")
  tcorrelator.param("OutputRootFileName","TXCorrelator.root")
  tcorrelator.param("TrackCollectionName","tracks")
  tcorrelator.param("InputHitCollectionNameVec","hit_m26  hit_tj2")
  alignert_path.add_processor(tcorrelator)

  # Finished with path for triplett correlator
  calpaths.append(alignert_path) 
  calpaths.append(prealigner_path)  
  calpaths.append(prealigner_path)  

  # Create path for alignment with tight cut track sample 
  aligner_path = Env.create_path('aligner_path')
  aligner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : maxRecordNrShort, 'LCIOInputFiles': "tmp.slcio" })
  aligner_path.add_processor(geo)
  aligner_path = add_hitmakers(aligner_path)
  
  trackfinder_tightcut = Processor(name="AlignTF_TC",proctype="FastTracker")
  trackfinder_tightcut.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_tj2")
  trackfinder_tightcut.param("ExcludeDetector", "")
  # HitQualitySelection==1: use cog hits and goehits for tracking
  # HitQualitySelection==0: only use goehits for tracking
  #trackfinder_tightcut.param("HitQualitySelection", 1)
  trackfinder_tightcut.param("MaxTrackChi2", 100)
  trackfinder_tightcut.param("MaximumGap", 1)
  trackfinder_tightcut.param("MinimumHits",6)
  trackfinder_tightcut.param("OutlierChi2Cut", 20)
  trackfinder_tightcut.param("ParticleCharge","-1")
  trackfinder_tightcut.param("ParticleMass",mass)
  trackfinder_tightcut.param("ParticleMomentum", energy)
  trackfinder_tightcut.param("SingleHitSeeding", "0")
  trackfinder_tightcut.param("MaxResidualU","0.4")
  trackfinder_tightcut.param("MaxResidualV","0.4")
  aligner_path.add_processor(trackfinder_tightcut)
   
  aligner = Processor(name="Aligner",proctype="KalmanAligner")
  aligner.param('ErrorsShiftX' , '0 10 10 10 10 10 10 0' )
  aligner.param('ErrorsShiftY' , '0 10 10 10 10 10 10 0')
  aligner.param('ErrorsShiftZ' , '0 10 10 10 10 10 10 0')
  aligner.param('ErrorsAlpha'  , '0 0 0 0.01 0 0 0 0')
  aligner.param('ErrorsBeta'   , '0 0 0 0.01 0 0 0 0')
  aligner.param('ErrorsGamma'  , '0 0.01 0.01 0.01 0.01 0.01 0.01 0')
  aligner_path.add_processor(aligner)    
  
  # Finished with path for aligner
  # Repeat this 3x
  calpaths.append(aligner_path)  
  calpaths.append(aligner_path)
  calpaths.append(aligner_path)
   
  # Creeate path for some track based dqm using current calibrations
  dqm_path = Env.create_path('dqm_path')
  dqm_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : maxRecordNrShort, 'LCIOInputFiles': "tmp.slcio" })
  dqm_path.add_processor(geo)
  dqm_path = add_hitmakers(dqm_path)
  dqm_path.add_processor(trackfinder_tightcut)

  teldqm = Processor(name="TelescopeDQM", proctype="TrackFitDQM") 
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
    preclustercal_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : maxRecordNrLong, 'LCIOInputFiles': "tmp.slcio" })
    preclustercal_path.add_processor(geo)
    preclustercal_path = add_hitmakers(preclustercal_path) 
    preclustercal_path.add_processor(trackfinder_tightcut)      
    preclustercal_path = add_clustercalibrators(preclustercal_path)
    
    # Finished with path for pre cluster calibration 
    calpaths.append(preclustercal_path)
    
    # Create path for alignment with tight cut track sample and cluster DB
    aligner_db_path = Env.create_path('aligner_db_path')
    aligner_db_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : maxRecordNrShort, 'LCIOInputFiles': "tmp.slcio" })
    aligner_db_path.add_processor(geo)
    aligner_db_path = add_hitmakersDB(aligner_db_path) 
    aligner_db_path.add_processor(trackfinder_tightcut) 
    aligner_db_path.add_processor(aligner)   
    
    # Finished with path for alignemnt with hits from pre clusterDB 
    # Repeat this 2x
    for i in range(2):
      calpaths.append(aligner_db_path) 
    
    # Creeate path for next iterations for computing clusterDBs for all sensors 
    clustercal_path = Env.create_path('clustercal_path')
    clustercal_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : maxRecordNrLong, 'LCIOInputFiles': "tmp.slcio" })
    clustercal_path.add_processor(geo)
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
    dqm_db_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : maxRecordNrShort, 'LCIOInputFiles': "tmp.slcio" })
    dqm_db_path.add_processor(geo)
    dqm_db_path = add_hitmakersDB(dqm_db_path)   
    dqm_db_path.add_processor(trackfinder_tightcut)  
    
    teldqm_db = Processor(name="TelescopeDQM_DB", proctype="TrackFitDQM") 
    teldqm_db.param("RootFileName","TelescopeDQM_DB.root")
    dqm_db_path.add_processor(teldqm_db)  
    
    # Finished with path for dqm with cluster calibration
    calpaths.append(dqm_db_path)
    
  return calpaths


def create_reco_path(Env, rawfile, gearfile, energy, useClusterDB):
  """
  Returns a list of tbsw path objects for reconstruciton of a test beam run 
  """
  
  linkpixel = dict()
  linkpixel["OF"] = "765:123 765:248 767:61 767:186"
  linkpixel["OB"] = "0:61 0:186 2:123 2:248"
  linkpixel["IF"] = "0:61 0:186 2:123 2:248"
  linkpixel["IB"] = "765:123 765:248 767:61 767:186"
  
  reco_path = Env.create_path('reco_path')
  reco_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : maxRecordNrLong })
  
  reco_path = add_rawinput(reco_path)

  geo = Processor(name="Geo",proctype="Geometry")
  geo.param("AlignmentDBFilePath", "localDB/alignmentDB.root")
  geo.param("ApplyAlignment", "true")
  geo.param("OverrideAlignment", "true")
  reco_path.add_processor(geo)  
  
  # Create path for all reconstruction up to hits
  reco_path = add_unpackers(reco_path)    
  reco_path = add_clusterizers(reco_path)    
   
  if useClusterDB: 
    reco_path = add_hitmakersDB(reco_path)   
  else: 
    reco_path = add_hitmakers(reco_path) 

  trackfinder = Processor(name="TrackFinder",proctype="FastTracker")
  trackfinder.param("InputHitCollectionNameVec","hit_m26 hit_fei4")
  trackfinder.param("ExcludeDetector", "3 8")
  # HitQualitySelection==1: use cog hits and goehits for tracking
  # HitQualitySelection==0: only use goehits for tracking
  #trackfinder.param("HitQualitySelection", 1) 
  trackfinder.param("MaxTrackChi2", "100")
  trackfinder.param("MaximumGap", "1")
  trackfinder.param("MinimumHits","6")
  trackfinder.param("OutlierChi2Cut", "20")
  trackfinder.param("ParticleCharge","-1")
  trackfinder.param("ParticleMass",mass)
  trackfinder.param("ParticleMomentum", energy)
  trackfinder.param("SingleHitSeeding", "0")
  trackfinder.param("MaxResidualU","0.4")
  trackfinder.param("MaxResidualV","0.4")
  reco_path.add_processor(trackfinder)  


  tj2_analyzer = Processor(name="TJ2Analyzer",proctype="PixelDUTAnalyzer")
  tj2_analyzer.param("NoiseDBFileName","localDB/NoiseDB-TJ2.root")  # for flagging hits at hot/dead channels
  tj2_analyzer.param("HitCollection","hit_tj2")  
  tj2_analyzer.param("DigitCollection","zsdata_tj2")
  tj2_analyzer.param("DUTPlane","3")
  tj2_analyzer.param("ReferencePlane","6")
  tj2_analyzer.param("MaxResidualU","0.2")
  tj2_analyzer.param("MaxResidualV","0.2")
  tj2_analyzer.param("RootFileName","Histos-TJ2.root")
  #tj2_analyzer.param("DUTTestPixels",linkpixel[mapping])
  reco_path.add_processor(tj2_analyzer)   
  
  return [ reco_path ]  
  
  
def calibrate(params):
  
  rawfile, steerfiles, gearfile, caltag = params
   
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
  CalObj = Calibration(steerfiles=steerfiles, name=os.path.splitext(os.path.basename(rawfile))[0] + '-cal')
  #CalObj.profile = profile
  # Create list of calibration paths
  calpaths = create_calibration_path(CalObj, rawfile, gearfile, energy, useClusterDB)
  # Run the calibration steps 
  CalObj.calibrate(paths=calpaths,ifile=rawfile,caltag=caltag)
   
  
def reconstruct(params):
  
  rawfile, steerfiles, gearfile, caltag = params
   
  # Reconsruct the rawfile using caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=os.path.splitext(os.path.basename(rawfile))[0] + '-reco' )
  #RecObj.profile = profile
  # Create reconstuction path
  recopath = create_reco_path(RecObj, rawfile, gearfile, energy, useClusterDB)  
  
  # Run the reconstuction  
  RecObj.reconstruct(paths=recopath,ifile=rawfile,caltag=caltag) 

def str2bool(v):
  if v.lower() in ('yes', 'true', 'on','t', 'y', '1'):
    return True
  elif v.lower() in ('no', 'false', 'off','f', 'n', '0'):
    return False
  else:
    raise argparse.ArgumentTypeError('Boolean value expected.')



if __name__ == '__main__':

  import argparse
  parser = argparse.ArgumentParser(description="Perform calibration and reconstruction of a test beam run")
  parser.add_argument('--rawfile', dest='rawfile', default='/home/benjamin/desy_tb/run388.txt', type=str, help='Full path to rawfile to process')
  parser.add_argument('--steerfiles', dest='steerfiles', default='steering-files/desy-tb-tj2/', type=str, help='Path to steerfiles')
  parser.add_argument('--gearfile', dest='gearfile', default='gear.xml', type=str, help='Name of gearfile inside steerfiles folder')
  parser.add_argument('--runno', dest='runno', type=int, help='Run number')
  args = parser.parse_args()

  steerfiles = args.steerfiles
  gearfile = args.gearfile  
  rawfile = args.rawfile
  runno = args.runno
  

  # Make sure that we have an absolute path
  rawfile = os.path.abspath(rawfile)
  
  # Tag for calibration data
  caltag = os.path.splitext(os.path.basename(rawfile))[0] + '-test'
  
  params = ( rawfile, steerfiles, gearfile, caltag )
  
  # Calibrate the telescope and reconstruct the rawfile 
  calibrate( params )
  
  # Reconstruct the rawfile 
  reconstruct( params )  
  
