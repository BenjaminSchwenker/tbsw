#!/usr/bin/env python
# coding: utf8

"""
Script for processing mini strip sensor data from test beam at DESY in 2016. Rawdata files 
are kept by Luise Poley. 

The input data consists in 'merged' lcio files with digits from two ministrip sensors 
and six m26 sensors forming a reference telescope. 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

import tbsw 
import os

# Currently not used, but can be added if needed
#import multiprocessing

maxRecordNrLong  = -1
maxRecordNrShort = 200000

def add_unpackers(path):
  """
  Adds unpackers to the path
  """  
  
  m26unpacker = tbsw.Processor(name="M26Unpacker",proctype="TelUnpacker")   
  m26unpacker.param("InputCollectionName", "original_zsdata")
  m26unpacker.param("OutputCollectionName","zsdata_m26")
  m26unpacker.param("FilterIDs", "0 1 2 3 4 5" )
  m26unpacker.param("Modulus", 4)
  path.add_processor(m26unpacker)
  
  stripunpacker = tbsw.Processor(name="H5Unpacker",proctype="DEPFETUnpacker")
  stripunpacker.param('InputCollectionName','original_zsdata')
  stripunpacker.param('OutputCollectionName','zsdata_strip')
  stripunpacker.param("FilterIDs", "6 7" )
  stripunpacker.param("Modulus", 4)
  path.add_processor(stripunpacker)   
  
  return path

def add_clusterizers(path):
  """
  Adds clusterizers to the path
  """  
    
  m26clust = tbsw.Processor(name="M26Clusterizer",proctype="PixelClusterizer")   
  m26clust.param("NoiseDBFileName","localDB/NoiseDB-M26.root")
  m26clust.param("SparseDataCollectionName","zsdata_m26")
  m26clust.param("ClusterCollectionName","zscluster_m26")
  m26clust.param("SparseClusterCut",0)
  m26clust.param("SparseSeedCut", 0)
  m26clust.param("SparseZSCut", 0)   
  path.add_processor(m26clust)  
  
  stripclust = tbsw.Processor(name="StripClusterizer",proctype="PixelClusterizer")   
  stripclust.param("NoiseDBFileName","localDB/NoiseDB-Strip.root")
  stripclust.param("SparseDataCollectionName","zsdata_strip")
  stripclust.param("ClusterCollectionName","zscluster_strip")
  stripclust.param("SparseClusterCut",0)
  stripclust.param("SparseSeedCut", 0)
  stripclust.param("SparseZSCut", 30)   
  path.add_processor(stripclust)  

  return path

def add_hitmakers(path):
  """
  Adds center of gravity hitmakers to the path
  """  

  m26hitmaker = tbsw.Processor(name="M26CogHitMaker",proctype="CogHitMaker")
  m26hitmaker.param("ClusterCollection","zscluster_m26")
  m26hitmaker.param("HitCollectionName","hit_m26")
  m26hitmaker.param("SigmaUCorrections", "0.698 0.31 0.315")  
  m26hitmaker.param("SigmaVCorrections", "0.698 0.31 0.315")
  path.add_processor(m26hitmaker)

  striphitmaker = tbsw.Processor(name="StripCogHitMaker",proctype="CogHitMaker")
  striphitmaker.param("ClusterCollection","zscluster_strip")
  striphitmaker.param("HitCollectionName","hit_strip")
  striphitmaker.param("SigmaUCorrections", "0.8 0.4 0.3")  
  striphitmaker.param("SigmaVCorrections", "1.0 1.0 1.0")
  path.add_processor(striphitmaker)   
  
  return path


def add_hitmakersDB(path):
  """
  Add cluster shape hitmakers to the path (requiring clusterDBs)
  """  
  
  m26goehitmaker = tbsw.Processor(name="M26GoeHitMaker",proctype="GoeHitMaker")   
  m26goehitmaker.param("ClusterCollection","zscluster_m26")
  m26goehitmaker.param("HitCollectionName","hit_m26")
  m26goehitmaker.param("ClusterDBFileName","localDB/clusterDB-M26.root")
  path.add_processor(m26goehitmaker)  
    
  stripid6goehitmaker = tbsw.Processor(name="StripID6GoeHitMaker",proctype="GoeHitMaker")   
  stripid6goehitmaker.param("ClusterCollection","zscluster_strip")
  stripid6goehitmaker.param("HitCollectionName","hit_stripid6")
  stripid6goehitmaker.param("ClusterDBFileName","localDB/clusterDB-StripID6.root")
  stripid6goehitmaker.param("UseCenterOfGravityFallback","true")
  stripid6goehitmaker.param("SigmaUCorrections", "0.8 0.4 0.3")  
  stripid6goehitmaker.param("SigmaVCorrections", "1.0 1.0 1.0")
  path.add_processor(stripid6goehitmaker)   

  stripid7goehitmaker = tbsw.Processor(name="StripID7GoeHitMaker",proctype="GoeHitMaker")   
  stripid7goehitmaker.param("ClusterCollection","zscluster_strip")
  stripid7goehitmaker.param("HitCollectionName","hit_stripid7")
  stripid7goehitmaker.param("ClusterDBFileName","localDB/clusterDB-StripID7.root")
  stripid7goehitmaker.param("UseCenterOfGravityFallback","true")
  stripid7goehitmaker.param("SigmaUCorrections", "0.8 0.4 0.3")  
  stripid7goehitmaker.param("SigmaVCorrections", "1.0 1.0 1.0")
  path.add_processor(stripid7goehitmaker)   

  return path

def add_clustercalibrators(path):
  """
  Add cluster calibration processors to create clusterDB's
  """
  
  m26clustdb = tbsw.Processor(name="M26ClusterCalibrator",proctype="GoeClusterCalibrator")   
  m26clustdb.param("ClusterDBFileName","localDB/clusterDB-M26.root")  
  m26clustdb.param("MinClusters","200")
  m26clustdb.param("IgnoreIDs","6 7 21")
  path.add_processor(m26clustdb)  
   
  stripid6clustdb = tbsw.Processor(name="StripID6ClusterCalibrator",proctype="GoeClusterCalibrator")   
  stripid6clustdb.param("ClusterDBFileName","localDB/clusterDB-StripID6.root")  
  stripid6clustdb.param("MinClusters","200")
  stripid6clustdb.param("IgnoreIDs","0 1 2 3 4 5 7")
  stripid6clustdb.param("MaxEtaBins","7")
  path.add_processor(stripid6clustdb)  

  stripid7clustdb = tbsw.Processor(name="StripID7ClusterCalibrator",proctype="GoeClusterCalibrator")   
  stripid7clustdb.param("ClusterDBFileName","localDB/clusterDB-StripID7.root")  
  stripid7clustdb.param("MinClusters","200")
  stripid7clustdb.param("IgnoreIDs","0 1 2 3 4 5 6")
  stripid7clustdb.param("MaxEtaBins","7")
  path.add_processor(stripid7clustdb)  
  
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
  mask_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : maxRecordNrLong, 'LCIOInputFiles': rawfile })
  
  geo = tbsw.Processor(name="Geo",proctype="Geometry")
  geo.param("AlignmentDBFilePath", "localDB/alignmentDB.root")
  geo.param("ApplyAlignment", "false")
  geo.param("OverrideAlignment", "true")
  mask_path.add_processor(geo)
  
  mask_path = add_unpackers(mask_path)
  
  m26hotpixelkiller = tbsw.Processor(name="M26HotPixelKiller",proctype="HotPixelKiller")
  m26hotpixelkiller.param("InputCollectionName", "zsdata_m26")
  m26hotpixelkiller.param("MaxOccupancy", 0.001)
  m26hotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-M26.root")
  m26hotpixelkiller.param("OfflineZSThreshold", 0)
  mask_path.add_processor(m26hotpixelkiller)

  striphotpixelkiller = tbsw.Processor(name="StripHotPixelKiller",proctype="HotPixelKiller")
  striphotpixelkiller.param("InputCollectionName", "zsdata_strip")
  striphotpixelkiller.param("MaxOccupancy", 0.5)
  striphotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-Strip.root")
  striphotpixelkiller.param("OfflineZSThreshold", 0)
  mask_path.add_processor(striphotpixelkiller)
  
  # Add path for masking
  calpaths.append(mask_path)   
  
  # Create path for detector level creation of clusters
  clusterizer_path = Env.create_path('clusterizer_path')
  clusterizer_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  maxRecordNrLong, 'LCIOInputFiles': rawfile})
  
  clusterizer_path.add_processor(geo)
  clusterizer_path = add_unpackers(clusterizer_path) 
  clusterizer_path = add_clusterizers(clusterizer_path)    
   
  lciooutput = tbsw.Processor(name="LCIOOutput",proctype="LCIOOutputProcessor")
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
  
  hitdqm = tbsw.Processor(name="RawDQM",proctype="RawHitDQM")
  hitdqm.param("InputHitCollectionNameVec","hit_m26 hit_strip")  
  hitdqm.param("RootFileName","RawDQM.root")
  correlator_path.add_processor(hitdqm)  
   
  correlator = tbsw.Processor(name="TelCorrelator", proctype="Correlator")
  correlator.param("InputHitCollectionNameVec","hit_m26 hit_strip")
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
  prealigner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : maxRecordNrShort, 'LCIOInputFiles': "tmp.slcio" })
  prealigner_path.add_processor(geo)
  prealigner_path = add_hitmakers(prealigner_path)
   
  trackfinder_loosecut = tbsw.Processor(name="AlignTF_LC",proctype="FastTracker")
  trackfinder_loosecut.param("InputHitCollectionNameVec","hit_m26 hit_strip")
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
  trackfinder_loosecut.param("MaxResidualV","0.5 0.5 0.5 -1 -1 0.5 0.5 0.5")
  prealigner_path.add_processor(trackfinder_loosecut)
 
  prealigner = tbsw.Processor(name="PreAligner",proctype="KalmanAligner")
  prealigner.param('ErrorsShiftX' , '0 10 10 10 10 10 10 0')
  prealigner.param('ErrorsShiftY' , '0 10 10  0  0 10 10 0')
  prealigner.param('ErrorsShiftZ' , '0 0 0 0 0 0 0 0')
  prealigner.param('ErrorsAlpha'  , '0 0 0 0 0 0 0 0')
  prealigner.param('ErrorsBeta'   , '0 0 0 0 0 0 0 0')
  prealigner.param('ErrorsGamma'  , '0 0.01 0.01 0.01 0.01 0.01 0.01 0')
  prealigner_path.add_processor(prealigner)  

  # Finished with path for prealigner
  calpaths.append(prealigner_path)  
  
  # Create path for alignment with tight cut track sample 
  aligner_path = Env.create_path('aligner_path')
  aligner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : maxRecordNrShort, 'LCIOInputFiles': "tmp.slcio" })
  aligner_path.add_processor(geo)
  aligner_path = add_hitmakers(aligner_path)
  
  trackfinder_tightcut = tbsw.Processor(name="AlignTF_TC",proctype="FastTracker")
  trackfinder_tightcut.param("InputHitCollectionNameVec","hit_m26 hit_strip")
  trackfinder_tightcut.param("ExcludeDetector", "")
  trackfinder_tightcut.param("MaxTrackChi2", 100)
  trackfinder_tightcut.param("MaximumGap", 1)
  trackfinder_tightcut.param("MinimumHits",7)
  trackfinder_tightcut.param("OutlierChi2Cut", 20)
  trackfinder_tightcut.param("ParticleCharge","-1")
  trackfinder_tightcut.param("ParticleMass","0.000511")
  trackfinder_tightcut.param("ParticleMomentum", energy)
  trackfinder_tightcut.param("SingleHitSeeding", "0")
  trackfinder_tightcut.param("MaxResidualU","0.5")
  trackfinder_tightcut.param("MaxResidualV","0.5 0.5 0.5 -1 -1 0.5 0.5 0.5")
  aligner_path.add_processor(trackfinder_tightcut)
   
  aligner = tbsw.Processor(name="Aligner",proctype="KalmanAligner")
  aligner.param('ErrorsShiftX' , '0 10 10 10 10 10 10 0')
  aligner.param('ErrorsShiftY' , '0 10 10  0  0 10 10 0')
  aligner.param('ErrorsShiftZ' , '0 10 10 10 10 10 10 0')
  aligner.param('ErrorsAlpha'  , '0 0 0 0 0 0 0 0')
  aligner.param('ErrorsBeta'   , '0 0 0 0 0 0 0 0')
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

  teldqm = tbsw.Processor(name="TelescopeDQM", proctype="TrackFitDQM") 
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
    
    teldqm_db = tbsw.Processor(name="TelescopeDQM_DB", proctype="TrackFitDQM") 
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
  reco_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : maxRecordNrLong, 'LCIOInputFiles': rawfile })
  
  geo = tbsw.Processor(name="Geo",proctype="Geometry")
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
  
  trackfinder = tbsw.Processor(name="TrackFinder",proctype="FastTracker")
  trackfinder.param("InputHitCollectionNameVec","hit_m26")
  trackfinder.param("ExcludeDetector", "3 4")
  trackfinder.param("MaxTrackChi2", "100")
  trackfinder.param("MaximumGap", "1")
  trackfinder.param("MinimumHits","5")
  trackfinder.param("OutlierChi2Cut", "20")
  trackfinder.param("ParticleCharge","-1")
  trackfinder.param("ParticleMass","0.000511")
  trackfinder.param("ParticleMomentum", energy)
  trackfinder.param("SingleHitSeeding", "0")
  trackfinder.param("MaxResidualU","0.4")
  trackfinder.param("MaxResidualV","0.4")
  reco_path.add_processor(trackfinder)  

  if useClusterDB: 
    hit_collection_id6 = "hit_stripid6"
    hit_collection_id7 = "hit_stripid7"
  else: 
    hit_collection_id6 = "hit_strip"
    hit_collection_id7 = "hit_strip" 
    

  stripid6_analyzer = tbsw.Processor(name="StripID6Analyzer",proctype="PixelDUTAnalyzer")
  stripid6_analyzer.param("HitCollection",hit_collection_id6)  
  stripid6_analyzer.param("DUTPlane","3")
  stripid6_analyzer.param("MaxResidualU","0.2")
  stripid6_analyzer.param("MaxResidualV","-1")
  stripid6_analyzer.param("RootFileName","Histos-ITK-ID6.root")
  reco_path.add_processor(stripid6_analyzer)    
  
  stripid7_analyzer = tbsw.Processor(name="StripID7Analyzer",proctype="PixelDUTAnalyzer")
  stripid7_analyzer.param("HitCollection",hit_collection_id7)  
  stripid7_analyzer.param("DUTPlane","4")
  stripid7_analyzer.param("MaxResidualU","0.2")
  stripid7_analyzer.param("MaxResidualV","-1")
  stripid7_analyzer.param("RootFileName","Histos-ITK-ID7.root")
  reco_path.add_processor(stripid7_analyzer)   
  
  return [ reco_path ]   

def calibrate(params,profile):
  
  rawfile, steerfiles, gearfile, energy, caltag, useClusterDB = params
   
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
  CalObj = tbsw.Calibration(steerfiles=steerfiles, name=os.path.splitext(os.path.basename(rawfile))[0] + '-cal')
  CalObj.profile = profile
  # Create list of calibration paths
  calpaths = create_calibration_path(CalObj, rawfile, gearfile, energy, useClusterDB)
  # Run the calibration steps 
  CalObj.calibrate(paths=calpaths,ifile=rawfile,caltag=caltag) 


def reconstruct(params,profile):
  
  rawfile, steerfiles, gearfile, energy, caltag, useClusterDB = params 
   
  # Reconsruct the rawfile using caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = tbsw.Reconstruction(steerfiles=steerfiles, name=os.path.splitext(os.path.basename(rawfile))[0] + '-reco' )
  RecObj.profile = profile
  # Create reconstuction path
  recopath = create_reco_path(RecObj, rawfile, gearfile, energy, useClusterDB)  
  
  # Run the reconstuction  
  RecObj.reconstruct(paths=recopath,ifile=rawfile,caltag=caltag) 
  

if __name__ == '__main__':
   
  import argparse
  parser = argparse.ArgumentParser(description="Perform calibration and reconstruction of a test beam run")
  parser.add_argument('--rawfile', dest='rawfile', default='/my/path/data_from_kristin/run000282-merger.slcio', type=str, help='Location of rawfile to process')
  parser.add_argument('--gearfile', dest='gearfile', default='gear_desy_goeid1.xml', type=str, help='Location of gearfile')
  parser.add_argument('--energy', dest='energy', default=4.6, type=float, help='Beam energy in GeV')
  parser.add_argument('--steerfiles', dest='steerfiles', default='steering-files/luise-tb/', type=str, help='Path to steerfiles')
  parser.add_argument('--caltag', dest='caltag', default='', type=str, help='Name of calibration tag to use')
  parser.add_argument('--useClusterDB', dest='use_cluster_db', default=True, type=bool, help="Use cluster database")
  parser.add_argument('--skipCalibration', dest='skip_calibration', default=False, type=bool, help="Skip creating a new calibration tag")
  parser.add_argument('--skipReconstruction', dest='skip_reco', default=False, type=bool, help="Skip reconstruction of run")
  parser.add_argument('--profile', dest='profile', action='store_true', help="profile execution time")

  args = parser.parse_args()
  
  if args.caltag=='':
    args.caltag = os.path.splitext(os.path.basename(args.rawfile))[0]
    
  if not args.skip_calibration: 
    print("Creating new calibration tag {} from run {}".format(args.caltag, args.rawfile))
    calibrate((args.rawfile, args.steerfiles, args.gearfile, args.energy, args.caltag, args.use_cluster_db, args.mapping),args.profile)
  
  if not args.skip_reco: 
    print("Reconstruct run {} using caltag {}".format(args.rawfile, args.caltag))
    reconstruct((args.rawfile, args.steerfiles, args.gearfile, args.energy, args.caltag, args.use_cluster_db, args.mapping),args.profile)

