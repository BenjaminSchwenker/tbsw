"""
Some helper code to define processor paths

:author: ulf.stolzenberg@phys.uni-goettinge.de  
"""

from tbsw import *




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


def add_M26coghitmaker(path):
  """
  Adds M26 Center of gravity hitmaker to the path
  """  
    
  m26coghitmaker = Processor(name="M26CogHitMaker",proctype="CogHitMaker")
  m26coghitmaker.param("ClusterCollection","zscluster_m26")
  m26coghitmaker.param("HitCollectionName","hit_m26")
  m26coghitmaker.param("SigmaU1",0.0033)
  m26coghitmaker.param("SigmaU2",0.0030)
  m26coghitmaker.param("SigmaU3",0.0050)
  m26coghitmaker.param("SigmaV1",0.0033)
  m26coghitmaker.param("SigmaV2",0.0030)
  m26coghitmaker.param("SigmaV3",0.0050)
  path.add_processor(m26coghitmaker)  

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
  trackfinder_loosecut.param("MaximumGap", 1)
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
  trackfinder_loosecut.param("MaxResidualU","0.5")
  trackfinder_loosecut.param("MaxResidualV","0.5")
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
  trackfinder_tightcut.param("MaximumGap", 1)
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



def add_M26goehitmaker(path):
  """
  Adds M26 Goe hitmaker to the path
  """  

  m26goehitmaker = Processor(name="M26GoeHitMaker",proctype="GoeHitMaker")
  m26goehitmaker.param("ClusterCollection","zscluster_m26")
  m26goehitmaker.param("HitCollectionName","hit_m26")
  m26goehitmaker.param("ClusterDBFileName","localDB/clusterDB-M26.root")
  path.add_processor(m26goehitmaker)  

  return path 



def add_clustercalibrator(path):
  """
  Adds M26 Goe hitmaker to the path
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



# Processor settings and sequence during angle reconstruction
def create_x0reco_path(Env, rawfile, gearfile, numberofevents, usesinglehitseeding):
  """
  Returns a list of tbsw path objects to reconstruct a test beam run 
  """

  reco = Env.create_path('reco')
  reco.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : numberofevents}) 
  reco.add_processor(name="RawInputProcessor", params={'FileName': rawfile})
  reco.add_processor(name="M26Unpacker")
  reco.add_processor(name="M26Clusterizer")
  reco.add_processor(name="M26GoeHitMaker")

  if usesinglehitseeding:
    reco.add_processor(name="DownstreamFinder", params={'ForwardPass_FirstPlane': -1, 'ForwardPass_SecondPlane': -1, 'SingleHitSeeding': 6 })
    reco.add_processor(name="UpstreamFinder", params={'ForwardPass_FirstPlane': -1, 'ForwardPass_SecondPlane': -1, 'SingleHitSeeding': 0 })

  else:
    reco.add_processor(name="DownstreamFinder")
    reco.add_processor(name="UpstreamFinder")
 
  reco.add_processor(name="X0Imager")

  return [ reco ]



def create_mc_x0reco_path(Env, rawfile, gearfile, numberofevents, usesinglehitseeding, beamenergy):
  """
  Returns a list of tbsw path objects to reconstruct a test beam run 
  """


  reco_path = Env.create_path('reco')
  reco_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : numberofevents, 'LCIOInputFiles': rawfile }) 

  reco_path=add_M26clusterizer(reco_path);
  reco_path=add_M26goehitmaker(reco_path)

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
def create_x0calibration_path(Env, rawfile, gearfile, nevents_cali):
  """
  Returns a list of tbsw path objects to calibrate the tracking telescope
  """
  
  hotpixelkiller = Env.create_path('hotpixelkiller')
  hotpixelkiller.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000})  
  hotpixelkiller.add_processor(name="RawInputProcessor",params={'FileName': rawfile})
  hotpixelkiller.add_processor(name="M26Unpacker")
  hotpixelkiller.add_processor(name="M26HotPixelKiller")
  
  clusterizer = Env.create_path('clusterizer')
  clusterizer.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali}) 
  clusterizer.add_processor(name="RawInputProcessor",params={'FileName': rawfile})
  clusterizer.add_processor(name="M26Unpacker")
  clusterizer.add_processor(name="M26Clusterizer")
  clusterizer.add_processor(name="LCIOOutput")

  correlator = Env.create_path('correlator')
  correlator.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" }) 
  correlator.add_processor(name="M26CogHitMaker")
  correlator.add_processor(name="RawDQM") 
  correlator.add_processor(name="TelCorrelator")
  
  kalman_aligner_1 = Env.create_path('kalman_aligner_1')
  kalman_aligner_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_1.add_processor(name="M26CogHitMaker")
  kalman_aligner_1.add_processor(name="AlignTF_LC")
  kalman_aligner_1.add_processor(name="PreAligner")
  
  kalman_aligner_2 = Env.create_path('kalman_aligner_2')
  kalman_aligner_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_2.add_processor(name="M26CogHitMaker")
  kalman_aligner_2.add_processor(name="AlignTF_TC")
  kalman_aligner_2.add_processor(name="TelAligner")
  
  telescope_dqm = Env.create_path('telescope_dqm')
  telescope_dqm.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm.add_processor(name="M26CogHitMaker")
  telescope_dqm.add_processor(name="AlignTF_TC")
  telescope_dqm.add_processor(name="TelescopeDQM")
  
  cluster_calibration_1 = Env.create_path('cluster_calibration_1')
  cluster_calibration_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': "tmp.slcio" }) 
  cluster_calibration_1.add_processor(name="M26CogHitMaker") 
  cluster_calibration_1.add_processor(name="AlignTF_TC")
  cluster_calibration_1.add_processor(name="M26ClusterCalibrator")
  
  kalman_aligner_3 = Env.create_path('kalman_aligner_3')
  kalman_aligner_3.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_3.add_processor(name="M26GoeHitMaker")
  kalman_aligner_3.add_processor(name="AlignTF_TC")
  kalman_aligner_3.add_processor(name="TelAligner")
  
  cluster_calibration_2 = Env.create_path('cluster_calibration_2')
  cluster_calibration_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': "tmp.slcio" })  
  cluster_calibration_2.add_processor(name="M26GoeHitMaker") 
  cluster_calibration_2.add_processor(name="AlignTF_TC")
  cluster_calibration_2.add_processor(name="M26ClusterCalibrator")

  correlator2 = Env.create_path('correlator2')
  correlator2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })
  correlator2.add_processor(name="M26GoeHitMaker")  
  correlator2.add_processor(name="RawDQM", params={'RootFileName': 'RawDQM2.root'})
  correlator2.add_processor(name="TelCorrelator", params={'OutputRootFileName': 'XCorrelator2.root'})
  
  kalman_aligner_4 = Env.create_path('kalman_aligner_4')
  kalman_aligner_4.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" }) 
  kalman_aligner_4.add_processor(name="M26GoeHitMaker")   
  kalman_aligner_4.add_processor(name="AlignTF_LC")
  kalman_aligner_4.add_processor(name="PreAligner")
  
  kalman_aligner_5 = Env.create_path('kalman_aligner_5')
  kalman_aligner_5.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" }) 
  kalman_aligner_5.add_processor(name="M26GoeHitMaker")  
  kalman_aligner_5.add_processor(name="AlignTF_TC")
  kalman_aligner_5.add_processor(name="TelAligner") 
  
  telescope_dqm2 = Env.create_path('telescope_dqm2')
  telescope_dqm2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': "tmp.slcio" }) 
  telescope_dqm2.add_processor(name="M26GoeHitMaker")  
  telescope_dqm2.add_processor(name="AlignTF_TC")
  telescope_dqm2.add_processor(name="TelescopeDQM", params={'RootFileName' : 'TelescopeDQM2.root'})
  
  
  # create sequence of calibration paths 
  calpath= [ hotpixelkiller ,
			 clusterizer, 
             correlator, 
             kalman_aligner_1, 
             kalman_aligner_2, 
             kalman_aligner_2, 
             kalman_aligner_2, 
             telescope_dqm, 
             cluster_calibration_1,
             kalman_aligner_3, 
             kalman_aligner_3, 
             kalman_aligner_3,  
             cluster_calibration_2, 
             cluster_calibration_2, 
             cluster_calibration_2, 
             cluster_calibration_2, 
             cluster_calibration_2, 
             cluster_calibration_2, 
             cluster_calibration_2, 
             cluster_calibration_2, 
             correlator2, 
             kalman_aligner_4, 
             kalman_aligner_5, 
             kalman_aligner_5, 
             kalman_aligner_5, 
             telescope_dqm2, 
           ]
  
  return calpath



# Processor settings and sequence during telescope calibration
def create_x0calibration_longtelescope_path(Env, rawfile, gearfile, nevents_cali, useclusterdb):
  """
  Returns a list of tbsw path objects to calibrate the tracking telescope
  """
  
  hotpixelkiller = Env.create_path('hotpixelkiller')
  hotpixelkiller.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 1000000})  
  hotpixelkiller.add_processor(name="RawInputProcessor",params={'FileName': rawfile})
  hotpixelkiller.add_processor(name="M26Unpacker")
  hotpixelkiller.add_processor(name="M26HotPixelKiller")
  
  clusterizer = Env.create_path('clusterizer')
  clusterizer.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali}) 
  clusterizer.add_processor(name="RawInputProcessor",params={'FileName': rawfile})
  clusterizer.add_processor(name="M26Unpacker")
  clusterizer.add_processor(name="M26Clusterizer")
  clusterizer.add_processor(name="LCIOOutput")

  correlator = Env.create_path('correlator')
  correlator.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 1000000, 'LCIOInputFiles': "tmp.slcio" }) 
  correlator.add_processor(name="M26CogHitMaker")
  correlator.add_processor(name="RawDQM") 
  correlator.add_processor(name="TelCorrelator")
  
  kalman_aligner_triplet_0 = Env.create_path('kalman_aligner_triplet_0')
  kalman_aligner_triplet_0.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_triplet_0.add_processor(name="M26CogHitMaker")
  kalman_aligner_triplet_0.add_processor(name="AlignTF_LC", params={'ExcludeDetector' : '3 4 5 6', 'MinimumHits' : '3' })     
  kalman_aligner_triplet_0.add_processor(name="PreAligner", params={'ErrorsShiftX' : '0 10 0 0 0 0 0', 
                                                                     'ErrorsShiftY' : '0 10 0 0 0 0 0', 
                                                                     'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                     'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                     'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                     'ErrorsGamma'  : '0 0.01 0 0 0 0 0'})
  
  kalman_aligner_triplet_1 = Env.create_path('kalman_aligner_triplet_1')
  kalman_aligner_triplet_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_triplet_1.add_processor(name="M26CogHitMaker")
  kalman_aligner_triplet_1.add_processor(name="AlignTF_TC", params={'ExcludeDetector' : '3 4 5 6', 'MinimumHits' : '3' })    
  kalman_aligner_triplet_1.add_processor(name="TelAligner", params={'ErrorsShiftX' : '0 10 0 0 0 0 0', 
                                                                     'ErrorsShiftY' : '0 10 0 0 0 0 0', 
                                                                     'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                     'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                     'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                     'ErrorsGamma'  : '0 0.01 0 0 0 0 0'})
  
  triplet_dqm = Env.create_path('triplet_dqm')
  triplet_dqm.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  triplet_dqm.add_processor(name="M26CogHitMaker")
  triplet_dqm.add_processor(name="AlignTF_TC",params={'ExcludeDetector': "3 4 5 6", 'MinimumHits': 3})
  triplet_dqm.add_processor(name="TelescopeDQM", params={'RootFileName' : 'TelescopeDQM_triplet.root'})

  tripletcorrelator = Env.create_path('tripletcorrelator')
  tripletcorrelator.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" }) 
  tripletcorrelator.add_processor(name="M26CogHitMaker")
  tripletcorrelator.add_processor(name="AlignTF_TC",params={'ExcludeDetector': "3 4 5 6", 'MinimumHits': 3})
  tripletcorrelator.add_processor(name="TriplettCorrelator", params={'OutputRootFileName':'TripletCorrelator.root'})

  kalman_aligner_quadruplet_0 = Env.create_path('kalman_aligner_quadruplet_0')
  kalman_aligner_quadruplet_0.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_quadruplet_0.add_processor(name="M26CogHitMaker")
  kalman_aligner_quadruplet_0.add_processor(name="AlignTF_LC", params={'ExcludeDetector' : '3 5 6', 'MinimumHits' : '4' })     
  kalman_aligner_quadruplet_0.add_processor(name="PreAligner", params={'ErrorsShiftX' : '0 10 10 0 0 0 0', 
                                                                       'ErrorsShiftY' : '0 10 10 0 0 0 0', 
                                                                       'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                       'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                       'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                       'ErrorsGamma'  : '0 0.01 0.01 0 0 0 0'})

  telescope_dqm_LC_quadruplet = Env.create_path('telescope_dqm_LC_quadruplet')
  telescope_dqm_LC_quadruplet.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm_LC_quadruplet.add_processor(name="M26CogHitMaker")
  telescope_dqm_LC_quadruplet.add_processor(name="AlignTF_LC",params={'ExcludeDetector': "3 5 6", 'MinimumHits': 4})
  telescope_dqm_LC_quadruplet.add_processor(name="TelescopeDQM", params={'RootFileName' : 'TelescopeDQM_LC_quadruplet.root'})
  
  kalman_aligner_quadruplet_1 = Env.create_path('kalman_aligner_quadruplet_1')
  kalman_aligner_quadruplet_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_quadruplet_1.add_processor(name="M26CogHitMaker")
  kalman_aligner_quadruplet_1.add_processor(name="AlignTF_TC", params={'ExcludeDetector' : '3 5 6', 'MinimumHits' : '4' })    
  kalman_aligner_quadruplet_1.add_processor(name="TelAligner", params={'ErrorsShiftX' : '0 10 10 0 0 0 0', 
                                                                       'ErrorsShiftY' : '0 10 10 0 0 0 0', 
                                                                       'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                       'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                       'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                       'ErrorsGamma'  : '0 0.01 0.01 0 0 0 0'})

  telescope_dqm_TC_quadruplet = Env.create_path('telescope_dqm_TC_quadruplet')
  telescope_dqm_TC_quadruplet.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm_TC_quadruplet.add_processor(name="M26CogHitMaker")
  telescope_dqm_TC_quadruplet.add_processor(name="AlignTF_TC",params={'ExcludeDetector': "3 5 6", 'MinimumHits': 4})
  telescope_dqm_TC_quadruplet.add_processor(name="TelescopeDQM", params={'RootFileName' : 'TelescopeDQM_TC_quadruplet.root'})

  quadrupletcorrelator = Env.create_path('quadrupletcorrelator')
  quadrupletcorrelator.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" }) 
  quadrupletcorrelator.add_processor(name="M26CogHitMaker")
  quadrupletcorrelator.add_processor(name="AlignTF_TC",params={'ExcludeDetector': "3 5 6", 'MinimumHits': 4})
  quadrupletcorrelator.add_processor(name="TelescopeDQM", params={'RootFileName' : 'TelescopeDQM_TC_sensor4_precorrelator.root'})
  quadrupletcorrelator.add_processor(name="TriplettCorrelator", params={'OutputRootFileName' : 'QuadrupletCorrelator.root'})

  kalman_aligner_quintet_0 = Env.create_path('kalman_aligner_quintet_0')
  kalman_aligner_quintet_0.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_quintet_0.add_processor(name="M26CogHitMaker")
  kalman_aligner_quintet_0.add_processor(name="AlignTF_LC", params={'ExcludeDetector' : '3 6', 'MinimumHits' : '5' })     
  kalman_aligner_quintet_0.add_processor(name="PreAligner", params={'ErrorsShiftX' : '0 10 10 0 10 0 0', 
                                                                    'ErrorsShiftY' : '0 10 10 0 10 0 0', 
                                                                    'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                    'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                    'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                    'ErrorsGamma'  : '0 0.01 0.01 0 0.01 0 0'})

  telescope_dqm_LC_quintet = Env.create_path('telescope_dqm_LC_quintet')
  telescope_dqm_LC_quintet.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm_LC_quintet.add_processor(name="M26CogHitMaker")
  telescope_dqm_LC_quintet.add_processor(name="AlignTF_LC",params={'ExcludeDetector': "3 6", 'MinimumHits': 5})
  telescope_dqm_LC_quintet.add_processor(name="TelescopeDQM", params={'RootFileName' : 'TelescopeDQM_LC_quintet.root'})
  
  kalman_aligner_quintet_1 = Env.create_path('kalman_aligner_quintet_1')
  kalman_aligner_quintet_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_quintet_1.add_processor(name="M26CogHitMaker")
  kalman_aligner_quintet_1.add_processor(name="AlignTF_TC", params={'ExcludeDetector' : '3 6', 'MinimumHits' : '5' })    
  kalman_aligner_quintet_1.add_processor(name="TelAligner", params={ 'ErrorsShiftX' : '0 10 10 0 10 0 0', 
                                                                     'ErrorsShiftY' : '0 10 10 0 10 0 0', 
                                                                     'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                     'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                     'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                     'ErrorsGamma'  : '0 0.01 0.01 0 0.01 0 0'})

  telescope_dqm_TC_quintet = Env.create_path('telescope_dqm_TC_quintet')
  telescope_dqm_TC_quintet.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm_TC_quintet.add_processor(name="M26CogHitMaker")
  telescope_dqm_TC_quintet.add_processor(name="AlignTF_TC",params={'ExcludeDetector': "3 6", 'MinimumHits': 5})
  telescope_dqm_TC_quintet.add_processor(name="TelescopeDQM", params={'RootFileName' : 'TelescopeDQM_TC_quintet.root'})

  quintetcorrelator = Env.create_path('quintetcorrelator')
  quintetcorrelator.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" }) 
  quintetcorrelator.add_processor(name="M26CogHitMaker")
  quintetcorrelator.add_processor(name="AlignTF_TC",params={'ExcludeDetector': "3 6", 'MinimumHits': 5})
  quintetcorrelator.add_processor(name="TriplettCorrelator", params={'OutputRootFileName' : 'QuintetCorrelator.root'})

  kalman_aligner_1 = Env.create_path('kalman_aligner_1')
  kalman_aligner_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_1.add_processor(name="M26CogHitMaker")
  kalman_aligner_1.add_processor(name="AlignTF_LC")
  kalman_aligner_1.add_processor(name="PreAligner")

  telescope_dqm_loose = Env.create_path('telescope_dqm_loose')
  telescope_dqm_loose.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm_loose.add_processor(name="M26CogHitMaker")
  telescope_dqm_loose.add_processor(name="AlignTF_LC")
  telescope_dqm_loose.add_processor(name="TelescopeDQM", params={'RootFileName' : 'TelescopeDQM_loose.root'})
  
  kalman_aligner_2 = Env.create_path('kalman_aligner_2')
  kalman_aligner_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 40000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_2.add_processor(name="M26CogHitMaker")
  kalman_aligner_2.add_processor(name="AlignTF_TC")
  kalman_aligner_2.add_processor(name="TelAligner")
  
  telescope_dqm = Env.create_path('telescope_dqm')
  telescope_dqm.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm.add_processor(name="M26CogHitMaker")
  telescope_dqm.add_processor(name="AlignTF_TC")
  telescope_dqm.add_processor(name="TelescopeDQM")

  cluster_calibration_1 = Env.create_path('cluster_calibration_1')
  cluster_calibration_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': "tmp.slcio" }) 
  cluster_calibration_1.add_processor(name="M26CogHitMaker") 
  cluster_calibration_1.add_processor(name="AlignTF_TC")
  cluster_calibration_1.add_processor(name="M26ClusterCalibrator")
  
  kalman_aligner_3 = Env.create_path('kalman_aligner_3')
  kalman_aligner_3.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_3.add_processor(name="M26GoeHitMaker")
  kalman_aligner_3.add_processor(name="AlignTF_TC")
  kalman_aligner_3.add_processor(name="TelAligner")
  
  cluster_calibration_2 = Env.create_path('cluster_calibration_2')
  cluster_calibration_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': "tmp.slcio" })  
  cluster_calibration_2.add_processor(name="M26GoeHitMaker") 
  cluster_calibration_2.add_processor(name="AlignTF_TC")
  cluster_calibration_2.add_processor(name="M26ClusterCalibrator")

  correlator2 = Env.create_path('correlator2')
  correlator2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })
  correlator2.add_processor(name="M26GoeHitMaker")  
  correlator2.add_processor(name="RawDQM", params={'RootFileName': 'RawDQM2.root'})
  correlator2.add_processor(name="TelCorrelator", params={'OutputRootFileName': 'XCorrelator2.root'})

  kalman_aligner_triplet_0_goehits = Env.create_path('kalman_aligner_triplet_0_goehits')
  kalman_aligner_triplet_0_goehits.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_triplet_0_goehits.add_processor(name="M26GoeHitMaker")
  kalman_aligner_triplet_0_goehits.add_processor(name="AlignTF_LC", params={'ExcludeDetector' : '3 4 5 6', 'MinimumHits' : '3' })     
  kalman_aligner_triplet_0_goehits.add_processor(name="PreAligner", params={'ErrorsShiftX' : '0 10 0 0 0 0 0', 
                                                                     'ErrorsShiftY' : '0 10 0 0 0 0 0', 
                                                                     'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                     'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                     'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                     'ErrorsGamma'  : '0 0.01 0 0 0 0 0'})
  
  kalman_aligner_triplet_1_goehits = Env.create_path('kalman_aligner_triplet_1_goehits')
  kalman_aligner_triplet_1_goehits.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_triplet_1_goehits.add_processor(name="M26GoeHitMaker")
  kalman_aligner_triplet_1_goehits.add_processor(name="AlignTF_TC", params={'ExcludeDetector' : '3 4 5 6', 'MinimumHits' : '3' })    
  kalman_aligner_triplet_1_goehits.add_processor(name="TelAligner", params={'ErrorsShiftX' : '0 10 0 0 0 0 0', 
                                                                     'ErrorsShiftY' : '0 10 0 0 0 0 0', 
                                                                     'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                     'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                     'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                     'ErrorsGamma'  : '0 0.01 0 0 0 0 0'})
  
  triplet_dqm_goehits = Env.create_path('triplet_dqm_goehits')
  triplet_dqm_goehits.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  triplet_dqm_goehits.add_processor(name="M26GoeHitMaker")
  triplet_dqm_goehits.add_processor(name="AlignTF_TC",params={'ExcludeDetector': "3 4 5 6", 'MinimumHits': 3})
  triplet_dqm_goehits.add_processor(name="TelescopeDQM", params={'RootFileName' : 'TelescopeDQM_triplet.root'})

  tripletcorrelator_goehits = Env.create_path('tripletcorrelator_goehits')
  tripletcorrelator_goehits.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" }) 
  tripletcorrelator_goehits.add_processor(name="M26GoeHitMaker")
  tripletcorrelator_goehits.add_processor(name="AlignTF_TC",params={'ExcludeDetector': "3 4 5 6", 'MinimumHits': 3})
  tripletcorrelator_goehits.add_processor(name="TriplettCorrelator", params={'OutputRootFileName':'TripletCorrelator.root'})

  kalman_aligner_quadruplet_0_goehits = Env.create_path('kalman_aligner_quadruplet_0_goehits')
  kalman_aligner_quadruplet_0_goehits.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_quadruplet_0_goehits.add_processor(name="M26GoeHitMaker")
  kalman_aligner_quadruplet_0_goehits.add_processor(name="AlignTF_LC", params={'ExcludeDetector' : '3 5 6', 'MinimumHits' : '4' })     
  kalman_aligner_quadruplet_0_goehits.add_processor(name="PreAligner", params={'ErrorsShiftX' : '0 10 10 0 0 0 0', 
                                                                       'ErrorsShiftY' : '0 10 10 0 0 0 0', 
                                                                       'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                       'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                       'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                       'ErrorsGamma'  : '0 0.01 0.01 0 0 0 0'})

  telescope_dqm_LC_quadruplet_goehits = Env.create_path('telescope_dqm_LC_quadruplet_goehits')
  telescope_dqm_LC_quadruplet_goehits.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm_LC_quadruplet_goehits.add_processor(name="M26GoeHitMaker")
  telescope_dqm_LC_quadruplet_goehits.add_processor(name="AlignTF_LC",params={'ExcludeDetector': "3 5 6", 'MinimumHits': 4})
  telescope_dqm_LC_quadruplet_goehits.add_processor(name="TelescopeDQM", params={'RootFileName' : 'TelescopeDQM_LC_quadruplet.root'})
  
  kalman_aligner_quadruplet_1_goehits = Env.create_path('kalman_aligner_quadruplet_1_goehits')
  kalman_aligner_quadruplet_1_goehits.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_quadruplet_1_goehits.add_processor(name="M26GoeHitMaker")
  kalman_aligner_quadruplet_1_goehits.add_processor(name="AlignTF_TC", params={'ExcludeDetector' : '3 5 6', 'MinimumHits' : '4' })    
  kalman_aligner_quadruplet_1_goehits.add_processor(name="TelAligner", params={'ErrorsShiftX' : '0 10 10 0 0 0 0', 
                                                                       'ErrorsShiftY' : '0 10 10 0 0 0 0', 
                                                                       'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                       'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                       'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                       'ErrorsGamma'  : '0 0.01 0.01 0 0 0 0'})

  telescope_dqm_TC_quadruplet_goehits = Env.create_path('telescope_dqm_TC_quadruplet_goehits')
  telescope_dqm_TC_quadruplet_goehits.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm_TC_quadruplet_goehits.add_processor(name="M26GoeHitMaker")
  telescope_dqm_TC_quadruplet_goehits.add_processor(name="AlignTF_TC",params={'ExcludeDetector': "3 5 6", 'MinimumHits': 4})
  telescope_dqm_TC_quadruplet_goehits.add_processor(name="TelescopeDQM", params={'RootFileName' : 'TelescopeDQM_TC_quadruplet.root'})

  quadrupletcorrelator_goehits = Env.create_path('quadrupletcorrelator_goehits')
  quadrupletcorrelator_goehits.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" }) 
  quadrupletcorrelator_goehits.add_processor(name="M26GoeHitMaker")
  quadrupletcorrelator_goehits.add_processor(name="AlignTF_TC",params={'ExcludeDetector': "3 5 6", 'MinimumHits': 4})
  quadrupletcorrelator_goehits.add_processor(name="TelescopeDQM", params={'RootFileName' : 'TelescopeDQM_TC_sensor4_precorrelator.root'})
  quadrupletcorrelator_goehits.add_processor(name="TriplettCorrelator", params={'OutputRootFileName' : 'QuadrupletCorrelator.root'})

  kalman_aligner_quintet_0_goehits = Env.create_path('kalman_aligner_quintet_0_goehits')
  kalman_aligner_quintet_0_goehits.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_quintet_0_goehits.add_processor(name="M26GoeHitMaker")
  kalman_aligner_quintet_0_goehits.add_processor(name="AlignTF_LC", params={'ExcludeDetector' : '3 6', 'MinimumHits' : '5' })     
  kalman_aligner_quintet_0_goehits.add_processor(name="PreAligner", params={'ErrorsShiftX' : '0 10 10 0 10 0 0', 
                                                                    'ErrorsShiftY' : '0 10 10 0 10 0 0', 
                                                                    'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                    'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                    'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                    'ErrorsGamma'  : '0 0.01 0.01 0 0.01 0 0'})

  telescope_dqm_LC_quintet_goehits = Env.create_path('telescope_dqm_LC_quintet_goehits')
  telescope_dqm_LC_quintet_goehits.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm_LC_quintet_goehits.add_processor(name="M26GoeHitMaker")
  telescope_dqm_LC_quintet_goehits.add_processor(name="AlignTF_LC",params={'ExcludeDetector': "3 6", 'MinimumHits': 5})
  telescope_dqm_LC_quintet_goehits.add_processor(name="TelescopeDQM", params={'RootFileName' : 'TelescopeDQM_LC_quintet.root'})
  
  kalman_aligner_quintet_1_goehits = Env.create_path('kalman_aligner_quintet_1_goehits')
  kalman_aligner_quintet_1_goehits.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_quintet_1_goehits.add_processor(name="M26GoeHitMaker")
  kalman_aligner_quintet_1_goehits.add_processor(name="AlignTF_TC", params={'ExcludeDetector' : '3 6', 'MinimumHits' : '5' })    
  kalman_aligner_quintet_1_goehits.add_processor(name="TelAligner", params={ 'ErrorsShiftX' : '0 10 10 0 10 0 0', 
                                                                     'ErrorsShiftY' : '0 10 10 0 10 0 0', 
                                                                     'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                     'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                     'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                     'ErrorsGamma'  : '0 0.01 0.01 0 0.01 0 0'})

  telescope_dqm_TC_quintet_goehits = Env.create_path('telescope_dqm_TC_quintet_goehits')
  telescope_dqm_TC_quintet_goehits.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm_TC_quintet_goehits.add_processor(name="M26GoeHitMaker")
  telescope_dqm_TC_quintet_goehits.add_processor(name="AlignTF_TC",params={'ExcludeDetector': "3 6", 'MinimumHits': 5})
  telescope_dqm_TC_quintet_goehits.add_processor(name="TelescopeDQM", params={'RootFileName' : 'TelescopeDQM_TC_quintet.root'})

  quintetcorrelator_goehits = Env.create_path('quintetcorrelator_goehits')
  quintetcorrelator_goehits.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 400000, 'LCIOInputFiles': "tmp.slcio" }) 
  quintetcorrelator_goehits.add_processor(name="M26GoeHitMaker")
  quintetcorrelator_goehits.add_processor(name="AlignTF_TC",params={'ExcludeDetector': "3 6", 'MinimumHits': 5})
  quintetcorrelator_goehits.add_processor(name="TriplettCorrelator", params={'OutputRootFileName' : 'QuintetCorrelator.root'})
  
  kalman_aligner_4 = Env.create_path('kalman_aligner_4')
  kalman_aligner_4.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" }) 
  kalman_aligner_4.add_processor(name="M26GoeHitMaker")   
  kalman_aligner_4.add_processor(name="AlignTF_LC")
  kalman_aligner_4.add_processor(name="PreAligner")
  
  kalman_aligner_5 = Env.create_path('kalman_aligner_5')
  kalman_aligner_5.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" }) 
  kalman_aligner_5.add_processor(name="M26GoeHitMaker")  
  kalman_aligner_5.add_processor(name="AlignTF_TC")
  kalman_aligner_5.add_processor(name="TelAligner") 
  
  telescope_dqm2 = Env.create_path('telescope_dqm2')
  telescope_dqm2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': "tmp.slcio" }) 
  telescope_dqm2.add_processor(name="M26GoeHitMaker")  
  telescope_dqm2.add_processor(name="AlignTF_TC")
  telescope_dqm2.add_processor(name="TelescopeDQM", params={'RootFileName' : 'TelescopeDQM2.root'})
  
  
  # create sequence of calibration paths 
  calpath= [ hotpixelkiller ,
			 clusterizer, 
			 correlator, 
			 kalman_aligner_triplet_0,
			 kalman_aligner_triplet_1,
			 triplet_dqm, 
			 tripletcorrelator,
			 kalman_aligner_quadruplet_0,
			 telescope_dqm_LC_quadruplet,
			 kalman_aligner_quadruplet_1,
			 telescope_dqm_TC_quadruplet,
			 quadrupletcorrelator,
			 kalman_aligner_quintet_0,
			 telescope_dqm_LC_quintet,
			 kalman_aligner_quintet_1,
			 telescope_dqm_TC_quintet,
			 quintetcorrelator,
			 kalman_aligner_1, 
			 kalman_aligner_1,
			 kalman_aligner_1,
			 kalman_aligner_1,
			 telescope_dqm_loose, 
			 kalman_aligner_2, 
			 kalman_aligner_2, 
			 kalman_aligner_2, 
			 kalman_aligner_2, 
			 kalman_aligner_2, 
			 telescope_dqm,
           ]

  calpath2=[ cluster_calibration_1,
			 kalman_aligner_3, 
			 kalman_aligner_3, 
			 kalman_aligner_3,  
			 cluster_calibration_2, 
			 cluster_calibration_2, 
			 cluster_calibration_2, 
			 cluster_calibration_2, 
			 cluster_calibration_2, 
			 cluster_calibration_2, 
			 cluster_calibration_2, 
			 cluster_calibration_2, 
			 correlator2, 
			 kalman_aligner_triplet_0_goehits,
			 kalman_aligner_triplet_1_goehits,
			 triplet_dqm_goehits, 
			 tripletcorrelator_goehits,
			 kalman_aligner_quadruplet_0_goehits,
			 telescope_dqm_LC_quadruplet_goehits,
			 kalman_aligner_quadruplet_1_goehits,
			 telescope_dqm_TC_quadruplet_goehits,
			 quadrupletcorrelator_goehits,
			 kalman_aligner_quintet_0_goehits,
			 telescope_dqm_LC_quintet_goehits,
			 kalman_aligner_quintet_1_goehits,
			 telescope_dqm_TC_quintet_goehits,
			 quintetcorrelator_goehits,
			 kalman_aligner_4, 
			 kalman_aligner_5, 
			 kalman_aligner_5, 
			 kalman_aligner_5, 
			 telescope_dqm2, 
           ]

  if useclusterdb:
    calpath.extend(calpath2)
  
  return calpath




def create_mc_x0calibration_path(Env, rawfile, gearfile, nevents_cali, beamenergy):
  """
  Returns a list of tbsw path objects to calibrate the tracking telescope
  """

  # Calibrations are organized in a sequence of calibration paths. 
  # The calibration paths are collected in a list for later execution
  calpath = []

  # Create path for noisy pixel masking
  mask_path = Env.create_path('mask_path')
  mask_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': rawfile }) 

  m26hotpixelkiller = Processor(name="M26HotPixelKiller",proctype="HotPixelKiller",)
  m26hotpixelkiller.param("InputCollectionName", "zsdata_m26")
  m26hotpixelkiller.param("MaxOccupancy", 0.001)
  m26hotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-M26.root")
  m26hotpixelkiller.param("OfflineZSThreshold", 0)
  mask_path.add_processor(m26hotpixelkiller)

  # Finished with path for noisy pixel masking
  calpath.append(mask_path)

  # Create path for detector level creation of clusters
  clusterizer_path = Env.create_path('clusterizer_path')
  clusterizer_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': rawfile }) 
  clusterizer_path=add_M26clusterizer(clusterizer_path)

  lciooutput = Processor(name="LCIOOutput",proctype="LCIOOutputProcessor")
  lciooutput.param("LCIOOutputFile","tmp.slcio")
  lciooutput.param("LCIOWriteMode","WRITE_NEW")
  clusterizer_path.add_processor(lciooutput)  

  # Finished with path for clusters
  calpath.append(clusterizer_path) 


  # Create path for m26 mc cluster calibration
  mc_clustercal_path = Env.create_path('mc_clustercal_path')
  mc_clustercal_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': "tmp.slcio" })  

  mc_clustercal_path=add_M26coghitmaker(mc_clustercal_path) 

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

  correlator_path=add_M26coghitmaker(correlator_path) 
  
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
  
  prealigner_path=add_M26coghitmaker(prealigner_path)
  prealigner_path=add_trackfinder_loosecut(prealigner_path, beamenergy)
  prealigner_path=add_prealigner(prealigner_path)
  
  # Finished with path for prealigner
  calpath.append(prealigner_path)  

  # Create path for alignment with tight cut track sample 
  aligner_path = Env.create_path('aligner_path')
  aligner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  
  aligner_path=add_M26coghitmaker(aligner_path)
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

  dqm_path=add_M26coghitmaker(dqm_path)
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

  cluster_calibration1_path=add_M26coghitmaker(cluster_calibration1_path)
  cluster_calibration1_path=add_trackfinder_tightcut(cluster_calibration1_path, beamenergy)
  cluster_calibration1_path=add_clustercalibrator(cluster_calibration1_path)

  # Finished with first part of cluster calibration
  calpath.append(cluster_calibration1_path)


  # Create path for pre alignment2 with loose cut track sample 
  prealigner2_path = Env.create_path('prealigner2_path')
  prealigner2_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" }) 

  prealigner2_path=add_M26coghitmaker(prealigner2_path)
  prealigner2_path=add_trackfinder_loosecut(prealigner2_path, beamenergy)
  prealigner2_path=add_prealigner(prealigner2_path)
  
  # Finished with path for second prealigner
  calpath.append(prealigner2_path)   


  # Create path for second part of cluster calibration
  cluster_calibration2_path = Env.create_path('cluster_calibration2_path')
  cluster_calibration2_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': "tmp.slcio" }) 

  cluster_calibration2_path=add_M26goehitmaker(cluster_calibration2_path)
  cluster_calibration2_path=add_trackfinder_tightcut(cluster_calibration2_path, beamenergy)
  cluster_calibration2_path=add_clustercalibrator(cluster_calibration2_path)

  # Finished with second part of cluster calibration
  # Repeat this 7x
  calpath.append(cluster_calibration2_path)
  calpath.append(cluster_calibration2_path)
  calpath.append(cluster_calibration2_path)
  calpath.append(cluster_calibration2_path)
  calpath.append(cluster_calibration2_path)
  calpath.append(cluster_calibration2_path)
  calpath.append(cluster_calibration2_path)


  # Create path for second correlator
  correlator2_path = Env.create_path('correlator2_path')
  correlator2_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  nevents_cali, 'LCIOInputFiles': "tmp.slcio"  }) 

  correlator2_path=add_M26goehitmaker(correlator2_path)
  
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
  calpath.append(correlator2_path) 


  # Create path for pre alignment with loose cut track sample 
  prealigner3_path = Env.create_path('prealigner3_path')
  prealigner3_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  

  prealigner3_path=add_M26goehitmaker(prealigner3_path)
  prealigner3_path=add_trackfinder_loosecut(prealigner3_path, beamenergy)
  prealigner3_path=add_prealigner(prealigner3_path)
  
  # Finished with path for prealigner
  calpath.append(prealigner3_path)  

  # Create second path for alignment with tight cut track sample 
  aligner2_path = Env.create_path('aligner2_path')
  aligner2_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  

  aligner2_path=add_M26goehitmaker(aligner2_path)
  aligner2_path=add_trackfinder_tightcut(aligner2_path, beamenergy)
  aligner2_path=add_aligner(aligner2_path)   
  
  # Finished with path for aligner
  # Repeat this 3x
  calpath.append(aligner2_path)  
  calpath.append(aligner2_path)
  calpath.append(aligner2_path)

  # Creeate second path for some track based dqm using current calibrations
  dqm2_path = Env.create_path('dqm2_path')
  dqm2_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  

  dqm2_path=add_M26goehitmaker(dqm2_path)
  dqm2_path=add_trackfinder_tightcut(dqm2_path, beamenergy)
  
  teldqm = Processor(name="TelescopeDQM", proctype="TrackFitDQM") 
  teldqm.param("AlignmentDBFileName","localDB/alignmentDB.root")
  teldqm.param("RootFileName","TelescopeDQM2.root")
  dqm2_path.add_processor(teldqm)  
  
  # Finished with path for teldqm
  calpath.append(dqm2_path)

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
  gun.param("BeamIntensity","2000")
  gun.param("BeamMomentum", str(beamenergy))
  gun.param("BeamVertexX","0")
  gun.param("BeamVertexY","0")
  gun.param("BeamVertexZ","-10")  
  gun.param("BeamVertexXSigma","10")
  gun.param("BeamVertexYSigma","10")  
  gun.param("BeamSlopeXSigma","0.0035") 
  gun.param("BeamSlopeXSigma","0.0035")  
  gun.param("PDG","11")
  gun.param("ParticleCharge","-1")
  gun.param("ParticleMass","0.000511")
  sim_air.add_processor(gun)
  sim_alu.add_processor(gun)

  fastsim = Processor(name="FastSim",proctype="FastSimulation")
  fastsim.param("AlignmentDBFileName","localDB/alignmentDB.root")
  fastsim.param("ScatterModel","0")
  fastsim.param("DoEnergyLossStraggling","true")
  fastsim.param("DoFractionalBetheHeitlerEnergyLoss","true")
  sim_air.add_processor(fastsim)
  sim_alu.add_processor(fastsim)

  m26digi = Processor(name="M26Digitizer",proctype="SiPixDigitizer")
  m26digi.param("DigitCollectionName","zsdata_m26")  
  m26digi.param("FrontEndType","1") 
  m26digi.param("ComparatorThrehold","1100")
  m26digi.param("ElectronicNoise","300")
  m26digi.param("FilterIDs","0 1 2 3 4 5")
  m26digi.param("IntegrationWindow","true")
  m26digi.param("StartIntegration","0")
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

