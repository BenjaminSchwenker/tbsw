"""
Some helper code to define processor paths

:author: ulf.stolzenberg@phys.uni-goettinge.de  
"""

from tbsw import *

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

def create_mc_x0reco_path(Env, rawfile, gearfile, numberofevents, usesinglehitseeding):
  """
  Returns a list of tbsw path objects to reconstruct a test beam run 
  """
  
  reco = Env.create_path('reco')
  reco.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : numberofevents, 'LCIOInputFiles': rawfile }) 
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


def create_mc_x0calibration_path(Env, rawfile, gearfile, nevents_cali):
  """
  Returns a list of tbsw path objects to calibrate the tracking telescope
  """
  
  hotpixelkiller = Env.create_path('hotpixelkiller')
  hotpixelkiller.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': rawfile })  
  hotpixelkiller.add_processor(name="M26HotPixelKiller")

  cluster_calibrator_mc = Env.create_path('cluster_calibrator_mc')
  cluster_calibrator_mc.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': rawfile })  
  cluster_calibrator_mc.add_processor(name="M26Clusterizer")
  cluster_calibrator_mc.add_processor(name="M26CogHitMaker")
  cluster_calibrator_mc.add_processor(name="M26ClusterCalibrationFromMC")

  clusterizer = Env.create_path('clusterizer')
  clusterizer.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_cali, 'LCIOInputFiles': rawfile }) 
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
             cluster_calibrator_mc,
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


def create_x0sim_path_air(Env, rawfile_air, rawfile_alu, gearfile, nevents_air, nevents_alu):
  """
  Returns a list of tbsw path objects to simulate a test beam run without any DUT
  """
  
  sim_air = Env.create_path('sim')
  sim_air.set_globals(params={'GearXMLFile': 'gear_air.xml' , 'MaxRecordNumber' : nevents_air})   
  sim_air.add_processor(name="InfoSetter")
  sim_air.add_processor(name="ParticleGun")
  sim_air.add_processor(name="FastSim")
  sim_air.add_processor(name="M26Digitizer")
  sim_air.add_processor(name="LCIOOutput",params={"LCIOOutputFile" : rawfile_air })

  sim_alu = Env.create_path('sim2')
  sim_alu.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_alu})   
  sim_alu.add_processor(name="InfoSetter")
  sim_alu.add_processor(name="ParticleGun")
  sim_alu.add_processor(name="FastSim")
  sim_alu.add_processor(name="M26Digitizer")
  sim_alu.add_processor(name="LCIOOutput",params={"LCIOOutputFile" : rawfile_alu })
    
  simpath = [ sim_air,
              sim_alu, 
            ]

  return simpath

