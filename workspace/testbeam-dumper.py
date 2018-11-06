"""
Script for processing mini 3D diamond test beam data for Helge. 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *
import multiprocessing



def create_calibration_path(Env, rawfile, gearfile, energy):
  """
  Returns a list of tbsw path objects needed to calibrate the tracking telescope
  """
  
  hotpixelkiller = Env.create_path('hotpixelkiller')
  hotpixelkiller.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1 })  
  hotpixelkiller.add_processor(name="RawInputProcessor",params={'FileName': rawfile})
  hotpixelkiller.add_processor(name="M26Unpacker")
  hotpixelkiller.add_processor(name="FEI4Unpacker")
  hotpixelkiller.add_processor(name="DEPBIGUnpacker",params={'InputCollectionName': 'DEPFET', 'OutputCollectionName':'zsdata_dep_big'})
  hotpixelkiller.add_processor(name="DEPH5Unpacker",params={'InputCollectionName': 'DEPFE5', 'OutputCollectionName':'zsdata_dep_h5',"Mapping":"Hybrid5"})
  hotpixelkiller.add_processor(name="M26HotPixelKiller")
  hotpixelkiller.add_processor(name="FEI4HotPixelKiller")
  hotpixelkiller.add_processor(name="DEPBIGHotPixelKiller",params={"InputCollectionName":"zsdata_dep_big",
                                                             "MaxOccupancy":"0.001",
                                                             "NoiseDBFileName":"localDB/NoiseDB-DEP-BIG.root",
                                                             "OfflineZSThreshold":"0"})
  hotpixelkiller.add_processor(name="DEPH5HotPixelKiller",params={"InputCollectionName":"zsdata_dep_h5",
                                                             "MaxOccupancy":"0.001",
                                                             "NoiseDBFileName":"localDB/NoiseDB-DEP-H5.root",
                                                             "OfflineZSThreshold":"0"})
 
  
  hitmaker = Env.create_path('hitmaker')
  hitmaker.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 200000 }) 
  hitmaker.add_processor(name="RawInputProcessor",params={'FileName': rawfile})
  hitmaker.add_processor(name="M26Unpacker")
  hitmaker.add_processor(name="M26Clusterizer")
  hitmaker.add_processor(name="M26CogHitMaker")
  hitmaker.add_processor(name="FEI4Unpacker")
  hitmaker.add_processor(name="FEI4Clusterizer")
  hitmaker.add_processor(name="FEI4CogHitMaker")
  hitmaker.add_processor(name="DEPBIGUnpacker",params={'InputCollectionName': 'DEPFET', 'OutputCollectionName':'zsdata_dep_big'})
  hitmaker.add_processor(name="DEPBIGClusterizer",params={"SparseDataCollectionName":"zsdata_dep_big",
                                                    "NoiseDBFileName":"localDB/NoiseDB-DEP-BIG.root",
                                                    "ClusterCollectionName":"zscluster_dep_big",
                                                    "SparseClusterCut":"5",
                                                    "SparseSeedCut":"5",
                                                    "SparseZSCut":"5",})
  hitmaker.add_processor(name="DEPBIGCogHitMaker",params={"ClusterCollection":"zscluster_dep_big",
                                                    "HitCollectionName":"hit_dep_big",
                                                    "SigmaU1":"0.0134",
                                                    "SigmaU2":"0.0077",
                                                    "SigmaU3":"0.0077",
                                                    "SigmaV1":"0.024",
                                                    "SigmaV2":"0.014",
                                                    "SigmaV3":"0.014",})  
   
  hitmaker.add_processor(name="DEPH5Unpacker",params={'InputCollectionName': 'DEPFE5', 'OutputCollectionName':'zsdata_dep_h5', "Mapping":"Hybrid5"})
  hitmaker.add_processor(name="DEPH5Clusterizer",params={"SparseDataCollectionName":"zsdata_dep_h5",
                                                    "NoiseDBFileName":"localDB/NoiseDB-DEP-H5.root",
                                                    "ClusterCollectionName":"zscluster_dep_h5",
                                                    "SparseClusterCut":"5",
                                                    "SparseSeedCut":"5",
                                                    "SparseZSCut":"5",})
  hitmaker.add_processor(name="DEPH5CogHitMaker",params={"ClusterCollection":"zscluster_dep_h5",
                                                    "HitCollectionName":"hit_dep_h5",
                                                    "SigmaU1":"0.0134",
                                                    "SigmaU2":"0.0077",
                                                    "SigmaU3":"0.0077",
                                                    "SigmaV1":"0.0134",
                                                    "SigmaV2":"0.0077",
                                                    "SigmaV3":"0.0077",})  
  hitmaker.add_processor(name="RawDQM",params={"InputHitCollectionNameVec":"hit_m26 hit_fei4 hit_dep_big hit_dep_h5", "RootFileName":"RawDQM.root"})
  hitmaker.add_processor(name="TelCorrelator",params={"InputHitCollectionNameVec":"hit_m26 hit_fei4 hit_dep_big hit_dep_h5",
                                                      "AlignmentDBFileName":"localDB/alignmentDB.root",
                                                      "OutputRootFileName":"XCorrelator.root",
                                                      "ParticleMomentum": energy })
  hitmaker.add_processor(name="LCIOOutput")

  kalman_aligner_1 = Env.create_path('kalman_aligner_1')
  kalman_aligner_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 200000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_1.add_processor(name="AlignTF_LC",params={"InputHitCollectionNameVec":"hit_m26 hit_fei4 hit_dep_big hit_dep_h5",
                                                           "ExcludeDetector": "",
                                                           "MaxTrackChi2": "10000000",
                                                           "MaximumGap": "1",
                                                           "MinimumHits":"7",
                                                           "OutlierChi2Cut": "100000000",
                                                           "ParticleMomentum": energy,
                                                           "SingleHitSeeding":"0 1"  
                                                            })
  kalman_aligner_1.add_processor(name="PreAligner", params={'ErrorsShiftX' : '0 10 10 10 10 10 0 10 10', 
                                                            'ErrorsShiftY' : '0 10 10 10 10 10 0 10 10', 
                                                            'ErrorsShiftZ' : '0 0 0 0 0 0 0 0 0', 
                                                            'ErrorsAlpha'  : '0 0 0 0 0 0 0 0 0',
                                                            'ErrorsBeta'   : '0 0 0 0 0 0 0 0 0', 
                                                            'ErrorsGamma'  : '0 0.01 0.01 0.01 0.01 0.01 0 0.01 0.01'})

  
     
   
  kalman_aligner_2 = Env.create_path('kalman_aligner_2')
  kalman_aligner_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 200000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_2.add_processor(name="AlignTF_TC",params={"InputHitCollectionNameVec":"hit_m26 hit_fei4 hit_dep_big hit_dep_h5",
                                                           "ExcludeDetector": "",
                                                           "MaxTrackChi2": "100",
                                                           "MaximumGap": "1",
                                                           "MinimumHits":"7",
                                                           "OutlierChi2Cut": "20",
                                                           "ParticleMomentum": energy ,
                                                           "SingleHitSeeding":"0 1"  
                                                            })
  kalman_aligner_2.add_processor(name="TelAligner", params={'ErrorsShiftX' : '0 10 10 10 10 10 0 10 10', 
                                                            'ErrorsShiftY' : '0 10 10 10 10 10 0 10 10', 
                                                            'ErrorsShiftZ' : '0 10 10 10 10 10 0 10 10', 
                                                            'ErrorsAlpha'  : '0 0 0 0.01 0 0 0 0.01 0.01',
                                                            'ErrorsBeta'   : '0 0 0 0.01 0 0 0 0.01 0.01', 
                                                            'ErrorsGamma'  : '0 0.01 0.01 0.01 0.01 0.01 0 0.01 0.01'})

  telescope_dqm = Env.create_path('telescope_dqm')
  telescope_dqm.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 200000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm.add_processor(name="AlignTF_TC",params={"InputHitCollectionNameVec":"hit_m26 hit_fei4 hit_dep_big hit_dep_h5",
                                                           "ExcludeDetector": "",
                                                           "MaxTrackChi2": "100",
                                                           "MaximumGap": "1",
                                                           "MinimumHits":"7",
                                                           "OutlierChi2Cut": "20",
                                                           "ParticleMomentum": energy,
                                                           "SingleHitSeeding":"0 1"  
                                                            })
  telescope_dqm.add_processor(name="TelescopeDQM",params={"RootFileName":"TelescopeDQM.root"}) 
  
  # create sequence of calibration paths 
  telescope_path = [ hotpixelkiller, 
                     hitmaker, 
                     kalman_aligner_1, 
                     kalman_aligner_2, 
                     kalman_aligner_2, 
                     kalman_aligner_2, 
                     telescope_dqm, 
                   ]
    
  calpath = []
  calpath.extend(telescope_path) 
  
  
  return calpath


def create_reco_path(Env, rawfile, gearfile, energy):
  """
  Returns a list of tbsw path objects for reconstruciton of a test beam run 
  """
  
  reco = Env.create_path('reco')
  reco.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1 }) 
  reco.add_processor(name="RawInputProcessor",params={'FileName': rawfile})
  reco.add_processor(name="M26Unpacker")
  reco.add_processor(name="M26Clusterizer")
  reco.add_processor(name="M26CogHitMaker")
  reco.add_processor(name="FEI4Unpacker")
  reco.add_processor(name="FEI4Clusterizer")
  reco.add_processor(name="FEI4CogHitMaker")
  reco.add_processor(name="DEPBIGUnpacker",params={'InputCollectionName': 'DEPFET', 'OutputCollectionName':'zsdata_dep_big'})
  reco.add_processor(name="DEPBIGClusterizer",params={"SparseDataCollectionName":"zsdata_dep_big",
                                                    "NoiseDBFileName":"localDB/NoiseDB-DEP-BIG.root",
                                                    "ClusterCollectionName":"zscluster_dep_big",
                                                    "SparseClusterCut":"5",
                                                    "SparseSeedCut":"5",
                                                    "SparseZSCut":"5",})
  reco.add_processor(name="DEPBIGCogHitMaker",params={"ClusterCollection":"zscluster_dep_big",
                                                    "HitCollectionName":"hit_dep_big",
                                                    "SigmaU1":"0.0134",
                                                    "SigmaU2":"0.0077",
                                                    "SigmaU3":"0.0077",
                                                    "SigmaV1":"0.024",
                                                    "SigmaV2":"0.014",
                                                    "SigmaV3":"0.014",})    
  
  reco.add_processor(name="DEPH5Unpacker",params={'InputCollectionName': 'DEPFE5', 'OutputCollectionName':'zsdata_dep_h5', "Mapping":"Hybrid5"})
  reco.add_processor(name="DEPH5Clusterizer",params={"SparseDataCollectionName":"zsdata_dep_h5",
                                                    "NoiseDBFileName":"localDB/NoiseDB-DEP-H5.root",
                                                    "ClusterCollectionName":"zscluster_dep_h5",
                                                    "SparseClusterCut":"5",
                                                    "SparseSeedCut":"5",
                                                    "SparseZSCut":"5",})
  reco.add_processor(name="DEPH5CogHitMaker",params={"ClusterCollection":"zscluster_dep_h5",
                                                    "HitCollectionName":"hit_dep_h5",
                                                    "SigmaU1":"0.0134",
                                                    "SigmaU2":"0.0077",
                                                    "SigmaU3":"0.0077",
                                                    "SigmaV1":"0.0134",
                                                    "SigmaV2":"0.0077",
                                                    "SigmaV3":"0.0077",})  
  
  reco.add_processor(name="RecoTF",params={"InputHitCollectionNameVec":"hit_m26 hit_fei4",
                                                           "ExcludeDetector": "3 8",
                                                           "MaxTrackChi2": "100",
                                                           "MaximumGap": "1",
                                                           "MinimumHits":"6",
                                                           "OutlierChi2Cut": "20",
                                                           "ParticleMomentum": energy,
                                                           "SingleHitSeeding":"0 1"  
                                                            })
   
  reco.add_processor(name="DEPH5Analyzer",params={"HitCollection":"hit_dep_h5",
                                                  "DigitCollection":"zsdata_dep_h5",
                                                  "DUTPlane":8,
                                                  "ReferencePlane":"7",
                                                  "MaxResidualU":0.2,
                                                  "MaxResidualV":0.2,
                                                  "RootFileName":"Histos-DEPH5.root"})

  reco.add_processor(name="DEPBIGAnalyzer",params={"HitCollection":"hit_dep_big",
                                                   "DigitCollection":"zsdata_dep_big",
                                                   "DUTPlane":3,
                                                   "ReferencePlane":"7",
                                                   "MaxResidualU":0.2,
                                                   "MaxResidualV":0.2,
                                                   "RootFileName":"Histos-DEPBIG.root"})
 
  reco.add_processor(name="FEI4Analyzer",params={"HitCollection":"hit_fei4",
                                                 "DigitCollection":"zsdata_fei4",
                                                 "DUTPlane":7,
                                                 "ReferencePlane":"7",
                                                 "MaxResidualU":0.2,
                                                 "MaxResidualV":0.2,
                                                 "RootFileName":"Histos-FEI4.root"})  
   
  return [ reco ]  

  
def calibrate(params):
  
  rawfile, steerfiles, gearfile, energy, caltag = params
   
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
  CalObj = Calibration(steerfiles=steerfiles, name=caltag + '-cal') 
  
  # Create list of calibration paths
  calpath = create_calibration_path(CalObj, rawfile, gearfile, energy)
  
  # Run the calibration steps 
  CalObj.calibrate(path=calpath,ifile=rawfile,caltag=caltag)  
   
 
def reconstruct(params):
  
  rawfile, steerfiles, gearfile, energy, caltag = params 
   
  # Reconsruct the rawfile using caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=caltag + '-reco' )
  
  # Create reconstuction path
  recopath = create_reco_path(RecObj, rawfile, gearfile, energy)  
  
  # Run the reconstuction  
  RecObj.reconstruct(path=recopath,ifile=rawfile,caltag=caltag) 

if __name__ == '__main__':
  
  import argparse
  parser = argparse.ArgumentParser(description="Perform reconstruction of a test beam run")
  parser.add_argument('--rawfile', dest='rawfile', default='/home/benjamin/Desktop/run000020.raw', type=str, help='Location of rawfile to process')
  parser.add_argument('--gearfile', dest='gearfile', default='gear_desy_v1.xml', type=str, help='Location of gearfile')
  parser.add_argument('--energy', dest='energy', default=5.0, type=float, help='Beam energy in GeV')
  parser.add_argument('--steerfiles', dest='steerfiles', default='steering-files/depfet-tb/', type=str, help='Path to steerfiles')
  parser.add_argument('--caltag', dest='caltag', default='', type=str, help='Name of calibration tag to use')
  args = parser.parse_args()
  
  if args.caltag=='':
    print("Compute a new calibration tag directly from the rawfile {}".format(args.rawfile))
    args.caltag = os.path.splitext(os.path.basename(args.rawfile))[0] 
    calibrate( (args.rawfile, args.steerfiles, args.gearfile, args.energy, args.caltag) )   
  else: 
    print("Use existing caltag {}".format(args.caltag))
  
  reconstruct( (args.rawfile, args.steerfiles, args.gearfile, args.energy, args.caltag) ) 
  



