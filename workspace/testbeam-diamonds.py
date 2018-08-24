"""
Script for processing mini 3D diamond test beam data for Helge. 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *
import multiprocessing

# Path to steering files 
steerfiles = 'steering-files/diamonds-tb/'

# Path to gear file 
gearfile = 'gear_desy_dia_feonly.xml'

# Select diamond pixel type: 0 not on diamond, 1 FE only, 2 rectangular, 3 hexagonal
diamond_pixeltype = 1


# List of runs to be processed
runlist = [
            '/home/benjamin/Desktop/workspace_3dDiamond/3DDiamond/rawdata/run002261.raw',
          ]

# Number of events to be processed (put -1 to process all)
nevents = -1

# Use cluster_db flag
use_clusterDB = False


def create_calibration_path(Env, rawfile, gearfile):
  """
  Returns a list of tbsw path objects needed to calibrate the tracking telescope
  """
  

  hotpixelkiller = Env.create_path('hotpixelkiller')
  hotpixelkiller.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 300000 })  
  hotpixelkiller.add_processor(name="RawInputProcessor",params={'FileName': rawfile})
  hotpixelkiller.add_processor(name="M26Unpacker")
  hotpixelkiller.add_processor(name="FEI4Unpacker")
  hotpixelkiller.add_processor(name="DiamondUnpacker")
  hotpixelkiller.add_processor(name="DiamondRawHitSorter", params={'PixelType': diamond_pixeltype, })
  hotpixelkiller.add_processor(name="M26HotPixelKiller")
  hotpixelkiller.add_processor(name="FEI4HotPixelKiller")
  hotpixelkiller.add_processor(name="DiamondHotPixelKiller")
  
  correlator = Env.create_path('correlator')
  correlator.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 300000 }) 
  correlator.add_processor(name="RawInputProcessor",params={'FileName': rawfile})
  correlator.add_processor(name="M26Unpacker")
  correlator.add_processor(name="M26Clusterizer")
  correlator.add_processor(name="M26CogHitMaker")
  correlator.add_processor(name="DiamondUnpacker")
  correlator.add_processor(name="DiamondRawHitSorter", params={'PixelType': diamond_pixeltype, })
  correlator.add_processor(name="DiamondClusterizer")
  correlator.add_processor(name="DiamondCogHitMaker")
  correlator.add_processor(name="FEI4Unpacker")
  correlator.add_processor(name="FEI4Clusterizer")
  correlator.add_processor(name="FEI4CogHitMaker")
  correlator.add_processor(name="RawDQM")
  correlator.add_processor(name="TelCorrelator")
  correlator.add_processor(name="LCIOOutput")
  
  kalman_aligner_1 = Env.create_path('kalman_aligner_1')
  kalman_aligner_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_1.add_processor(name="M26CogHitMaker")
  kalman_aligner_1.add_processor(name="DEPCogHitMaker")
  kalman_aligner_1.add_processor(name="AlignTF_LC")
  kalman_aligner_1.add_processor(name="PreAligner", params={'ErrorsShiftX' : '0 10 10 10 10 10 10 0', 
                                                            'ErrorsShiftY' : '0 10 10 10 10 10 10 0', 
                                                            'ErrorsShiftZ' : '0 0 0 0 0 0 0 0', 
                                                            'ErrorsAlpha'  : '0 0 0 0 0 0 0 0',
                                                            'ErrorsBeta'   : '0 0 0 0 0 0 0 0', 
                                                            'ErrorsGamma'  : '0 0.01 0.01 0.01 0.01 0.01 0.01 0'})
  
  kalman_aligner_2 = Env.create_path('kalman_aligner_2')
  kalman_aligner_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_2.add_processor(name="M26CogHitMaker")
  kalman_aligner_2.add_processor(name="DEPCogHitMaker")
  kalman_aligner_2.add_processor(name="AlignTF_TC")
  kalman_aligner_2.add_processor(name="TelAligner", params={'ErrorsShiftX' : '0 10 10 10 10 10 10 0', 
                                                            'ErrorsShiftY' : '0 10 10 10 10 10 10 0', 
                                                            'ErrorsShiftZ' : '0 10 10 10 10 10 10 0', 
                                                            'ErrorsAlpha'  : '0 0 0 0 0 0 0 0',
                                                            'ErrorsBeta'   : '0 0 0 0 0 0 0 0', 
                                                            'ErrorsGamma'  : '0 0.01 0.01 0.01 0.01 0.01 0.01 0'})


  kalman_aligner_dia_1 = Env.create_path('kalman_aligner_dia_1')
  kalman_aligner_dia_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 300000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_dia_1.add_processor(name="AlignTF_Dia_LC")
  kalman_aligner_dia_1.add_processor(name="PreAligner", params={'ErrorsShiftX' : '0 0 0 10 0 0 0 0', 
                                                                    'ErrorsShiftY' : '0 0 0 10 0 0 0 0', 
                                                                    'ErrorsShiftZ' : '0 0 0 0 0 0 0 0', 
                                                                    'ErrorsAlpha'  : '0 0 0 0 0 0 0 0',
                                                                    'ErrorsBeta'   : '0 0 0 0 0 0 0 0', 
                                                                    'ErrorsGamma'  : '0 0 0 0.01 0 0 0 0'})
  
  kalman_aligner_dia_2 = Env.create_path('kalman_aligner_dia_2')
  kalman_aligner_dia_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 300000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_dia_2.add_processor(name="AlignTF_Dia_TC")
  kalman_aligner_dia_2.add_processor(name="TelAligner", params={'ErrorsShiftX' : '0 0 0 10 0 0 0 0', 
                                                                'ErrorsShiftY' : '0 0 0 10 0 0 0 0', 
                                                                'ErrorsShiftZ' : '0 0 0 10 0 0 0 0', 
                                                                'ErrorsAlpha'  : '0 0 0 0 0 0 0 0',
                                                                'ErrorsBeta'   : '0 0 0 0 0 0 0 0', 
                                                                'ErrorsGamma'  : '0 0 0 0.01 0 0 0 0'})

  telescope_dqm = Env.create_path('telescope_dqm')
  telescope_dqm.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 300000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm.add_processor(name="AlignTF_Dia_TC")
  telescope_dqm.add_processor(name="TelescopeDQM")
  
  cluster_calibration_1 = Env.create_path('cluster_calibration_1')
  cluster_calibration_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents, 'LCIOInputFiles': "tmp.slcio" })  
  cluster_calibration_1.add_processor(name="AlignTF_TC")
  cluster_calibration_1.add_processor(name="FEI4ClusterCalibrator")
  cluster_calibration_1.add_processor(name="DiamondClusterCalibrator")
  cluster_calibration_1.add_processor(name="M26ClusterCalibrator")

  kalman_aligner_3 = Env.create_path('kalman_aligner_3')
  kalman_aligner_3.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_3.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  kalman_aligner_3.add_processor(name="DiamondGoeHitMaker", params={'HitCollectionName' : 'goehit_dia' })
  kalman_aligner_3.add_processor(name="FEI4GoeHitMaker", params={'HitCollectionName' : 'goehit_fei' })
  kalman_aligner_3.add_processor(name="AlignTF_TC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_fei goehit_dia'})
  kalman_aligner_3.add_processor(name="TelAligner")
  
  cluster_calibration_2 = Env.create_path('cluster_calibration_2')
  cluster_calibration_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents, 'LCIOInputFiles': "tmp.slcio" })  
  cluster_calibration_2.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  cluster_calibration_2.add_processor(name="DiamondGoeHitMaker", params={'HitCollectionName' : 'goehit_dia' })
  cluster_calibration_2.add_processor(name="FEI4GoeHitMaker", params={'HitCollectionName' : 'goehit_fei' })
  cluster_calibration_2.add_processor(name="AlignTF_TC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_fei goehit_dia'})
  cluster_calibration_2.add_processor(name="DiamondClusterCalibrator")
  cluster_calibration_2.add_processor(name="FEI4ClusterCalibrator")
  cluster_calibration_2.add_processor(name="M26ClusterCalibrator")
    
  # create sequence of calibration paths 
  calpath= [ hotpixelkiller , 
             correlator, 
             kalman_aligner_1, 
             kalman_aligner_2, 
             kalman_aligner_2, 
             kalman_aligner_2, 
             kalman_aligner_dia_1,
             kalman_aligner_dia_2,
             kalman_aligner_dia_2,
             kalman_aligner_dia_2,    
             telescope_dqm,]
    
  if use_clusterDB: 
    calpath_db = [ cluster_calibration_1, 
                   kalman_aligner_3, 
                   kalman_aligner_3, 
                   kalman_aligner_3, 
                   cluster_calibration_2, 
                   cluster_calibration_2, 
                   cluster_calibration_2, 
                   cluster_calibration_2, 
                   cluster_calibration_2, 
                   cluster_calibration_2, 
                   kalman_aligner_3, 
                   kalman_aligner_3, 
                   kalman_aligner_3, ] 
    calpath.extend(calpath_db)
  
  return calpath


def create_reco_path(Env, rawfile, gearfile):
  """
  Returns a list of tbsw path objects for reconstruciton of a test beam run 
  """
  
  reco = Env.create_path('reco')
  reco.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents }) 
  reco.add_processor(name="RawInputProcessor",params={'FileName': rawfile})
  reco.add_processor(name="M26Unpacker")
  reco.add_processor(name="M26Clusterizer")
  reco.add_processor(name="M26CogHitMaker")
  reco.add_processor(name="DiamondUnpacker")
  reco.add_processor(name="DiamondRawHitSorter", params={'PixelType': diamond_pixeltype, })
  reco.add_processor(name="DiamondClusterizer")
  reco.add_processor(name="DiamondCogHitMaker")
  reco.add_processor(name="FEI4Unpacker")
  reco.add_processor(name="FEI4Clusterizer")
  reco.add_processor(name="FEI4CogHitMaker")
  reco.add_processor(name="RecoTF")
  reco.add_processor(name="DiamondAnalyzer")
  reco.add_processor(name="FEI4Analyzer")

  return [ reco ]

  
def calibrate_and_reconstruct(params):
  
  rawfile, gearfile = params
  
 
  # Tag for calibration data
  caltag = os.path.splitext(os.path.basename(rawfile))[0] 
  
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
  CalObj = Calibration(steerfiles=steerfiles, name=caltag + '-cal') 
  
  # Create list of calibration paths
  calpath = create_calibration_path(CalObj, rawfile, gearfile)
  
  # Run the calibration steps 
  CalObj.calibrate(path=calpath,ifile=rawfile,caltag=caltag)  
  
  # Reconsruct the rawfile using caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=caltag + '-reco' )

  # Create reconstuction path
  recopath = create_reco_path(RecObj, rawfile, gearfile)  

  # Run the reconstuction  
  RecObj.reconstruct(path=recopath,ifile=rawfile,caltag=caltag) 
  
  
if __name__ == '__main__':
  count = 1 # multiprocessing.cpu_count()
  pool = multiprocessing.Pool(processes=count)
   
  params = [ (rawfile, gearfile ) for rawfile in runlist ]
  pool.map(calibrate_and_reconstruct, params)

  
  



