"""
Script for processing mini 3D diamond test beam data for Helge. 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *
import multiprocessing

# Path to steering files 
steerfiles = 'steering-files/diamonds-tb/'

# List of runs to be processed
runlist = [
            '/home/benjamin/Desktop/3DDiamond/rawdata/run002261.raw',
          ]

def create_calibration_path(Env, rawfile, gearfile, diamond_pixeltype):
  """
  Returns a list of tbsw path objects needed to calibrate the tracking telescope
  """

  if diamond_pixeltype == 1: 
    uCellPeriod = 2
    vCellPeriod = 1
  elif diamond_pixeltype == 2:
    uCellPeriod = 5
    vCellPeriod = 1
  elif diamond_pixeltype == 3:
    uCellPeriod = 4
    vCellPeriod = 2
  
  hotpixelkiller = Env.create_path('hotpixelkiller')
  hotpixelkiller.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1 })  
  hotpixelkiller.add_processor(name="RawInputProcessor",params={'FileName': rawfile})
  hotpixelkiller.add_processor(name="M26Unpacker")
  hotpixelkiller.add_processor(name="FEI4Unpacker")
  hotpixelkiller.add_processor(name="DiamondUnpacker")
  hotpixelkiller.add_processor(name="DiamondRawHitSorter", params={'PixelType': diamond_pixeltype, })
  hotpixelkiller.add_processor(name="M26HotPixelKiller")
  hotpixelkiller.add_processor(name="FEI4HotPixelKiller")
  hotpixelkiller.add_processor(name="DiamondHotPixelKiller")
  
  hitmaker = Env.create_path('hitmaker')
  hitmaker.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1 }) 
  hitmaker.add_processor(name="RawInputProcessor",params={'FileName': rawfile})
  hitmaker.add_processor(name="M26Unpacker")
  hitmaker.add_processor(name="M26Clusterizer")
  hitmaker.add_processor(name="M26CogHitMaker")
  hitmaker.add_processor(name="DiamondUnpacker")
  hitmaker.add_processor(name="DiamondRawHitSorter", params={'PixelType': diamond_pixeltype, })
  hitmaker.add_processor(name="DiamondClusterizer")
  hitmaker.add_processor(name="DiamondCogHitMaker")
  hitmaker.add_processor(name="FEI4Unpacker")
  hitmaker.add_processor(name="FEI4Clusterizer")
  hitmaker.add_processor(name="FEI4CogHitMaker")
  hitmaker.add_processor(name="RawDQM")
  hitmaker.add_processor(name="TelCorrelator")
  hitmaker.add_processor(name="LCIOOutput")
  
  kalman_aligner_1 = Env.create_path('kalman_aligner_1')
  kalman_aligner_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_1.add_processor(name="AlignTF_LC")
  kalman_aligner_1.add_processor(name="PreAligner", params={'ErrorsShiftX' : '0 10 10 10 10 10 10 0', 
                                                            'ErrorsShiftY' : '0 10 10 10 10 10 10 0', 
                                                            'ErrorsShiftZ' : '0 0 0 0 0 0 0 0', 
                                                            'ErrorsAlpha'  : '0 0 0 0 0 0 0 0',
                                                            'ErrorsBeta'   : '0 0 0 0 0 0 0 0', 
                                                            'ErrorsGamma'  : '0 0.01 0.01 0.01 0.01 0.01 0.01 0'})
  
  kalman_aligner_2 = Env.create_path('kalman_aligner_2')
  kalman_aligner_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_2.add_processor(name="AlignTF_TC")
  kalman_aligner_2.add_processor(name="TelAligner", params={'ErrorsShiftX' : '0 10 10 10 10 10 10 0', 
                                                            'ErrorsShiftY' : '0 10 10 10 10 10 10 0', 
                                                            'ErrorsShiftZ' : '0 10 10 10 10 10 10 0', 
                                                            'ErrorsAlpha'  : '0 0 0 0 0 0 0 0',
                                                            'ErrorsBeta'   : '0 0 0 0 0 0 0 0', 
                                                            'ErrorsGamma'  : '0 0.01 0.01 0.01 0.01 0.01 0.01 0'})
   
  
  telescope_dqm_1 = Env.create_path('telescope_dqm_1')
  telescope_dqm_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm_1.add_processor(name="AlignTF_TC")
  telescope_dqm_1.add_processor(name="TelescopeDQM", params={'RootFileName':'TelescopeDQM_basic.root'})
 

  cluster_calibration_1 = Env.create_path('cluster_calibration_1')
  cluster_calibration_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1, 'LCIOInputFiles': "tmp.slcio" })  
  cluster_calibration_1.add_processor(name="AlignTF_TC")
  cluster_calibration_1.add_processor(name="FEI4ClusterCalibrator")
  cluster_calibration_1.add_processor(name="M26ClusterCalibrator")

  kalman_aligner_3 = Env.create_path('kalman_aligner_3')
  kalman_aligner_3.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_3.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  kalman_aligner_3.add_processor(name="FEI4GoeHitMaker", params={'HitCollectionName' : 'goehit_fei' })
  kalman_aligner_3.add_processor(name="AlignTF_TC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_fei'})
  kalman_aligner_3.add_processor(name="TelAligner", params={'ErrorsShiftX' : '0 10 10 10 10 10 10 0', 
                                                            'ErrorsShiftY' : '0 10 10 10 10 10 10 0', 
                                                            'ErrorsShiftZ' : '0 10 10 10 10 10 10 0', 
                                                            'ErrorsAlpha'  : '0 0 0 0 0 0 0 0',
                                                            'ErrorsBeta'   : '0 0 0 0 0 0 0 0', 
                                                            'ErrorsGamma'  : '0 0.01 0.01 0.01 0.01 0.01 0.01 0'})
  
  cluster_calibration_2 = Env.create_path('cluster_calibration_2')
  cluster_calibration_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1, 'LCIOInputFiles': "tmp.slcio" })  
  cluster_calibration_2.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  cluster_calibration_2.add_processor(name="FEI4GoeHitMaker", params={'HitCollectionName' : 'goehit_fei' })
  cluster_calibration_2.add_processor(name="AlignTF_TC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_fei'})
  cluster_calibration_2.add_processor(name="FEI4ClusterCalibrator")
  cluster_calibration_2.add_processor(name="M26ClusterCalibrator") 
  
  telescope_dqm_2 = Env.create_path('telescope_dqm_2')
  telescope_dqm_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm_2.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  telescope_dqm_2.add_processor(name="FEI4GoeHitMaker", params={'HitCollectionName' : 'goehit_fei' })
  telescope_dqm_2.add_processor(name="AlignTF_TC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_fei'})
  telescope_dqm_2.add_processor(name="TelescopeDQM", params={'RootFileName':'TelescopeDQM_advanced.root'})
  
  
  kalman_aligner_dia_1 = Env.create_path('kalman_aligner_dia_1')
  kalman_aligner_dia_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :-1, 'LCIOInputFiles': "tmp.slcio" }) 
  kalman_aligner_dia_1.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  kalman_aligner_dia_1.add_processor(name="FEI4GoeHitMaker", params={'HitCollectionName' : 'goehit_fei' }) 
  kalman_aligner_dia_1.add_processor(name="AlignTF_Dia_LC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_fei hit_dia'})
  kalman_aligner_dia_1.add_processor(name="PreAligner", params={'ErrorsShiftX' : '0 0 0 10 0 0 0 0', 
                                                                'ErrorsShiftY' : '0 0 0 10 0 0 0 0', 
                                                                'ErrorsShiftZ' : '0 0 0 0 0 0 0 0', 
                                                                'ErrorsAlpha'  : '0 0 0 0 0 0 0 0',
                                                                'ErrorsBeta'   : '0 0 0 0 0 0 0 0', 
                                                                'ErrorsGamma'  : '0 0 0 0.01 0 0 0 0'})
  
  kalman_aligner_dia_2 = Env.create_path('kalman_aligner_dia_2')
  kalman_aligner_dia_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_dia_2.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  kalman_aligner_dia_2.add_processor(name="FEI4GoeHitMaker", params={'HitCollectionName' : 'goehit_fei' }) 
  kalman_aligner_dia_2.add_processor(name="AlignTF_Dia_LC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_fei hit_dia'})
  kalman_aligner_dia_2.add_processor(name="TelAligner", params={'ErrorsShiftX' : '0 0 0 10 0 0 0 0', 
                                                                'ErrorsShiftY' : '0 0 0 10 0 0 0 0', 
                                                                'ErrorsShiftZ' : '0 0 0 10 0 0 0 0', 
                                                                'ErrorsAlpha'  : '0 0 0 0 0 0 0 0',
                                                                'ErrorsBeta'   : '0 0 0 0 0 0 0 0', 
                                                                'ErrorsGamma'  : '0 0 0 0.01 0 0 0 0'})
  
  telescope_dqm_3 = Env.create_path('telescope_dqm_3')
  telescope_dqm_3.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm_3.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  telescope_dqm_3.add_processor(name="FEI4GoeHitMaker", params={'HitCollectionName' : 'goehit_fei' }) 
  telescope_dqm_3.add_processor(name="AlignTF_Dia_LC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_fei hit_dia'})
  telescope_dqm_3.add_processor(name="TelescopeDQM", params={'RootFileName':'DiamondDQM_basic.root'})
   
  diamond_calibration_1 = Env.create_path('diamond_calibration_1')
  diamond_calibration_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1, 'LCIOInputFiles': "tmp.slcio" })  
  diamond_calibration_1.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  diamond_calibration_1.add_processor(name="FEI4GoeHitMaker", params={'HitCollectionName' : 'goehit_fei' }) 
  diamond_calibration_1.add_processor(name="AlignTF_Dia_LC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_fei hit_dia'})
  diamond_calibration_1.add_processor(name="DiamondClusterCalibrator", params={'uCellPeriod': uCellPeriod, 'vCellPeriod': vCellPeriod})

  kalman_aligner_dia_3 = Env.create_path('kalman_aligner_dia_3')
  kalman_aligner_dia_3.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_dia_3.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  kalman_aligner_dia_3.add_processor(name="DiamondGoeHitMaker", params={'HitCollectionName' : 'goehit_dia' })
  kalman_aligner_dia_3.add_processor(name="FEI4GoeHitMaker", params={'HitCollectionName' : 'goehit_fei' })
  kalman_aligner_dia_3.add_processor(name="AlignTF_Dia_LC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_fei goehit_dia'})
  kalman_aligner_dia_3.add_processor(name="TelAligner", params={'ErrorsShiftX' : '0 0 0 10 0 0 0 0', 
                                                                'ErrorsShiftY' : '0 0 0 10 0 0 0 0', 
                                                                'ErrorsShiftZ' : '0 0 0 10 0 0 0 0', 
                                                                'ErrorsAlpha'  : '0 0 0 0 0 0 0 0',
                                                                'ErrorsBeta'   : '0 0 0 0 0 0 0 0', 
                                                                'ErrorsGamma'  : '0 0 0 0.01 0 0 0 0'})
  
  diamond_calibration_2 = Env.create_path('diamond_calibration_2')
  diamond_calibration_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1, 'LCIOInputFiles': "tmp.slcio" })  
  diamond_calibration_2.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  diamond_calibration_2.add_processor(name="DiamondGoeHitMaker", params={'HitCollectionName' : 'goehit_dia' })
  diamond_calibration_2.add_processor(name="FEI4GoeHitMaker", params={'HitCollectionName' : 'goehit_fei' })
  diamond_calibration_2.add_processor(name="DiamondGoeHitMaker", params={'HitCollectionName' : 'goehit_dia' })
  diamond_calibration_2.add_processor(name="AlignTF_Dia_LC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_fei goehit_dia'})
  diamond_calibration_2.add_processor(name="DiamondClusterCalibrator", params={'uCellPeriod': uCellPeriod, 'vCellPeriod': vCellPeriod})
  
  telescope_dqm_4 = Env.create_path('telescope_dqm_4')
  telescope_dqm_4.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm_4.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  telescope_dqm_4.add_processor(name="FEI4GoeHitMaker", params={'HitCollectionName' : 'goehit_fei' }) 
  telescope_dqm_4.add_processor(name="DiamondGoeHitMaker", params={'HitCollectionName' : 'goehit_dia' })
  telescope_dqm_4.add_processor(name="AlignTF_Dia_LC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_fei goehit_dia'})
  telescope_dqm_4.add_processor(name="TelescopeDQM", params={'RootFileName':'DiamondDQM_advanced.root'})
  
  # create sequence of calibration paths 
  telescope_path = [ hotpixelkiller , 
                    hitmaker, 
                    kalman_aligner_1, 
                    kalman_aligner_2, 
                    kalman_aligner_2, 
                    kalman_aligner_2, 
                    telescope_dqm_1,
                  ]
     
  telescope_db_path = [ cluster_calibration_1, 
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
                        kalman_aligner_3,
                        telescope_dqm_2,  
                      ] 
   
  diamond_path = [ kalman_aligner_dia_1,
                   kalman_aligner_dia_2,
                   kalman_aligner_dia_2,
                   kalman_aligner_dia_2,
                   telescope_dqm_3, 
                 ]   
  
  diamond_db_path = [ diamond_calibration_1, 
                      kalman_aligner_dia_3, 
                      kalman_aligner_dia_3, 
                      kalman_aligner_dia_3, 
                      diamond_calibration_2, 
                      diamond_calibration_2, 
                      diamond_calibration_2, 
                      #diamond_calibration_2, 
                      #diamond_calibration_2, 
                      #diamond_calibration_2, 
                      kalman_aligner_dia_3, 
                      kalman_aligner_dia_3, 
                      kalman_aligner_dia_3,
                      telescope_dqm_4,  
                    ] 

  calpath = []
  calpath.extend(telescope_path) 
  calpath.extend(telescope_db_path) 
  calpath.extend(diamond_path)
  calpath.extend(diamond_db_path)  
  
  return calpath


def create_reco_path(Env, rawfile, gearfile, diamond_pixeltype):
  """
  Returns a list of tbsw path objects for reconstruciton of a test beam run 
  """
  
  reco = Env.create_path('reco')
  reco.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1 }) 
  reco.add_processor(name="RawInputProcessor",params={'FileName': rawfile})
  reco.add_processor(name="M26Unpacker")
  reco.add_processor(name="M26Clusterizer")
  #reco.add_processor(name="M26CogHitMaker")
  reco.add_processor(name="M26GoeHitMaker")
  reco.add_processor(name="DiamondUnpacker")
  reco.add_processor(name="DiamondRawHitSorter", params={'PixelType': diamond_pixeltype, })
  reco.add_processor(name="DiamondClusterizer")
  #reco.add_processor(name="DiamondCogHitMaker")
  reco.add_processor(name="DiamondGoeHitMaker")
  reco.add_processor(name="FEI4Unpacker")
  reco.add_processor(name="FEI4Clusterizer")
  #reco.add_processor(name="FEI4CogHitMaker")
  reco.add_processor(name="FEI4GoeHitMaker")
  reco.add_processor(name="RecoTF")
  reco.add_processor(name="DiamondAnalyzer", params={'RootFileName':'Histos-DIA.root'})
  reco.add_processor(name="FEI4Analyzer")
  
  return [ reco ]

  
def calibrate_and_reconstruct(params):
  
  rawfile, gearfile, diamond_pixeltype = params
  
  # Tag for calibration data
  caltag = os.path.splitext(os.path.basename(rawfile))[0] + '-type-{}'.format(diamond_pixeltype) 
  
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
  CalObj = Calibration(steerfiles=steerfiles, name=caltag + '-cal') 
  
  # Create list of calibration paths
  calpath = create_calibration_path(CalObj, rawfile, gearfile, diamond_pixeltype)
  
  # Run the calibration steps 
  CalObj.calibrate(path=calpath,ifile=rawfile,caltag=caltag)  
   
  # Reconsruct the rawfile using caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=caltag + '-reco' )
  
  # Create reconstuction path
  recopath = create_reco_path(RecObj, rawfile, gearfile, diamond_pixeltype)  
  
  # Run the reconstuction  
  RecObj.reconstruct(path=recopath,ifile=rawfile,caltag=caltag) 
  
  
if __name__ == '__main__':
  
  # Select diamond pixel type: 0 not on diamond, 1 FE only, 2 rectangular, 3 hexagonal
  diamond_pixeltype = 1
  gearfile = 'gear_desy_dia_type_{}.xml'.format(diamond_pixeltype)
  calibrate_and_reconstruct( (runlist[0], gearfile, diamond_pixeltype) ) 
  
  diamond_pixeltype = 2
  gearfile = 'gear_desy_dia_type_{}.xml'.format(diamond_pixeltype)
  calibrate_and_reconstruct( (runlist[0], gearfile, diamond_pixeltype) ) 
  
  diamond_pixeltype = 3
  gearfile = 'gear_desy_dia_type_{}.xml'.format(diamond_pixeltype)
  calibrate_and_reconstruct( (runlist[0], gearfile, diamond_pixeltype) ) 
  
