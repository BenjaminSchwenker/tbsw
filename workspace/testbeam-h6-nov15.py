"""
Script for processing DEPFET Hyrbid6 data from test beam at DESY in november 2015. Rawdata files 
are kept by Florian Luetticke in Bonn. 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *
import multiprocessing

# Path to steering files 
steerfiles = 'steering-files/depfet-H6-tb/'



# List of runs to be processed using gear file 'gear_geoid2.xml'
runlist_partone = [
                    '/home/benjamin/Desktop/rawdata-nov15/run000065.raw',
                  ]

# List of runs to be processed using gear file 'gear_geoid16.xml'
runlist_parttwo = [
                 
                  ]


def create_calibration_path(Env, rawfile, gearfile):
  """
  Returns a list of tbsw path objects needed to calibrate the tracking telescope
  """
  
  hotpixelkiller = Env.create_path('hotpixelkiller')
  hotpixelkiller.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000 })  
  hotpixelkiller.add_processor(name="RawInputProcessor",params={'FileName': rawfile})
  hotpixelkiller.add_processor(name="M26Unpacker")
  hotpixelkiller.add_processor(name="DEPHybrid6Unpacker")
  hotpixelkiller.add_processor(name="M26HotPixelKiller")
  hotpixelkiller.add_processor(name="DEPHotPixelKiller")
  hotpixelkiller.add_processor(name="DEPHotPixelKiller_2")

  correlator = Env.create_path('correlator')
  correlator.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 1000000 }) 
  correlator.add_processor(name="RawInputProcessor",params={'FileName': rawfile})
  correlator.add_processor(name="M26Unpacker")
  correlator.add_processor(name="M26Clusterizer")
  correlator.add_processor(name="M26CogHitMaker")
  correlator.add_processor(name="DEPHybrid6Unpacker")
  correlator.add_processor(name="DEPClusterizer_2")
  correlator.add_processor(name="DEPCogHitMaker_2")
  correlator.add_processor(name="RawDQM")
  correlator.add_processor(name="TelCorrelator")
  correlator.add_processor(name="LCIOOutput")
  
  kalman_aligner_1 = Env.create_path('kalman_aligner_1')
  kalman_aligner_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_1.add_processor(name="AlignTF_LC")
  kalman_aligner_1.add_processor(name="PreAligner")
  
  kalman_aligner_2 = Env.create_path('kalman_aligner_2')
  kalman_aligner_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_2.add_processor(name="AlignTF_TC")
  kalman_aligner_2.add_processor(name="TelAligner")
  
  telescope_dqm = Env.create_path('telescope_dqm')
  telescope_dqm.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm.add_processor(name="AlignTF_TC")
  telescope_dqm.add_processor(name="TelescopeDQM")
  
  cluster_calibration_1 = Env.create_path('cluster_calibration_1')
  cluster_calibration_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 4000000, 'LCIOInputFiles': "tmp.slcio" })  
  cluster_calibration_1.add_processor(name="AlignTF_TC")
  cluster_calibration_1.add_processor(name="DEPClusterCalibrator")
  cluster_calibration_1.add_processor(name="M26ClusterCalibrator")

  kalman_aligner_3 = Env.create_path('kalman_aligner_3')
  kalman_aligner_3.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_3.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  kalman_aligner_3.add_processor(name="DEPGoeHitMaker_2", params={'HitCollectionName' : 'goehit_dep2' })
  kalman_aligner_3.add_processor(name="AlignTF_TC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_dep2'})
  kalman_aligner_3.add_processor(name="TelAligner")

  cluster_calibration_2 = Env.create_path('cluster_calibration_2')
  cluster_calibration_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 4000000, 'LCIOInputFiles': "tmp.slcio" })  
  cluster_calibration_2.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  cluster_calibration_2.add_processor(name="DEPGoeHitMaker_2", params={'HitCollectionName' : 'goehit_dep2' })
  cluster_calibration_2.add_processor(name="AlignTF_TC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_dep2'})
  cluster_calibration_2.add_processor(name="DEPClusterCalibrator")
  cluster_calibration_2.add_processor(name="M26ClusterCalibrator")
  
  
  # create sequence of calibration paths 
  calpath= [ hotpixelkiller , 
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
             kalman_aligner_3, 
             kalman_aligner_3, 
             kalman_aligner_3, 
           ]
  
  return calpath


def create_reco_path(Env, rawfile, gearfile):
  """
  Returns a list of tbsw path objects for reconstruciton of a test beam run 
  """
  
  reco = Env.create_path('reco')
  reco.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 1000000})   
  reco.add_processor(name="RawInputProcessor", params={'FileName': rawfile})
  reco.add_processor(name="M26Unpacker")
  reco.add_processor(name="M26Clusterizer")
  reco.add_processor(name="M26GoeHitMaker")
  reco.add_processor(name="DEPHybrid6Unpacker")
  reco.add_processor(name="DEPClusterizer")
  reco.add_processor(name="DEPGoeHitMaker")
  reco.add_processor(name="DEPClusterizer_2")
  reco.add_processor(name="DEPGoeHitMaker_2")
  reco.add_processor(name="RecoTF", params={'OutlierChi2Cut': 20})
  reco.add_processor(name="DEPFETAnalyzer_2")
  return [ reco ]

  
def calibrate(params):
  
  rawfile, gearfile = params
  
  # Tag for calibration data
  caltag = os.path.splitext(os.path.basename(rawfile))[0] + '-test'
  
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
  CalObj = Calibration(steerfiles=steerfiles, name=caltag + '-cal') 
  
  # Create list of calibration paths
  calpath = create_calibration_path(CalObj, rawfile, gearfile)
  
  # Run the calibration steps 
  CalObj.calibrate(path=calpath,ifile=rawfile,caltag=caltag)  
  

def reconstruct(params):
  
  rawfile, gearfile = params
  
  # Tag for calibration data
  caltag = os.path.splitext(os.path.basename(rawfile))[0] + '-test'
  
  # Reconsruct the rawfile using caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=caltag + '-reco2' )

  # Create reconstuction path
  recopath = create_reco_path(RecObj, rawfile, gearfile)  

  # Run the reconstuction  
  RecObj.reconstruct(path=recopath,ifile=rawfile,caltag=caltag) 
  
  
if __name__ == '__main__':
  count = multiprocessing.cpu_count()
  pool = multiprocessing.Pool(processes=count)
   
  try:
    import psutil
    parent = psutil.Process()
    parent.nice(18)
    for child in parent.children():
      child.nice(18)
    print("All %d processes set to minimum priority"%(count))
  except:
    print("Could not change process priority. Please install psutil")
    pass
  
  
  # Runs from first part must be processed  with gearfile geoid 2
  params = [ (rawfile, 'gear_geoid2.xml' ) for rawfile in runlist_partone ]
  pool.map(calibrate, params)

  # Runs from first part must be processed  with gearfile geoid 2
  params = [ (rawfile, 'gear_geoid2.xml' ) for rawfile in runlist_partone ]
  pool.map(reconstruct, params)


  # Runs from second part must be processed with gearfile geoid 16
  params = [ (rawfile, 'gear_geoid16.xml') for rawfile in runlist_parttwo ]
  pool.map(calibrate, params)
  
  # Runs from second part must be processed with gearfile geoid 16
  params = [ (rawfile, 'gear_geoid16.xml') for rawfile in runlist_parttwo ]
  pool.map(reconstruct, params)

  
  # Create some efficiency and residual plots from root files after 
  # reconstruction. 
  trackfiles =  glob.glob('root-files/Histos-H6-*.root')  
  
  for trackfile in trackfiles: 
        
    ofile = 'Residuals-' + os.path.basename(trackfile)
    residuals.plot(inputfilename = trackfile, histofilename = ofile, basecut = "hasTrack==0")
      
    ofile = 'Efficiency-' + os.path.basename(trackfile)  
    efficiency.plot(inputfilename = trackfile, histofilename = ofile, basecut = "nTelTracks == 1 && nDutHits > 0 && cellU_fit >= 0 && cellU_fit < 64", matchcut="hasHit == 0", ucells=64, vcells=480)



