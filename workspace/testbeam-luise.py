"""
Script for processing mini strip sensor data from test beam at DESY in 2016. Rawdata files 
are kept by Luise Poley. 

The input data consists in 'merged' lcio files with digits from two ministrip sensors 
and six m26 sensors forming a reference telescope. 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *
import multiprocessing

# Path to steering files 
steerfiles = 'steering-files/luise-tb/'


# List of runs to be processed
runlist = [
            '/home/benjamin/hep/testbeam/workspace_desy_itk_luise_test/data_from_kristin/run000282-merger.slcio',
          ]

# Number of events to be processed (put -1 to process all)
nevents = -1


def create_calibration_path(Env, rawfile, gearfile):
  """
  Returns a list of tbsw path objects needed to calibrate the tracking telescope
  """
  
  hotpixelkiller = Env.create_path('hotpixelkiller')
  hotpixelkiller.set_globals(params={'GearXMLFile': gearfile, 'MaxRecordNumber' : 100000, 'LCIOInputFiles': rawfile })  
  hotpixelkiller.add_processor(name="M26Unpacker")
  hotpixelkiller.add_processor(name="StripUnpacker")
  hotpixelkiller.add_processor(name="M26HotPixelKiller")
  hotpixelkiller.add_processor(name="StripHotPixelKiller")

  correlator = Env.create_path('correlator')
  correlator.set_globals(params={'GearXMLFile': gearfile, 'MaxRecordNumber' : nevents, 'LCIOInputFiles': rawfile }) 
  correlator.add_processor(name="M26Unpacker")
  correlator.add_processor(name="M26Clusterizer")
  correlator.add_processor(name="M26CogHitMaker")
  correlator.add_processor(name="StripUnpacker")
  correlator.add_processor(name="StripClusterizer", params={'SparseZSCut': 30})
  correlator.add_processor(name="StripCogHitMaker")
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
  cluster_calibration_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents, 'LCIOInputFiles': "tmp.slcio" })  
  cluster_calibration_1.add_processor(name="AlignTF_TC")
  cluster_calibration_1.add_processor(name="StripID6ClusterCalibrator")
  cluster_calibration_1.add_processor(name="StripID7ClusterCalibrator")
  cluster_calibration_1.add_processor(name="M26ClusterCalibrator")

  kalman_aligner_3 = Env.create_path('kalman_aligner_3')
  kalman_aligner_3.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_3.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  kalman_aligner_3.add_processor(name="StripID6GoeHitMaker", params={'HitCollectionName' : 'goehit_stripid6' })
  kalman_aligner_3.add_processor(name="StripID7GoeHitMaker", params={'HitCollectionName' : 'goehit_stripid7' })
  kalman_aligner_3.add_processor(name="AlignTF_TC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_stripid6 goehit_stripid7'})
  kalman_aligner_3.add_processor(name="TelAligner")
  
  cluster_calibration_2 = Env.create_path('cluster_calibration_2')
  cluster_calibration_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents, 'LCIOInputFiles': "tmp.slcio" })  
  cluster_calibration_2.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  cluster_calibration_2.add_processor(name="StripID6GoeHitMaker", params={'HitCollectionName' : 'goehit_stripid6' })
  cluster_calibration_2.add_processor(name="StripID7GoeHitMaker", params={'HitCollectionName' : 'goehit_stripid7' })
  cluster_calibration_2.add_processor(name="AlignTF_TC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_stripid6 goehit_stripid7'})
  cluster_calibration_2.add_processor(name="StripID6ClusterCalibrator")
  cluster_calibration_2.add_processor(name="StripID7ClusterCalibrator")
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
  reco.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents, 'LCIOInputFiles': rawfile})   
  reco.add_processor(name="M26Unpacker")
  reco.add_processor(name="M26Clusterizer")
  reco.add_processor(name="M26GoeHitMaker")
  reco.add_processor(name="StripUnpacker")
  reco.add_processor(name="StripClusterizer", params={'SparseZSCut': 30})
  reco.add_processor(name="StripID6GoeHitMaker", params={'HitCollectionName' : 'goehit_stripid6' })
  reco.add_processor(name="StripID7GoeHitMaker", params={'HitCollectionName' : 'goehit_stripid7' })
  reco.add_processor(name="RecoTF")
  reco.add_processor(name="StripID6Analyzer", params={'HitCollection' : 'goehit_stripid6' })
  reco.add_processor(name="StripID7Analyzer", params={'HitCollection' : 'goehit_stripid7' })

  return [ reco ]

  
def calibrate_and_reconstruct(params):
  
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
  
  # Reconsruct the rawfile using caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=caltag + '-reco2' )

  # Create reconstuction path
  recopath = create_reco_path(RecObj, rawfile, gearfile)  

  # Run the reconstuction  
  RecObj.reconstruct(path=recopath,ifile=rawfile,caltag=caltag) 
  
  
if __name__ == '__main__':
  count = 1 # multiprocessing.cpu_count()
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
  
  
  params = [ (rawfile, 'gear_desy_goeid1.xml' ) for rawfile in runlist ]
  pool.map(calibrate_and_reconstruct, params)

  
  # Create some efficiency and residual plots from root files after 
  # reconstruction. 
  trackfiles =  glob.glob('root-files/Histos-ITK*.root')  
  
  for trackfile in trackfiles: 
        
    ofile = 'Residuals-' + os.path.basename(trackfile)
    residuals.plot(inputfilename = trackfile, histofilename = ofile, basecut = "hasTrack==0")
    
    ofile = 'Occupancy-' + os.path.basename(trackfile)
    occupancy.plot(inputfilename = trackfile, histofilename = ofile, ucells=128, vcells=1)
    
    ofile = 'Efficiency-' + os.path.basename(trackfile)  
    efficiency.plot(inputfilename = trackfile, histofilename = ofile, ucells=128, vcells=1)
    
    ofile = 'InpixEfficiency-' + os.path.basename(trackfile)  
    efficiency.plot_unit_inpix(inputfilename = trackfile, histofilename = ofile, basecut = "nTelTracks == 1", matchcut="hasHit == 0", upitch=0.0745, vpitch=0.5, ubins=20, vbins=10)
    
    ofile = 'UnitInpixEfficiency-' + os.path.basename(trackfile)
    efficiency.plot_inpix(inputfilename = trackfile, histofilename = ofile, basecut = "nTelTracks == 1", matchcut="hasHit == 0", usize=1.0, vsize=1.0, ubins=200, vbins=1)
    
    ofile = 'InpixSignal-' + os.path.basename(trackfile)
    inpixel.plot(inputfilename = trackfile, histofilename = ofile, usize=1.0, vsize=1.0, ubins=100, vbins=10)
    
    ofile = 'UnitInpixSignal-' + os.path.basename(trackfile)
    inpixel.plot_unit(inputfilename = trackfile, histofilename = ofile, upitch=0.0745, vpitch=0.5, ubins=20, vbins=10)



