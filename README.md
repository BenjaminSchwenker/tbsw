
# Introduction to the test beam software framework (tbsw) 

The test beam software framework (tbsw) deals with the reconstruction and analysis of test beam data obtained with the 
EUDET/AIDA reference telescope. Test beams are all about learning about a device under test, or DUT, installed in the 
centre of the reference telescope. The device under test is position sensitive detector, typically a pixel or strip 
silicon sensor. The idea is to probe the device under test with charged tracks crossing the DUT at known times and 
precisely estimated kinematic parameters. 

The tbsw framework is designed to run on local machines (Desktop / Laptop) to enable detector groups doing their own data 
analysis right during the beam test (or shortly afterwards). Some effort was invested to make this as smooth and simple
as possible. The installation of the framework is covered in INSTALL. 


# For the impatient reader: The text below is a step by step instruction how to run the code on raw data


After installation, the root directory of tbsw contains a 'workspace' folder. The first step to analyse data from a new beam test
is to make a working copy of this folder. For example: 

```
$ cp -r workspace ~/workspace-my-beam-test
```

The idea of the folder 'workspace-my-beam-test' is to collect all steering files and analysis results from a specific beam test. 
Please replace the 'my-beam-test' part with something more specific like 'desy_october_2016' or similar. For simplicity, we will
call this folder just workspace in the following. 

It is a good idea to copy the workspace to some place outside of the tbsw installation. In case you want to reprocess your data with 
a new (different) version of tbsw, you typically just have to make sure the environment variables in 'init_tbsw.sh' point to the
desired installation.

Raw data files (.raw) should not be stored directly in your workspace folder. Again, this simplifies switching versions of tbsw. It is
advised collect all raw data files in a seperate folder. 

For example: 

```
$ mkdir path-to-big-disk/rawdata-my-beam-test 
```

Raw data from a one week beam test is typically ~1TB and fits on most Desktop PCs or laptops.  

Next, go to your new workspace folder and setup all environment variables. 

```
$ . init_tbsw.sh
```

A test beam with the EUDET/AIDA telescope results in a number of .raw files. Each raw file corresponds to a particular run having 
a specific test beam geometry, beam energy and DUT configuration.  Before we can start with the reconstruction and analysis of data, 
we need to prepare a set of steering files to describe this situation. 

In order to have an introductory example, we will perform a simulation experiment. In this way, we do not need access to real test beam
data but generate the raw data for one run by our own. The added benefit of this approach is that we can compare the reconstructed data 
with the simulated truth data. The steerfiles folder 'workspace/steering-files/depfet-H5-tb/' contains: 


- gear_desynov15_geoID15.xml         : Telescope geometry file, describes a small depfet sensor installed in the EUDET telescope 
- processors.xml                     : XML file containing default parameters for all Marlin processors needed for this example 
- align-config-iteration-1.cfg       : Alignment configuration file, describes which sensor parameters should be corrected 
- align-config-iteration-2.cfg       : Alignment configuration file, describes which sensor parameters should be corrected 


All these files area readable by a text editor and heavily commented. The Python script tbsw_example.py runs us through all steps of the 
example. The source code of the tbsw_example.py is: 

```
#!python

"""
This is an example script to demonstrate how tbsw can be used to analyze test beam 
data using Python scripts.

The script below simulates a test beam experiment where charged tracks cross a misaligned
pixel telescope containing six Mimosa 26 detector planes forming a reference telescope.
The reference telescope has two arms with three sensors. A small DEPFET sensor with 32x64
pixels is installed in the center of the reference telescope as the device under test.

The script creates a lcio file containing simulated digits, performs a full calibration 
of the misaligned telescope and reconstructs the data. Root files containing tracks and 
hits at the device under test are created in the folder root-files. 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *
from gear import *
import tbsw_tools

# Path to steering files 
# The folder steerfiles contains one or more gear file describing the 
# nominal telescope geometry. 
steerfiles = 'steering-files/depfet-H5-tb/'
# Select the name of a gearfile to use from the steerfiles folder  
gearfile = 'gear_desynov15_geoID15.xml'
# Select filename for the simulated test beam run  
rawfile = os.getcwd() + '/simrun.slcio'
# Number of events to simulate 
nevents = 5000000

def create_sim_path(Env):
  """
  Returns a list of tbsw path objects to simulate a test beam run 
  """
  
  sim = Env.create_path('sim')
  sim.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents})   
  sim.add_processor(name="InfoSetter")
  sim.add_processor(name="ParticleGun")
  sim.add_processor(name="FastSim")
  sim.add_processor(name="TriggerGenerator")
  sim.add_processor(name="M26Digitizer")
  sim.add_processor(name="DEPFETDigitizer")
  sim.add_processor(name="LCIOOutput",params={"LCIOOutputFile" : rawfile })

  cluster_calibrator_mc = Env.create_path('cluster_calibrator_mc')
  cluster_calibrator_mc.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1, 'LCIOInputFiles': rawfile })  
  cluster_calibrator_mc.add_processor(name="M26Clusterizer")
  cluster_calibrator_mc.add_processor(name="M26CogHitMaker")
  cluster_calibrator_mc.add_processor(name="DEPClusterizer")
  cluster_calibrator_mc.add_processor(name="DEPCogHitMaker")
  cluster_calibrator_mc.add_processor(name="M26ClusterCalibrationFromMC")
  cluster_calibrator_mc.add_processor(name="DEPClusterCalibrationFromMC")
  
  return [ sim , cluster_calibrator_mc]

def create_calibration_path(Env):
  """
  Returns a list of tbsw path objects to calibrate the tracking telescope
  """
  
  hotpixelkiller = Env.create_path('hotpixelkiller')
  hotpixelkiller.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': rawfile })  
  hotpixelkiller.add_processor(name="M26HotPixelKiller")
  hotpixelkiller.add_processor(name="DEPHotPixelKiller")
  
  correlator = Env.create_path('correlator')
  correlator.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents, 'LCIOInputFiles': rawfile }) 
  correlator.add_processor(name="M26Clusterizer")
  correlator.add_processor(name="M26CogHitMaker")
  correlator.add_processor(name="DEPClusterizer")
  correlator.add_processor(name="DEPCogHitMaker")
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
  cluster_calibration_1.add_processor(name="M26ClusterCalibrator")
  cluster_calibration_1.add_processor(name="DEPClusterCalibrator")

  kalman_aligner_3 = Env.create_path('kalman_aligner_3')
  kalman_aligner_3.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_3.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  kalman_aligner_3.add_processor(name="DEPGoeHitMaker", params={'HitCollectionName' : 'goehit_dep' })
  kalman_aligner_3.add_processor(name="AlignTF_TC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_dep'})
  kalman_aligner_3.add_processor(name="TelAligner")
  
  cluster_calibration_2 = Env.create_path('cluster_calibration_2')
  cluster_calibration_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 4000000, 'LCIOInputFiles': "tmp.slcio" })  
  cluster_calibration_2.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' }) 
  cluster_calibration_2.add_processor(name="DEPGoeHitMaker", params={'HitCollectionName' : 'goehit_dep' })
  cluster_calibration_2.add_processor(name="AlignTF_TC", params={'InputHitCollectionNameVec': 'goehit_m26 goehit_dep'})
  cluster_calibration_2.add_processor(name="M26ClusterCalibrator")
  cluster_calibration_2.add_processor(name="DEPClusterCalibrator")
  cluster_calibration_2.add_processor(name="TelescopeDQM", params={'RootFileName': 'ClusterDQM.root'})
  
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


def create_reco_path(Env):
  """
  Returns a list of tbsw path objects to reconstruct a test beam run 
  """
  
  reco = Env.create_path('reco')
  reco.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents, 'LCIOInputFiles': rawfile }) 
  reco.add_processor(name="M26Clusterizer")
  reco.add_processor(name="M26GoeHitMaker")
  reco.add_processor(name="DEPClusterizer")
  reco.add_processor(name="DEPGoeHitMaker")
  reco.add_processor(name="RecoTF")         
  reco.add_processor(name="DEPFETAnalyzer")      
  reco.add_processor(name="TelescopeDQM")                              
    
  return [ reco ]

  
def simulate(params): 
  """
  Simulates a rawfile from a simulated test beam experiment
  Creates a folder tmp-runs/name-sim/ and populates it with 
  Marlin steering and logfiles.  
  """ 
  
  rawfile, steerfiles, gearfile = params
  
  # Create tmpdir to hold all steerfiles and log files 
  SimObj = Simulation(steerfiles=steerfiles, name=os.path.splitext(os.path.basename(rawfile))[0] + '-sim' )

  # Create steerfiles for processing
  simpath = create_sim_path(SimObj)

  # Randomize sensor positions in gear file to create misalignment
  randomize_telescope(gearfile=SimObj.get_filename(gearfile), mean_pos=0, sigma_pos=0.1, mean_rot=0, sigma_rot=0.1)
   
  # Run simulation to create rawfile with simulated digits 
  SimObj.simulate(path=simpath)  

  # Export clusterDB created from truth hits
  SimObj.export_caltag(caltag='simulation')
   

def calibrate_and_reconstruct(params):
  """
  Calibrates an misaligned tracking telescope from run data. 
  Creates a folder localDB/caltag in workspace containing 
  calibration results. 
  Creates a folder tmp-runs/name-sim/ and populates it with 
  Marlin steering and logfiles.  
  """ 
  
  rawfile, steerfiles, gearfile = params
  
  # Tag for calibration data
  caltag = os.path.splitext(os.path.basename(rawfile))[0] + '-test'
  
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
  CalObj = Calibration(steerfiles=steerfiles, name=caltag + '-cal') 
  
  # Create list of calibration steps 
  calpath = create_calibration_path(CalObj)
  
  # Run the calibration steps 
  CalObj.calibrate(path=calpath,ifile=rawfile,caltag=caltag)  
  
  # Reconsruct the rawfile using the caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=caltag + '-reco2' )

  # Create reconstuction path
  recopath = create_reco_path(RecObj)  

  # Run the reconstuction  
  RecObj.reconstruct(path=recopath,ifile=rawfile,caltag=caltag) 
  
  
if __name__ == '__main__':
  
  params = ( rawfile, steerfiles, gearfile )

  # Create a simulated rawfile 
  simulate( params )

  # Calibrate the telescope and reconstruct the rawfile 
  calibrate_and_reconstruct( params )
  
  # Make a list of root files containing reconstructed trees 
  # for tracks / hits / events
  trackfile = root-files/Histos-H5-simrun-test-reco2.root'  
  
  # Plot DUT residuals and cluster signal histograms from the 'Hit'
  # tree in the workspace. 
  ofile = 'Example-Residuals.root'
  tbsw_tools.residuals.plot(inputfilename = trackfile, histofilename = ofile, basecut = "hasTrack==0")
      
  # Plot DUT hit efficiency histograms from the 'Track' tree 
  # in the workspace. 
  ofile = 'Example-Efficiency.root' 
  selection = "nTelTracks == 1 && cellU_fit >= 0 && cellU_fit < 64 && cellV_fit >= 0 && cellV_fit < 64"
  tbsw_tools.efficiency.plot(inputfilename = trackfile, histofilename = ofile, basecut = selection, ucells=64, vcells=64)
    
```

You can run the script by typing:  

```
$ python tbsw_example.py
```

In case the script is working, your console should display messages like: 

```
[INFO] Starting to simulate a test beam run ...
[INFO] Marlin sim.xml is done
[INFO] Marlin cluster_calibrator_mc.xml is done
('[INFO] Created new caltag ', 'simulation')
[INFO] Starting to calibrate file /home/benjamin/hep/testbeam/workspace_july16_itk/simrun.slcio ...
[INFO] Marlin hotpixelkiller.xml is done
[INFO] Marlin correlator.xml is done
[INFO] Marlin kalman_aligner_1.xml is done
[INFO] Marlin kalman_aligner_2.xml is done
[INFO] Marlin kalman_aligner_2.xml is done
[INFO] Marlin kalman_aligner_2.xml is done
[INFO] Marlin telescope_dqm.xml is done
[INFO] Marlin cluster_calibration_1.xml is done
[INFO] Marlin kalman_aligner_3.xml is done
[INFO] Marlin kalman_aligner_3.xml is done
[INFO] Marlin kalman_aligner_3.xml is done
[INFO] Marlin cluster_calibration_2.xml is done
[INFO] Marlin cluster_calibration_2.xml is done
[INFO] Marlin cluster_calibration_2.xml is done
[INFO] Marlin cluster_calibration_2.xml is done
[INFO] Marlin cluster_calibration_2.xml is done
[INFO] Marlin cluster_calibration_2.xml is done
[INFO] Marlin kalman_aligner_3.xml is done
[INFO] Marlin kalman_aligner_3.xml is done
[INFO] Marlin kalman_aligner_3.xml is done
('[INFO] Created new caltag ', 'simrun-test')
[INFO] Done processing file /home/benjamin/hep/testbeam/workspace_july16_itk/simrun.slcio
[INFO] Starting to reconstruct file /home/benjamin/hep/testbeam/workspace_july16_itk/simrun.slcio ...
[INFO] Marlin reco.xml is done
[INFO] Done processing file /home/benjamin/hep/testbeam/workspace_july16_itk/simrun.slcio
```

All calibration data is collected in the folder 'workspace/localDB' and root files with reconstructed tracks and hits are collected in 
the folder 'workspace/root-files'. The structure of the root trees for tracks hits and evenst is explained in the next section. The example
script creates two root files with simple residual histograms ('Example-Residuals.root') and DUT hit efficiency plots ('Example-Efficiency.root')
to show what can be done using the hit and track trees.   


Have fun with test beams ;)  

benjamin.schwenker@phys.uni-goettingen.de
