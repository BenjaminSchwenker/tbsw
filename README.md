
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
  residuals.plot(inputfilename = trackfile, histofilename = ofile, basecut = "hasTrack==0")
      
  # Plot DUT hit efficiency histograms from the 'Track' tree 
  # in the workspace. 
  ofile = 'Example-Efficiency.root' 
  selection = "nTelTracks == 1 && cellU_fit >= 0 && cellU_fit < 64 && cellV_fit >= 0 && cellV_fit < 64"
  efficiency.plot(inputfilename = trackfile, histofilename = ofile, basecut = selection, ucells=64, vcells=64)
    
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


# What is the output data format of tbsw?

The tbsw framework uses Marlin (http://ilcsoft.desy.de/portal/software_packages/marlin/index_eng.html) to organize the calibration and reconstruction of 
test beam data into a ```Path``` of small steps handled by so called ```Processors```. The reconstruction path from the example script can illustrate the
situation:

```
#!python
def create_reco_path(Env):
  """
  Returns a list of tbsw path objects to reconstruct a test beam run 
  """
  
  reco = Env.create_path('reco')
  # Tells Marlin to read input data from lcio file with name rawfile
  reco.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents, 'LCIOInputFiles': rawfile }) 
  # Adds cluster of Mimosa26 digits to the path 
  reco.add_processor(name="M26Clusterizer")
  # Adds hit position reconstruction for Mimosa26 clusters to the path 
  reco.add_processor(name="M26GoeHitMaker")
  # Like above, but for the small DEPFET senor used as device under test
  reco.add_processor(name="DEPClusterizer")
  reco.add_processor(name="DEPGoeHitMaker")
  # Adds track finding and fitting to the path, track finder does not require hit on plane 3(DUT)
  reco.add_processor(name="RecoTF", params={'ExcludeDetector'=3})
  # Adds PixelDUTAnalyzer processor to path, this processor creats output root trees         
  reco.add_processor(name="DEPFETAnalyzer") 
  # Adds TrackFitDQM processor to path, creates histograms for track fit quality       
  reco.add_processor(name="TelescopeDQM")                              
    
  return [ reco ] 
```

The last two processors in this reconstruction path create the output data useful for studying the performance of the device under test. The 
default parameters for the DEPFETAnalyzer processor are taken from the file ```processors.xml``` living in the steerfiles folder: 

```xml
<processor name="DEPFETAnalyzer" type="PixelDUTAnalyzer">
  <!--PixelDUTAnalyzer: DUT Analysis Processor-->
  <!--Name of DUT hit collection-->
  <parameter name="HitCollection" type="string" lcioInType="TrackerHit" value="hit_dep"/>
  <!--Name of telescope track collection-->
  <parameter name="TrackCollection" type="string" lcioInType="Track" value="tracks"/>
  <!--This is the name of the LCIO file with the alignment constants (add .slcio)-->
  <parameter name="AlignmentDBFileName" type="string" value="localDB/eudet-alignmentDB.slcio"/>
  <!--Plane number of DUT along the beam line-->
  <parameter name="DUTPlane" type="int" value="3"/>
  <!--Maximum u residual for matching DUT hits to telescope track [mm]. Put -1 to deactivate cut.-->
  <parameter name="MaxResidualU" type="double" value="0.2"/>
  <!--Maximum v residual for matching DUT hits to telescope track [mm]. Put -1 to deactivate cut.-->
  <parameter name="MaxResidualV" type="double" value="0.2"/>
  <!--Output root file name-->
  <parameter name="RootFileName" type="string" value="Histos-H5.root"/>
</processor>
```

The root file ```Histos-H5.root``` will contain three TTree's called Hit, Track and Event. These trees will be filled with reference tracks 
interpolated to the third plane along the beam line, i.e. the position of the DEPFET sensor in our example. The main advantage of this choice 
is that we can always work in the local uvw coordinates of our device under test. 

- Origin of local coordinate system is center of the sensitive volume
- u,v,w positions are always in millimeters (mm)
- u axis points in direction of increasing ucells
- v axis points in direction of increasing vcells 
- w axis completes a right handed cartesian coordinate system
- The interpolated track state does not include information from the DUT hit

The structure of the root trees is always the same: 

```c
int _rootEventNumber;             // Event number from lcio file
int _rootRunNumber;               // Run number from lcio file 
int _rootSensorID;                // SensorID from lcio file (this is typically NOT the plane number!!)
int _rootNTelTracks;              // Number of tracks in reference telescope in same event as hit
int _rootNDUTHits;                // Number of DUT hits in the same event as hit   
int _rootHitQuality;              // GoodCluster == 0, BadCluster != 0
double _rootHitU;                 // Hit coordinate u reconstructed from DUT cluster in mm, in local DUT uvw coordinates       
double _rootHitV;                 // Hit coordinate v reconstructed from DUT cluster in mm, in local DUT uvw coordinates     
double _rootHitClusterCharge;     // Sum over all charges in the cluster 
double _rootHitSeedCharge;        // Highest charge in cluster
int _rootHitSize;                 // Number of hit cells (pixels/strips) in cluster
int _rootHitSizeU;                // Number of hit cells along u direction in cluster
int _rootHitSizeV;                // Number of hit cells along v direction in cluster
int _rootHitCellU;                // Hit u coordinate lies on this u cell
int _rootHitCellV;                // Hit v coordinate lies on this v cell
int _rootHitHasTrack;             // Hit can be matched to track (== 0)     
double _rootHitFitMomentum;       // Estimated track momentum from fit, only filled in case HasTrack==0            
double _rootHitFitU;              // Estimated track intersection u coordimate in mm, in local DUT uvw coordinates        
double _rootHitFitV;              // Estimated track intersection v coordimate in mm, in local DUT uvw coordinates                  
double _rootHitFitdUdW;           // Estimated track slope du/dw in radians, in local DUT uvw coordinates       
double _rootHitFitdVdW;           // Estimated track slope dv/dw in radians, in local DUT uvw coordinates    
double _rootHitFitErrorU;         // Estimated 1x sigma uncertainty for track intersection u coordinate
double _rootHitFitErrorV;         // Estimated 1x sigma uncertainty for track intersection v coordinate
double _rootHitPullResidualU;     // Standardized residual in u direction, should have mean = 0 and rms = 1
double _rootHitPullResidualV;     // Standardized residual in v direction, should have mean = 0 and rms = 1        
int _rootHitFitCellU;             // Estimated track intersection u coordinate lies on this u cell      
int _rootHitFitCellV;             // Estimated track intersection v coordinate lies on this v cell         
double _rootHitFitCellUCenter;    // Central coordinate of cell 'FitCellU' in mm        
double _rootHitFitCellVCenter;    // Central coordinate of cell 'FitCellV' in mm 
double _rootHitTrackChi2 ;        // Chi2 value from fit of reference track
double _rootHitLocalChi2;         // Chi2 value from hit-track residual on device under test 
int _rootHitTrackNDF;             // Number of degrees of freedom of track fit
int _rootHitTrackNHits;           // Number of telescope hits used for track fitting 
int _rootTrackHasHit;             // Track can be matched to a DUT hit (== 0) 
double _rootTrackFitMomentum;     // Estimated track momentum from fit    
int _rootTrackNDF;                // Number of degrees of freedom of track fit
double _rootTrackChi2;            // Chi2 value from fit of reference track
double _rootTrackLocalChi2;       // Chi2 value from hit-track residual on device under test 
double _rootTrackFitU ;           // Estimated track intersection u coordimate in mm, in local DUT uvw coordinates 
double _rootTrackFitV ;           // Estimated track intersection v coordimate in mm, in local DUT uvw coordinates 
double _rootTrackFitdUdW;         // Estimated track slope du/dw in radians, in local DUT uvw coordinates     
double _rootTrackFitdVdW;         // Estimated track slope dv/dw in radians, in local DUT uvw coordinates     
int _rootTrackFitCellU;           // Estimated track intersection u coordinate lies on this u cell  
int _rootTrackFitCellV;           // Estimated track intersection v coordinate lies on this v cell   
int _rootTrackNHits;              // Number of telescope hits used for track fitting 
double _rootTrackFitCellUCenter;  // Central coordinate of cell 'FitCellU' in mm 
double _rootTrackFitCellVCenter;  // Central coordinate of cell 'FitCellV' in mm 
double _rootTrackSeedCharge;      // Highest charge in cluster, only filled if cluster matched

// Hit Tree, filled once per DUT hit
_rootHitTree = new TTree("Hit","Hit info");
_rootHitTree->Branch("iRun"            ,&_rootRunNumber           ,"iRun/I");
_rootHitTree->Branch("iEvt"            ,&_rootEventNumber         ,"iEvt/I");
_rootHitTree->Branch("sensorID"        ,&_rootSensorID            ,"sensorID/I");
_rootHitTree->Branch("DEPFETGoodEvent" ,&_rootDEPFETGoodEvent     ,"DEPFETGoodEvent/I");
_rootHitTree->Branch("DEPFETStartgate" ,&_rootDEPFETStartGate     ,"DEPFETStartgate/I");       
_rootHitTree->Branch("nTelTracks"      ,&_rootNTelTracks          ,"nTelTracks/I"); 
_rootHitTree->Branch("nDutHits"        ,&_rootNDUTHits            ,"nDutHits/I");
_rootHitTree->Branch("clusterQuality"  ,&_rootHitQuality          ,"clusterQuality/I");
_rootHitTree->Branch("u_hit"           ,&_rootHitU                ,"u_hit/D");
_rootHitTree->Branch("v_hit"           ,&_rootHitV                ,"v_hit/D");     
_rootHitTree->Branch("clusterCharge"   ,&_rootHitClusterCharge    ,"clusterCharge/D");
_rootHitTree->Branch("seedCharge"      ,&_rootHitSeedCharge       ,"seedCharge/D");
_rootHitTree->Branch("sizeU"           ,&_rootHitSizeU            ,"sizeU/I");
_rootHitTree->Branch("sizeV"           ,&_rootHitSizeV            ,"sizeV/I");
_rootHitTree->Branch("size"            ,&_rootHitSize             ,"size/I");
_rootHitTree->Branch("hasTrack"        ,&_rootHitHasTrack         ,"hasTrack/I");   
_rootHitTree->Branch("u_fit"           ,&_rootHitFitU             ,"u_fit/D");
_rootHitTree->Branch("v_fit"           ,&_rootHitFitV             ,"v_fit/D"); 
_rootHitTree->Branch("dudw_fit"        ,&_rootHitFitdUdW          ,"dudw_fit/D");
_rootHitTree->Branch("dvdw_fit"        ,&_rootHitFitdVdW          ,"dvdw_fit/D");    
_rootHitTree->Branch("u_fiterr"        ,&_rootHitFitErrorU        ,"u_fiterr/D");
_rootHitTree->Branch("v_fiterr"        ,&_rootHitFitErrorV        ,"v_fiterr/D");   
_rootHitTree->Branch("pull_resu"       ,&_rootHitPullResidualU    ,"pull_resu/D");
_rootHitTree->Branch("pull_resv"       ,&_rootHitPullResidualV    ,"pull_resv/D");  
_rootHitTree->Branch("cellU_fit"       ,&_rootHitFitCellU         ,"cellU_fit/I");
_rootHitTree->Branch("cellV_fit"       ,&_rootHitFitCellV         ,"cellV_fit/I");
_rootHitTree->Branch("cellU_hit"       ,&_rootHitCellU            ,"cellU_hit/I");
_rootHitTree->Branch("cellV_hit"       ,&_rootHitCellV            ,"cellV_hit/I");
_rootHitTree->Branch("cellUCenter_fit" ,&_rootHitFitCellUCenter   ,"cellUCenter_fit/D");
_rootHitTree->Branch("cellVCenter_fit" ,&_rootHitFitCellVCenter   ,"cellVCenter_fit/D");                                      
_rootHitTree->Branch("trackChi2"       ,&_rootHitTrackChi2        ,"trackChi2/D");
_rootHitTree->Branch("trackNdof"       ,&_rootHitTrackNDF         ,"trackNdof/I");
_rootHitTree->Branch("trackNHits"      ,&_rootHitTrackNHits       ,"trackNHits/I");  
_rootHitTree->Branch("momentum"        ,&_rootHitFitMomentum      ,"momentum/D");    
_rootHitTree->Branch("localChi2"       ,&_rootHitLocalChi2        ,"localChi2/D"); 
   
// Track tree, filled once per reference track 
_rootTrackTree = new TTree("Track","Track info");
_rootTrackTree->Branch("iRun"            ,&_rootRunNumber           ,"iRun/I");
_rootTrackTree->Branch("iEvt"            ,&_rootEventNumber         ,"iEvt/I");
_rootTrackTree->Branch("sensorID"        ,&_rootSensorID            ,"sensorID/I");
_rootTrackTree->Branch("DEPFETGoodEvent" ,&_rootDEPFETGoodEvent     ,"DEPFETGoodEvent/I");
_rootTrackTree->Branch("DEPFETStartgate" ,&_rootDEPFETStartGate     ,"DEPFETStartgate/I");    
_rootTrackTree->Branch("nTelTracks"      ,&_rootNTelTracks          ,"nTelTracks/I"); 
_rootTrackTree->Branch("nDutHits"        ,&_rootNDUTHits            ,"nDutHits/I");
_rootTrackTree->Branch("hasHit"          ,&_rootTrackHasHit         ,"hasHit/I");
_rootTrackTree->Branch("momentum"        ,&_rootTrackFitMomentum    ,"momentum/D");                                                           
_rootTrackTree->Branch("u_fit"           ,&_rootTrackFitU           ,"u_fit/D");
_rootTrackTree->Branch("v_fit"           ,&_rootTrackFitV           ,"v_fit/D");
_rootTrackTree->Branch("dudw_fit"        ,&_rootTrackFitdUdW        ,"dudw_fit/D");
_rootTrackTree->Branch("dvdw_fit"        ,&_rootTrackFitdVdW        ,"dvdw_fit/D");
_rootTrackTree->Branch("cellU_fit"       ,&_rootTrackFitCellU       ,"cellU_fit/I");
_rootTrackTree->Branch("cellV_fit"       ,&_rootTrackFitCellV       ,"cellV_fit/I");
_rootTrackTree->Branch("cellUCenter_fit" ,&_rootTrackFitCellUCenter ,"cellUCenter_fit/D");
_rootTrackTree->Branch("cellVCenter_fit" ,&_rootTrackFitCellVCenter ,"cellVCenter_fit/D");
_rootTrackTree->Branch("trackChi2"       ,&_rootTrackChi2           ,"trackChi2/D");
_rootTrackTree->Branch("trackNdof"       ,&_rootTrackNDF            ,"trackNdof/I");
_rootTrackTree->Branch("trackNHits"      ,&_rootTrackNHits          ,"trackNHits/I");  
_rootTrackTree->Branch("seedCharge"      ,&_rootTrackSeedCharge     ,"seedCharge/D");  
_rootTrackTree->Branch("localChi2"       ,&_rootTrackLocalChi2      ,"localChi2/D"); 
      
// Event tree, filled once per event 
_rootEventTree = new TTree("Event","Event info");
_rootEventTree->Branch("iRun"            ,&_rootRunNumber       ,"iRun/I");
_rootEventTree->Branch("iEvt"            ,&_rootEventNumber     ,"iEvt/I");
_rootEventTree->Branch("sensorID"        ,&_rootSensorID        ,"sensorID/I");   
_rootEventTree->Branch("DEPFETGoodEvent" ,&_rootDEPFETGoodEvent ,"DEPFETGoodEvent/I");
_rootEventTree->Branch("DEPFETStartgate" ,&_rootDEPFETStartGate ,"DEPFETStartgate/I");     
_rootEventTree->Branch("nTelTracks"      ,&_rootNTelTracks      ,"nTelTracks/I"); 
_rootEventTree->Branch("nDutHits"        ,&_rootNDUTHits        ,"nDutHits/I");

```

The plotting of root trees is well known in the HEP community. The script tbsw_example.py contains some instructive examples
based on PyRoot showing how to get residuals and efficiencies plots from the root trees. 
 

Have fun with test beams ;)  

benjamin.schwenker@phys.uni-goettingen.de
