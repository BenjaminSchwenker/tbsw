
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
with the simulated truth data. The steerfiles folder 'workspace/steering-files/depfet-tb/' contains: 


- gear_desy_W11OF2.xml: Telescope geometry file, describes a Belle II PXD sensor installed in the EUDET telescope (6x M26 planes and a FEI4 timing plane)  


The Python script example.py runs us through all steps of the example which includes the generation of simulated data, the calibration of 
telescope data and the full reconstruction with final calibration constants. 

You can run the script by typing:  

```
$ python example.py
```

In case the script is working, your console should display messages like: 

```
benjamin@benjamin-ThinkPad-T470s:~/Desktop/TBSW/tbsw/workspace$ python example.py 
[INFO] Misalign position of plane with ID 0
[INFO] Misalign position of plane with ID 1
[INFO] Misalign position of plane with ID 2
[INFO] Misalign position of plane with ID 3
[INFO] Misalign position of plane with ID 4
[INFO] Misalign position of plane with ID 5
[INFO] Misalign position of plane with ID 21
[INFO] Misalign position of plane with ID 7
[INFO] Misalign position of plane with ID 6
[INFO] Starting to simulate a test beam run ...
[INFO] Marlin sim.xml is done
('[INFO] Created new caltag ', 'default')
('[INFO] Created new caltag ', 'simulation')
[INFO] Starting to calibrate file /home/benjamin/Desktop/TBSW/tbsw/workspace/simrun.slcio ...
[INFO] Marlin mask_path.xml is done
[INFO] Marlin clusterizer_path.xml is done
[INFO] Marlin correlator_path.xml is done
[INFO] Marlin prealigner_path.xml is done
[INFO] Marlin aligner_path.xml is done
[INFO] Marlin aligner_path.xml is done
[INFO] Marlin aligner_path.xml is done
[INFO] Marlin dqm_path.xml is done
[INFO] Marlin preclustercal_path.xml is done
[INFO] Marlin aligner_db_path.xml is done
[INFO] Marlin aligner_db_path.xml is done
[INFO] Marlin clustercal_path.xml is done
[INFO] Marlin clustercal_path.xml is done
[INFO] Marlin clustercal_path.xml is done
[INFO] Marlin clustercal_path.xml is done
[INFO] Marlin clustercal_path.xml is done
[INFO] Marlin clustercal_path.xml is done
[INFO] Marlin aligner_db_path.xml is done
[INFO] Marlin aligner_db_path.xml is done
[INFO] Marlin dqm_db_path.xml is done
('[INFO] Created new caltag ', 'simrun-test')
[INFO] Done processing file /home/benjamin/Desktop/TBSW/tbsw/workspace/simrun.slcio
[INFO] Starting to reconstruct file /home/benjamin/Desktop/TBSW/tbsw/workspace/simrun.slcio ...
[INFO] Marlin reco_path.xml is done
[INFO] Done processing file /home/benjamin/Desktop/TBSW/tbsw/workspace/simrun.slcio

```

All calibration data is collected in the folder 'workspace/localDB' and root files with reconstructed tracks and hits are collected in 
the folder 'workspace/root-files'. The structure of the root trees for tracks hits and evenst is explained in the next section. The example
script creates two root files with simple residual histograms ('Example-Residuals.root') and DUT hit efficiency plots ('Example-Efficiency.root')
to show what can be done using the hit and track trees.   


# What is the output data format of tbsw?

The tbsw framework uses Marlin[1] and LCIO[2] to organize the calibration and reconstruction of test beam data into a ```Path``` of small steps 
handled by so called ```Processors```. In order to read in raw files from EUDAQ a part of the EUDAQ 1 code base [3] is hosted to read files and 
convert event raw data into LCIO collections for further processing. The reconstruction path from the example script can illustrate the situation:

```
#!python

#Create a path for full reconstruction of telescope data
reco_path = Env.create_path('reco_path')

# Read input data from simulated raw file in lcio format 
# The file contains zero suppressed digits for all sensor
reco_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents, 'LCIOInputFiles': rawfile  }) 

# Create a precessor for clustering of telescope data   
m26clust = Processor(name="M26Clusterizer",proctype="PixelClusterizer")   
  
# Now you can set various parameters of the processor 
# You can obtain a full list of all processor types with  
# all their parameters and default values by executing 
# the following command in a shell: 
#  
# $MARLIN/bin/Marlin -x 
#
# All parameters can be overriding from the python script

# Set the path to the hotpixel mask (generated in the calibration before)
m26clust.param("NoiseDBFileName","localDB/NoiseDB-M26.root")

# Set the name of the input collection containing the m26 digits
m26clust.param("SparseDataCollectionName","zsdata_m26")

# Set the name of the output collection with with clusters
m26clust.param("ClusterCollectionName","zscluster_m26")

# Add the processor to the path   
reco_path.add_processor(m26clust)  

# Create a processor for clustering pxd data 
pxdclust = Processor(name="PXDClusterizer",proctype="PixelClusterizer")   
pxdclust.param("NoiseDBFileName","localDB/NoiseDB-PXD.root")
pxdclust.param("SparseDataCollectionName","zsdata_pxd")
pxdclust.param("ClusterCollectionName","zscluster_pxd")
reco_path.add_processor(pxdclust)  

# Create center of gravity hits from m26 clusters
m26hitmaker = Processor(name="M26CogHitMaker",proctype="CogHitMaker")
m26hitmaker.param("ClusterCollection","zscluster_m26")
m26hitmaker.param("HitCollectionName","hit_m26")
reco_path.add_processor(m26hitmaker)

# Create center of gravity hits for pxd clusters
pxdhitmaker = Processor(name="PXDCogHitMaker",proctype="CogHitMaker")
pxdhitmaker.param("ClusterCollection","zscluster_pxd")
pxdhitmaker.param("HitCollectionName","hit_pxd")
reco_path.add_processor(pxdhitmaker) 

# Make similar steps for the data from the FEI4 timing plane.
# You can easily extend to include more sensors. 

# Create (find/fit)tracks from M26 and FEI4 and extrapolate to PXD.
trackfinder = Processor(name="TrackFinder",proctype="FastTracker")

# Use hits from the following collections. You can add additional collections 
# as needed.  
trackfinder.param("InputHitCollectionNameVec","hit_m26 hit_fei4")

# Set the file from which to take the alignment. Computed during calibration.
trackfinder.param("AlignmentDBFileName","localDB/alignmentDB.root")

# PXD hits (plane 3 in geometry file) are explicitely from 
# the trackfinding to avoid bias in measuring PXD efficiency 
# resolution.  
# This is not strictly necessary here, because the PXD hits are not 
# in the input hit collection list. But one can use this to measure  
# a specific M26 plane, for example
trackfinder.param("ExcludeDetector", "3")
# Add to path 
reco_path.add_processor(trackfinder)  

# Create a processor to compute out root file for a specific device under
# test (DUT). Here, the DUT is the Belle II sensor at plane 3. 
pxd_analyzer = Processor(name="PXDAnalyzer",proctype="PixelDUTAnalyzer")
pxd_analyzer.param("HitCollection","hit_pxd")  
pxd_analyzer.param("DigitCollection","zsdata_pxd")
pxd_analyzer.param("AlignmentDBFileName","localDB/alignmentDB.root")
pxd_analyzer.param("DUTPlane","3")
reco_path.add_processor(pxd_analyzer)   

```

The last processor in this reconstruction path create the output root data for studying the performance of the device under test. 
The root file Histos-PXD.root in workspace/root-files will contain three TTree's called Hit, Track and Event. These trees will 
be filled with reference tracks interpolated to the third plane along the beam line, i.e. the position of the DUT sensor in our 
example. The main advantage of this choice is that we can always work in the local uvw coordinates of our device under test. 

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
int nDutDigits;                   // Number of DUT digits in the same event as hit   
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

// Hit Tree, filled once per DUT cluster 
_rootHitTree = new TTree("Hit","Hit info");
_rootHitTree->Branch("iRun"            ,&_rootRunNumber        ,"iRun/I");
_rootHitTree->Branch("iEvt"            ,&_rootEventNumber      ,"iEvt/I");
_rootHitTree->Branch("sensorID"        ,&_rootSensorID       ,"sensorID/I");
_rootHitTree->Branch("nTelTracks"      ,&_rootNTelTracks       ,"nTelTracks/I"); 
_rootHitTree->Branch("nDutDigits"        ,&_rootNDUTDigits          ,"nDutDigits/I");
_rootHitTree->Branch("clusterQuality"  ,&_rootHitQuality   ,"clusterQuality/I");
_rootHitTree->Branch("u_hit"           ,&_rootHitU             ,"u_hit/D");
_rootHitTree->Branch("v_hit"           ,&_rootHitV             ,"v_hit/D");     
_rootHitTree->Branch("clusterCharge"   ,&_rootHitClusterCharge ,"clusterCharge/D");
_rootHitTree->Branch("seedCharge"      ,&_rootHitSeedCharge       ,"seedCharge/D");
_rootHitTree->Branch("sizeU"           ,&_rootHitSizeU        ,"sizeU/I");
_rootHitTree->Branch("sizeV"           ,&_rootHitSizeV        ,"sizeV/I");
_rootHitTree->Branch("size"            ,&_rootHitSize         ,"size/I");
_rootHitTree->Branch("hasTrack"        ,&_rootHitHasTrack         ,"hasTrack/I");  
_rootHitTree->Branch("hasTrackWithRefHit", &_rootHitHasTrackWithRefHit,"hasTrackWithRefHit/I"); 
_rootHitTree->Branch("u_fit"           ,&_rootHitFitU             ,"u_fit/D");
_rootHitTree->Branch("v_fit"           ,&_rootHitFitV             ,"v_fit/D"); 
_rootHitTree->Branch("dudw_fit"        ,&_rootHitFitdUdW          ,"dudw_fit/D");
_rootHitTree->Branch("dvdw_fit"        ,&_rootHitFitdVdW          ,"dvdw_fit/D");    
_rootHitTree->Branch("u_fiterr"        ,&_rootHitFitErrorU        ,"u_fiterr/D");
_rootHitTree->Branch("v_fiterr"        ,&_rootHitFitErrorV        ,"v_fiterr/D");   
_rootHitTree->Branch("pull_resu"       ,&_rootHitPullResidualU    ,"pull_resu/D");
_rootHitTree->Branch("pull_resv"       ,&_rootHitPullResidualV    ,"pull_resv/D");  
_rootHitTree->Branch("cellU_fit"       ,&_rootHitFitCellU           ,"cellU_fit/I");
_rootHitTree->Branch("cellV_fit"       ,&_rootHitFitCellV           ,"cellV_fit/I");
_rootHitTree->Branch("cellU_hit"       ,&_rootHitCellU           ,"cellU_hit/I");
_rootHitTree->Branch("cellV_hit"       ,&_rootHitCellV           ,"cellV_hit/I");
_rootHitTree->Branch("startCellU_cluster"       ,&_rootClusterStartCellU           ,"startCellU_cluster/I");
_rootHitTree->Branch("srartCellV_cluster"       ,&_rootClusterStartCellV           ,"startCellV_cluster/I");
_rootHitTree->Branch("cellUCenter_fit" ,&_rootHitFitCellUCenter  ,"cellUCenter_fit/D");
_rootHitTree->Branch("cellVCenter_fit" ,&_rootHitFitCellVCenter  ,"cellVCenter_fit/D");                                      
_rootHitTree->Branch("trackChi2"       ,&_rootHitTrackChi2      ,"trackChi2/D");
_rootHitTree->Branch("trackNdof"       ,&_rootHitTrackNDF       ,"trackNdof/I");
_rootHitTree->Branch("trackNHits"      ,&_rootHitTrackNHits     ,"trackNHits/I");  
_rootHitTree->Branch("momentum"        ,&_rootHitFitMomentum      ,"momentum/D");    
_rootHitTree->Branch("localChi2"       ,&_rootHitLocalChi2       ,"localChi2/D"); 
_rootHitTree->Branch("pixeltype"       ,&_rootHitPixelType     ,"pixeltype/I");  
   
// Track tree, filled once per reference track 
_rootTrackTree = new TTree("Track","Track info");
_rootTrackTree->Branch("iRun"            ,&_rootRunNumber        ,"iRun/I");
_rootTrackTree->Branch("iEvt"            ,&_rootEventNumber      ,"iEvt/I");
_rootTrackTree->Branch("sensorID"        ,&_rootSensorID         ,"sensorID/I");
_rootTrackTree->Branch("nTelTracks"      ,&_rootNTelTracks       ,"nTelTracks/I"); 
_rootTrackTree->Branch("nDutDigits"      ,&_rootNDUTDigits       ,"nDutDigits/I");
_rootTrackTree->Branch("hasHit"          ,&_rootTrackHasHit         ,"hasHit/I");
_rootTrackTree->Branch("hasRefHit"       ,&_rootTrackWithRefHit     ,"hasRefHit/I");
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
_rootTrackTree->Branch("pixeltype"       ,&_rootTrackPixelType      ,"pixeltype/I");     
      
// Event tree, filled once per event 
_rootEventTree = new TTree("Event","Event info");
_rootEventTree->Branch("iRun"            ,&_rootRunNumber      ,"iRun/I");
_rootEventTree->Branch("iEvt"            ,&_rootEventNumber    ,"iEvt/I");
_rootEventTree->Branch("sensorID"        ,&_rootSensorID       ,"sensorID/I");   
_rootEventTree->Branch("nTelTracks"      ,&_rootNTelTracks     ,"nTelTracks/I"); 
_rootEventTree->Branch("nDutDigits"      ,&_rootNDUTDigits     ,"nDutDigits/I");

```

The plotting of root trees is well known in the HEP community. The script tbsw_example.py contains some instructive examples
based on PyRoot showing how to get residuals and efficiencies plots from the root trees. 
 

Have fun with test beams ;)  

benjamin.schwenker@phys.uni-goettingen.de


References: 

[1] F. Gaede, Marlin and LCCD: Software tools for the ILC, Nucl.Instrum.Meth. A559 (2006) 177–180.

[2] S. Aplin, J. Engels, F. Gaede, N. A. Graf, T. Johnson, et al., LCIO: A Persistency Framework and Event Data Model for HEP

[3] https://eudaq.github.io/
