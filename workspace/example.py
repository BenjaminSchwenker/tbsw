"""
This is an example script to demonstrate how tbsw can be used to analyze test beam 
data using Python scripts.

The script below simulates a test beam experiment where charged tracks cross a misaligned
pixel telescope containing six Mimosa 26 detector planes forming a reference  telescope.
The reference telescope has two arms with three sensors. A BELLE II PXD sensor is installed in 
the center of the reference telescope as device under test. A reference FEI4 plane is used
as a timing plane. 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *

# Path to steering files 
# The folder steerfiles contains one or more gear file describing the 
# nominal telescope geometry. 
steerfiles = 'steering-files/depfet-tb/'
# Select the name of a gearfile to use from the steerfiles folder  
gearfile = 'gear_desy_W11OF2.xml'
# Select filename for the simulated test beam run  
rawfile = os.getcwd() + '/simrun.slcio'
# Number of events to simulate 
nevents = 3000000
# Beam energy in GeV
energy = 4

#Parameters x,y,z,alpha,beta,gamma for simulation of misalignment
#Position parameters in mm
mean_list=[0.0,0.0,0.0,0.0,0.0,0.0] 
sigma_list=[0.1,0.1,0.1,0.1,0.1,0.1]

# List of sensor ids and modes, which are excluded during misalignment
sensorexception_list=[] 
modeexception_list=[]

def add_clusterizers(path):
  """
  Adds clusterizers to the path
  """  
    
  m26clust = Processor(name="M26Clusterizer",proctype="PixelClusterizer")   
  m26clust.param("NoiseDBFileName","localDB/NoiseDB-M26.root")
  m26clust.param("SparseDataCollectionName","zsdata_m26")
  m26clust.param("ClusterCollectionName","zscluster_m26")
  m26clust.param("SparseClusterCut",0)
  m26clust.param("SparseSeedCut", 0)
  m26clust.param("SparseZSCut", 0)   
  path.add_processor(m26clust)  

  fei4clust = Processor(name="FEI4Clusterizer",proctype="PixelClusterizer")   
  fei4clust.param("NoiseDBFileName","localDB/NoiseDB-FEI4.root")
  fei4clust.param("SparseDataCollectionName","zsdata_fei4")
  fei4clust.param("ClusterCollectionName","zscluster_fei4")
  fei4clust.param("SparseClusterCut",0)
  fei4clust.param("SparseSeedCut", 0)
  fei4clust.param("SparseZSCut", 0)   
  path.add_processor(fei4clust)  

  pxdclust = Processor(name="PXDClusterizer",proctype="PixelClusterizer")   
  pxdclust.param("NoiseDBFileName","localDB/NoiseDB-PXD.root")
  pxdclust.param("SparseDataCollectionName","zsdata_pxd")
  pxdclust.param("ClusterCollectionName","zscluster_pxd")
  pxdclust.param("SparseClusterCut",4)
  pxdclust.param("SparseSeedCut", 4)
  pxdclust.param("SparseZSCut", 4)   
  path.add_processor(pxdclust)  
  
  return path

def add_hitmakers(path):
  """
  Adds hitmakers to the path
  """  

  m26hitmaker = Processor(name="M26CogHitMaker",proctype="CogHitMaker")
  m26hitmaker.param("ClusterCollection","zscluster_m26")
  m26hitmaker.param("HitCollectionName","hit_m26")
  m26hitmaker.param("SigmaU1",0.0037)
  m26hitmaker.param("SigmaU2",0.0033)
  m26hitmaker.param("SigmaU3",0.0050)
  m26hitmaker.param("SigmaV1",0.0037)
  m26hitmaker.param("SigmaV2",0.0033)
  m26hitmaker.param("SigmaV3",0.0050)
  path.add_processor(m26hitmaker)

  fei4hitmaker = Processor(name="FEI4CogHitMaker",proctype="CogHitMaker")
  fei4hitmaker.param("ClusterCollection","zscluster_fei4")
  fei4hitmaker.param("HitCollectionName","hit_fei4")
  fei4hitmaker.param("SigmaU1",0.072)
  fei4hitmaker.param("SigmaU2",0.072)
  fei4hitmaker.param("SigmaU3",0.072)
  fei4hitmaker.param("SigmaV1",0.0144)
  fei4hitmaker.param("SigmaV2",0.0144)
  fei4hitmaker.param("SigmaV3",0.0144)
  path.add_processor(fei4hitmaker)

  pxdhitmaker = Processor(name="PXDCogHitMaker",proctype="CogHitMaker")
  pxdhitmaker.param("ClusterCollection","zscluster_pxd")
  pxdhitmaker.param("HitCollectionName","hit_pxd")
  pxdhitmaker.param("SigmaU1",0.0134)
  pxdhitmaker.param("SigmaU2",0.0077)
  pxdhitmaker.param("SigmaU3",0.0077)
  pxdhitmaker.param("SigmaV1",0.024)
  pxdhitmaker.param("SigmaV2",0.014)
  pxdhitmaker.param("SigmaV3",0.014)
  path.add_processor(pxdhitmaker)
  
  return path

def create_sim_path(Env):
  """
  Returns a list of tbsw path objects to simulate a test beam run 
  """
  
  sim_path = Env.create_path('sim')
  sim_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents})   
  infosetter = Processor(name="InfoSetter", proctype='EventInfoSetter')
  infosetter.param("RunNumber","0")
  sim_path.add_processor(infosetter)
  
  gun = Processor(name="ParticleGun",proctype="ParticleGunGenerator")
  gun.param("BeamIntensity","2000")
  gun.param("BeamMomentum","4.0")
  gun.param("BeamVertexX","0")
  gun.param("BeamVertexY","0")
  gun.param("BeamVertexZ","-10")  
  gun.param("BeamVertexXSigma","7")
  gun.param("BeamVertexYSigma","7")   
  gun.param("PDG","11")
  gun.param("ParticleCharge","-1")
  gun.param("ParticleMass","0.000511")
  sim_path.add_processor(gun)
  
  fastsim = Processor(name="FastSim",proctype="FastSimulation")
  fastsim.param("AlignmentDBFileName","localDB/alignmentDB.root")
  fastsim.param("ScatterModel","0")
  fastsim.param("DoEnergyLossStraggling","true")
  fastsim.param("DoFractionalBetheHeitlerEnergyLoss","true")
  sim_path.add_processor(fastsim)
  
  tlu = Processor(name="TLU",proctype="TriggerGenerator")
  tlu.param("FakeTriggerPeriod","0")
  tlu.param("ScinitNo1", "0 -5 -5 5 5")
  tlu.param("ScinitNo2", "")
  tlu.param("ScinitNo3", "")
  tlu.param("ScinitNo4", "") 
  sim_path.add_processor(tlu)
   
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
  sim_path.add_processor(m26digi)
  
  pxddigi = Processor(name="DEPFETDigitizer",proctype="SiPixDigitizer")
  pxddigi.param("DigitCollectionName","zsdata_pxd")  
  pxddigi.param("FrontEndType","0")
  pxddigi.param("ADCBits","8")
  pxddigi.param("ADCRange","300000")
  pxddigi.param("ElectronicNoise","300")
  pxddigi.param("FilterIDs","6")
  pxddigi.param("IntegrationWindow","true")
  pxddigi.param("StartIntegration","0")
  pxddigi.param("StopIntegration","20000")
  pxddigi.param("uSideBorderLength","6")
  pxddigi.param("vSideBorderLength","6")
  pxddigi.param("ZSThreshold", "5")
  sim_path.add_processor(pxddigi)
  
  fei4digi = Processor(name="FEI4Digitizer",proctype="SiPixDigitizer")
  fei4digi.param("DigitCollectionName","zsdata_fei4")  
  fei4digi.param("FrontEndType","0") 
  fei4digi.param("ADCBits","2")
  fei4digi.param("ADCRange","300000")
  fei4digi.param("ElectronicNoise","300")
  fei4digi.param("FilterIDs","21")
  fei4digi.param("IntegrationWindow","true")
  fei4digi.param("StartIntegration","0")
  fei4digi.param("StopIntegration","20")
  fei4digi.param("uSideBorderLength","5")
  fei4digi.param("vSideBorderLength","5")
  fei4digi.param("ZSThreshold", "0")
  sim_path.add_processor(fei4digi)
   
  lciooutput = Processor(name="LCIOOutput",proctype="LCIOOutputProcessor")
  lciooutput.param("LCIOOutputFile",rawfile)
  lciooutput.param("LCIOWriteMode","WRITE_NEW")  
  sim_path.add_processor(lciooutput)
  
  return [ sim_path]

def create_calibration_path(Env):
  """
  Returns a list of tbsw path objects to calibrate the tracking telescope
  """
  
  # Calibrations are organized in a sequence of calibration paths. 
  # The calibration paths are collected in a list for later execution
  calpaths = []

  # Create path for detector level masking of hot channels 
  mask_path = Env.create_path('mask_path')
  mask_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : -1, 'LCIOInputFiles': rawfile }) 
   
  m26hotpixelkiller = Processor(name="M26HotPixelKiller",proctype="HotPixelKiller",)
  m26hotpixelkiller.param("InputCollectionName", "zsdata_m26")
  m26hotpixelkiller.param("MaxOccupancy", 0.001)
  m26hotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-M26.root")
  m26hotpixelkiller.param("OfflineZSThreshold", 0)
  mask_path.add_processor(m26hotpixelkiller)
   
  fei4hotpixelkiller = Processor(name="FEI4HotPixelKiller", proctype="HotPixelKiller")
  fei4hotpixelkiller.param("InputCollectionName", "zsdata_fei4")
  fei4hotpixelkiller.param("MaxOccupancy", 0.001)
  fei4hotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-FEI4.root")
  fei4hotpixelkiller.param("OfflineZSThreshold", 0)
  mask_path.add_processor(fei4hotpixelkiller)
   
  pxdhotpixelkiller = Processor(name="PXDHotPixelKiller", proctype="HotPixelKiller")
  pxdhotpixelkiller.param("InputCollectionName", "zsdata_pxd")
  pxdhotpixelkiller.param("MaxOccupancy", 0.001)
  pxdhotpixelkiller.param("NoiseDBFileName", "localDB/NoiseDB-PXD.root")
  pxdhotpixelkiller.param("OfflineZSThreshold", 0)
  mask_path.add_processor(pxdhotpixelkiller)  
  
  # Add path for masking
  calpaths.append(mask_path)  
  
  # Create path for detector level creation of hits
  hit_path = Env.create_path('hit_path')
  hit_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' :  nevents, 'LCIOInputFiles': rawfile  }) 
  
  hit_path = add_clusterizers(hit_path)    
  hit_path = add_hitmakers(hit_path) 
  
  hitdqm = Processor(name="RawDQM",proctype="RawHitDQM")
  hitdqm.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_pxd")  
  hitdqm.param("RootFileName","RawDQM.root")
  hit_path.add_processor(hitdqm)  
  
  correlator = Processor(name="TelCorrelator", proctype="Correlator")
  correlator.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_pxd hit_h5")
  correlator.param("AlignmentDBFileName", "localDB/alignmentDB.root")
  correlator.param("OutputRootFileName","XCorrelator.root")
  correlator.param("ReferencePlane","0")
  correlator.param("ParticleCharge","-1")
  correlator.param("ParticleMass","0.000511")
  correlator.param("ParticleMomentum", energy)
  hit_path.add_processor(correlator)  
  
  lciooutput = Processor(name="LCIOOutput",proctype="LCIOOutputProcessor")
  lciooutput.param("LCIOOutputFile","tmp.slcio")
  lciooutput.param("LCIOWriteMode","WRITE_NEW")
  hit_path.add_processor(lciooutput)  
  
  # Finished with path for hits
  calpaths.append(hit_path)  
    
  # Create path for pre alignment with loose cut track sample 
  prealigner_path = Env.create_path('prealigner_path')
  prealigner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 200000, 'LCIOInputFiles': "tmp.slcio" })  
  
  trackfinder_loosecut = Processor(name="AlignTF_LC",proctype="FastTracker")
  trackfinder_loosecut.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_pxd hit_h5")
  trackfinder_loosecut.param("AlignmentDBFileName","localDB/alignmentDB.root")
  trackfinder_loosecut.param("ExcludeDetector", "")
  trackfinder_loosecut.param("MaxTrackChi2", 10000000)
  trackfinder_loosecut.param("MaximumGap", 1)
  trackfinder_loosecut.param("MinimumHits",7)
  trackfinder_loosecut.param("OutlierChi2Cut", 100000000)
  trackfinder_loosecut.param("ParticleCharge","-1")
  trackfinder_loosecut.param("ParticleMass","0.000511")
  trackfinder_loosecut.param("ParticleMomentum", energy)
  trackfinder_loosecut.param("SingleHitSeeding", "0")
  trackfinder_loosecut.param("MaxResidualU","0.5")
  trackfinder_loosecut.param("MaxResidualV","0.5")
  prealigner_path.add_processor(trackfinder_loosecut)
  
  prealigner = Processor(name="PreAligner",proctype="KalmanAligner")
  prealigner.param("AlignmentDBFileName","localDB/alignmentDB.root")
  prealigner.param('ErrorsShiftX' , '0 10 10 10 10 10 0 10 10')
  prealigner.param('ErrorsShiftY' , '0 10 10 10 10 10 0 10 10')
  prealigner.param('ErrorsShiftZ' , '0 0 0 0 0 0 0 0 0')
  prealigner.param('ErrorsAlpha'  , '0 0 0 0 0 0 0 0 0')
  prealigner.param('ErrorsBeta'   , '0 0 0 0 0 0 0 0 0')
  prealigner.param('ErrorsGamma'  , '0 0.01 0.01 0.01 0.01 0.01 0 0.01 0.01')
  prealigner_path.add_processor(prealigner)  
  
  # Finished with path for prealigner
  calpaths.append(prealigner_path)  
  
  # Create path for alignment with tight cut track sample 
  aligner_path = Env.create_path('aligner_path')
  aligner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 200000, 'LCIOInputFiles': "tmp.slcio" })  
  
  trackfinder_tightcut = Processor(name="AlignTF_TC",proctype="FastTracker")
  trackfinder_tightcut.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_pxd hit_h5")
  trackfinder_tightcut.param("AlignmentDBFileName","localDB/alignmentDB.root")
  trackfinder_tightcut.param("ExcludeDetector", "")
  trackfinder_tightcut.param("MaxTrackChi2", 100)
  trackfinder_tightcut.param("MaximumGap", 1)
  trackfinder_tightcut.param("MinimumHits",7)
  trackfinder_tightcut.param("OutlierChi2Cut", 20)
  trackfinder_tightcut.param("ParticleCharge","-1")
  trackfinder_tightcut.param("ParticleMass","0.000511")
  trackfinder_tightcut.param("ParticleMomentum", energy)
  trackfinder_tightcut.param("SingleHitSeeding", "0")
  trackfinder_tightcut.param("MaxResidualU","0.4")
  trackfinder_tightcut.param("MaxResidualV","0.4")
  aligner_path.add_processor(trackfinder_tightcut)
   
  aligner = Processor(name="Aligner",proctype="KalmanAligner")
  aligner.param("AlignmentDBFileName","localDB/alignmentDB.root")
  aligner.param('ErrorsShiftX' , '0 10 10 10 10 10 0 10 10' )
  aligner.param('ErrorsShiftY' , '0 10 10 10 10 10 0 10 10')
  aligner.param('ErrorsShiftZ' , '0 10 10 10 10 10 0 10 10')
  aligner.param('ErrorsAlpha'  , '0 0 0 0.01 0 0 0 0.01 0.01')
  aligner.param('ErrorsBeta'   , '0 0 0 0.01 0 0 0 0.01 0.01')
  aligner.param('ErrorsGamma'  , '0 0.01 0.01 0.01 0.01 0.01 0 0.01 0.01')
  aligner_path.add_processor(aligner)    
  
  # Finished with path for aligner
  # Repeat this 3x
  calpaths.append(aligner_path)  
  calpaths.append(aligner_path)
  calpaths.append(aligner_path)
   
  # Creeate path for some track based dqm using current calibrations
  dqm_path = Env.create_path('dqm_path')
  dqm_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 200000, 'LCIOInputFiles': "tmp.slcio" })  
  dqm_path.add_processor(trackfinder_tightcut)  
  
  teldqm = Processor(name="TelescopeDQM", proctype="TrackFitDQM") 
  teldqm.param("AlignmentDBFileName","localDB/alignmentDB.root")
  teldqm.param("RootFileName","TelescopeDQM.root")
  dqm_path.add_processor(teldqm)  
  
  # Finished with path for teldqm
  calpaths.append(dqm_path)
  
  return calpaths


def create_reco_path(Env):
  """
  Returns a list of tbsw path objects to reconstruct a test beam run 
  """
   
  reco_path = Env.create_path('reco_path')
  reco_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents, 'LCIOInputFiles': rawfile  }) 
  
  # Create path for all reconstruction up to hits
  reco_path = add_clusterizers(reco_path)    
  reco_path = add_hitmakers(reco_path) 
  
  trackfinder = Processor(name="TrackFinder",proctype="FastTracker")
  trackfinder.param("InputHitCollectionNameVec","hit_m26 hit_fei4")
  trackfinder.param("AlignmentDBFileName","localDB/alignmentDB.root")
  trackfinder.param("ExcludeDetector", "3 8")
  trackfinder.param("MaxTrackChi2", "100")
  trackfinder.param("MaximumGap", "1")
  trackfinder.param("MinimumHits","6")
  trackfinder.param("OutlierChi2Cut", "20")
  trackfinder.param("ParticleCharge","-1")
  trackfinder.param("ParticleMass","0.000511")
  trackfinder.param("ParticleMomentum", energy)
  trackfinder.param("SingleHitSeeding", "0")
  trackfinder.param("MaxResidualU","0.4")
  trackfinder.param("MaxResidualV","0.4")
  reco_path.add_processor(trackfinder)  

  pxd_analyzer = Processor(name="PXDAnalyzer",proctype="PixelDUTAnalyzer")
  pxd_analyzer.param("HitCollection","hit_pxd")  
  pxd_analyzer.param("DigitCollection","zsdata_pxd")
  pxd_analyzer.param("AlignmentDBFileName","localDB/alignmentDB.root")
  pxd_analyzer.param("DUTPlane","3")
  pxd_analyzer.param("ReferencePlane","7")
  pxd_analyzer.param("MaxResidualU","0.2")
  pxd_analyzer.param("MaxResidualU","0.2")
  pxd_analyzer.param("RootFileName","Histos-PXD.root")
  reco_path.add_processor(pxd_analyzer)   
  
  return [ reco_path ]  
 
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

  # Misalign gear file
  randomize_telescope(gearfile=SimObj.get_filename(gearfile), mean_list=mean_list, sigma_list=sigma_list, sensorexception_list=sensorexception_list, modeexception_list=modeexception_list)
   
  # Run simulation to create rawfile with simulated digits 
  SimObj.simulate(path=simpath)  

  # Export clusterDB created from truth hits
  SimObj.export_caltag(caltag='simulation')
   

def calibrate(params):
  """
  Calibrates an misaligned tracking telescope from run data. 
  Creates a folder localDB/caltag in workspace containing 
  calibration results.   
  """ 
  
  rawfile, steerfiles, gearfile, caltag = params
  
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
  CalObj = Calibration(steerfiles=steerfiles, name=caltag + '-cal') 
  
  # Create list of calibration steps 
  calpath = create_calibration_path(CalObj)
  
  # Run the calibration steps 
  CalObj.calibrate(path=calpath,ifile=rawfile,caltag=caltag)  
  
  
def reconstruct(params):
  """
  Reconstruct raw data from a tracking telescope. 
  """ 
  
  rawfile, steerfiles, gearfile, caltag = params
   
  # Reconsruct the rawfile using the caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=caltag + '-reco' )

  # Create reconstuction path
  recopath = create_reco_path(RecObj)  
  
  # Run the reconstuction  
  RecObj.reconstruct(path=recopath,ifile=rawfile,caltag=caltag) 

if __name__ == '__main__':
  
  # Tag for calibration data
  caltag = os.path.splitext(os.path.basename(rawfile))[0] + '-test'
  
  params = ( rawfile, steerfiles, gearfile, caltag )
  
  # Create a simulated rawfile 
  simulate( params )
  
  # Calibrate the telescope and reconstruct the rawfile 
  calibrate( params )
  
  # Reconstruct the rawfile 
  reconstruct( params )
  
  # Make a list of root files containing reconstructed trees 
  # for tracks / hits / events
  trackfile = 'root-files/Histos-PXD-simrun-test-reco.root'  
  
  # Plot DUT residuals and cluster signal histograms from the 'Hit'
  # tree in the workspace. 
  ofile = 'Example-Residuals.root'
  residuals.plot(inputfilename=trackfile, histofilename=ofile, basecut = "hasTrack==0", nbins=201, urange=400, vrange=400)
      
  # Plot DUT hit efficiency histograms from the 'Track' tree 
  # in the workspace. 
  ofile = 'Example-Efficiency.root' 
  efficiency.plot(inputfilename=trackfile, histofilename=ofile, basecut="trackNHits==7 && nDutDigits>=0", matchcut="hasHit==0", uaxis=(250,0,250), vaxis=(768,0,768))
    
    
    
