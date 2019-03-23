"""
This is an example script to quickly check the correctness of the track fitting in a test beam
environment. 

The script below simulates a test beam experiment where charged tracks cross a misaligned
pixel telescope containing six Mimosa 26 detector planes forming a reference  telescope.
The reference telescope has two arms with three sensors. A BELLE II PXD sensor is installed in 
the center of the reference telescope as device under test. A reference FEI4 plane is used
as a timing plane. 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

import tbsw 
import os

# Path to steering files 
# The folder steerfiles contains one or more gear file describing the 
# nominal telescope geometry. 
steerfiles = 'steering-files/depfet-tb/'
# Select the name of a gearfile to use from the steerfiles folder  
gearfile = 'gear_desy_W11OF2_perp_geoid2.xml'
# Select filename for the simulated test beam run  
rawfile = os.getcwd() + '/simrun.slcio'
# Number of events to simulate 
nevents = 300000
# Beam energy in GeV
energy = 4
# Flag to skip calibration and use truth misalgnment 
useTruthMisalignment = True



#Parameters x,y,z,alpha,beta,gamma for simulation of misalignment
#Position parameters in mm and degree for rotations
mean_list=[0.0,0.0,0.0,0.0,0.0,0.0] 
sigma_list=[0.1,0.1,0.1,0.3,0.3,1.5]

# List of sensor ids and modes ('positionX','positionY','positionZ','alpha','beta','gamma'), which are excluded during misalignment
sensorexception_list=[] 
modeexception_list=[]


def create_sim_path(Env):
  """
  Returns a list of tbsw path objects to simulate a test beam run 
  """
  
  sim_path = Env.create_path('sim')
  sim_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents})   
  infosetter = tbsw.Processor(name="InfoSetter", proctype='EventInfoSetter')
  infosetter.param("RunNumber","0")
  infosetter.param("DetectorName","EUTelescope")
  sim_path.add_processor(infosetter)
   
  geo_noalign = tbsw.Processor(name="Geo",proctype="Geometry")
  geo_noalign.param("AlignmentDBFilePath", "localDB/alignmentDB.root")
  geo_noalign.param("ApplyAlignment", "false")
  geo_noalign.param("OverrideAlignment", "true")
  sim_path.add_processor(geo_noalign)
  
  gun = tbsw.Processor(name="ParticleGun",proctype="ParticleGunGenerator")
  gun.param("BeamIntensity","2000")
  gun.param("BeamMomentum",energy)
  gun.param("BeamVertexX","0")
  gun.param("BeamVertexY","0")
  gun.param("BeamVertexZ","-10")  
  gun.param("BeamVertexXSigma","7")
  gun.param("BeamVertexYSigma","7")   
  gun.param("PDG","11")
  gun.param("ParticleCharge","-1")
  gun.param("ParticleMass","0.000511")
  sim_path.add_processor(gun)
  
  fastsim = tbsw.Processor(name="FastSim",proctype="FastSimulation")
  fastsim.param("ScatterModel","0")
  fastsim.param("DoEnergyLossStraggling","false")
  fastsim.param("DoFractionalBetheHeitlerEnergyLoss","false")
  sim_path.add_processor(fastsim)
  
  tlu = tbsw.Processor(name="TLU",proctype="TriggerGenerator")
  tlu.param("FakeTriggerPeriod","0")
  tlu.param("ScinitNo1", "0 -5 -5 5 5")
  tlu.param("ScinitNo2", "")
  tlu.param("ScinitNo3", "")
  tlu.param("ScinitNo4", "") 
  sim_path.add_processor(tlu)
    
  m26digi = tbsw.Processor(name="M26SimpleDigit",proctype="SmearingDigitizer")
  m26digi.param("FilterIDs", "0 1 2 3 4 5")
  m26digi.param("HitCollectionName", "hit_m26")
  m26digi.param("IntegrationWindow","true")
  m26digi.param("StartIntegration","0")
  m26digi.param("StopIntegration","100000")
  m26digi.param("ClusterSigmaU", "0.0032" )
  m26digi.param("ClusterSigmaV", "0.0032" )
  sim_path.add_processor(m26digi)
  
  pxddigi = tbsw.Processor(name="PXDSimpleDigit",proctype="SmearingDigitizer")
  pxddigi.param("FilterIDs", "6")
  pxddigi.param("HitCollectionName", "hit_pxd")
  pxddigi.param("IntegrationWindow","true")
  pxddigi.param("StartIntegration","0")
  pxddigi.param("StopIntegration","20000")
  pxddigi.param("ClusterSigmaU", "0.006" )
  pxddigi.param("ClusterSigmaV", "0.010" )
  sim_path.add_processor(pxddigi)
  
  feidigi = tbsw.Processor(name="FEISimpleDigit",proctype="SmearingDigitizer")
  feidigi.param("FilterIDs", "21")
  feidigi.param("HitCollectionName", "hit_fei4")
  feidigi.param("IntegrationWindow","true")
  feidigi.param("StartIntegration","0")
  feidigi.param("StopIntegration","20")
  feidigi.param("ClusterSigmaU", "0.030" )
  feidigi.param("ClusterSigmaV", "0.070" )
  sim_path.add_processor(feidigi)
  
  lciooutput = tbsw.Processor(name="LCIOOutput",proctype="LCIOOutputProcessor")
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
  
  # Create path for pre alignmnet and dqm based on hits
  correlator_path = Env.create_path('correlator_path')
  correlator_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 200000, 'LCIOInputFiles': rawfile })  

  geo = tbsw.Processor(name="Geo",proctype="Geometry")
  geo.param("AlignmentDBFilePath", "localDB/alignmentDB.root")
  geo.param("ApplyAlignment", "true")
  geo.param("OverrideAlignment", "true")
  correlator_path.add_processor(geo)
  
  hitdqm = tbsw.Processor(name="RawDQM",proctype="RawHitDQM")
  hitdqm.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_pxd")  
  hitdqm.param("RootFileName","RawDQM.root")
  correlator_path.add_processor(hitdqm)  

  correlator = tbsw.Processor(name="TelCorrelator", proctype="Correlator")
  correlator.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_pxd hit_h5")
  correlator.param("OutputRootFileName","XCorrelator.root")
  correlator.param("ReferencePlane","0")
  correlator.param("ParticleCharge","-1")
  correlator.param("ParticleMass","0.000511")
  correlator.param("ParticleMomentum", energy)
  correlator_path.add_processor(correlator)  
  
  # Finished with path for hit based pre alignment
  calpaths.append(correlator_path)  
    
  # Create path for pre alignment with loose cut track sample 
  prealigner_path = Env.create_path('prealigner_path')
  prealigner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 200000, 'LCIOInputFiles': rawfile })  
  prealigner_path.add_processor(geo)
  
  trackfinder_loosecut = tbsw.Processor(name="AlignTF_LC",proctype="FastTracker")
  trackfinder_loosecut.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_pxd hit_h5")
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
  
  prealigner = tbsw.Processor(name="PreAligner",proctype="KalmanAligner")
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
  aligner_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 200000, 'LCIOInputFiles': rawfile })  
  aligner_path.add_processor(geo)
 
  trackfinder_tightcut = tbsw.Processor(name="AlignTF_TC",proctype="FastTracker")
  trackfinder_tightcut.param("InputHitCollectionNameVec","hit_m26 hit_fei4 hit_pxd hit_h5")
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
   
  aligner = tbsw.Processor(name="Aligner",proctype="KalmanAligner")
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
   
  # Creeate path for some track based dqm using center of gravity hits
  dqm_path = Env.create_path('dqm_path')
  dqm_path.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 200000, 'LCIOInputFiles': rawfile })
  dqm_path.add_processor(geo)  
  dqm_path.add_processor(trackfinder_tightcut)  
   
  teldqm = tbsw.Processor(name="TelescopeDQM", proctype="TrackFitDQM") 
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
  
  geo = tbsw.Processor(name="Geo",proctype="Geometry")
  geo.param("AlignmentDBFilePath", "localDB/alignmentDB.root")
  geo.param("ApplyAlignment", "true")
  geo.param("OverrideAlignment", "true")
  reco_path.add_processor(geo)
   
  trackfinder = tbsw.Processor(name="TrackFinder",proctype="FastTracker")
  trackfinder.param("InputHitCollectionNameVec","hit_m26 hit_fei4")
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

  pxd_analyzer = tbsw.Processor(name="PXDAnalyzer",proctype="PixelDUTAnalyzer")
  pxd_analyzer.param("HitCollection","hit_pxd")  
  pxd_analyzer.param("DigitCollection","zsdata_pxd")
  pxd_analyzer.param("DUTPlane","3")
  pxd_analyzer.param("ReferencePlane","7")
  pxd_analyzer.param("MaxResidualU","0.2")
  pxd_analyzer.param("MaxResidualU","0.2")
  pxd_analyzer.param("RootFileName","Histos-PXD.root")
  reco_path.add_processor(pxd_analyzer)   
  
  fittester = tbsw.Processor(name="fittester", proctype="TrackFitValidation")
  fittester.param("RootFileName", "Validation.root")  
  reco_path.add_processor(fittester)
  
  return [ reco_path ]  
 
def simulate(params): 
  """
  Simulates a rawfile from a simulated test beam experiment
  Creates a folder tmp-runs/name-sim/ and populates it with 
  Marlin steering and logfiles.  
  """ 
  
  rawfile, steerfiles, gearfile, caltag = params
  
  # Create tmpdir to hold all steerfiles and log files 
  SimObj = tbsw.Simulation(steerfiles=steerfiles, name=os.path.splitext(os.path.basename(rawfile))[0] + '-sim' )
  
  # Create steerfiles for processing
  simpath = create_sim_path(SimObj)
  
  # Misalign gear file
  tbsw.gear.randomize_telescope(gearfile=SimObj.get_filename(gearfile), mean_list=mean_list, sigma_list=sigma_list, sensorexception_list=sensorexception_list, modeexception_list=modeexception_list)
   
  # Run simulation to create rawfile with simulated digits 
  SimObj.simulate(paths=simpath, caltag=caltag)  
    
def calibrate(params):
  """
  Calibrates an misaligned tracking telescope from run data. 
  Creates a folder localDB/caltag in workspace containing 
  calibration results.   
  """ 
  
  rawfile, steerfiles, gearfile, caltag = params
  
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
  CalObj = tbsw.Calibration(steerfiles=steerfiles, name=caltag + '-cal') 
  
  # Create list of calibration steps 
  calpaths = create_calibration_path(CalObj)
  
  # Run the calibration steps 
  CalObj.calibrate(paths=calpaths,ifile=rawfile,caltag=caltag)  
  
  
def reconstruct(params):
  """
  Reconstruct raw data from a tracking telescope. 
  """ 
  
  rawfile, steerfiles, gearfile, caltag = params
   
  # Reconsruct the rawfile using the caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = tbsw.Reconstruction(steerfiles=steerfiles, name=caltag + '-reco' )

  # Create reconstuction path
  recopath = create_reco_path(RecObj)  
  
  # Run the reconstuction  
  RecObj.reconstruct(paths=recopath,ifile=rawfile,caltag=caltag) 

if __name__ == '__main__':
  
  # Tag for calibration data
  caltag = os.path.splitext(os.path.basename(rawfile))[0] + '-test'
  
  params = ( rawfile, steerfiles, gearfile, caltag )
  
  # Create a simulated rawfile 
  simulate( params )
  
  if not useTruthMisalignment:
    # Align the telescope  
    calibrate( params )
  else: 
    print("Skipping alignment and using truth misaligned geometry instead")
  
  # Reconstruct the rawfile 
  reconstruct( params )
  
  # Make a list of root files containing reconstructed trees 
  # for tracks / hits / events
  trackfile = 'root-files/Histos-PXD-simrun-test-reco.root'  
  
  # Plot DUT residuals and cluster signal histograms from the 'Hit'
  # tree in the workspace. 
  ofile = 'Example-Residuals.root'
  tbsw.residuals.plot(inputfilename=trackfile, histofilename=ofile, basecut = "hasTrack==0", nbins=501, urange=400, vrange=400)
      
  # Plot DUT hit efficiency histograms from the 'Track' tree 
  # in the workspace. 
  ofile = 'Example-Efficiency.root' 
  tbsw.efficiency.plot(inputfilename=trackfile, histofilename=ofile, basecut="trackNHits==7 && nDutDigits>=0", matchcut="hasHit==0", uaxis=(250,0,250), vaxis=(768,0,768))
    
    
    
