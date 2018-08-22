"""
This is an example script to demonstrate how to use the triplett correlator to reconstruct
telescope+DUT data for cases when the telescope arms are far away and the normal correlator 
does not yield a acceptable pre-alignment. 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *

# Path to steering files 
# The folder steerfiles contains one or more gear file describing the 
# nominal telescope geometry. 
steerfiles = 'steering-files/depfet-H5-tb/'
# Select the name of a gearfile to use from the steerfiles folder  
gearfile = 'gear_desynov15_very_long.xml'
# Select filename for the simulated test beam run  
rawfile = os.getcwd() + '/simrun.slcio'
# Number of events to simulate 
nevents = 100000

#Parameters for simulation of misalignment
#Position parameters in mm
mean_list=[0.0,0.0,0.0,0.0,0.0,0.0] 
sigma_list=[0.1,0.1,0.1,0.1,0.1,0.1]

# List of sensor ids and modes, which are excluded during misalignment
#sensorexception_list=[11,5,0] 
#modeexception_list=['positionZ']

sensorexception_list=[] 
modeexception_list=[]

def create_sim_path(Env):
  """
  Returns a list of tbsw path objects to simulate a test beam run 
  """
  
  sim = Env.create_path('sim')
  sim.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents})   
  sim.add_processor(name="InfoSetter")
  sim.add_processor(name="ParticleGun")
  sim.add_processor(name="FastSim")
  sim.add_processor(name="TriggerGeneratorLargeSct")
  sim.add_processor(name="M26Digitizer")
  sim.add_processor(name="DEPFETDigitizer")
  sim.add_processor(name="LCIOOutput",params={"LCIOOutputFile" : rawfile })

  return [ sim ]

def create_calibration_path(Env):
  """
  Returns a list of tbsw path objects to calibrate the tracking telescope
  """
  
  hotpixelkiller = Env.create_path('hotpixelkiller')
  hotpixelkiller.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': rawfile })  
  hotpixelkiller.add_processor(name="M26HotPixelKiller")
  hotpixelkiller.add_processor(name="DEPHotPixelKiller")
  
  clusterizer = Env.create_path('clusterizer')
  clusterizer.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents, 'LCIOInputFiles': rawfile }) 
  clusterizer.add_processor(name="M26Clusterizer")
  clusterizer.add_processor(name="DEPClusterizer")
  clusterizer.add_processor(name="LCIOOutput")
  
  correlator = Env.create_path('correlator')
  correlator.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" }) 
  correlator.add_processor(name="M26CogHitMaker")
  correlator.add_processor(name="DEPCogHitMaker")
  correlator.add_processor(name="RawDQM") 
  correlator.add_processor(name="TelCorrelator")
  
  
  
  # Deal with first arm in one pass 
  kalman_aligner_triplet_0 = Env.create_path('kalman_aligner_triplet_0')
  kalman_aligner_triplet_0.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_triplet_0.add_processor(name="M26CogHitMaker")
  kalman_aligner_triplet_0.add_processor(name="AlignTF_LC", params={'ExcludeDetector' : '3 4 5 6', 'MinimumHits' : '3' })     
  kalman_aligner_triplet_0.add_processor(name="PreAligner", params={'ErrorsShiftX' : '0 10 0 0 0 0 0', 
                                                                     'ErrorsShiftY' : '0 10 0 0 0 0 0', 
                                                                     'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                     'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                     'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                     'ErrorsGamma'  : '0 0.01 0 0 0 0 0'})
  
  kalman_aligner_triplet_1 = Env.create_path('kalman_aligner_triplet_1')
  kalman_aligner_triplet_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_triplet_1.add_processor(name="M26CogHitMaker")
  kalman_aligner_triplet_1.add_processor(name="AlignTF_TC", params={'ExcludeDetector' : '3 4 5 6', 'MinimumHits' : '3' })    
  kalman_aligner_triplet_1.add_processor(name="TelAligner", params={'ErrorsShiftX' : '0 10 0 0 0 0 0', 
                                                                     'ErrorsShiftY' : '0 10 0 0 0 0 0', 
                                                                     'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                     'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                     'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                     'ErrorsGamma'  : '0 0.01 0 0 0 0 0'})
  
  correlator_triplet = Env.create_path('correlator_triplet')
  correlator_triplet.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  correlator_triplet.add_processor(name="M26CogHitMaker")
  correlator_triplet.add_processor(name="DEPCogHitMaker")
  correlator_triplet.add_processor(name="AlignTF_TC", params={'ExcludeDetector' : '3 4 5 6', 'MinimumHits' : '3' })    
  correlator_triplet.add_processor(name="TriplettCorrelator", params={'OutputRootFileName':'TripletCorrelator.root'})

  # Add fourth m26
  kalman_aligner_quadruplet_0 = Env.create_path('kalman_aligner_quadruplet_0')
  kalman_aligner_quadruplet_0.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_quadruplet_0.add_processor(name="M26CogHitMaker")
  kalman_aligner_quadruplet_0.add_processor(name="AlignTF_LC", params={'ExcludeDetector' : '3 5 6', 'MinimumHits' : '4' })     
  kalman_aligner_quadruplet_0.add_processor(name="PreAligner", params={'ErrorsShiftX' : '0 10 10 0 0 0 0', 
                                                                       'ErrorsShiftY' : '0 10 10 0 0 0 0', 
                                                                       'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                       'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                       'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                       'ErrorsGamma'  : '0 0.01 0.01 0 0 0 0'})
  
  kalman_aligner_quadruplet_1 = Env.create_path('kalman_aligner_quadruplet_1')
  kalman_aligner_quadruplet_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_quadruplet_1.add_processor(name="M26CogHitMaker")
  kalman_aligner_quadruplet_1.add_processor(name="AlignTF_TC", params={'ExcludeDetector' : '3 5 6', 'MinimumHits' : '4' })    
  kalman_aligner_quadruplet_1.add_processor(name="TelAligner", params={'ErrorsShiftX' : '0 10 10 0 0 0 0', 
                                                                       'ErrorsShiftY' : '0 10 10 0 0 0 0', 
                                                                       'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                       'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                       'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                       'ErrorsGamma'  : '0 0.01 0.01 0 0 0 0'})
  
  correlator_quadruplet = Env.create_path('correlator_quadruplet')
  correlator_quadruplet.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  correlator_quadruplet.add_processor(name="M26CogHitMaker")
  correlator_quadruplet.add_processor(name="DEPCogHitMaker")
  correlator_quadruplet.add_processor(name="AlignTF_TC", params={'ExcludeDetector' : '3 5 6', 'MinimumHits' : '4' })    
  correlator_quadruplet.add_processor(name="TriplettCorrelator", params={'OutputRootFileName':'QuadrupletCorrelator.root'})

  # Add fifth m26
  kalman_aligner_quintet_0 = Env.create_path('kalman_aligner_quintet_0')
  kalman_aligner_quintet_0.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_quintet_0.add_processor(name="M26CogHitMaker")
  kalman_aligner_quintet_0.add_processor(name="AlignTF_LC", params={'ExcludeDetector' : '3 6', 'MinimumHits' : '5' })     
  kalman_aligner_quintet_0.add_processor(name="PreAligner", params={'ErrorsShiftX' : '0 10 10 0 10 0 0', 
                                                                    'ErrorsShiftY' : '0 10 10 0 10 0 0', 
                                                                    'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                    'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                    'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                    'ErrorsGamma'  : '0 0.01 0.01 0 0.01 0 0'})
  
  kalman_aligner_quintet_1 = Env.create_path('kalman_aligner_quintet_1')
  kalman_aligner_quintet_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_quintet_1.add_processor(name="M26CogHitMaker")
  kalman_aligner_quintet_1.add_processor(name="AlignTF_TC", params={'ExcludeDetector' : '3 6', 'MinimumHits' : '5' })    
  kalman_aligner_quintet_1.add_processor(name="TelAligner", params={ 'ErrorsShiftX' : '0 10 10 0 10 0 0', 
                                                                     'ErrorsShiftY' : '0 10 10 0 10 0 0', 
                                                                     'ErrorsShiftZ' : '0 0 0 0 0 0 0', 
                                                                     'ErrorsAlpha'  : '0 0 0 0 0 0 0',
                                                                     'ErrorsBeta'   : '0 0 0 0 0 0 0', 
                                                                     'ErrorsGamma'  : '0 0.01 0.01 0 0.01 0 0'})
  
  correlator_quintet = Env.create_path('correlator_quintet')
  correlator_quintet.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  correlator_quintet.add_processor(name="M26CogHitMaker")
  correlator_quintet.add_processor(name="DEPCogHitMaker")
  correlator_quintet.add_processor(name="AlignTF_TC", params={'ExcludeDetector' : '3 6', 'MinimumHits' : '5' })    
  correlator_quintet.add_processor(name="TriplettCorrelator", params={'OutputRootFileName':'QuintetCorrelator.root'})
  
  kalman_aligner_1 = Env.create_path('kalman_aligner_1')
  kalman_aligner_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_1.add_processor(name="M26CogHitMaker")
  kalman_aligner_1.add_processor(name="DEPCogHitMaker")
  kalman_aligner_1.add_processor(name="AlignTF_LC")
  kalman_aligner_1.add_processor(name="PreAligner")
  
  kalman_aligner_2 = Env.create_path('kalman_aligner_2')
  kalman_aligner_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_2.add_processor(name="M26CogHitMaker")
  kalman_aligner_2.add_processor(name="DEPCogHitMaker")
  kalman_aligner_2.add_processor(name="AlignTF_TC")
  kalman_aligner_2.add_processor(name="TelAligner")
  
  telescope_dqm = Env.create_path('telescope_dqm')
  telescope_dqm.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm.add_processor(name="M26CogHitMaker")
  telescope_dqm.add_processor(name="DEPCogHitMaker")
  telescope_dqm.add_processor(name="AlignTF_TC")
  telescope_dqm.add_processor(name="TelescopeDQM")
  
  # create sequence of calibration paths 
  calpath= [ hotpixelkiller , 
             clusterizer, 
             correlator, 
             kalman_aligner_triplet_0,
             kalman_aligner_triplet_1,
             correlator_triplet,
             kalman_aligner_quadruplet_0,
             kalman_aligner_quadruplet_1,
             correlator_quadruplet,
             kalman_aligner_quintet_0,
             kalman_aligner_quintet_1,
             correlator_quintet,
             kalman_aligner_1, 
             kalman_aligner_2, 
             kalman_aligner_2, 
             kalman_aligner_2, 
             telescope_dqm, 
           ]
  
  return calpath

  
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
  
  
if __name__ == '__main__':
  
  params = ( rawfile, steerfiles, gearfile )
  
  # Create a simulated rawfile 
  simulate( params )
  
  # Calibrate the telescope and reconstruct the rawfile 
  calibrate( params )
  
    
    
    
