"""
This is an example script to demonstrate how TBSW can be used to create an X0 image from a simulated
test beam experiment.

The script below simulates a test beam experiment where charged tracks from a monoenergetic beam 
cross a a misaligned pixel telescope containing six Mimosa 26 detector planes. Two data sets are 
simulated. A first 'air' run is done with no additional scatterer between the telescope arms. The 
'air' run is used to calibrate the telescope. In a second 'aluminium' run, a aluminium plate with 
a well known thickness profile is inserted in between the telescope arms. This second run is used 
to compute a X0 image from the reconstructed scattering angles. The known comparison between the 
reconstructed X0 image and the a priori known image is used to calibrate the beam energy and the 
angular resolution of the telescope. This second step completes the calibration of the telescope
for X0 imaging. 

Author: Ulf Stolzenberg <ulf.stolzenberg@phys.uni-goettingen.de>  
"""

from tbsw import *

# Path to steering files 
# Steeringfiles are xml files and define details of the simulation like how many events are produced
# or how M26 sensors are digitized. XML parameters can be adjusted using any test editor
steerfiles = 'steering-files/x0-tb/'
#cal tag
air_caltag='mc-air-test'
# Gearfile for runs 
gearfile = 'gear.xml'
# File name for raw data 
rawfile_air = os.getcwd()+'/mc-air.slcio'
rawfile_alu = os.getcwd()+'/mc-alu.slcio'
# Number of events to simulate 
nevents_air = 700000
nevents_alu = 5000000

#Parameters for simulation of misalignment
#Position parameters in mm
mean_pos=0.0 
sigma_pos=.1
#Rotation parameters in degrees
mean_rot=0.0 
sigma_rot=0.1

def create_sim_path_air(Env):
  """
  Returns a list of tbsw path objects to simulate a test beam run without any DUT
  """
  
  sim_air = Env.create_path('sim')
  sim_air.set_globals(params={'GearXMLFile': 'gear_air.xml' , 'MaxRecordNumber' : nevents_air})   
  sim_air.add_processor(name="InfoSetter")
  sim_air.add_processor(name="ParticleGun")
  sim_air.add_processor(name="FastSim")
  sim_air.add_processor(name="M26Digitizer")
  sim_air.add_processor(name="LCIOOutput",params={"LCIOOutputFile" : rawfile_air })

  sim_alu = Env.create_path('sim2')
  sim_alu.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_alu})   
  sim_alu.add_processor(name="InfoSetter")
  sim_alu.add_processor(name="ParticleGun")
  sim_alu.add_processor(name="FastSim")
  sim_alu.add_processor(name="M26Digitizer")
  sim_alu.add_processor(name="LCIOOutput",params={"LCIOOutputFile" : rawfile_alu })
    
  simpath = [ sim_air,
              sim_alu, 
            ]

  return simpath


def create_calibration_path(Env):
  """
  Returns a list of tbsw path objects to calibrate the tracking telescope
  """
  
  hotpixelkiller = Env.create_path('hotpixelkiller')
  hotpixelkiller.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': rawfile_air })  
  hotpixelkiller.add_processor(name="M26HotPixelKiller")

  cluster_calibrator_mc = Env.create_path('cluster_calibrator_mc')
  cluster_calibrator_mc.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': rawfile_air })  
  cluster_calibrator_mc.add_processor(name="M26Clusterizer")
  cluster_calibrator_mc.add_processor(name="M26CogHitMaker")
  cluster_calibrator_mc.add_processor(name="M26ClusterCalibrationFromMC")

  coghitmaker = Env.create_path('coghitmaker')
  coghitmaker.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_air, 'LCIOInputFiles': rawfile_air }) 
  coghitmaker.add_processor(name="M26Clusterizer")
  coghitmaker.add_processor(name="M26CogHitMaker")
  coghitmaker.add_processor(name="RawDQM")
  coghitmaker.add_processor(name="LCIOOutput")
  
  correlator = Env.create_path('correlator')
  correlator.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_air, 'LCIOInputFiles': "tmp.slcio" }) 
  correlator.add_processor(name="TelCorrelator")
  
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
  cluster_calibration_1.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_air, 'LCIOInputFiles': "tmp.slcio" })  
  cluster_calibration_1.add_processor(name="AlignTF_TC")
  cluster_calibration_1.add_processor(name="M26ClusterCalibrator")

  kalman_aligner_3 = Env.create_path('kalman_aligner_3')
  kalman_aligner_3.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_3.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : 'goehit_m26' })
  kalman_aligner_3.add_processor(name="AlignTF_TC")
  kalman_aligner_3.add_processor(name="TelAligner")
  
  cluster_calibration_2 = Env.create_path('cluster_calibration_2')
  cluster_calibration_2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_air, 'LCIOInputFiles': "tmp.slcio" })  
  cluster_calibration_2.add_processor(name="M26GoeHitMaker", params={'HitCollectionName' : "goehit_m26" }) 
  cluster_calibration_2.add_processor(name="AlignTF_TC", params={'InputHitCollectionNameVec': "goehit_m26"})
  cluster_calibration_2.add_processor(name="M26ClusterCalibrator")

  goehitmaker = Env.create_path('goehitmaker')
  goehitmaker.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_air, 'LCIOInputFiles': rawfile_air }) 
  goehitmaker.add_processor(name="M26Clusterizer")
  goehitmaker.add_processor(name="M26GoeHitMaker")
  goehitmaker.add_processor(name="RawDQM", params={'RootFileName': 'RawDQM2.root'})
  goehitmaker.add_processor(name="LCIOOutput", params={'LCIOWriteMode': "WRITE_NEW"})

  correlator2 = Env.create_path('correlator2')
  correlator2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_air, 'LCIOInputFiles': rawfile_air }) 
  correlator2.add_processor(name="TelCorrelator", params={'OutputRootFileName': 'XCorrelator2.root')
  
  kalman_aligner_4 = Env.create_path('kalman_aligner_4')
  kalman_aligner_4.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_4.add_processor(name="AlignTF_LC")
  kalman_aligner_4.add_processor(name="PreAligner", params={'RootFileName': 'KalmanAlign-iteration-3.root'})
  
  kalman_aligner_5 = Env.create_path('kalman_aligner_5')
  kalman_aligner_5.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : 100000, 'LCIOInputFiles': "tmp.slcio" })  
  kalman_aligner_5.add_processor(name="AlignTF_TC")
  kalman_aligner_5.add_processor(name="TelAligner", params={'RootFileName': 'KalmanAlign-final2.root'}) 
  
  telescope_dqm2 = Env.create_path('telescope_dqm2')
  telescope_dqm2.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_air, 'LCIOInputFiles': "tmp.slcio" })  
  telescope_dqm2.add_processor(name="AlignTF_TC")
  telescope_dqm2.add_processor(name="TelescopeDQM", params={'RootFileName' : 'TelescopeDQM2.root'})
  
  
  # create sequence of calibration paths 
  calpath= [ hotpixelkiller , 
             cluster_calibrator_mc,
             coghitmaker,
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
             cluster_calibration_2, 
             cluster_calibration_2,
			 goehitmaker, 
             correlator2, 
             kalman_aligner_4, 
             kalman_aligner_5, 
             kalman_aligner_5, 
             kalman_aligner_5, 
             telescope_dqm2,   
           ]
  
  return calpath


def create_reco_path(Env):
  """
  Returns a list of tbsw path objects to reconstruct a test beam run 
  """
  
  reco = Env.create_path('reco')
  reco.set_globals(params={'GearXMLFile': gearfile , 'MaxRecordNumber' : nevents_alu, 'LCIOInputFiles': rawfile_alu }) 
  reco.add_processor(name="M26Clusterizer")
  reco.add_processor(name="M26GoeHitMaker")
  reco.add_processor(name="DownstreamFinder")
  reco.add_processor(name="UpstreamFinder")
  reco.add_processor(name="X0Imager")
    
  return [ reco ]

def simulate(): 
  """
  Simulates a rawfile from a simulated test beam experiment
  Creates a folder tmp-runs/name-sim/ and populates it with 
  Marlin steering and logfiles.  
  """ 
  
  # Create tmpdir to hold all steerfiles and log files 
  SimObj = Simulation(steerfiles=steerfiles, name=os.path.splitext(os.path.basename(rawfile_alu))[0] + '-sim' )

  # Create steerfiles for processing
  simpath = create_sim_path_air(SimObj)

  # Get gearfile
  localgearfile = SimObj.get_filename('gear.xml')

  # Misalign gear file
  randomize_telescope(gearfile=localgearfile, mean_pos=mean_pos, sigma_pos=sigma_pos, mean_rot=mean_rot, sigma_rot=sigma_rot)

  # Copy gearfile
  gearfile_air=SimObj.tmpdir+'/'+'gear_air.xml'
  shutil.copy2(localgearfile,gearfile_air)

  # Change DUT in copied gearfile
  set_parameter(gearfile=gearfile_air, sensorID=11, parametername='thickness', value=0.0001)
  set_parameter(gearfile=gearfile_air, sensorID=11, parametername='radLength', value=304.0)

  # Run simulation to create rawfile with simulated digits 
  SimObj.simulate(path=simpath)  


def calibrate():
  """
  Calibrates an misaligned tracking telescope from run data. 
  Creates a folder localDB/caltag in workspace containing 
  calibration results. 
  Creates a folder tmp-runs/name-sim/ and populates it with 
  Marlin steering and logfiles.  
  """ 
  
  # Tag for calibration data
  localcaltag = os.path.splitext(os.path.basename(rawfile_air))[0] + '-test'
  
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
  CalObj = Calibration(steerfiles=steerfiles, name=localcaltag + '-cal') 

  # Get gearfile and set air as DUT material
  localgearfile = CalObj.get_filename('gear.xml')
  set_parameter(gearfile=localgearfile, sensorID=11, parametername='thickness', value=0.0001)
  set_parameter(gearfile=localgearfile, sensorID=11, parametername='radLength', value=304.0)
  
  # Create list of calibration steps 
  calpath = create_calibration_path(CalObj)
  
  # Run the calibration steps 
  CalObj.calibrate(path=calpath,ifile=rawfile_air,caltag=localcaltag)  

def reconstruct():

  # Reconsruct the rawfile using the caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=air_caltag + '-reco' )

  # Create reconstuction path
  recopath = create_reco_path(RecObj)  

  # Run the reconstuction  
  RecObj.reconstruct(path=recopath,ifile=rawfile_alu,caltag=air_caltag) 
  


# Function which starts the imaging script
def x0imaging(filename,caltag,deletetag):

  flags='./tbsw/x0imaging/GenerateImage.py -i '+filename+' -f '+imagecfgfilename+' -c '+caltag+' -d '+`deletetag`
  print('Starting X0 imaging')
  print(flags)
  subprocess.call(flags, shell=True)


# Function which starts the x0 calibration script
def x0calibration(filename,imagefilename,caltag):

  flags='./tbsw/x0imaging/X0Calibration.py -i '+filename+' -f '+calibrationcfgfilename+' -m '+imagefilename+' -c '+caltag
  print('Starting X0 calibration')
  print(flags)
  subprocess.call(flags, shell=True)

  
if __name__ == '__main__':


  # Get current directory
  cwdir = os.getcwd()

  # Create a simulated rawfile with air data 
  simulate( )

  # Calibrate the telescope 
  calibrate( )

  # Reconstruct the alu rawfile 
  reconstruct( )

  imagefilename='root-files/X0-mc-air-test-reco-Uncalibrated-X0image.root'
  imagecfgfilename='steering-files/x0-tb/image.cfg'
  calibrationcfgfilename='steering-files/x0-tb/x0calibration.cfg'
  deletetag=1

  # Base filename of the X0 root file
  basefilename='X0-mc-air-test-reco'

  # Total path of X0 root file
  filename='root-files/'+basefilename+'.root'

  # Total path if the different kinds of X0 image files
  imagefile=cwdir+'/root-files/'+basefilename+'-X0image.root'
  uncalibratedimagefile=cwdir+'/root-files/'+basefilename+'-Uncalibrated-X0image.root'
  calibratedimagefile=cwdir+'/root-files/'+basefilename+'-Calibrated-X0image.root'

  # Total path to the different temporary work directories,
  # in which the images are generated
  tmpdir = cwdir+'/tmp-runs/'+basefilename+'-X0image'
  uncaltmpdir = cwdir+'/tmp-runs/'+basefilename+'-UncalibratedX0image'
  caltmpdir = cwdir+'/tmp-runs/'+basefilename+'-CalibratedX0image'

  # Generate a uncalibrated X/X0 image
  x0imaging(filename,'',deletetag)

  # Rename the image file
  if os.path.isfile(uncalibratedimagefile):
     os.remove(uncalibratedimagefile) 
  if os.path.isfile(imagefile):
     shutil.move(imagefile,uncalibratedimagefile)

  # Rename the directory in tmp-runs, in which the image was generated
  if os.path.isdir(uncaltmpdir):
     shutil.rmtree(uncaltmpdir) 
  if os.path.isdir(tmpdir):
     shutil.copytree(tmpdir, uncaltmpdir)

  # Do a calibration of the angle resolution
  x0calibration(filename,imagefilename,air_caltag)

  # Generate a calibrated X/X0 image
  x0imaging(filename,air_caltag,deletetag)

  # Rename the image file
  if os.path.isfile(calibratedimagefile):
     os.remove(calibratedimagefile) 
  if os.path.isfile(imagefile):
     shutil.move(imagefile,calibratedimagefile)

  # Rename the directory in tmp-runs, in which the image was generated
  if os.path.isdir(caltmpdir):
     shutil.rmtree(caltmpdir) 
  if os.path.isdir(tmpdir):
     shutil.copytree(tmpdir, caltmpdir)

