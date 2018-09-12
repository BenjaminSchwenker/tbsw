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

import tbsw.x0imaging.X0Calibration

# Path to steering files 
# Steeringfiles are xml files and define details of the simulation like how many events are produced
# or how M26 sensors are digitized. XML parameters can be adjusted using any test editor
steerfiles = 'steering-files/x0-tb/'

# Gearfile for runs 
gearfile = 'gear.xml'

# Tag for x0 calibration
x0tag='alutarget'

# Name of the truth db file
truthdb_filename='alignmentDB_truth.root'
alignmentdb_filename='alignmentDB.root'

# File name for raw data 
rawfile_air = os.getcwd()+'/mc-air.slcio'
rawfile_alu = os.getcwd()+'/mc-alu.slcio'

# Tag for calibration data
caltag = os.path.splitext(os.path.basename(rawfile_air))[0] + '-test'

# Number of events to simulate 
nevents_air = 1000000
nevents_TA = 1000000
nevents_alu = 6000000

#Parameters for simulation of misalignment
#Position parameters in mm
mean_list=[0.0,0.0,0.0,0.0,0.0,0.0] 
sigma_list=[0.1,0.1,1.0,0.1,0.1,0.1]

# List of sensor ids and modes, which are excluded during misalignment
sensorexception_list=[5,0] 
modeexception_list=['']

# Number of iterations during target alignment
# Set to 0 or negative integer to disable target alignment
targetalignment_iterations=0

# Nominal Beam energy
beamenergy=2.0

#Delete partial images after generating image
deletetag='1'

# Use Single Hit seeding to speed up track finding?
Use_SingleHitSeeding=False

def simulate(): 
  """
  Simulates a rawfile from a simulated test beam experiment
  Creates a folder tmp-runs/name-sim/ and populates it with 
  Marlin steering and logfiles.  
  """ 
  
  # Create tmpdir to hold all steerfiles and log files 
  SimObj = Simulation(steerfiles=steerfiles, name=os.path.splitext(os.path.basename(rawfile_alu))[0] + '-sim' )

  # Set Beam energy
  SimObj.set_beam_momentum(beamenergy)

  # Create steerfiles for processing
  simpath = create_x0sim_path_air(SimObj, rawfile_air, rawfile_alu, gearfile, nevents_air, nevents_alu)

  # Get gearfile
  localgearfile = SimObj.get_filename('gear.xml')

  # Misalign gear file
  randomize_telescope(gearfile=localgearfile, mean_list=mean_list, sigma_list=sigma_list, sensorexception_list=sensorexception_list, modeexception_list=modeexception_list)

  localtruthdb_filename=SimObj.create_dbfilename(truthdb_filename)

  # Convert gear file to alignmentDB root file, which will be stored in the sim folder
  Create_AlignmentDBFile_From_Gear(gearfile=SimObj.get_filename('gear.xml'), truthdbfilename=localtruthdb_filename)

  # Copy gearfile
  SimObj.copy_file('gear.xml','gear_air.xml')

  # Get air gearfile
  gearfile_air = SimObj.get_filename('gear_air.xml')

  # Change DUT in copied gearfile
  set_parameter(gearfile=gearfile_air, sensorID=11, parametername='thickness', value=0.0001)
  set_parameter(gearfile=gearfile_air, sensorID=11, parametername='radLength', value=304000.0)


  # Create caltag for the truthdb
  caltag = os.path.splitext(os.path.basename(rawfile_air))[0] + '-test'
  simcaltag=caltag+ '-truthdb'

  # Run simulation to create rawfile with simulated digits 
  SimObj.simulate(path=simpath,caltag=simcaltag)  


def calibrate():
  """
  Calibrates an misaligned tracking telescope from run data. 
  Creates a folder localDB/caltag in workspace containing 
  calibration results. 
  Creates a folder tmp-runs/name-sim/ and populates it with 
  Marlin steering and logfiles.  
  """ 
  
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
  CalObj = Calibration(steerfiles=steerfiles, name=caltag + '-cal') 

  # Set Beam energy
  CalObj.set_beam_momentum(beamenergy)

  # Get gearfile and set air as DUT material
  localgearfile = CalObj.get_filename('gear.xml')
  set_parameter(gearfile=localgearfile, sensorID=11, parametername='thickness', value=0.0001)
  set_parameter(gearfile=localgearfile, sensorID=11, parametername='radLength', value=304000.0)
  
  # Create list of calibration steps 
  calpath = create_mc_x0calibration_path(CalObj, rawfile_air, gearfile, nevents_air)

  # Run the calibration steps 
  CalObj.calibrate(path=calpath,ifile=rawfile_air,caltag=caltag)  

def reconstruct():

  # Reconsruct the rawfile using the caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=caltag + '-reco' )

  # Set Beam energy
  RecObj.set_beam_momentum(beamenergy)

  # Create reconstuction path
  recopath = create_mc_x0reco_path(RecObj, rawfile_alu, gearfile, nevents_alu, Use_SingleHitSeeding)  

  # Use caltag of the last target alignment iteration
  iteration_string='-target-alignment-it'+str(targetalignment_iterations-1)
  localcaltag=caltag+iteration_string 

  if targetalignment_iterations < 1:
    localcaltag=caltag

  # Run the reconstuction  
  RecObj.reconstruct(path=recopath,ifile=rawfile_alu,caltag=localcaltag)   


def targetalignment(params):
  """
  Starts the scattering angle reconstruction and vertex fit on the central target
  plane. Afterwards the mean vertex z position is set as the new target z position in
  the aligment DB file and the calibration tag is exported
    :@params:       consists of rawfile, steerfiles, gearfile, caltag, iteration
    :@rawfile:      Input file for the reconstruction 
    :@BE            Nominal beam energy of the run
    :@nevents       Number of events
    :@steerfiles:   Directory with the steering files for the reconstruction
    :@gearfile:     Name of the gear file
    :@caltag:       calibration tag for the reconstruction
    :@iteration:    Target alignment iteration counter  
    :author: benjamin.schwenker@phys.uni-goettinge.de  
  """ 

  rawfile, steerfiles, calibrationtag, iteration = params

  if rawfile == None:
    return None

  if iteration == None:
    return None

  prev_iteration_string='-target-alignment-it'+str(iteration-1)
  curr_iteration_string='-target-alignment-it'+str(iteration)

  localcaltag=caltag+prev_iteration_string

  if iteration == 0:
    localcaltag=calibrationtag

  newcaltag=calibrationtag+curr_iteration_string

  # Reconsruct the rawfile using the caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=newcaltag )

  # Set Beam energy
  RecObj.set_beam_momentum(beamenergy)
  
  # Create reconstuction path
  recopath = create_mc_x0reco_path(RecObj, rawfile_alu, gearfile, nevents_TA, Use_SingleHitSeeding)  

  # Run the reconstuction  
  RecObj.reconstruct(path=recopath,ifile=rawfile,caltag=localcaltag)  

  # Read the vertex position and save it in the alignmentDB
  dbname=RecObj.create_dbfilename("alignmentDB.root")
  treename=RecObj.get_rootfilename('X0')
  save_targetpos(treename,dbname)
  RecObj.export_caltag(newcaltag)

  
if __name__ == '__main__':


  # Get current directory
  cwdir = os.getcwd()

  # Create a simulated rawfile with air data 
  simulate( )

  # Calibrate the telescope 
  calibrate( )


  for it in range(0,targetalignment_iterations):
    params_TA = (rawfile_alu, steerfiles, caltag, it)
    print "The parameters for the target alignment are: " 
    print params_TA
    targetalignment(params_TA)

  # Reconstruct the alu rawfile 
  reconstruct( )

  # Base filename of the X0 root file
  basefilename='X0-mc-air-test-reco'

  # Total path of X0 root file
  filename='root-files/'+basefilename+'.root'

  # Generate a uncalibrated X/X0 image
  tbsw.x0imaging.X0Calibration.x0imaging(filename=filename,caltag='',deletetag=deletetag,steerfiles=steerfiles,nametag='Uncalibrated')

  # Path to uncalibrated X0 image file
  imagefilename='/root-files/'+basefilename+'-UncalibratedX0image.root'

  # Do a calibration of the angle resolution
  tbsw.x0imaging.X0Calibration.x0calibration(filename=filename,imagefilename=imagefilename,caltag=x0tag,steerfiles=steerfiles)

  # Generate a calibrated X/X0 image
  tbsw.x0imaging.X0Calibration.x0imaging(filename=filename,caltag=x0tag,deletetag=deletetag,steerfiles=steerfiles,nametag='Calibrated')

