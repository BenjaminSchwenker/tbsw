"""
This is an example script to demonstrate how TBSW can be used to create an X0 image from a simulated
test beam experiment.

The script below simulates a test beam experiment where charged tracks from a monoenergetic beam 
cross a misaligned pixel telescope containing six Mimosa 26 detector planes. Two data sets are 
simulated. A first 'air' run is done with no additional scatterer between the telescope arms. The 
'air' run is used to calibrate the telescope. In a second 'aluminium' run, a aluminium plate with 
a well known thickness profile is inserted in between the telescope arms. This second run is used 
to compute a X0 image from the reconstructed scattering angles. The known X/X0 of the Al plate is 
used to calibrate the beam energy and the angular resolution of the telescope. This second step 
completes the calibration of the telescope for X0 imaging. 

Author: Ulf Stolzenberg <ulf.stolzenberg@phys.uni-goettingen.de>  
"""

import tbsw 
import os
import multiprocessing

# Determine maximum number of processes
nprocesses=2
count = min(nprocesses,multiprocessing.cpu_count())
pool = multiprocessing.Pool(processes=count)

# Path to steering files 
# Folder contains a gear file detailing the detector geometry and a config file
# for x0 calibration.
steerfiles = 'steering-files/x0-tb/'

# cal tags
# telescope calibration cal tag (typically named after telescope setup, beam energy etc.)
caltag='x0-sim'

# Tag for x0 calibration
x0caltag='alutarget'

# Name of gearfile
# This file describes the nominal geometry of a telescope 
# with two arms having three M26 planes each.  
gearfile = 'gear.xml'

# Use Single Hit seeding to speed up track finding?
Use_SingleHitSeeding=False

# Long telescopes may require a sensor by sensor alignment approach
Use_LongTelescopeCali=True

# Determine cluster resolution and store in cluster DB?
Use_clusterDB=True

# Flag to indicate that mc data is used (slcio format)
mcdata=True

# Script purpose option:
# 0: Script only processes imaging part
# 1: Script processes x0 calibration and imaging part
# Everything else: Script processes the whole chain: Simulation, Telescope calibration, angle reconstruction, x0 calibration and x0 imaging 
Script_purpose_option=2

# Number of iterations during target alignment
# Set to 0 or negative integer to disable target alignment
targetalignment_iterations=0

# File name for raw data 
rawfile_air = os.getcwd()+'/mc-air.slcio'
rawfile_alu_list = []
for nruns in range(0,nprocesses):
 rawfile_alu_list.append(os.getcwd()+'/mc-alu-run{:d}.slcio'.format(nruns+1))

# Set the name of this image
name_image1='alutarget-image'

# Number of events to simulate 
nevents_air = 1000000
nevents_TA = 1000000
nevents_alu = 6000000

#Parameters for simulation of misalignment
#Position parameters in mm and degree for rotations
mean_list=[0.0,0.0,0.0,0.0,0.0,0.0] 
sigma_list=[0.5,0.5,1.0,0.3,0.3,1.5]

# List of sensor ids and modes, which are excluded during misalignment
sensorexception_list=[5,0,11] 
modeexception_list=['']

# Nominal Beam energy
beamenergy=2.0

def simulate(): 
  """
  Simulates a rawfiles from a simulated test beam experiment
  Creates a folder in tmp-runs/ and populates it with 
  Marlin steering and logfiles.  
  """ 
  
  # Create tmpdir to hold all steerfiles and log files 
  SimObj = tbsw.Simulation(steerfiles=steerfiles, name='x0-sim' )
  
  # Get path to gearfile in simulation environment
  localgearfile = SimObj.get_filename(gearfile)
  
  # Misalign the telescope by applying random shifts to constants in local gearfile
  tbsw.gear.randomize_telescope(gearfile=localgearfile, mean_list=mean_list, sigma_list=sigma_list, sensorexception_list=sensorexception_list, modeexception_list=modeexception_list)
  
  # Create a gearfile for the run with an Al plate 
  # Only the Al plate is added in between the telescope arms
  # The misaligned positions of the telescope planes remain
  # the same. 
  gearfile_alu = "gear_alu.xml"
  localgearfile_alu = SimObj.copy_file(gearfile, gearfile_alu)
  
  # Replace the air layer  by a Al plate with 0.5mm thickness 
  tbsw.gear.set_parameter(gearfile=localgearfile_alu, sensorID=11, parametername='thickness', value=0.50)
  tbsw.gear.set_parameter(gearfile=localgearfile_alu, sensorID=11, parametername='radLength', value=89.0)
  tbsw.gear.set_parameter(gearfile=localgearfile_alu, sensorID=11, parametername='atomicNumber', value=13)
  tbsw.gear.set_parameter(gearfile=localgearfile_alu, sensorID=11, parametername='atomicMass', value=27)  
  
  # List to populate with simulation pathes (=jobs)
  simpaths = []
   
  # Create path to simulate air run (just air between telescope arms) 
  simpath_air = tbsw.path_utils.create_x0sim_path(SimObj, 'sim_air', rawfile_air, gearfile, nevents_air,  beamenergy)
  simpaths.append(simpath_air)
  
  # Create path to simulate alu run (Al plate between telescope arms)
  for nruns in range(0,nprocesses):
    simname='sim_alu_run{:d}'.format(nruns+1)
    simpath_alu = tbsw.path_utils.create_x0sim_path(SimObj, simname, rawfile_alu_list[nruns], gearfile_alu, nevents_alu/nprocesses, beamenergy)
    simpaths.append(simpath_alu)   
  
  # Run simulation to create rawfile with simulated digits 
  SimObj.simulate(paths=simpaths,caltag='x0-sim-truthdb')  
  

# Perform the telescope calibration
def calibrate(params):
  """
  Calibrates an misaligned tracking telescope from run data. 
  Creates a folder localDB/caltag in workspace containing 
  calibration results. 
  Creates a folder tmp-runs/name-sim/ and populates it with 
  Marlin steering and logfiles.  
  """ 
  
  rawfile, steerfiles, caltag = params
  
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
  CalObj = tbsw.Calibration(steerfiles=steerfiles, name=caltag + '-cal') 
  
  # Create list of calibration steps 
  calpaths = tbsw.path_utils.create_x0analysis_calibration_paths(CalObj, rawfile, gearfile, nevents_air, Use_clusterDB, beamenergy, mcdata, Use_LongTelescopeCali)
  
  # Run the calibration steps 
  CalObj.calibrate(paths=calpaths,ifile=rawfile,caltag=caltag)

  # Create DQM plots 
  tbsw.DQMplots.calibration_DQMPlots(name=caltag + '-cal')

# Perform the angle reconstruction of a single run
def reconstruct(params):

  rawfile, steerfiles, caltag = params

  # Set cal tag that includes run name
  name = os.path.splitext(os.path.basename(rawfile))[0] + '-' + caltag
  
  # Reconsruct the rawfile using the caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = tbsw.Reconstruction(steerfiles=steerfiles, name=name )
  
  # Create reconstuction path
  recopath = tbsw.path_utils.create_anglereco_path(RecObj, rawfile, gearfile, nevents_alu, Use_SingleHitSeeding, Use_clusterDB, beamenergy, mcdata)  
  
  # Run the reconstuction  
  RecObj.reconstruct(paths=recopath,ifile=rawfile,caltag=caltag) 


# Create angle reconstruction DQM plots
def reconstruction_DQM(params):

  rawfile, caltag = params

  # Set cal tag that includes run name
  name = os.path.splitext(os.path.basename(rawfile))[0] + '-' + caltag
  
  # Create DQM plots
  tbsw.DQMplots.anglereco_DQMPlots(filepath='root-files/X0-{}.root'.format(name)) 
  

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
   
  rawfile, steerfiles, caltag, iteration = params
  
  if iteration == 0:
    oldcaltag=caltag
  else: 
    oldcaltag=caltag+'-target-align-it-{:d}'.format(iteration-1)
     
  newcaltag=caltag+'-target-align-it-{:d}'.format(iteration)
  
  # Reconsruct the rawfile using the caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = tbsw.Reconstruction(steerfiles=steerfiles, name='x0-reco-targetalign-it-{:d}'.format(iteration) )
  
  # Create reconstuction path
  recopath = tbsw.path_utils.create_anglereco_path(RecObj, rawfile, gearfile, nevents_TA, Use_SingleHitSeeding, Use_clusterDB, beamenergy, mcdata)  
  
  # Run the reconstuction  
  RecObj.reconstruct(path=recopath,ifile=rawfile,caltag=oldcaltag)  
  
  # Read the vertex position and save it in the alignmentDB
  dbname=RecObj.create_dbfilename("alignmentDB.root")
  treename=RecObj.get_rootfilename('X0')
  tbsw.gear.save_targetpos(treename,dbname)
  RecObj.export_caltag(newcaltag)


# Perform x0 calibration
def xx0calibration(params):

  x0caltag, RunList, steerfiles, caltag = params

  # Create list with input root files from list of input raw files
  RootFileList_x0cali=[]
  tbsw.x0imaging.X0Calibration.CreateRootFileList(rawlist=RunList,rootlist=RootFileList_x0cali, caltag=caltag)

  # Generate a uncalibrated X/X0 image
  imagenametag='X0image-calitarget-Uncalibrated'
  tbsw.x0imaging.GenerateImage.x0imaging(rootfilelist=[RootFileList_x0cali[0]],caltag='',steerfiles=steerfiles,nametag=imagenametag)

  # Path to uncalibrated X0 image file
  imagefilename='/root-files/'+imagenametag+'.root'

  # Do a calibration of the angle resolution
  tbsw.x0imaging.X0Calibration.x0calibration(rootfilelist=RootFileList_x0cali,imagefile=imagefilename,caltag=x0caltag,steerfiles=steerfiles)

  nametag='X0Calibration-'+x0caltag
  tbsw.DQMplots.x0calibration_DQMPlots(nametag=nametag)


# Generate x0 image
def xx0image(params):

  x0caltag, RunList, steerfiles, caltag, listnametag = params

  RootFileList_x0image=[]
  tbsw.x0imaging.X0Calibration.CreateRootFileList(rawlist=RunList,rootlist=RootFileList_x0image, caltag=caltag)

  if listnametag=='':
    print("No image name found. Using default naming scheme!")
    listnametag='X0image'

  # Determine name of image
  if caltag=='':
    nametag=listnametag+'-Uncalibrated'
  else:
    nametag=listnametag+'-Calibrated-'+x0caltag

  # Do a calibration of the angle resolution
  tbsw.x0imaging.GenerateImage.x0imaging(rootfilelist=RootFileList_x0image,caltag=x0caltag,steerfiles=steerfiles,nametag=nametag)

  tbsw.DQMplots.x0image_Plots(nametag=nametag)
  
  
if __name__ == '__main__':
   
  # Get current directory
  cwdir = os.getcwd()
  
  # Create a simulated lcio for run with no target (air run) and 
  # multiple run s with a Al plate as scattering material
  simulate( )
  
  # Calibrate the telescope 
  # In case you already have all the DB files from another telescope calibration 
  # and want to reuse it, just switch to Script_purpose_option 0 or 1
  #
  # DQM plots like track p/chi2 values, residuals and other interesting parameters
  # from this telescope calibration step can be found as pdf files in 
  # workspace/results/telescopeDQM
  if Script_purpose_option !=0 and Script_purpose_option !=1:
    params_cali = ( rawfile_air, steerfiles, caltag)
    calibrate( params_cali )
   
  for it in range(0,targetalignment_iterations):
    params_TA = (rawfile_alu, steerfiles, caltag, it)
    print "The parameters for the target alignment are: " 
    print params_TA
    targetalignment(params_TA)
  
  # Angle reconstruction
  # In case you already have reconstructed the scattering angles for all
  # the runs you are interested in, just switch to Script_purpose_option 0 or 1
  #
  # The root files with the reconstructed angles and other parameters (see 
  # README_X0.md for a full list and some descriptions) can be found in 
  # workspace/root-files/X0-run*runnumber, etc*-reco.root
  # The histmap and angle resolution for every single run can be found in 
  # workspace/results/anglerecoDQM/
  if Script_purpose_option !=0 and Script_purpose_option !=1:
    params_reco=[(x, steerfiles, caltag) for x in rawfile_alu_list]
    print "The parameters for the reconstruction are: " 
    print params_reco

    count = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=count)
    pool.map(reconstruct, params_reco)

    for rawfile in rawfile_alu_list:
      params=(rawfile, caltag)
      reconstruction_DQM(params)

  # Start x0 calibration
  # In case you already have the x0 calibration DB file from a previous x0 calibration 
  # and want to reuse it, just switch to Script_purpose_option 0
  #
  # The fitted distributions and self-consistency plots in pdf format from this 
  # x0 calibration can be found in the workspace/tmp-runs/*X0Calibration/ directory
  if Script_purpose_option !=0:
    params_x0cali = ( x0caltag, rawfile_alu_list, steerfiles, caltag)
    xx0calibration(params_x0cali)

  # Generate a calibrated X/X0 image
  #
  # The calibrated radiation length image and other images, such as the beamspot
  # etc can be found in the workspace/root-files/*CalibratedX0Image.root
  params_x0image = ( x0caltag, rawfile_alu_list, steerfiles, caltag, name_image1)
  xx0image(params_x0image)

