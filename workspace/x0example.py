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

# Path to steering files 
# Folder contains a gear file detailing the detector geometry and a config file
# for x0 calibration.
steerfiles = 'steering-files/x0-tb/'

# Name of gearfile
# This file describes the nominal geometry of a telescope 
# with two arms having three M26 planes each.  
gearfile = 'gear.xml'

# Tag for x0 calibration
x0caltag='alutarget'

# File name for raw data 
rawfile_air = os.getcwd()+'/mc-air.slcio'
rawfile_alu = os.getcwd()+'/mc-alu.slcio'

# Tag for calibration data
caltag = 'airtarget'

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

# Number of iterations during target alignment
# Set to 0 or negative integer to disable target alignment
targetalignment_iterations=0

# Nominal Beam energy
beamenergy=2.0

# Use Single Hit seeding to speed up track finding?
Use_SingleHitSeeding=False

# Determine cluster resolution and store in cluster DB?
Use_clusterDB=True

# Flag to indicate that mc data is used (slcio format)
mcdata=True


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
  simpath_alu = tbsw.path_utils.create_x0sim_path(SimObj, 'sim_alu', rawfile_alu, gearfile_alu, nevents_alu, beamenergy)
  simpaths.append(simpath_alu)   
  
  # Run simulation to create rawfile with simulated digits 
  SimObj.simulate(paths=simpaths,caltag='x0-sim-truthdb')  
  

def calibrate():
  """
  Calibrates an misaligned tracking telescope from air run. 
  Creates a folder localDB/caltag in workspace containing 
  calibration results. 
  Creates a folder tmp-runs/ and populates it with 
  Marlin steering and logfiles.  
  """ 
  
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
  CalObj = tbsw.Calibration(steerfiles=steerfiles, name='x0-cal') 
  
  # Create list of calibration steps 
  calpaths = tbsw.path_utils.create_x0analysis_calibration_paths(CalObj, rawfile_air, gearfile, nevents_air, Use_clusterDB, beamenergy, mcdata)
  
  # Run the calibration steps 
  CalObj.calibrate(paths=calpaths,ifile=rawfile_air,caltag=caltag)  
  
  # Create DQM plots 
  tbsw.DQMplots.calibration_DQMPlots(name='x0-cal')
  
def reconstruct():
  
  # Reconsruct the rawfile using the caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = tbsw.Reconstruction(steerfiles=steerfiles, name='x0-reco')
  
  # Create reconstuction path
  recopath = tbsw.path_utils.create_anglereco_path(RecObj, rawfile_alu, gearfile, nevents_alu, Use_SingleHitSeeding, Use_clusterDB, beamenergy, mcdata)  
  
  # Run the reconstuction  
  RecObj.reconstruct(paths=recopath,ifile=rawfile_alu,caltag=caltag)  
  
  # Create DQM plots
  tbsw.DQMplots.anglereco_DQMPlots(filepath='root-files/X0-x0-reco.root') 
  

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
  
  
if __name__ == '__main__':
   
  # Get current directory
  cwdir = os.getcwd()
  
  # Create a simulated rawfiles for run with no target (air run) and Al plate 
  #simulate( )
  
  # Calibrate the telescope using air run  
  #calibrate( )
   
  for it in range(0,targetalignment_iterations):
    params_TA = (rawfile_alu, steerfiles, caltag, it)
    print "The parameters for the target alignment are: " 
    print params_TA
    targetalignment(params_TA)
  
  # Reconstruct the alu rawfile 
  #reconstruct( )
   
  # Create root file list as input for x0 analysis 
  rootlist=['root-files/X0-x0-reco.root']
  
  # Generate a uncalibrated X/X0 image
  nametag='X0image-Uncalibrated'
  tbsw.x0imaging.X0Calibration.x0imaging(filelist=rootlist,caltag='',steerfiles=steerfiles,nametag=nametag)
  tbsw.DQMplots.x0image_Plots(nametag=nametag)

  # Path to uncalibrated X0 image file
  imagefilename='/root-files/'+nametag+'.root'

  # Do a calibration of the angle resolution
  nametag='X0Calibration-'+x0caltag
  tbsw.x0imaging.X0Calibration.x0calibration(filelist=rootlist,imagefilename=imagefilename,caltag=x0caltag,steerfiles=steerfiles)
  tbsw.DQMplots.x0calibration_DQMPlots(nametag=nametag)

  # Generate a calibrated X/X0 image
  nametag='X0image-Calibrated-'+x0caltag
  tbsw.x0imaging.X0Calibration.x0imaging(filelist=rootlist,caltag=x0caltag,steerfiles=steerfiles,nametag=nametag)
  tbsw.DQMplots.x0image_Plots(nametag=nametag)

