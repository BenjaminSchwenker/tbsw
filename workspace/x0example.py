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
# 1: Script only processes x0 calibration part
# 2: Script processes x0 calibration and imaging part
# 3 and larger: Process x0 calibration and imaging part and angle reconstruction
# 4 and larger: Process x0 calibration and imaging part, angle reconstruction, telescope calibration and target alignment
# 5 and larger: Process everything
Script_purpose_option=5

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


  
if __name__ == '__main__':
   
  # Get current directory
  cwdir = os.getcwd()
  
  # Create a simulated lcio for run with no target (air run) and 
  # multiple run s with a Al plate as scattering material
  if Script_purpose_option > 4:
    tbsw.x0script_functions.simulate(rawfile_air, rawfile_alu_list, steerfiles, gearfile, nevents_air, nevents_alu, nprocesses, mean_list, sigma_list, sensorexception_list, modeexception_list, beamenergy)
  
  
  # Calibrate the telescope 
  # In case you already have all the DB files from another telescope calibration 
  # and want to reuse it, just switch to Script_purpose_option 0 or 1
  #
  # DQM plots like track p/chi2 values, residuals and other interesting parameters
  # from this telescope calibration step can be found as pdf files in 
  # workspace/results/telescopeDQM
  if Script_purpose_option > 3:
    tbsw.x0script_functions.calibrate( rawfile_air, steerfiles, caltag, gearfile, nevents_air, Use_clusterDB, beamenergy, mcdata, Use_LongTelescopeCali)
   
    # Target alignment
    for it in range(0,targetalignment_iterations):
      tbsw.x0script_functions.targetalignment(rawfile_alu_list[0], steerfiles, it, caltag, gearfile, nevents_TA, Use_SingleHitSeeding, Use_clusterDB, beamenergy, mcdata)
  
  # Angle reconstruction
  # In case you already have reconstructed the scattering angles for all
  # the runs you are interested in, just switch to Script_purpose_option 0 or 1
  #
  # The root files with the reconstructed angles and other parameters (see 
  # README_X0.md for a full list and some descriptions) can be found in 
  # workspace/root-files/X0-run*runnumber, etc*-reco.root
  # The histmap and angle resolution for every single run can be found in 
  # workspace/results/anglerecoDQM/
  if Script_purpose_option > 2:
    params_reco=[(x, steerfiles, caltag, gearfile, nevents_alu, Use_SingleHitSeeding, Use_clusterDB, beamenergy, mcdata) for x in rawfile_alu_list]
    print "The parameters for the reconstruction are: " 
    print params_reco

    count = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=count)
    pool.map(tbsw.x0script_functions.reconstruct, params_reco)

    for rawfile in rawfile_alu_list:
      tbsw.x0script_functions.reconstruction_DQM(rawfile, caltag)

  if Script_purpose_option > 1:
    # Start x0 calibration
    # In case you already have the x0 calibration DB file from a previous x0 calibration 
    # and want to reuse it, just switch to Script_purpose_option 0
    #
    # The fitted distributions and self-consistency plots in pdf format from this 
    # x0 calibration can be found in the workspace/tmp-runs/*X0Calibration/ directory
    tbsw.x0script_functions.xx0calibration(rawfile_alu_list, steerfiles, x0caltag, caltag)

    # Generate a calibrated X/X0 image
    #
    # The calibrated radiation length image and other images, such as the beamspot
    # etc can be found in the workspace/root-files/*CalibratedX0Image.root
    tbsw.x0script_functions.xx0image(rawfile_alu_list, steerfiles, x0caltag, caltag, name_image1)

  if Script_purpose_option == 1:
    tbsw.x0script_functions.xx0calibration(rawfile_alu_list, steerfiles, x0caltag, caltag)

  if Script_purpose_option == 0:
    tbsw.x0script_functions.xx0image(rawfile_alu_list, steerfiles, x0caltag, caltag, name_image1)

