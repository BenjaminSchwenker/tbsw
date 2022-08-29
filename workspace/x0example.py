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

Usage: 

python x0example.py


Author: Ulf Stolzenberg <ulf.stolzenberg@phys.uni-goettingen.de>  
"""

from tbsw import x0script_functions
import os
import multiprocessing
import argparse

# Script purpose option: Determines which steps should be 
# processed by the script. All steps with associated integer values
# between the start and stop parameter are conducted. The following list 
# provides the integer values of the individual reconstruction steps:

# 0: Test beam simulation
# 1: Telescope calibration (and target alignment if enabled)
# 2: Angle reconstruction
# 3: X0 calibration
# 4: X0 imaging

# This default setting means that all reconstruction steps are performed
parser = argparse.ArgumentParser(description="Perform calibration and reconstruction of a test beam run")
parser.add_argument('--startStep', dest='startStep', default=0, type=int, help='Start processing at this step number. Steps are 0) Test beam simulation, 1) Telescope calibration, 2) Angle reconstruction, 3) X0 calibration, 4) X0 imaging')
parser.add_argument('--stopStep', dest='stopStep', default=4, type=int, help='Stop processing at this step number. Steps are 0) Test beam simulation, 1) Telescope calibration, 2) Angle reconstruction, 3) X0 calibration, 4) X0 imaging')
args = parser.parse_args()

# Determine maximum number of processes
nprocesses=2



# Path to steering files 
# Folder contains a gear file detailing the detector geometry and a config file
# for x0 calibration.
steerfiles = 'steering-files/x0-tb/'

# Nominal Beam energy
beamenergy=2.0

# Definition of the calibration tag. It is typically named after telescope setup, beam energy, x0calibration target etc.
# The caltag is used to generate a directory under localDB/*caltag* where all calibration parameters are stored
# in local DB files. Additionally DQM plots of the calibration steps will be stored under results/ to cross check the
# calibration results. The telescope calibration step (Step 1 in the enumeration above) will generate a hotpixel mask (NoiseDB-M26.root),
# a files with alignment information (alignmentDB.root) and a data base containing cluster resolutions (clusterDB-M26.root). DQM plots
# of the cluster calibration are stored under results/clusterDB-M26/*caltag*. Other track based DQM plots such as track p values,
# pulls and the mean number of tracks per event are stored for example in results/TelescopeDQM2/*caltag*.
# During the radiation length calibration step (Step 3) the beam energy, the beam energy gradients and a global offset of the telescope
# angle resolution will be determined and stored in a text file (x0cal_result.cfg). The DQM plots such as a selfconsistency diagram and
# angle distributions with their associated fits can be found in results/x0calibrationDQM/*caltag*.
caltag='x0-sim'

# Name of gearfile
# This file describes the nominal geometry of a telescope 
# with two arms having three M26 planes each.  
gearfile = 'gear.xml'

# By default,the track finder constructs track seeds using hit pairs from two 
# planes. With SingleHitSeeding, seed tracks are constructed from a single hit 
# and extrapolated parallel to the z axis. This can safe time but risks missing hits. 
Use_SingleHitSeeding=False

# Finding correlations between first and last sensor can be difficult
# for low momentum tracks and/or large distances between sensors.
# By default, a robust method is used that correlates sensors step by
# step. Only deactivate when you are really certain you do not need
# this feature (expert decision).  
Use_LongTelescopeCali=True

# Switch to use clusters on outer planes to calculate cluster resolution
# The track resolution is expected to be worse on the outer planes, using them may 
# have a negative impact on the determined cluster resolutions
UseOuterPlanesForClusterDB=False

# Determine cluster resolution and store in cluster DB?
Use_clusterDB=True

# Flag to indicate that mc data is used (slcio format)
mcdata=True

# By default, the z position of the X0 target is defined by an mechanical survey
# measurement. For sufficiently thick targets, we can correct the z position of 
# the X0 target by forming a vertex from the upstream and downstream track. 
# This option is deactivated by default. If you know what you are doing, you can 
# enable it by using a positive number of iterations (expert decision)  
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
  
if __name__ == '__main__':


   
  # Get current directory
  cwdir = os.getcwd()
  
  # Create a simulated lcio for run with no target (air run) and 
  # multiple run s with a Al plate as scattering material
  if args.startStep < 1 and args.stopStep >= 0:
    x0script_functions.simulate(rawfile_air, rawfile_alu_list, steerfiles, gearfile, nevents_air, nevents_alu, nprocesses, mean_list, sigma_list, sensorexception_list, modeexception_list, beamenergy)
  
  
  # Calibrate the telescope 
  # In case you already have all the DB files from another telescope calibration 
  # and want to reuse it, just switch to Script_purpose_option 0 or 1
  #
  # DQM plots like track p/chi2 values, residuals and other interesting parameters
  # from this telescope calibration step can be found as pdf files in 
  # workspace/results/telescopeDQM
  if args.startStep < 2 and args.stopStep >= 1:
    x0script_functions.calibrate( rawfile_air, steerfiles, caltag, gearfile, nevents_air, Use_clusterDB, beamenergy, mcdata, Use_LongTelescopeCali, UseOuterPlanesForClusterDB)
   
    # Target alignment
    for it in range(0,targetalignment_iterations):
      x0script_functions.targetalignment(rawfile_alu_list[0], steerfiles, it, caltag, gearfile, nevents_TA, Use_SingleHitSeeding, Use_clusterDB, beamenergy, mcdata)
  
  # Angle reconstruction
  # In case you already have reconstructed the scattering angles for all
  # the runs you are interested in, just switch to Script_purpose_option 0 or 1
  #
  # The root files with the reconstructed angles and other parameters (see 
  # README_X0.md for a full list and some descriptions) can be found in 
  # workspace/root-files/X0-run*runnumber, etc*-reco.root
  # The histmap and angle resolution for every single run can be found in 
  # workspace/results/anglerecoDQM/
  if args.startStep < 3 and args.stopStep >= 2:
    params_reco=[(x, steerfiles, caltag, gearfile, nevents_alu, Use_SingleHitSeeding, Use_clusterDB, beamenergy, mcdata) for x in rawfile_alu_list]
    print("The parameters for the reconstruction are: ") 
    print(params_reco)
    
    def work(params):
      """ Multiprocessing work
      """
      rawfile, steerfiles, caltag, gearfile, nevents, Use_SingleHitSeeding, Use_clusterDB, beamenergy, mcdata = params
      x0script_functions.reconstruct(rawfile, steerfiles, caltag, gearfile, nevents, Use_SingleHitSeeding, Use_clusterDB, beamenergy, mcdata)
                  
    count = min(nprocesses,multiprocessing.cpu_count())
    pool = multiprocessing.Pool(processes=count)
    results = pool.map_async(work, params_reco)
    pool.close()
    pool.join()
    
    for rawfile in rawfile_alu_list:
      x0script_functions.reconstruction_DQM(rawfile, caltag)
    
  if args.startStep < 4 and args.stopStep >= 3:
    # Start x0 calibration
    # In case you already have the x0 calibration DB file from a previous x0 calibration 
    # and want to reuse it, just switch to Script_purpose_option 0
    #
    # The fitted distributions and self-consistency plots in pdf format from this 
    # x0 calibration can be found in the workspace/tmp-runs/*X0Calibration/ directory
    x0script_functions.xx0calibration(rawfile_alu_list, steerfiles, caltag)

    # Generate a calibrated X/X0 image
    #
    # The calibrated radiation length image and other images, such as the beamspot
    # etc can be found in the workspace/root-files/*CalibratedX0Image.root
  if args.startStep < 5 and args.stopStep >= 4:
    x0script_functions.xx0image(rawfile_alu_list, steerfiles, caltag, name_image1)

