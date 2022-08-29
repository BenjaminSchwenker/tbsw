"""
This is an example script to demonstrate how TBSW can be used to create an X0 image from a 
test beam experiment.

The script assumes data was taken with the EUDET/AIDA tracking telescope in a monoenergetic beam.
In order to calibrate the tracking telescope, it is necessary to take data with no scattering object
placed in between the telescope arm. This situation is called an air run. The next step is to 
calibrate the X0 measurements by putting a series of Al plates with known thicknesses in between 
the telesocpe arms. After the x0 calibration data was taken, the scattering object to be studied can be 
placed in the telescope. Neither the beam energy nor the telescope planes should be touched in between 
these runs. 

The script has four main steps: 

1) Telescope calibration: Hot pixel masking, telescope alignment and cluster calibration (only from air run)
2) Angle reconstruction: Reconstruction of kink angles on the scattering object (from all runs: air, Al, DUT)   
3) X0 calibration: Obtain X0 calibration constants using X0 calibraiton runs (all runs with known scatterer) 
4) X0 imaging: Obtain calibrated X0 images of unknown scattering objects   

The script must be modified by most users to their specific use case. Places where such modifications 
are needed are commented in the text below. Examples for such modifactions are the pathes to the rawdata
files, settings for beam energy and the definition of objects used for the X0 calibration. 

Example .raw files from a X0 test beam can be downloaded. First create a directory (for example at 
'./rawdata') where they are stored. Afterwards download the files:

mkdir rawdata
wget -O ./rawdata/run006958.raw  https://owncloud.gwdg.de/index.php/s/NKdExF0pgz4G3UA/download
wget -O ./rawdata/run006965.raw  https://owncloud.gwdg.de/index.php/s/7a3SXRqHGVTOQnn/download
wget -O ./rawdata/run006973.raw  https://owncloud.gwdg.de/index.php/s/DkiJTIaHSJWgXgf/download

Take into account that the total size of these files is ~3 GB and the download may take a while.

The runs have three different targets: Air (run006973.raw), 0.5mm (run006958.raw) and 1mm (run006965.raw) 
of aluminium. These three runs can be directly used to perform a X/X0 analysis with this script.

The steering files for in order to process this test beam is shipped with tbsw and 
can be found at steering-files/x0-tb-oct16. 

Usage: 

python x0-reco.py --startStep 1 --stopStep 4

Have fun with X0 imaging. 

Author: Ulf Stolzenberg <ulf.stolzenberg@phys.uni-goettingen.de>  
Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import x0script_functions
import os
import multiprocessing
import argparse

# Script purpose option: Determines which steps should be 
# processed by the script. All steps with associated integer values
# between the start and stop parameter are conducted. The following list 
# provides the integer values of the individual reconstruction steps:

# 1: Telescope calibration (and target alignment if enabled)
# 2: Angle reconstruction
# 3: X0 calibration
# 4: X0 imaging

# This default setting means that all reconstruction steps are performed
parser = argparse.ArgumentParser(description="Perform calibration and reconstruction of a test beam run")
parser.add_argument('--startStep', dest='startStep', default=0, type=int, help='Start processing at this step number. Steps are 1) Telescope calibration, 2) Angle reconstruction, 3) X0 calibration, 4) X0 imaging')
parser.add_argument('--stopStep', dest='stopStep', default=4, type=int, help='Stop processing at this step number. Steps are 1) Telescope calibration, 2) Angle reconstruction, 3) X0 calibration, 4) X0 imaging')
args = parser.parse_args()

# Path to steering files 
# Folder contains a gear file detailing the detector geometry and a config file
# for x0 calibration. Users will likely want to rename this folder. 
steerfiles = 'steering-files/x0-tb-oct16/'

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
caltag='tboct16-2GeV'

# Name of the gearfile, which describes the telescope setup 
# Must be placed in the steerfiles folder
gearfile = 'gear.xml'

# Determine cluster resolution and store in cluster DB?
Use_clusterDB=True

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

# Flag to indicate that real EUTelescope data is used (raw format)
mcdata=False

# By default, the z position of the X0 target is defined by an mechanical survey
# measurement. For sufficiently thick targets, we can correct the z position of 
# the X0 target by forming a vertex from the upstream and downstream track. 
# This option is deactivated by default. If you know what you are doing, you can 
# enable it by using a positive number of iterations (expert decision)  
targetalignment_iterations=0

# File names and lists of filenames for the different steps 

# global path to raw files
rawfile_path = os.path.abspath('./rawdata') 

# Eudaq raw files used during telescope calibration. Telescope calibration 
# includes the alignment of the reference telescope and the calibration
# of its spatial resolution. 
# Best use a run without a scattering target in between the telescope arms.
# The calibration has to be done for every telescope setup, beam energy and m26 threshold settings
cali_run='run006973.raw'
rawfile_cali = os.path.join(rawfile_path, cali_run)

# raw file used for target alignment (only useful with a thick (X/X0 > 5 %) scattering target)
TA_run='run006958.raw'
rawfile_TA = os.path.join(rawfile_path, TA_run)  

# List of runs, which are used as input for the scattering angle reconstruction
# The angle reconstruction step is essential and every run, that will be used later during the x0 calibration or x0 imaging steps, must be listed
RunList_reco = [
		    'run006973.raw', #air
		    'run006965.raw', #0.5 mm Alu
		    'run006958.raw', #1 mm Alu
          ]

RawfileList_reco = [os.path.join(rawfile_path, x)  for x in RunList_reco]

# List of runs, which are input for the x0 calibration
# Typically runs with various different materials and thicknesses have to be used to achieve a sensible calibration.
# Good results can for example be achieved with air (no material between the telescope arms) and two different 
# aluminium thicknesses.
#
# In most cases two 1-2 runs with approximately 1 million events per thickness/material are sufficient for a good x0 calibration
#
# The different measurement regions and other options have to be set in the x0.cfg file in the steer files directory
RunList_x0cali = [
		    'run006973.raw', #air
		    'run006965.raw', #0.5 mm Alu
		    'run006958.raw', #1 mm Alu
          ]

RawfileList_x0cali = [os.path.join(rawfile_path, x) for x in RunList_x0cali]

# List of runs, which are input for the first x0 image
# Use only runs, with exactly the same target material and positioning
RunList_x0image = [
		    'run006958.raw', #1mm alu
          ]

# Set the name of this image
name_image1='1mm-alu'

RawfileList_x0image = [os.path.join(rawfile_path, x) for x in RunList_x0image]


# List of runs, which are input for the second x0 image
# Remove comment in case you want to produce more than one image
#RunList_x0image2 = [
#		    'run000209.raw', #other material
#		    'run000210.raw', #other material
#          ]

# Set the name of this image
name_image2='image2'

#RawfileList_x0image2 = [rawfile_path+x for x in RunList_x0image2]

# Number of events ...
# for telescope calibration
nevents_cali = 700000

# for target alignment
nevents_TA = 1000000

# for angle reconstruction (-1 use all available events)
nevents_reco = -1   

# Format of telescope data
csvdata = False
  
if __name__ == '__main__':


  # Calibrate the telescope 
  # In case you already have all the DB files from another telescope calibration 
  # and want to reuse it, just set startStep > 1
  #
  # DQM plots like track p/chi2 values, residuals and other interesting parameters
  # from this telescope calibration step can be found as pdf files in 
  # workspace/results/telescopeDQM

  if args.startStep < 2 and args.stopStep >= 1:
    x0script_functions.calibrate( rawfile_cali, steerfiles, caltag, gearfile, nevents_cali, Use_clusterDB, beamenergy, mcdata, Use_LongTelescopeCali, UseOuterPlanesForClusterDB, csvdata=csvdata)

    # Target alignment
    for it in range(0,targetalignment_iterations):
      x0script_functions.targetalignment(rawfile_TA, steerfiles, it, caltag, gearfile, nevents_TA, Use_clusterDB, beamenergy, mcdata, csvdata=csvdata)

  # Angle reconstruction
  # In case you already have reconstructed the scattering angles for all
  # the runs you are interested in, just set startStep > 2
  #
  # The root files with the reconstructed angles and other parameters (see 
  # README_X0.md for a full list and some descriptions) can be found in 
  # workspace/root-files/X0-run*runnumber, etc*-reco.root
  # The histmap and angle resolution for every single run can be found in 
  # workspace/results/anglerecoDQM/
  if args.startStep < 3 and args.stopStep >= 2:
    params_reco=[(x, steerfiles, caltag, gearfile, nevents_reco, Use_SingleHitSeeding, Use_clusterDB, beamenergy, mcdata) for x in RawfileList_reco]
    print("The parameters for the reconstruction are: ") 
    print(params_reco)
    
    def work(params):
      """ Multiprocessing work
      """
      rawfile, steerfiles, caltag, gearfile, nevents, Use_SingleHitSeeding, Use_clusterDB, beamenergy, mcdata = params
      x0script_functions.reconstruct(rawfile, steerfiles, caltag, gearfile, nevents, Use_SingleHitSeeding, Use_clusterDB, beamenergy, mcdata, csvdata=csvdata)
                  
    count = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=count)
    results = pool.map_async(work, params_reco)
    pool.close()
    pool.join()

    for rawfile in RawfileList_reco:
      params=(rawfile, caltag)
      x0script_functions.reconstruction_DQM(rawfile, caltag)

  # Start x0 calibration
  # In case you already have the x0 calibration DB file from a previous x0 calibration 
  # and want to reuse it, just set startStep > 3
  #
  # The fitted distributions and self-consistency plots in pdf format from this 
  # x0 calibration can be found in the workspace/tmp-runs/*X0Calibration/ directory
  if args.startStep < 4 and args.stopStep >= 3:
    x0script_functions.xx0calibration(RawfileList_x0cali, steerfiles, caltag)

  # Generate a calibrated X/X0 image
  #
  # The calibrated radiation length image and other images, such as the beamspot
  # etc can be found in the workspace/root-files/*CalibratedX0Image.root
  if args.startStep < 5 and args.stopStep >= 4:
    x0script_functions.xx0image(RawfileList_x0image, steerfiles, caltag, name_image1)

  # Generate another calibrated X/X0 image
  # The X/X0 image step can be repeated multiple times to generate a set of images
  # Just remove the comment and add a run list for each image
    #x0script_functions.xx0image(RawfileList_x0image2, steerfiles, caltag, name_image2)

