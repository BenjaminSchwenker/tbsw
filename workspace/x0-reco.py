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

Have fun with X0 imaging. 

Author: Ulf Stolzenberg <ulf.stolzenberg@phys.uni-goettingen.de>  
Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

import tbsw 
import os
import multiprocessing

# Path to steering files 
# Folder contains a gear file detailing the detector geometry and a config file
# for x0 calibration. Users will likely want to rename this folder. 
steerfiles = 'steering-files/x0-tb-june17/'

# Nominal Beam energy
beamenergy=2.0

# cal tags
# telescope calibration cal tag (typically named after telescope setup, beam energy etc.)
caltag='tb2017_PB'

# x0 calibration cal tag
x0caltag='air-alu-2GeV'

# Name of the gearfile, which describes the telescope setup 
# Must be placed in the steerfiles folder
gearfile = 'gear.xml'

# Determine cluster resolution and store in cluster DB?
Use_clusterDB=True

# Use Single Hit seeding to speed up track finding?
Use_SingleHitSeeding=False

# Long telescopes may require a sensor by sensor alignment approach
Use_LongTelescopeCali=True

# Switch to use clusters on outer planes on outer planes to calculate cluster resolution
# The track resolution is expected to be worse on the outer planes, using them may 
# have a negative impact on the determined cluster resolutions
Use_OuterPlanes=False

# Flag to indicate that real EUTelescope data is used (raw format)
mcdata=False

# Script purpose option:
# 0: Script only processes imaging part
# 1: Script processes x0 calibration and imaging part
# Everything else: Script processes the whole chain: Telescope calibration, angle reconstruction, x0 calibration and x0 imaging 
Script_purpose_option=2

# Number of iterations during target alignment
# Set to 0 or negative integer to disable target alignment
targetalignment_iterations=0

# File names and lists of filenames for the different steps 

# global path to raw files
rawfile_path='/work1/rawdata/luise/'

# Eudaq raw files used during telescope calibration. Telescope calibration 
# includes the alignment of the reference telescope and the calibration
# of its spatial resolution. 
# Best use a run without a scattering target in between the telescope arms.
# The calibration has to be done for every telescope setup, beam energy and m26 threshold settings
cali_run='run000210.raw'
rawfile_cali = rawfile_path + cali_run

# raw file used for target alignment (only useful with a thick (X/X0 > 5 %) scattering target)
#TA_run='run006958.raw'
#rawfile_TA = rawfile_path + TA_run

# List of runs, which are used as input for the scattering angle reconstruction
# The angle reconstruction step is essential and every run, that will be used later during the x0 calibration or x0 imaging steps, must be listed
RunList_reco = [
		    'run000210.raw', #air
		    'run000154.raw', #4 mm Alu
		    'run000161.raw', #6 mm Alu
		    'run000171.raw', #PB 1
          ]

RawfileList_reco = [rawfile_path+x for x in RunList_reco]

# List of runs, which are input for the x0 calibration
# Typically runs with various different materials and thicknesses have to be used to achieve a sensible calibration.
# Good results can for example be achieved with air (no material between the telescope arms) and two different 
# aluminium thicknesses.
#
# In most cases two 1-2 runs with approximately 1 million events per thickness/material are sufficient for a good x0 calibration
#
# The different measurement regions and other options have to be set in the x0.cfg file in the steer files directory
RunList_x0cali = [
		    'run000210.raw', #air
		    'run000154.raw', #4 mm Alu
		    'run000161.raw', #6 mm Alu
          ]

RawfileList_x0cali = [rawfile_path+x for x in RunList_x0cali]

# List of runs, which are input for the first x0 image
# Use only runs, with exactly the same target material and positioning
RunList_x0image = [
		    'run000171.raw', #PB 1
          ]

# Set the name of this image
name_image1='image1'

RawfileList_x0image = [rawfile_path+x for x in RunList_x0image]


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
nevents_cali = 50000

# for target alignment
nevents_TA = 1000000

# for angle reconstruction (-1 use all available events)
nevents_reco = -1


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
  calpaths = tbsw.path_utils.create_x0analysis_calibration_paths(CalObj, rawfile, gearfile, nevents_cali, Use_clusterDB, beamenergy, mcdata, Use_LongTelescopeCali, Use_OuterPlanes)
  
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
  recopath = tbsw.path_utils.create_anglereco_path(RecObj, rawfile, gearfile, nevents_reco, Use_SingleHitSeeding, Use_clusterDB, beamenergy, mcdata)  
  
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
  RecObj.reconstruct(paths=recopath,ifile=rawfile,caltag=oldcaltag)  
  
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


  # Calibrate the telescope 
  # In case you already have all the DB files from another telescope calibration 
  # and want to reuse it, just switch to Script_purpose_option 0 or 1
  #
  # DQM plots like track p/chi2 values, residuals and other interesting parameters
  # from this telescope calibration step can be found as pdf files in 
  # workspace/results/telescopeDQM
  
  if Script_purpose_option !=0 and Script_purpose_option !=1:
    params_cali = ( rawfile_cali, steerfiles, caltag)
    calibrate( params_cali )

  # Target alignment
  if Script_purpose_option !=0 and Script_purpose_option !=1:
    for it in range(0,targetalignment_iterations):
      params_TA = (rawfile_TA, steerfiles, caltag, it)
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
    params_reco=[(x, steerfiles, caltag) for x in RawfileList_reco]
    print "The parameters for the reconstruction are: " 
    print params_reco

    count = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=count)
    pool.map(reconstruct, params_reco)

    for rawfile in RawfileList_reco:
      params=(rawfile, caltag)
      reconstruction_DQM(params)

  # Start x0 calibration
  # In case you already have the x0 calibration DB file from a previous x0 calibration 
  # and want to reuse it, just switch to Script_purpose_option 0
  #
  # The fitted distributions and self-consistency plots in pdf format from this 
  # x0 calibration can be found in the workspace/tmp-runs/*X0Calibration/ directory
  if Script_purpose_option !=0:
    params_x0cali = ( x0caltag, RawfileList_x0cali, steerfiles, caltag)
    xx0calibration(params_x0cali)

  # Generate a calibrated X/X0 image
  #
  # The calibrated radiation length image and other images, such as the beamspot
  # etc can be found in the workspace/root-files/*CalibratedX0Image.root
  params_x0image = ( x0caltag, RawfileList_x0image, steerfiles, caltag, name_image1)
  xx0image(params_x0image)

  # Generate another calibrated X/X0 image
  # The X/X0 image step can be repeated multiple times to generate a set of images
  # Just remove the comment and add a run list for each image
  #params_x0image = ( x0caltag, RawfileList_x0image2, steerfiles, caltag, name_image2)
  #xx0image(params_x0image)

