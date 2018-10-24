"""
This is an example script to demonstrate how TBSW can be used to create an X0 image from a simulated
test beam experiment.

The script below simulates a test beam experiment where charged tracks from a monoenergetic beam 
cross a a misaligned pixel telescope containing six Mimosa 26 detector planes. This script is 
specificially setup to handle test beam experiments with a very long telescope. In this case the 
default script might not work, because hit correlations between the first M26 sensor and the sensors
in the second telescope arm may be very weak. The solution implemented here is to employ a track based
correlation algorithm and iteratively add hits from the second telescope arm to the tracks.

Two data sets are simulated. A first 'air' run is done with no additional scatterer between the 
telescope arms. The 'air' run is used to calibrate the telescope. In a second 'aluminium' run, a 
aluminium plate with a well known thickness profile is inserted in between the telescope arms. 
This second run is used to compute a X0 image from the reconstructed scattering angles. The known 
comparison between the reconstructed X0 image and the a priori known image is used to calibrate
the beam energy and the angular resolution of the telescope. This second step completes the 
calibration of the telescope for X0 imaging. 

Author: Ulf Stolzenberg <ulf.stolzenberg@phys.uni-goettingen.de>  
"""

from tbsw import *
import multiprocessing

import tbsw.x0imaging.X0Calibration

# Path to steering files 
# Steeringfiles are xml files and define details of the simulation like how many events are produced
# or how M26 sensors are digitized. XML parameters can be adjusted using any test editor

# Steerfiles for the telescope calibration
steerfiles_cali = 'steering-files/x0-tb17/'

# Steerfiles for the angle reconstruction (can be the same directory as telescope calibration steerfiles)
steerfiles_reco = 'steering-files/x0-tb17/'

# Steerfiles for the x0calibration/x0imaging (can be the same directory as telescope calibration steerfiles)
steerfiles_x0 = 'steering-files/x0-tb17/'

# Nominal Beam energy
beamenergy=2.0

# cal tags
# telescope calibration cal tag (typically named after telescope setup, beam energy etc.)
caltag='tb2017_PB'

# x0 calibration cal tag
x0tag='air-alu-2GeV'

# Name of the gearfile, which describes the telescope setup 
gearfile = 'gear.xml'
gearfile_longtelescope = 'gear_long_telescope.xml'

# Alignment DB file name
alignmentdb_filename='alignmentDB.root'

# Determine cluster resolution and store in cluster DB?
Use_clusterDB=True

# Use Single Hit seeding to speed up track finding?
Use_SingleHitSeeding=True

# Use Single Hit seeding to speed up track finding?
Use_LongTelescopeCali=True

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

# raw file used during telescope calibration (best use data with scattering target)
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

RawfileList_x0image = [rawfile_path+x for x in RunList_x0image]


# List of runs, which are input for the second x0 image
# Remove comment in case you want to produce more than one image
#RunList_x0image2 = [
#		    'run000209.raw', #other material
#		    'run000210.raw', #other material
#          ]

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
  CalObj = Calibration(steerfiles=steerfiles, name=caltag + '-cal') 
  CalObj.set_beam_momentum(beamenergy)
  
  # Create list of calibration steps 
  if Use_LongTelescopeCali:
    calpath = create_x0calibration_longtelescope_path(CalObj, rawfile, gearfile_longtelescope, nevents_cali, Use_clusterDB)
  else:
    calpath = create_x0calibration_path(CalObj, rawfile, gearfile, nevents_cali, Use_clusterDB)
  
  # Run the calibration steps 
  CalObj.calibrate(path=calpath,ifile=rawfile,caltag=caltag)

  DQMplots.calibration_DQMPlots(caltag, Use_clusterDB)


# Perform the angle reconstruction of a single run
def reconstruct(params):

  rawfile, steerfiles, caltag = params

  # Set cal tag that includes run name
  name = os.path.splitext(os.path.basename(rawfile))[0] + '-' + caltag
  
  # Reconsruct the rawfile using the caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=name )
  RecObj.set_beam_momentum(beamenergy)

  # Create reconstuction path
  if Use_LongTelescopeCali:
    recopath = create_x0reco_path(RecObj, rawfile, gearfile_longtelescope, nevents_reco, Use_SingleHitSeeding)  

  else:
    recopath = create_x0reco_path(RecObj, rawfile, gearfile, nevents_reco, Use_SingleHitSeeding)  

  # Use caltag of the last target alignment iteration
  iteration_string='-target-alignment-it'+str(targetalignment_iterations-1)
  localcaltag=caltag+iteration_string
 
  if targetalignment_iterations < 1:
    localcaltag=caltag

  # Run the reconstuction  
  RecObj.reconstruct(path=recopath,ifile=rawfile,caltag=localcaltag) 

# Perform target alignment
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
  recopath = create_x0reco_path(RecObj, rawfile, gearfile, nevents_TA, Use_SingleHitSeeding)  

  # Run the reconstuction  
  RecObj.reconstruct(path=reco,ifile=rawfile,caltag=localcaltag) 

  # Read the vertex position and save it in the alignmentDB
  dbname=RecObj.create_dbfilename("alignmentDB.root")
  treename=RecObj.get_rootfilename('X0')
  save_targetpos(treename,dbname)
  RecObj.export_caltag(newcaltag)

# Perform x0 calibration
def xx0calibration(params):

  x0tag, RunList, steerfiles, calibrationtag= params

  # Generate a uncalibrated X/X0 image
  imagenametag='X0image-calitarget-Uncalibrated'
  tbsw.x0imaging.X0Calibration.x0imaging(filelist=RunList,caltag='',steerfiles=steerfiles,nametag=imagenametag)

  # Path to uncalibrated X0 image file
  imagefilename='/root-files/'+imagenametag

  # Do a calibration of the angle resolution
  tbsw.x0imaging.X0Calibration.x0calibration(filelist=RunList,imagefilename=imagefilename,caltag=x0tag,steerfiles=steerfiles)

  nametag='X0Calibration-'+x0tag
  DQMplots.x0calibration_DQMPlots(nametag=nametag)


# Generate x0 image
def xx0image(params):

  x0caltag, RunList, steerfiles, calibrationtag, listnametag = params

  if listnametag=='':
    print("No image name found. Using default naming scheme!")
    listnametag='X0image'

  # Determine name of image
  if caltag=='':
    nametag=listnametag+'-Uncalibrated'
  else:
    nametag=listnametag+'-Calibrated-'+x0tag

  # Do a calibration of the angle resolution
  tbsw.x0imaging.X0Calibration.x0imaging(filelist=RunList,caltag=x0caltag,steerfiles=steerfiles_reco,nametag=nametag)

  DQMplots.x0image_Plots(nametag=nametag)
    
  
if __name__ == '__main__':


  # Calibrate the telescope 
  # In case you already have all the DB files from another telescope calibration 
  # and want to reuse it, just switch to Script_purpose_option 0 or 1
  #
  # DQM plots like track p/chi2 values, residuals and other interesting parameters
  # from this telescope calibration step can be found as pdf files in 
  # workspace/results/telescopeDQM
  
  if Script_purpose_option !=0 and Script_purpose_option !=1:
    params_cali = ( rawfile_cali, steerfiles_cali, caltag)
    calibrate( params_cali )

  # Target alignment
  if Script_purpose_option !=0 and Script_purpose_option !=1:
    for it in range(0,targetalignment_iterations):
      params_TA = (rawfile_TA, steerfiles_reco, caltag, it)
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
    params_reco=[(x, steerfiles_reco, caltag) for x in RawfileList_reco]
    print "The parameters for the reconstruction are: " 
    print params_reco

    count = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=count)
    pool.map(reconstruct, params_reco)

    for rawfile in RawfileList_reco:
      runspec = os.path.splitext(os.path.basename(rawfile))[0] + '-'
      DQMplots.anglereco_DQMPlots(runspec,caltag)

  # start x0 calibration
  # In case you already have the x0 calibration DB file from a previous x0 calibration 
  # and want to reuse it, just switch to Script_purpose_option 0
  #
  # The fitted distributions and self-consistency plots in pdf format from this 
  # x0 calibration can be found in the workspace/tmp-runs/*X0Calibration/ directory
  if Script_purpose_option !=0:
    RootFileList_x0cali=[]
    tbsw.x0imaging.X0Calibration.CreateRootFileList(rawlist=RawfileList_x0cali,rootlist=RootFileList_x0cali, caltag=caltag)
    params_x0cali = ( x0tag, RootFileList_x0cali, steerfiles_x0, caltag)
    xx0calibration(params_x0cali)

  # Generate a calibrated X/X0 image
  #
  # The calibrated radiation length image and other images, such as the beamspot
  # etc can be found in the workspace/root-files/*CalibratedX0Image.root
  listnametag='image1'
  RootFileList_x0image=[]
  tbsw.x0imaging.X0Calibration.CreateRootFileList(rawlist=RawfileList_x0image,rootlist=RootFileList_x0image, caltag=caltag)
  params_x0image = ( x0tag, RootFileList_x0image, steerfiles_x0, caltag, listnametag)
  xx0image(params_x0image)

  # Generate another calibrated X/X0 image
  # The X/X0 image step can be repeated multiple times to generate a set of images
  # Just remove the comment and add a run list for each image
  #nametag='image2'
  #params_x0image = ( x0tag, RawfileList_x0image2, steerfiles_x0, caltag, deletetag, nametag)
  #xx0image(params_x0image)

