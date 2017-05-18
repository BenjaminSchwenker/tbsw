
# Introduction to X/X0 measurements with the test beam analysis software (tbsw)

This README is a step-by-step explanation of how to generate a calibrated X0 image with the test beam software framework. The 
data reconstruction will be explained on the simple example script workspace/tbsw_x0.py. Its source code is:


```
#!python

"""
This is an example script to demonstrate how TBSW can be used to analyze test beam 
data using Python scripts.

The script below simulates a test beam experiment where charged tracks cross a misaligned
pixel telescope containing six Mimosa 26 detector planes and a 05mm aluminium plate,
centered in the telescope. Afterwards, the simulated raw data is calibrated and
reconstucted. Afterwards an X/X0 of the aluminium plate is generated and a
calibration of the angle resolution of the telescope is performed

Author: Ulf Stolzenberg <ulf.stolzenberg@phys.uni-goettingen.de>  
"""

from tbsw import *

# Path to steering files 
steerfiles = 'steering-files/x0-sim/'
# Tag for calibration data 
caltag = 'default'
# File name for raw data 
runname='mc' 
rawfile = runname+'.slcio'

# Defines the sequence of calibration steps. 
# XML steer files are taken from steerfiles. 
calpath = [ 
           'hotpixelkiller.xml' ,              
           'cluster-calibration-mc.xml',     # creates clusterDB, but it will not be used
           'correlator-iteration-1.xml' ,
           'kalmanalign-iteration-1.xml',
           'kalmanalign-iteration-2.xml',
           'kalmanalign-iteration-2.xml',
           'kalmanalign-iteration-2.xml',
           'telescope-dqm-iteration-1.xml',
           'cluster-calibration-tb-iteration-1.xml',
           'cluster-calibration-tb-iteration-2.xml',
           'cluster-calibration-tb-iteration-2.xml',
           'cluster-calibration-tb-iteration-2.xml',
           'cluster-calibration-tb-iteration-2.xml',
           'cluster-calibration-tb-iteration-2.xml',
           'cluster-calibration-tb-iteration-2.xml',
           'correlator-iteration-2.xml',
           'kalmanalign-iteration-3.xml',
           'kalmanalign-iteration-4.xml',
           'kalmanalign-iteration-4.xml',
           'kalmanalign-iteration-4.xml',
           'telescope-dqm-iteration-2.xml',
         ]

# Base name for temporary folder created in tmp-runs/ 
name = os.path.splitext(os.path.basename(rawfile))[0] + '-' + caltag  

# Simulate a rawfile from a test beam experiment
# SimObj creates folder tmp-runs/name-sim/ and populates it with 
# copies of steering files. After processing, the folder contains
# also logfiles. 
#SimObj = Simulation(steerfiles=steerfiles, name=name + '-sim' )
#SimObj.simulate(path=['simulation.xml'], ofile=rawfile, caltag=None)  
   
# Calibrate the telescope using the rawfile. Creates a folder caltag 
# containing all calibrations. 
#CalObj = Calibration(steerfiles=steerfiles, name=name + '-cal') 
#CalObj.calibrate(path=calpath,ifile=rawfile,caltag=caltag)  
   
# Reconsruct the rawfile using caltag. Resulting root files are 
# written to folder root-files/
#RecObj = Reconstruction(steerfiles=steerfiles, name=name + '-reco' )
#RecObj.reconstruct(path=['reco.xml'],ifile=rawfile,caltag=caltag) 

filename='root-files/cal-tag-default/X0-mc-default-reco.root'
imagefilename='root-files/cal-tag-default/uncalibrated_x0image/X0-completeimage.root'
deletetag=1

# Function which starts the imaging script
def x0imaging(filename,caltag,deletetag):

  flag='./root-scripts/x0imaging/GenerateImage.py -i '+filename+' -c '+caltag+' -d '+`deletetag`
  print(flag)
  subprocess.call(flag, shell=True)

  return None

# Function which starts the x0 calibration script
def x0calibration(filename,imagefilename,caltag):

  flag='./root-scripts/x0imaging/X0Calibration.py -i '+filename+'-m '+imagefilename+' -c '+caltag
  print(flag)
  subprocess.call(flag, shell=True)

  return None


# Create a new folder at workspace/root-files/cal-tag-default/ and store the X/X0 results there
cwdir = os.getcwd()
rootpath=cwdir+'/root-files/cal-tag-'+caltag+'/' 
if not os.path.isdir(rootpath):
   os.mkdir(rootpath)

# Copy root file
rootfile=cwdir+'/root-files/X0-'+name+'-'+'reco.root'
shutil.copy(cwdir+'/root-files/X0-mc-default-reco.root', rootpath)  
  
# Generate a uncalibrated X/X0 image
x0imaging(filename,'',deletetag)

# Rename the directory in which the imaging results are stored
if os.path.isdir(rootpath+'x0image'):
   shutil.move(rootpath+'x0image', rootpath+'uncalibrated_x0image')  

# Do a calibration of the angle resolution
x0calibration(filename,imagefilename,caltag)

# Generate a calibrated X/X0 image
x0imaging(filename,caltag,deletetag)

# Rename the directory in which the imaging results are stored
if os.path.isdir(rootpath+'x0image'):
   shutil.move(rootpath+'x0image', rootpath+'calibrated_x0image')  


```

There are a number of different reconstruction steps, which are necessary for
generating a calibrated X/X0 image. All of these steps are included in the
example script. The different steps are described in the following bullet points:

1) Telescope Calibration: 

   In this first step a noisy pixel mask, a cluster calibration and telescope alignment is performed.
   The steering files, which are used during this process are described in README.md and can be
   found in workspace/steering-files/x0-sim/ .The calibration results can be found in 
   workspace/cal-files/default/ .

2) Scattering angle Reconstruction:

   Using the calibration results from the previous step the scattering angles on the DUT are
   reconstructed. During the angle reconstruction step the workspace/steering-files/x0-sim/reco.xml
   steering file is used. The reco.xml for the X0 analysis contains the following processors:

     2.1) M26Clusterizer:			Input: NoiseDB (must be present in cal-files/default) and M26 digit collection, 
                         	    	Output: M26 clusters
     2.2) GoeHitMaker:				Input: clusterDB (must be present in cal-files/default) and M26 clusters, 
                        	    	Output: M26 hits 
	 2.3) Fastracker (downstream):  Input: alignmentDB (must be in cal-files/default) and M26 hits (only downstream hits), 
                                    Output: downstream tracks 
     2.4) Fastracker (upstream):    Input: alignmentDB (must be in cal-files/cal-tag) and M26 hits (only upstream hits), 
                                    Output: upstream tracks
     2.7) X0ImageProducer:          Input: alignmentDB (must be present in cal-files/cal-tag),
                                    get kink calculates angles from tracks, Output: X0 root file 

   The results of the angle reconstruction can be found at 
   workspace/root-files/cal-tag-default/X0-mc-reco-cal-tag-default.root . The results file contains a 
   root tree named MSCTree. The following variables are stored in the tree:

     *  iRun            : Run number of the track
     *  iEvt            : Event number of the track
     *  prob_up         : p value of the upstream track
     *  prob_down       : p value of the downstream track
     *  prob_combo      : p value of the combined track
     *  du/dw           : mean of the two track slopes in u-w plane (rad)
     *  dv/dw           : mean of the two track slopes in v-w plane (rad)
     *  u_in            : u intersection of upstream track on the target (in local coordinates (mm))
     *  v_in            : v intersection of upstream track on the target (in local coordinates (mm))
     *  u_out           : u intersection of downstream track on the target (in local coordinates (mm))
     *  v_out           : v intersection of downstream track on the target (in local coordinates (mm))
     *  u               : mean of u_in and u_out (in local coordinates (mm))
     *  v               : mean of v_in and v_out (in local coordinates (mm))
     *  u_var           : variance of u (mm^2)
     *  v_var           : variance of v (mm^2)
     *  theta1          : Projected kink angle in u-w plane (rad)
     *  theta2          : Projected kink angle inv-w plane (rad)
     *  theta1_var      : Theta1 variance calculated via error propagation from Kalman filter (rad^2)
     *  theta2_var      : Theta2 variance calculated via error propagation from Kalman filter (rad^2)
     *  momentum        : Track momentum in GeV
     *  vertex_chi2     : Chi2 of the determined vertex of the up- and downstream track
     *  momentum        : p value of the determined vertex of the up- and downstream track

3) X/X0 Imaging (uncalibrated):

     The tree, which was generated in the last step will now be used to generate a uncalibrated
     radiation length image of the DUT. The imaging procedure employs a cfg file 
     (see "workspace/root-scripts/x0imaging/image.cfg"). The parameters in the config file are

     * lambda			        : The calibration factor of the angle resolution sigma, during this first
                                  uncalibrated imaging, lambda should be 1.0
     * momentumoffset           : mean value of the momentum/beam energy in the center of the image (u=0,v=0)
     * momentumugradient        : momentum gradient of the beam in u direction
     * momentumvgradient        : momentum gradient of the beam in u direction
     * resultsfilename          : name of the results file, which is produced by the GenerateImage.py script
     * u_length/v_length        : side lengths of the total image, should be 20 x 10 mm² for an image of the whole beam spot 
     * umin/vmax                : Position of the upper left corner of the total image, should be (umin,vmax)=(-10mm,5mm) for 
                                  an image of the whole beam spot
     * u/v pixelsize            : Pixel size of the image in µm, should be chosen in such away, that at least 1000 particles 
                                  pass through the area of most pixels

     After the imaging process the uncalibrated radiation length image is located in 
     workspace/root-files/cal-tag-default/uncalibrated_x0image/X0-completeimage.root . The file contains many 2D histograms, the
     most important ones are:

     * x0_image                         : Radiation length (X/X0) image 
     * x0err_image                      : Image of the statistical errors of the image (due to the fit)
     * x0relerr_image                   : Image of relative X/X0 errors
     * fit(1/2/sum)chi2ndof_image       : Images of the fit chi2 values
     * fit(1/2/sum)prob_image          	: p values of the fit of the projected angle distributions and the combined distribution
     * theta(1/2)mean_image             : Image of the mean values of the projected kink angle distributions (before correction)
     * correctedtheta(1/2)mean_image    : Image of the mean values of the projected kink angle distributions (after correction)
     * (u/v)residualmean_image    		: Image of the u and v residuals of down and upstream track
     * beamspot                         : Image of number of tracks 
     * BE_image                         : Image of the particle momentum/beam energy

3) Calibration of BE and telescope angle resolution:

     The relevant calibration constants ares: 

     * lambda:                  calibration factor of the angle resolution sigma on the central target plane
     * momentumugradient:       Momentum gradient (in u direction) of beam particles in GeV/mm
     * momentumvgradient:       Momentum gradient (in v direction) of beam particles in GeV/mm
     * momentumoffset:          Momentum value of beam particles at (u,v)=(0,0) in GeV

     The average beam momentum at DESY is parametrized by a linear function 

     momentum = momentumoffset+momentumugradient*u+momentumvgradient*v

     The calibration parameters are measured by fitting a X0 image from a calibration target. The calibration target is
     expected to be a planar object with areas of well defined  X/X0. For the calibration software, the calibration target
     must be described in a cfg file. This file is located at workspace/root-scripts/x0imaging/x0calibration.cfg. It is 
     used to to set starting parameters, measurement areas (MA) and fit options. The measurement areas can be added
     as simple rectangular shapes (MA in the cfg file) with a center position, length, thickness and material parameters. 
     Alternatively a line can be defined, which constructs multiple measurement areas at the same time. More detailed 
     descriptions of the options can be found in the cfg file itself. 

     In the example script a 0.5mm thick aluminium plate is used as the DUT. The cfg file we are using here employs
     two measurement areas and a line in order to fit the alibration parameters.

     The results of the calibration, including pictures of the fits, can be found in 
     workspace/root-files/cal-tag-default/x0calibration/ . The results of the calibration are also stored as a cfg file
     in workspace/cal-files/default/x0cal_result.cfg.

3) X/X0 Imaging (calibrated):

     In the last step an calibrated X/X0 image is produced, which used the cfg file from the previous calibration step.
     The results can be found in workspace/root-files/cal-tag-default/calibrated_x0image/X0-completeimage.root


If you are performing an analysis of actual test beam data, it is probably a good idea to start with the example script,
which was explained here and change it accordingly in order to generate radiation length images of your DUTs.


Ulf Stolzenberg

Göttingen 2017 

ulf.stolzenberg@phys.uni-goettingen.de 



 

















