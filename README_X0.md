
#Introduction to X/X0 measurements with the test beam software (tbsw)

This README is a step-by-step explanation of how to generate a calibrated X0 image with the test beam software framework. The 
data reconstruction will be explained on the simple example script workspace/tbsw_x0.py. The script simulates a test beam 
experiment where charged tracks cross a misaligned pixel telescope containing six Mimosa 26 detector planes and a centered DUT.
Afterwards, the simulated raw data is calibrated and reconstucted. 

__ Due to the large track sample needed for the X/X0 imaging running this script will take several hours!__

There are a number of different reconstruction steps, which are necessary for
generating a calibrated X/X0 image. All of these steps are included in the
example script. The different steps are described in the following bullet points:


##0. Simulation of telescope digits:

The example script generates two runs with raw digits stored in two slcio files. The first run, called air run,
is an simulation of the telescope without any additional material within the telescope arms. It is used
for the calibration of the telescope, i.e. alignment, hot pixel masking and cluster calibration. The second run, called 
aluminium run, includes an aluminium plate as a calibration object in the center of the telescope. The data from 
the aluminium run  will be reconstructed (using the calibration from the air run) and a calibrated image of the the 
material distribution will be generated. 

##1. Telescope calibration:

During the telescope calibration, noisy pixel masks is produced and cluster calibration and telescope alignment
is carried out. The steering files, which are used during this process are described in README.md and can
be found in workspace/steering-files/x0-sim/ .The calibration results can be found in the folder  
workspace/cal-files/default/. 

It is preferred to use an air run for the telescope calibration. Afterwards, scattering objects must be installed 
w/o moving the telescope arms. The reason for this approach is that unkown, or wrongly modelled, material inside
the telescope may bias the telescope alignment. The size of possible bias effects can be estimated using toy simulations 
within tbsw. 

##2. Scattering angle reconstruction:

Using the calibration results from the previous step the scattering angles on the aluminium DUT are
reconstructed. During the angle reconstruction step the workspace/steering-files/x0-sim/reco.xml
steering file is employed. The reco.xml for the X0 analysis contains the following processors:

   _1. M26Clusterizer:_			  Input: NoiseDB (must be present in cal-files/default) and M26 digit collectiion, Output: M26 clusters
   
   _2. GoeHitMaker:_			  Input: clusterDB (must be present in cal-files/default) and M26 clusters, Output: M26 hits 
   
   _3. Fastracker (downstream):_  Input: alignmentDB (must be in cal-files/default) and M26 hits (only downstream hits), Output: downstream tracks 
   
   _4. Fastracker (upstream):_    Input: alignmentDB (must be in cal-files/cal-tag) and M26 hits (only upstream hits), Output: upstream tracks
   
   _5. X0ImageProducer:_          Input: alignmentDB (must be present in cal-files/cal-tag), Get kink calculates angles from tracks, Output: X0 root file 

The results of the angle reconstruction can be found at workspace/root-files/X0-mc-reco-cal-tag-default.root. The results file contains a root 
tree named MSCTree. The following variables are stored in the tree:

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
   

##3. X/X0 imaging (uncalibrated):

The tree, which was generated in the last step will now be used to generate a uncalibrated
radiation length image of the aluminium DUT. The imaging procedure employs a cfg file 
(see "workspace/steering-files/x0-sim/image.cfg"). The parameters in the config file are

  * lambda			         : The calibration factor of the angle resolution sigma, during this first
                              uncalibrated imaging, lambda should be 1.0
  * momentumoffset           : mean value of the momentum/beam energy in the center of the image (u=0,v=0)
  * momentumugradient        : momentum gradient of the beam in u direction
  * momentumvgradient        : momentum gradient of the beam in u direction
  * resultsfilename          : name of the results file, which is produced by the GenerateImage.py script
  * u_length/v_length        : side lengths of the total image, should be 20 x 10 mm^2 for an image of the whole beam spot 
  * umin/vmax                : Position of the upper left corner of the total image, should be (umin,vmax)=(-10mm,5mm) for 
                               an image of the whole beam spot
  * u/v pixelsize            : Pixel size of the image in microns, should be chosen in such away, that at least 1000 particles 
                               pass through the area of most pixels

After the imaging process the uncalibrated radiation length image is located in 
workspace/root-files/X0-mc-alu-default-reco-Uncalibrated-X0image.root . The file contains many 2D histograms, the
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

## Calibration of beam energy and telescope angle resolution:

The relevant calibration constants ares: 

   * lambda:                  calibration factor of the angle resolution sigma on the central target plane
   * momentumugradient:       Momentum gradient (in u direction) of beam particles in GeV/mm
   * momentumvgradient:       Momentum gradient (in v direction) of beam particles in GeV/mm
   * momentumoffset:          Momentum value of beam particles at (u,v)=(0,0) in GeV

The average beam momentum at DESY is parametrized by a linear function 

momentum = momentumoffset+momentumugradient*u+momentumvgradient*v

The calibration parameters are measured by fitting a X0 image from a calibration target. The calibration target is
expected to be a planar object with areas of well defined  X/X0. For the calibration software, the calibration target
must be described in a cfg file. This file is located at workspace/steering-files/x0-sim/x0calibration.cfg. It is 
used to to set starting parameters, measurement areas (MA) and fit options. The measurement areas can be added
as simple rectangular shapes (MA in the cfg file) with a center position, length, thickness and material parameters. 
Alternatively a line can be defined, which constructs multiple measurement areas at the same time. More detailed 
descriptions of the options can be found in the cfg file itself. 

In the example script a aluminium plate with a material step is used as the DUT. The cfg file we are using here employs
several measurement areas and a line in order to fit the alibration parameters. The measurement areas are marked on the
following image of the uncalibrated aluminium plate:

[picture](workspace/tbsw_tools/validation/X0image_Boxes.png)

There is one measurement area in the 0.5mm thick aluminium (1), 6 measurement areas in the 1mm aluminium area (2,3 and 9-12)
the remaining measurement areas lie in the air region.

The results of the calibration, including pictures of the fits of the individual measurement areas, can be found in 
workspace/tmp-runs/X0-mc-alu-default-reco-X0Calibration/ . The results of the calibration are also stored as a cfg file
in workspace/cal-files/default/x0cal_result.cfg.

## X/X0 imaging (calibrated):

In the last step an calibrated X/X0 image is produced, which used the cfg file from the previous calibration step. The
calibrated X/X0 image should look like this:

[picture](workspace/tbsw_tools/validation/X0image.png)

The different aluminium thicknesses can be seen. The area with small X/X0 values surrounding the aluminium plate is the surrounding air.

A X/X0 profile cut can be used for X/X0 measurements. In this example case the steps between thick and thin aluminium, as well as air is clearly visible:

[picture](workspace/tbsw_tools/validation/X0profile.png)

The complete results can be found in workspace/root-files/X0-mc-alu-default-reco-Calibrated-X0image.root

## Analysis of test beam data:

If you are performing an analysis of actual test beam data, it is probably a good idea to start with the example script,
which was explained here and change it accordingly in order to generate radiation length images of your DUTs. The most
important changes are the following

   * You will want to use .raw files recorded during the test beam instead of simulated .slcio files. This change is not
     requiring any manual changes, you just have to use another steering-files folder (steering-files/x0-tb)
   * The gear files in x0-tb, describing the telescope geometry, should be changed according to the real setup 
   * Some changes in the x0 cfg files might be necessary (changing the measurement areas during the X0 calibration etc.)

An example script for analysing test beam data can be found in workspace/testbeam_x0.py.


Ulf Stolzenberg

Goettingen 2017

ulf.stolzenberg@phys.uni-goettingen.de 



 

















