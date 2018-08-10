
#Introduction to X/X0 measurements with the test beam software (tbsw)

This README is a step-by-step explanation of how to generate a calibrated X0 image with the test beam software framework. The 
data processing is explained on the heavily commented example script tbsw_x0example.py. The script simulates
a test beam experiment where charged tracks cross a misaligned pixel telescope containing six Mimosa 26 detector planes and a 
centered device under test (DUT). Afterwards, the simulated raw data is calibrated and reconstucted. Finally, a calibrated X0 
image of the DUT is computed.


```
$ cp -r <path-to-tbsw-installdir>/workspace ~/workspace-x0-example
$ cd ~/workspace-x0-example
$ source init_tbsw.sh
$ python tbsw_x0example.py
```  

__ Due to the large track sample needed for the X0 imaging, running this script can take several hours!__

There are a number of different reconstruction steps, which are necessary for generating a calibrated X0 image. All of these steps 
are included in the example script. The different steps are described in the following bullet points:


##0. Simulation of telescope digits:

The example script generates two runs with raw digits stored in two slcio files. The first run, called air run,
is an simulation of the telescope without any additional material within the telescope arms. It is used
for the calibration of the telescope, i.e. alignment, hot pixel masking and cluster calibration. The second run, called 
aluminium run, includes an aluminium plate as a calibration object in the center of the telescope. The data from 
the aluminium run  will be reconstructed (using the calibration from the air run) and a calibrated image of the the 
material distribution will be generated. 

##1. Telescope calibration:

During the telescope calibration, noisy pixel masks is produced and cluster calibration and telescope alignment
is carried out. The MARLIN Processors, which are used during this process are described in README.md and can
be found in workspace/steering-files/x0-tb/processors.xml .The calibration results can be found in the folder  
workspace/localDB/mc-air-test/. 

It is preferred to use an air run for the telescope calibration. Afterwards, scattering objects must be installed 
w/o moving the telescope arms. The reason for this approach is that unkown, or wrongly modelled, material inside
the telescope may bias the telescope alignment. The size of possible bias effects can be estimated using toy simulations 
within tbsw. 

##2. Scattering angle reconstruction:

Using the calibration results from the previous step the scattering angles on the aluminium DUT are
reconstructed. During the angle reconstruction step the workspace/steering-files/x0-tb/processors.xml
steering file is employed. The following processors are used in the X0 analysis:

   _1. M26Clusterizer:_			  Input: NoiseDB (must be present in cal-files/default) and M26 digit collectiion, Output: M26 clusters
   
   _2. GoeHitMaker:_			  Input: clusterDB (must be present in cal-files/default) and M26 clusters, Output: M26 hits 
   
   _3. Fastracker (downstream):_  Input: alignmentDB (must be in cal-files/default) and M26 hits (only downstream hits), Output: downstream tracks 
   
   _4. Fastracker (upstream):_    Input: alignmentDB (must be in cal-files/cal-tag) and M26 hits (only upstream hits), Output: upstream tracks
   
   _5. X0ImageProducer:_          Input: alignmentDB (must be present in cal-files/cal-tag), Get kink calculates angles from tracks, Output: X0 root file 

The results of the angle reconstruction can be found at workspace/root-files/X0-mc-air-test-reco.root. The results file contains a root 
tree named MSCTree. The following variables are stored in the tree:

   *  iRun            			: Run number of the track
   *  iEvt            			: Event number of the track
   *  prob_up         			: p value of the upstream track
   *  prob_down       			: p value of the downstream track
   *  prob_combo      			: p value of the combined track
   *  u_in            			: u intersection of upstream track on the target (in local coordinates (mm))
   *  v_in            			: v intersection of upstream track on the target (in local coordinates (mm))
   *  u_out           			: u intersection of downstream track on the target (in local coordinates (mm))
   *  v_out           			: v intersection of downstream track on the target (in local coordinates (mm))
   *  u               			: mean of u_in and u_out (in local coordinates (mm))
   *  v               			: mean of v_in and v_out (in local coordinates (mm))
   *  u_var           			: variance of u (mm^2)
   *  v_var           			: variance of v (mm^2)
   *  theta1          			: Projected kink angle in u-w plane (rad)
   *  theta2          			: Projected kink angle inv-w plane (rad)
   *  theta1_var      			: Theta1 variance calculated via error propagation from Kalman filter (rad^2)
   *  theta2_var      			: Theta2 variance calculated via error propagation from Kalman filter (rad^2)
   *  momentum        			: Track momentum in GeV
   *  chi2            			: Chi2 calculated from the up- and downstream track residuals on the central plane
   *  prob            			: p value calculated from the up- and downstream track residuals on the central plane
   *  vertex_u        			: Vertex u position in local coordinates of the target plane (mm)
   *  vertex_v        			: Vertex v position in local coordinates of the target plane (mm)
   *  vertex_w        			: Vertex w position in local coordinates of the target plane (mm)
   *  vertex_u_var    			: Vertex u position variance in local coordinates of the target plane (mm^2)
   *  vertex_v_var    			: Vertex v position variance in local coordinates of the target plane (mm^2)
   *  vertex_w_var    			: Vertex w position variance in local coordinates of the target plane (mm^2)
   *  vertex_x        			: Vertex x position in global coordinates (mm)
   *  vertex_y        			: Vertex y position in global coordinates (mm)
   *  vertex_z        			: Vertex z position in global coordinates (mm)
   *  vertex_x_var    			: Vertex x position variance in global coordinates (mm^2)
   *  vertex_y_var    			: Vertex y position variance in global coordinates (mm^2)
   *  vertex_z_var    			: Vertex z position variance in global coordinates (mm^2)
   *  vertex_chi2     			: Chi2 value of the vertex fit
   *  vertex_prob     			: p value of the vertex fit
   *  vertex_u_res    			: Vertex u residual (mm)
   *  vertex_v_res    			: Vertex v residual (mm)
   *  vertex_multiplicity	    : Vertex multiplicity (= number of downstream track propagating from this vertex)
   *  vertex_id	   				: Vertex ID (determined by the upstream track)
   

##3. X/X0 imaging (uncalibrated):

The tree, which was generated in the last step will now be used to generate a uncalibrated
radiation length image of the aluminium DUT. The imaging procedure employs a cfg file 
(see "workspace/steering-files/x0-tb/image.cfg"). The parameters in the config file are

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
workspace/root-files/X0-mc-air-test-reco-Uncalibrated-X0image.root . The file contains many 2D histograms, the
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
  * vertex_w_mean_image              : Image of the mean vertex w position , this might be useful, when trying to align the target plane
  * vertex_w_rms_image              : Image of the vertex w position RMS, this parameter can be used as a complementary thickness measurement
  * vertex_chi2_image              : Image of the vertex fit chi2 value
  * vertex_multiplicity_image              : Image of the vertex multiplicity, this parameter also has a thickness/X/X0 dependency

##4. Calibration of beam energy and telescope angle resolution:

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

In the example script, an aluminium plate with a material step is used as the DUT. The cfg file we are using here employs
several measurement areas and a line in order to fit the calibration parameters. The measurement areas are marked on the
following image of the uncalibrated aluminium plate:

[picture](workspace/tbsw/validation/X0image_Boxes.png)

There is one measurement area in the 0.5mm thick aluminium (1), 6 measurement areas in the 1mm aluminium area (2,3 and 9-12)
the remaining measurement areas lie in the air region.

The results of the calibration, including pictures of the fits of the individual measurement areas, can be found in 
workspace/tmp-runs/X0-mc-air-test-reco-X0Calibration/ . The results of the calibration are also stored as a cfg file
in workspace/localDB/mc-air-test/x0cal_result.cfg.

##5. X/X0 imaging (calibrated):

In the last step an calibrated X/X0 image is produced, which used the cfg file from the previous calibration step. The
calibrated X/X0 image should look like this:

[picture](workspace/tbsw/validation/X0image.png)

The area of the aluminium plate is 8x8mm^2 and has a thickness of 1mm. The central part has a reduced thickness of 0.5mm. The aluminium plate is freely 
hanging in air. 

A X/X0 profile cut can be used for X/X0 measurements. In this example case the steps between thick and thin aluminium, as well as air is clearly visible.

[picture](workspace/tbsw/validation/X0profile.png)

The complete results can be found in workspace/root-files/X0-mc-alu-default-reco-Calibrated-X0image.root

## Reconstruction of test beam data:

If you are performing an analysis of real test beam data, the basic processing steps do not change. The example script for processing test beam datat is 
testbeam_x0.py. The most important changes are the following:

   * You do not simulate the air and aluminium runs. You have to actually record the runs in a mono-energetic particle beam using a high resolution tracking 
     telescope. 
   * For test beams at DESY using the EUDET/AIDA telescope: We recommend to use 40mm as a spacing between telescope sensors and a beam energy of 2GeV. The air run can 
     be short (~30min) while the aluminium runs should last around two hours. 
   * The raw data format from EUDAQ will be .raw files. The reading of .raw and .slcio files is slightly different. The required changes are already done in 
     the steering files 'x0-tb' instead of 'x0-sim' used for simulations. 
   * The file gear.xml in 'x0-tb' needs to be adjusted to your telescope geometry. In particular, the z positions of Mimosa 26 sensors along the beam axis must be carefully 
     measured. 
   * The track fitting w/o magnetic field needs a mean beam momentum as input. Change the momentum in the testbeam_x0.py to the beam momentum value 
     selected for your runs. 
   * Some changes in the x0 cfg files will be necessary as well. In particular, the measurement areas with known thicknessess for the X0 calibration will be different. 
     The selection of the measurement areas can be done by looking at the uncalibrated X0 image.
   * The uncalibrated X0 image is also useful to check that your scattering object is actually inside the beam. It is highly recommended to check this whenever you switch 
     or move your object. Do it during the data taking, not weeks after you are home ;)

In the following paragraph a guideline of how you have to change the testbeam_x0.py script and the x0.cfg file in order to be able to analyse your test beam data.

#0. First steps

Change the beam energy to the correct value, which was used during your beam test experiment

```
 # Nominal Beam energy
 beamenergy=2.0
``` 

The next step is to copy the default x0 steering files (workspace/steering-files/x0-tb) to workspace/steering-files/x0-mytb. Then the gear file (workspace/steering-files/x0-mytb/gear.xml) has to be edited to
match the telescope setup you used during you beam test experiment. Especially the M26 spacings, target position and mean thickness will have to be modified.

#1. Telescope calibration 

The telescope calibration works with a single run raw file 
```
# global path to raw files
rawfile_path='/work1/rawdata/tboct16/'

# raw file used during telescope calibration (best use data with scattering target)
# The calibration has to be done for every telescope setup, beam energy and m26 threshold settings
cali_run='run006973.raw'
rawfile_cali = rawfile_path + cali_run
```
The path to the run raw files and the name of the run raw file has to be modified to match your data. The calibration run should have 500k to 1mio tracks without any additional material in the telescope. The calibration procedure 

```
 # Calibrate the telescope
 params_cali = ( rawfile_cali, steerfiles_cali, gearfile, caltag)
 calibrate( params_cali )
```

produces a calibration tag directory, that is written to workspace/local-DB/caltag. Choose an appropriate name for the caltag,

```
# telescope calibration cal tag (typically named after telescope setup, beam energy etc.)
caltag='40mm-spacing-2GeV'
```

so that you always know, which telescope setup, beam energy etc is covered by this tag. The telescope calibration has to be done only once for each beam energy, telescope setup and M26 threshold setting. Once you have the caltag directory you can comment out/remove the telescope calibration procedure from the testbeam_x0.py script.


#2. Angle reconstruction

The angle reconstruction is performed for a list of runs. The list should contain all runs, which were recorded with the same beam energy and telescope settings and which are needed either for the x0 calibration or the x0 imaging. 
The target position, material and thickness can be different for the individual runs in the list, because for each run a seperate root file with the reconstructed scattering angles is created.


```
 RunList_reco = [
             'run006958.raw',
             'run006959.raw',
             'run006960.raw',
             'run006961.raw',
           ]
```

Of course, this list also has to be modified to match the runs you want to analyse. The angle reconstruction is started by the following code snippet

```
  # Angle reconstruction
  params_reco=[(x, steerfiles_reco, gearfile, caltag) for x in RawfileList_reco]
  print "The parameters for the reconstruction are: " 
  print params_reco

  count = multiprocessing.cpu_count()
  pool = multiprocessing.Pool(processes=count)
  pool.map(reconstruct, params_reco)
```

As you can see the angle reconstruction uses multiprocessing. In case you don't want to use all system ressources feel free to set the variable count to some number smaller than the number of processor cores of your maschine.

#3. X/X0 calibration

Now the x0 calibration is performed. Once again a list of runs has to be modified:

```
 RunList_x0cali = [
             'run006958.raw',
             'run006961.raw',
           ]
```

The list should contain runs with many different well known calibration targets. Typically different aluminium layers with various thicknesses and an air measurement without a scattering target. 

Before starting the x0 calibration, the workspace/steering-files/x0-mytb/x0.cfg file will have to be modified. You have to define several measurement areas with different X/X0 values and positions on the target plane.
For each measurement area a range of runnumbers and u and v position limits must be specified. The angle distribution found for these specifications, will then be used in the x0 calibration measurements. It is important
to use different X/X0 values during the calibration and a run without a scattering target should always be included, as the calibration of the angle reconstruction error depends on it. 

The x0 calibration process is started via:

```
  # start x0 calibration
  params_x0cali = ( x0caltag, RawfileList_x0cali, steerfiles_x0, caltag, deletetag)
  xx0calibration(params_x0cali)
```

The xx0calibration() functions performs the following steps:

The root files with reconstructed angle information, that were reconstructed in the previous step are merged

```
   # Merge the root trees in the root files directory
   tbsw.x0imaging.X0Calibration.merge_rootfile(filename=filename,RunList=RunList,caltag=caltag)
```
A uncalibrated image is generated from the merged root file. This image will be used to visualize the position of the measurement areas, which are employed during the x0 calibration.

```
   # Generate a uncalibrated X/X0 image
   tbsw.x0imaging.X0Calibration.x0imaging(filename=filename,caltag='',deletetag=deletetag,steerfiles=steerfiles,nametag='Uncalibrated')
```

Afterwards the x0 calibration itself is performed

```
   # Do a calibration of the angle resolution
   tbsw.x0imaging.X0Calibration.x0calibration(filename=filename,imagefilename=imagefilename,caltag=x0tag,steerfiles=steerfiles)
```

The x0 calibration creates a new caltag at localDB/x0tag. Just like the telescope calibration the x0 calibration step has to be done only once. Before proceeding with the x0imaging step one should check the quality of the calibration measurement. The fitted distributions can be found in workspace/tmp-runs/...-X0Calibration/X0calibration_results.root.

#5. X/X0 Images 

Now calibrated images can be generated. The imaging process is started via:

```
  # Generate a calibrated X/X0 image
  nametag='x0image-list'
  params_x0image = ( x0caltag, RawfileList_x0image, steerfiles_x0, caltag, deletetag,nametag)
  xx0image(params_x0image)
```

The runs used for the generation of the image are defined in a list:

```
 RunList_image1 = [
             'run006958.raw',
             'run006959.raw',
             'run006960.raw',
           ]
```

In this step it is important that exactly the same target material and target position are used. The Imaging process can now be repeated for every image you want to create. The image settings can be changed in the
config file workspace/steering-files/x0-mytb/x0.cfg.

## Simple example analysis

The following section will demonstrate how to modify the testbeam_x0.py script in order to do a x0 measurement. The beam energy is 4GeV and the telescope setup is assumed to be as follows:


			   ID0     ID1     ID2        target	   ID3     ID4     ID5
beam -->		| 4.5cm | 4.5cm |   5.0 cm  |   5.0 cm  | 4.5cm | 4.5cm |
		

The measurement is based on 4 data sets with different scattering materials. For each new target the telescope setup/spacing itself is unchanged, only the thickness and radiation of the scattering target is
different. In particular the target position is unchanged (at least on the scale of a few mm). The 4 datasets are:

target: air 				 run000001 - run000005

target: 0.5mm alu			 run000006 - run000010

target: 1.0mm alu  			 run000011 - run000015

target: unknown material   	 run000011 - run000015

Each step of the required modifications in the gearfiles, testbeam_x0.py and x0.cfg files to perform the telescope calibration, angle reconstruction, x/x0 calibration and x/x0 imaging will be documented in the
bullet points below. Every unnecessary modification such as target alignment etc will just not be mentioned. The bullet point list only contains necessary modifications:

   * Create a new steeringfiles folder from the default steering files folder:

	 ```
	  cp -r steering-files/x0-tb steering-files/x0-tb2018-xymeasurement
	 ```

   * Edit the gear file in the newly created directory (steering-files/x0-tb2018-xymeasurement/gear.xml):

		* Enter the correct positions of the M26 sensors (ID0 to ID5) and the target (ID11) in the gearfile.
     	* Set the thickness of the target (ID11) to 0.0001 and the radLength to 304000.0 (corresponding to the radiation length constant of air)

     This gear file can be used for all target thicknesses as the angle reconstruction is completely independent of the nominal target thickness documented in the gear file.

   * Edit the testbeam_x0.py script :

		* Edit the steerfing file path(lines 27-34):

			```
			# Steerfiles for the telescope calibration
			steerfiles_cali = 'steering-files/x0-tb2018-xymeasurement/'

			# Steerfiles for the angle reconstruction (can be the same directory as telescope calibration steerfiles)
			steerfiles_reco = 'steering-files/x0-tb2018-xymeasurement/'

			# Steerfiles for the x0calibration/x0imaging (can be the same directory as telescope calibration steerfiles)
			steerfiles_x0 = 'steering-files/x0-tb2018-xymeasurement/'
			```

		* Set the nominal beam energy (line 37):

			```
			# Nominal Beam energy
			beamenergy=4.0
			```

		* Edit the calibration tag of the telescope calibration, typically combination of testbeam, beam energy and telescope setup (line 41):

			```
			# cal tags
			# telescope calibration cal tag (typically named after telescope setup, beam energy etc.)
			caltag='45mm-spacing-4GeV'
			```

		* Edit the calibration tag of the x/x0 calibration, typically combination of materials used and the beamenergy (line 44):

			```
			# x0 calibration cal tag
			x0tag='air-1mmalu-4GeV'
			```

		* Edit the directory, where your raw files are stored (line 59):

			```
			# global path to raw files
			rawfile_path='/work1/rawdata/tb2018/'
			```

		* Choose one raw file, which will be used during the telescope calibration, any air run (run000001 - run000005) will do (line 63):

			```
			# raw file used during telescope calibration (best use data with scattering target)
			# The calibration has to be done for every telescope setup, beam energy and m26 threshold settings
			cali_run='run000001.raw'
			```

		* Select all rawfiles, for which the angle reconstruction has to be done, in this case every available raw file (run000001 - run000020, lines 72-77 in the original file):

			```
			# List of runs, which are used as input for the scattering angle reconstruction
			# The angle reconstruction step is essential and every run, that will be used later during the x0 calibration or x0 imaging steps, must be listed
			RunList_reco = [
						'run000001.raw',
						'run000002.raw',
						'run000003.raw',
						'run000004.raw',
						'run000005.raw',
						'run000006.raw',
						'run000007.raw',
						'run000008.raw',
						'run000009.raw',
						'run000010.raw',
						'run000011.raw',
						'run000012.raw',
						'run000013.raw',
						'run000014.raw',
						'run000015.raw',
						'run000016.raw',
						'run000017.raw',
						'run000018.raw',
						'run000019.raw',
						'run000020.raw',
					  ]
			```

		* Select all rawfiles, which will be used in the x/x0 calibration (run000001 - run000015, lines 84-87 in the original file):

			```
			# List of runs, which are input for the x0 calibration
			# Typically runs with various different materials and thicknesses have to be used to achieve a sensible calibration
			# The different measurement regions and other options have to be set in the x0.cfg file in the steer files directory
			RunList_x0cali = [
						'run000001.raw',
						'run000002.raw',
						'run000003.raw',
						'run000004.raw',
						'run000005.raw',
						'run000006.raw',
						'run000007.raw',
						'run000008.raw',
						'run000009.raw',
						'run000010.raw',
						'run000011.raw',
						'run000012.raw',
						'run000013.raw',
						'run000014.raw',
						'run000015.raw',
					  ]
			```

		* Select all rawfiles, which will be used in the x/x0 imaging process (run000016 - run000020, lines 93-97 in the original file):

			```
			# List of runs, which are input for the first x0 image
			# Use only runs, with exactly the same target material and positioning
			RunList_x0image = [
						'run000016.raw',
						'run000017.raw',
						'run000018.raw',
						'run000019.raw',
						'run000020.raw',
					  ]
			```
     These are all necessary changes in the testbeam_x0.py script.


   * Edit the x0.cfg file in the newly created directory (steering-files/x0-tb2018-xymeasurement/x0.cfg):

		* Set the beam energy start value to the nominal beam energy (line 17):

			```
			# Choose a momentumoffset start value [GeV]
			momentumoffset:4.0
			```

		* Set the side length of the x0 image (lines 43-45):

			```
			# u and v length of complete X0 image  in mm
			u_length : 30.0
			v_length : 30.0
			```

		* Define the position of the x0 image (lines 47-49):

			```
			# umin and vmax of the complete X0 image in mm
			umin : -15
			vmax : +15
			```
			The given values represent the upper left corner of the x/x0 image.

		* Define the pixel pitches of the x0 image (lines 51-53):

			```
			# Pixel sizes of the image in µm
			u_pixel_size : 400.0
			v_pixel_size : 400.0
			```

		* Define the measurement areas during the x/x0 calibration (from line 116 onward):

		    First measurement area (air, center of the beamspot), here only the thickness entry and the min/max runnumber entries have to be edited.
			The maximum number of scattering angles per distribution can be limited by setting the corresponding variable (maxanglenumber)

			```
			# ---Measurement area settings---
			# Measurement areas are single rectangular areas on the target u-v plane with known material properties 

			# Use area in fit: 0(no), 1(yes)
			MA1.exist:1

			# Center position in mm
			MA1.ucenter:0.0

			# Center position in mm
			MA1.vcenter:0.0

			# Side length in mm
			MA1.ulength:3.0

			# Side length in mm
			MA1.vlength:3.0

			# Thickness in mm
			MA1.thickness:0.0

			# Atomic number Z
			MA1.atomicnumber:13.0

			# Atomic mass A
			MA1.atomicmass:27.0

			# Density in g/cm³
			MA1.density:2.7

			# Smallest run number to be used (-1 to disable)
			MA1.minrunnumber:1

			# Largest run number to be used (-1 to disable)
			MA1.maxrunnumber:5

			# Limit number of angles in distribution to thsi value (-1 use all available angles)
			MA1.maxanglenumber:10000
			```

			Create 4 more air measurement areas at other positions, which positions you choose depends on the shape of the beamspot:

			```
			MA2.exist:1            	
			MA2.ucenter:3.25
			MA2.vcenter:2.0
			MA2.ulength:0.8
			MA2.vlength:1.0
			MA2.thickness:0.0
			MA2.atomicnumber:13.0
			MA2.atomicmass:27.0
			MA2.density:2.7
			MA2.minrunnumber:1		    
			MA2.maxrunnumber:5	
			MA2.maxanglenumber:10000	    

			MA3.exist:1            	
			MA3.ucenter:-3.25
			MA3.vcenter:2.0
			MA3.ulength:0.8
			MA3.vlength:1.0
			MA3.thickness:0.0
			MA3.atomicnumber:13.0
			MA3.atomicmass:27.0
			MA3.density:2.7
			MA3.minrunnumber:1		    
			MA3.maxrunnumber:5	
			MA3.maxanglenumber:10000

			MA4.exist:1            	
			MA4.ucenter:3.25
			MA4.vcenter:-2.0
			MA4.ulength:0.8
			MA4.vlength:1.0
			MA4.thickness:0.0
			MA4.atomicnumber:13.0
			MA4.atomicmass:27.0
			MA4.density:2.7
			MA4.minrunnumber:1		    
			MA4.maxrunnumber:5	
			MA4.maxanglenumber:10000

			MA5.exist:1            	
			MA5.ucenter:-3.25
			MA5.vcenter:-2.0
			MA5.ulength:0.8
			MA5.vlength:1.0
			MA5.thickness:0.0
			MA5.atomicnumber:13.0
			MA5.atomicmass:27.0
			MA5.density:2.7
			MA5.minrunnumber:1		    
			MA5.maxrunnumber:5	
			MA5.maxanglenumber:10000

			```
			Afterwards define 5 additional measurement areas (MA6-MA10) for 0.5mm of aluminium and another 5 measurement areas (MA11-MA15) for 1.0mm of aluminium.



			```
			MA6.exist:1
			MA6.ucenter:0.0
			MA6.vcenter:0.0
			MA6.ulength:3.0
			MA6.vlength:3.0
			MA6.thickness:0.5
			MA6.atomicnumber:13.0
			MA6.atomicmass:27.0
			MA6.density:2.7
			MA6.minrunnumber:6
			MA6.maxrunnumber:10
			MA6.maxanglenumber:10000

			MA7.exist:1            	
			MA7.ucenter:3.25
			MA7.vcenter:2.0
			MA7.ulength:0.8
			MA7.vlength:1.0
			MA7.thickness:0.5
			MA7.atomicnumber:13.0
			MA7.atomicmass:27.0
			MA7.density:2.7
			MA7.minrunnumber:6		    
			MA7.maxrunnumber:10	
			MA7.maxanglenumber:10000	    

			MA8.exist:1            	
			MA8.ucenter:-3.25
			MA8.vcenter:2.0
			MA8.ulength:0.8
			MA8.vlength:1.0
			MA8.thickness:0.5
			MA8.atomicnumber:13.0
			MA8.atomicmass:27.0
			MA8.density:2.7
			MA8.minrunnumber:6		    
			MA8.maxrunnumber:10	
			MA8.maxanglenumber:10000

			MA9.exist:1            	
			MA9.ucenter:3.25
			MA9.vcenter:-2.0
			MA9.ulength:0.8
			MA9.vlength:1.0
			MA9.thickness:0.5
			MA9.atomicnumber:13.0
			MA9.atomicmass:27.0
			MA9.density:2.7
			MA9.minrunnumber:6		    
			MA9.maxrunnumber:10	
			MA9.maxanglenumber:10000

			MA10.exist:1            	
			MA10.ucenter:-3.25
			MA10.vcenter:-2.0
			MA10.ulength:0.8
			MA10.vlength:1.0
			MA10.thickness:0.5
			MA10.atomicnumber:13.0
			MA10.atomicmass:27.0
			MA10.density:2.7
			MA10.minrunnumber:6		    
			MA10.maxrunnumber:10	
			MA10.maxanglenumber:10000

			MA11.exist:1
			MA11.ucenter:0.0
			MA11.vcenter:0.0
			MA11.ulength:3.0
			MA11.vlength:3.0
			MA11.thickness:1.0
			MA11.atomicnumber:13.0
			MA11.atomicmass:27.0
			MA11.density:2.7
			MA11.minrunnumber:11
			MA11.maxrunnumber:15
			MA11.maxanglenumber:10000

			MA12.exist:1            	
			MA12.ucenter:3.25
			MA12.vcenter:2.0
			MA12.ulength:0.8
			MA12.vlength:1.0
			MA12.thickness:1.0
			MA12.atomicnumber:13.0
			MA12.atomicmass:27.0
			MA12.density:2.7
			MA12.minrunnumber:11		    
			MA12.maxrunnumber:15	
			MA12.maxanglenumber:10000	    

			MA13.exist:1            	
			MA13.ucenter:-3.25
			MA13.vcenter:2.0
			MA13.ulength:0.8
			MA13.vlength:1.0
			MA13.thickness:1.0
			MA13.atomicnumber:13.0
			MA13.atomicmass:27.0
			MA13.density:2.7
			MA13.minrunnumber:11		    
			MA13.maxrunnumber:15	
			MA13.maxanglenumber:10000

			MA14.exist:1            	
			MA14.ucenter:3.25
			MA14.vcenter:-2.0
			MA14.ulength:0.8
			MA14.vlength:1.0
			MA14.thickness:1.0
			MA14.atomicnumber:13.0
			MA14.atomicmass:27.0
			MA14.density:2.7
			MA14.minrunnumber:11		    
			MA14.maxrunnumber:15	
			MA14.maxanglenumber:10000

			MA15.exist:1            	
			MA15.ucenter:-3.25
			MA15.vcenter:-2.0
			MA15.ulength:0.8
			MA15.vlength:1.0
			MA15.thickness:1.0
			MA15.atomicnumber:13.0
			MA15.atomicmass:27.0
			MA15.density:2.7
			MA15.minrunnumber:11		    
			MA15.maxrunnumber:15	
			MA15.maxanglenumber:10000

			```

		* Remove the definition of the measurement area lines (from line 216 and 258 in the original file):

			```
			# Use line in fit: 0(no), 1(yes)
			line1.exist:0
			```

			```
			line2.exist:0
			```


     These are all necessary changes in the x0.cfg file.

Afterwards you can start the analysis via:

```
$ source init_tbsw.sh
$ python testbeam_x0.py
```




The technique and algorithm behind radiation length imaging was shown VCI 2016. The link to the proceedings paper is http://www.sciencedirect.com/science/article/pii/S0168900216306519.
Please you this reference for citing the method. 



Ulf Stolzenberg
Benjamin Schwenker

Goettingen 2018

ulf.stolzenberg@phys.uni-goettingen.de 
benjamin.schwenker@phys.uni-goettingen.de


 

















