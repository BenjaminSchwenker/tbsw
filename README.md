
# Introduction to the test beam analysis software (TBSW) 

The test beam software framework (tbsw) deals with the reconstruction and analysis of test beam data obtained with the 
EUDET/AIDA reference telescope. Test beams are all about learning about a device under test, or DUT, installed in the 
centre of the reference telescope. The idea is to probe the device under test with charged tracks crossing the DUT at 
known times and precisely estimated kinematic parameters. 

The tbsw framework is designed to run on local machines (Desktop / Laptop) to enable detector groups doing their own data 
analysis right during the beam test (or shortly afterwards). Some effort was invested to make this as smooth and simple
as possible. The installation of the framework is covered in INSTALL. 


For the impatient reader: The text below is a step by step instruction how to run the code on raw data.


After installation, the root directory of tbsw contains a 'workspace' folder. The first step to analyse data from a new beam test
is to make a working copy of this folder. For example: 

```
$ cp -r workspace ~/workspace-my-beam-test
```

The idea of the folder 'workspace-my-beam-test' is to collect all steering files and analysis results from a specific beam test. 
Please replace the 'my-beam-test' part with something more specific like 'desy_october_2016' or similar. For simplicity, we will
call this folder just workspace in the following. 

It is a good idea to copy the workspace to some place outside of the tbsw installation. In case you want to reprocess your data with 
a new (different) version of tbsw, you typically just have to make sure the environment variables in 'init_tbsw.sh' point to the
desired installation.

Raw data files (.raw) should not be stored directly in your workspace folder. Again, this simplifies switching version of tbsw. It is
advised to create a seperate folder and to collect all raw data files there.

For example: 

```
$ mkdir path-to-big-disk/eudaq-raw/my-beam-test 
```

Raw data from a one week beam test is typically ~1TB and fits on most Desktop PCs or laptops.  

Next, go to your new workspace folder and setup all environment variables. 

```
$ . init_tbsw.sh
```

A test beam with the EUDET/AIDA telescope results in a number of .raw files stored at 'eudaq-raw/my-beam-test'. Before we can start
with the reconstruction and analysis of data, we need to prepare a set of steering files. Steering files are simple text files and 
can be edited just using the text editor of your choice. Typically, a test beam period breaks down to into groups of runs where one 
DUT was studied under 'unchanged' conditions. The phrase unchanged conditions refers to all conditions relevant for data reconstruction
and include: 

- for reference telescope: distances between detector planes, beam energy,  TLU configuration, m26 configuration files (->thresholds).
- for the DUT: you name it, since it is your dut ;) You will probably want to put the name of the DUT and some important configurations. 

It is advised to assign a running number (GeoID) for each of these run groups and to write up all needed information together for the 
telescope and the DUT into your logbook.  

In order to provide an introductory example, the folder 'workspace/steering-files/x0-sim' contains a full set of steering files needed 
to simulate a beam run with a misaligned telescope and to calibrate and reconstruct the run afterwards. The folder has the content: 

- gear-misaligned.xml                : xml file to describe the misaligned telescope geometry (the one used for simulation)
- gear.xml                           : xml file to describe the nominal telescope geometry  
- simulation.xml                     : Marlin steer file, produces a lcio file with simulated digits and truth tracks 
- cluster-calibration-mc.xml         : Marlin steer file, produces ClusterDB calibration file based in MC truth tracks 
- hotpixelkiller.xml                 : Marlin steer file, produces NoiseDB calibration file(s) to mask hot pixels  
- correlator.xml                     : Marlin steer file, produces alignmentDB calibration file with pre-aligned XY positions 
- kalmanalign-iteration-1.xml        : Marlin steer file, produces alignmentDB calibration file using tracks (loose cuts)   
- align-config-iteration-1.cfg       : text file used to define alignment parameters 
- kalmanalign-iteration-2.xml        : Marlin steer file, produces alignmentDB calibration file using tracks (tight cuts) 
- align-config-iteration-2.cfg       : text file used to define alignment parameters 
- telescope-dqm.xml                  : Marlin steer file, produces DQM histos to monitor track fit quaility after alignment 
- reco.xml                           : Marlin steer file, runs the full reconstruction using calibration files  


All these files area readable by a text editor and heavily commented. The Python script tbsw_example.py demonstrates 
the basic use case. The source code of the tbsw_example.py is: 

```
#!python

"""
This is an example script to demonstrate how TBSW can be used to analyze test beam 
data using Python scripts.

The script below simulates a test beam experiment where charged tracks cross a misaligned
pixel telescope containing six Mimosa 26 detector planes. Afterwards, the simulated 
raw data is calibrated and reconstucted. Final results are prepared in form of root files 
that get copied to the folder root-files/.

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *
import getopt, sys

# Path to steering files 
steerfiles = 'steering-files/x0-sim/'
# Tag for calibration data 
caltag = 'default'
# File name for raw data  
rawfile = 'mc.slcio'

# Defines the sequence of calibration steps. 
# XML steer files are taken from steerfiles. 
calpath = [ 
           'hotpixelkiller.xml' ,              
           'cluster-calibration-mc.xml',
           'correlator.xml' ,
           'kalmanalign-iteration-1.xml',
           'kalmanalign-iteration-2.xml',
           'kalmanalign-iteration-2.xml',
           'kalmanalign-iteration-2.xml',
           'telescope-dqm.xml',
           #'cluster-calibration-tb.xml',
         ]

# Base name for temporary folder created in tmp-runs/ 
name = os.path.splitext(os.path.basename(rawfile))[0] + '-' + caltag  

# Simulate a rawfile from a test beam experiment
# SimObj creates folder tmp-runs/name-sim/ and populates it with 
# copies of steering files. After processing, the folder contains
# also logfiles. 
SimObj = Simulation(steerfiles=steerfiles, name=name + '-sim' )
SimObj.simulate(path=['simulation.xml'], ofile=rawfile, caltag=None)  
   
# Calibrate the telescope using the rawfile. Creates a folder caltag 
# containing all calibrations. 
CalObj = Calibration(steerfiles=steerfiles, name=name + '-cal') 

# The following lines show how to change parameters in copied 
# XML steer files managed by CalObj 
xmlfile = CalObj.get_filename('cluster-calibration-mc.xml')
override_xmlfile(xmlfile=xmlfile, procname='M26ClusterDBCreator', paramname='SoftScale', value=8) 

CalObj.calibrate(path=calpath,ifile=rawfile,caltag=caltag)  
   
# Reconsruct the rawfile using caltag. Resulting root files are 
# written to folder root-files/
RecObj = Reconstruction(steerfiles=steerfiles, name=name + '-reco' )
RecObj.reconstruct(path=['reco.xml'],ifile=rawfile,caltag=caltag) 
  

```

You can run the script by typing:  

```
$ python tbsw_example.py
```

In case the script is working, your console should display the following: 

```
[INFO] Starting to simulate mc.slcio ...
[INFO] Marlin simulation.xml is done
[INFO] Starting to calibrate file mc.slcio ...
[INFO] Marlin hotpixelkiller.xml is done
[INFO] Marlin cluster-calibration-mc.xml is done
[INFO] Marlin correlator.xml is done
[INFO] Marlin kalmanalign-iteration-1.xml is done
[INFO] Marlin kalmanalign-iteration-2.xml is done
[INFO] Marlin kalmanalign-iteration-2.xml is done
[INFO] Marlin kalmanalign-iteration-2.xml is done
[INFO] Marlin telescope-dqm.xml is done
[INFO] Created new caltag  default
[INFO] Done processing file mc.slcio
[INFO] Starting to reconstruct file mc.slcio ...
[INFO] Marlin reco.xml is done
[INFO] Done processing file mc.slcio
```



Have fun with test beams ;)  

benjamin.schwenker@phys.uni-goettingen.de
