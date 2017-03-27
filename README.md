Welcome to the world of test beam data analysis 


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

> cp -r workspace ~/workspace-my-beam-test

The idea of the folder 'workspace-my-beam-test' is to collect all steering files and analysis results from a specific beam test. 
Please replace the 'my-beam-test' part with something more specific like 'desy_october_2016' or similar. For simplicity, we will
call this folder just workspace in the following. 

It is a good idea to copy the workspace to some place outside of the tbsw installation. In case you want to reprocess your data with 
a new (different) version of tbsw, you typically just have to make sure the environment variables in 'init_tbsw.sh' point to the
desired installation.

Raw data files (.raw) should not be stored directly in your workspace folder. Again, this simplifies switching version of tbsw. It is
advised to create a seperate folder and to collect all raw data files there.

For example: 

>mkdir path-to-big-disk/eudaq-raw/my-beam-test 

Raw data from a one week beam test is typically ~1TB and fits on most Desktop PCs or laptops.  

Next, go to your new workspace folder and setup all environment variables. 

> . init_tbsw.sh

A test beam with the EUDET/AIDA telescope results in a number of .raw files stored at 'eudaq-raw/my-beam-test'. Before we can start
with the reconstruction and analysis of data, we need to prepare a set of steering files. Steering files are simple text files and 
can be edited just using the text editor of your choice. Typically, a test beam period breaks down to into groups of runs where one 
DUT was studied under 'unchanged' conditions. The phrase unchanged conditions refers to all conditions relevant for data reconstruction
and include: 

- for reference telescope: distances between detector planes, beam energy,  TLU configuration, m26 configuration files (->thresholds).
- for the DUT: you name it, since it is your dut ;) You will probably want to put the name of the DUT and some important configurations. 

It is advised to assign a running number (GeoID) for each of these run groups and to write up all needed information together for the 
telescope and the DUT into your logbook.  

How do steerign files look like and how to manage them? You should put your steering files in the folder 'workspace/steering-files' and 
create a subfolder for each GeoID containing all needed steering files. The folder 'workspace/steering-files/x0-tb' contains a full set 
of steering files to analyse a X0 test beam. The folder has the content: 

1) gear.xml                           : xml file to describe the nominal telescope geometry (before alignment corrections) 
2) hotpixelkiller.xml                 : Marlin steering file, part of calibration 
3) correlator.xml                     : Marlin steering file, part of calibration 
4) kalmanalign-iteration-1.xml        : Marlin steering file, part of calibration 
5) align-config-iteration-1.cfg       : text file used to define alignment parameters 
6) kalmanalign-iteration-2.xml        : Marlin steering file, part of calibration 
7) align-config-iteration-2.cfg       : text file used to define alignment parameters 
8) telescope-dqm.xml                  : Marlin steering file, part of calibration 
9) reco.xml                           : Marlin steering file, contains all parameters for reconstruction 
10)simulation.xml                     : Marlin steering file to produce lcio file with simulated digits
11)gear-misaligned.xml                : xml file to describe the misaligned telescope geometry (the one used for simulation)  

All these files area readable by a text editor and are heavily commented. The content of the x0-tb steer files are explained in more 
detail in README_X0Image. 

Once we have the a set of steer files for run, we can proceed to run the calibration script from the workspace:  

> ./calibrate.py -i path-to-big-disk/eudaq-raw/my-beam-test/run-XY.raw -x steering-files/my-x0-tb -c geoid

The -i option defines which run to use for calibration. The -x options defines the folder of steering files to be used. The option -c 
sets a calibration tag, a common name for all calibration data files that are created. For example, you can use the GeoID as a tag but any 
string can be used. The calibrate.py script makes starts a series of process to compute the calibration data. A typical execution example 
looks like this : 

> ./calibrate.py -i path-to-big-disk/eudaq-raw/my-beam-test/run-XY.raw -x steering-files/my-x0-tb -c geoid
[Print] HotPixelKiller done ...
[Print] Correlator done ...
[Print] Alignment first iteration done ...
[Print] Alignment second iteration done ...
[Print] Alignment third iteration done ...
[Print] Alignment fourth iteration done ...
[Print] TelescopeDQM done.
[Print] All calibration files written to ~/workspace_my_test_beam/cal-files/geoid

You can see that all needed calibration files aregrouped into a folder 'cal-files/geoid'. Also note that a folder 'tmp-runs/run-XY-geoid' gets 
created containing log files and DQM root files with many interesting plots to check the success of the calibration. The log files contain very
important information for trouble shooting and debugging. 

After the calibration is done, we can now run the full reconstruction on all runs taken under identical conditions (belonging to the same calibration 
tag). The ranges of run numbers taken under identical conditions and belonging to the same tag must(!) be written down in your logbook. 

For our example, we will just reconstruct the run also used to create the tag: 

> ./reco.py -i path-to-big-disk/eudaq-raw/my-beam-test/run-XY.raw -x steering-files/my-x0-tb/reco.xml -c geoid

The reco.py script will loop once over all events in the file input file and apply the reconstruction processors defined in the xml file 'reco.xml' 
to each event. The string '-c geoid' means that calibration results stored in 'cal-files/geoid' will be used. 

The reco.py script will populate the folder 'root-files' in your workspace. The root files can be browsed directly using ROOT TBrowser.  
Some useful root macros can be found in the 'root-scripts' folder in the workspace. For X0 studies, have a look in the 'root-scripts/x0imaging' 
sub folder. 


Have fun with test beams ;)  

benjamin.schwenker@phys.uni-goettingen.de
