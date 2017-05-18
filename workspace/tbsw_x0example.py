"""
This is an example script to demonstrate how TBSW can be used to analyze test beam 
data using Python scripts.

The script below simulates a test beam experiment where charged tracks cross a misaligned
pixel telescope containing six Mimosa 26 detector planes and a 05mm aluminium plate,
centered in the telescope. Afterwards, the simulated raw data is calibrated and
reconstucted. The last step is the creation of an X/X0 of the aluminium plate. Additionally
a calibration of the angle resolution of the telescope is performed.

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
SimObj = Simulation(steerfiles=steerfiles, name=name + '-sim' )
SimObj.simulate(path=['simulation.xml'], ofile=rawfile, caltag=None)  
   
# Calibrate the telescope using the rawfile. Creates a folder caltag 
# containing all calibrations. 
CalObj = Calibration(steerfiles=steerfiles, name=name + '-cal') 
CalObj.calibrate(path=calpath,ifile=rawfile,caltag=caltag)  
   
# Reconsruct the rawfile using caltag. Resulting root files are 
# written to folder root-files/
RecObj = Reconstruction(steerfiles=steerfiles, name=name + '-reco' )
RecObj.reconstruct(path=['reco.xml'],ifile=rawfile,caltag=caltag) 

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
x0imaging(filename,imagefilename,caltag)

# Generate a calibrated X/X0 image
x0imaging(filename,caltag,deletetag)

# Rename the directory in which the imaging results are stored
if os.path.isdir(rootpath+'x0image'):
   shutil.move(rootpath+'x0image', rootpath+'calibrated_x0image')  


