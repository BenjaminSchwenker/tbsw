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
  


