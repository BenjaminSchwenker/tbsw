"""
This is an example script to demonstrate how TBSW can be used to analyze test beam 
data using Python scripts.

The script below simulates a test beam experiment where charged tracks cross a misaligned
pixel telescope containing six Mimosa 26 detector planes. In order to have a realistic scenario
two runs are simulated. One with only air between the two telescope arms and one with an aluminium
plate, centered in the telescope. After the simulations, the simulated raw data without DUT is
calibrated. The calibration data is used in the next step to reconstuct the aluminium run. Then
an X/X0 of the aluminium plate is generated and a calibration of the angle resolution
of the telescope is performed

Author: Ulf Stolzenberg <ulf.stolzenberg@phys.uni-goettingen.de>  
"""

from tbsw import *

# Path to steering files 
# Steeringfiles are xml files and define details of the simulation like how many events are produced
# or how M26 sensors are digitized. XML parameters can be adjusted using any test editor
steerfiles = 'steering-files/x0-sim/'
# Tag for calibration data 
caltag = 'default'
# File name for raw data 
runname_air='mc-air' 
runname_alu='mc-alu'

rawfile_air = runname_air+'.slcio'
rawfile_alu = runname_alu+'.slcio'

# Defines the sequence of calibration steps. 
# XML steer files are taken from steerfiles. 
calpath = [ 
           'hotpixelkiller.xml' ,              
           'cluster-calibration-mc.xml',     # creates clusterDB using MC truth information, but will not be used for reconstruction
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

# Function which is required to change the global settings of a steer file
def override_xmlfileglobal(xmlfile=None, paramname=None, value=None):
  """
  Overrides proecessor parameters in Marlin XML steering file. 
    :@xmlfile:    Marlin steering file to be overwritten  
    :@procname:   name of global parameter   
    :@value:      value of the parameter 
 
    :author: benjamin.schwenker@phys.uni-goettinge.de  
    """   
  tree = xml.etree.ElementTree.parse(xmlfile)
  root = tree.getroot() 
  
  for glob in root.findall('global'):
      for param in glob.findall('parameter'):
        if param.get('name') ==  paramname:
          param.set('value', str(value)) 
  
  tree.write(xmlfile)  

# Base name for temporary folder created in tmp-runs/ 
name_air = os.path.splitext(os.path.basename(rawfile_air))[0] + '-' + caltag  

# Simulate a rawfile from a test beam experiment
# SimObj creates folder tmp-runs/name-sim/ and populates it with 
# copies of steering files. After processing, the folder contains
# also logfiles. 
# Air data for telescope calibration
SimObj_air = Simulation(steerfiles=steerfiles, name=name_air + '-sim' )

# The following lines show how to change parameters in copied 
# XML steer files managed by CalObj 
xmlfile = SimObj_air.get_filename('simulation.xml')
override_xmlfileglobal(xmlfile=xmlfile, paramname='GearXMLFile', value='gear-air.xml') 
override_xmlfileglobal(xmlfile=xmlfile, paramname='MaxRecordNumber', value=200000) 

SimObj_air.simulate(path=['simulation.xml'], ofile=rawfile_air, caltag=None)  


# Base name for temporary folder created in tmp-runs/ 
name_alu = os.path.splitext(os.path.basename(rawfile_alu))[0] + '-' + caltag  

# Simulate a rawfile from a test beam experiment
# SimObj creates folder tmp-runs/name-sim/ and populates it with 
# copies of steering files. After processing, the folder contains
# also logfiles. 
# Air data for telescope calibration
SimObj_alu = Simulation(steerfiles=steerfiles, name=name_alu + '-sim' )
SimObj_alu.simulate(path=['simulation.xml'], ofile=rawfile_alu, caltag=None)  
   
# Calibrate the telescope using the air rawfile. Creates a folder caltag 
# containing all calibrations. 
CalObj = Calibration(steerfiles=steerfiles, name=name_air + '-cal') 
CalObj.calibrate(path=calpath,ifile=rawfile_air,caltag=caltag)  
   
# Reconsruct the alu rawfile using caltag. Resulting root files are 
# written to folder root-files/
RecObj = Reconstruction(steerfiles=steerfiles, name=name_alu + '-reco' )
RecObj.reconstruct(path=['reco.xml'],ifile=rawfile_alu,caltag=caltag) 

imagefilename='root-files/cal-tag-default/uncalibrated_x0image/X0-completeimage.root'
deletetag=1

# Function which starts the imaging script
def x0imaging(filename,caltag,deletetag):

  flags='./root-scripts/x0imaging/GenerateImage.py -i '+filename+' -c '+caltag+' -d '+`deletetag`
  subprocess.call(flags, shell=True)

  return None

# Function which starts the x0 calibration script
def x0calibration(filename,imagefilename,caltag):

  flags='./root-scripts/x0imaging/X0Calibration.py -i '+filename+' -m '+imagefilename+' -c '+caltag
  subprocess.call(flags, shell=True)

  return None


# Create a new folder at workspace/root-files/cal-tag-default/ and store the X/X0 results there
cwdir = os.getcwd()
rootpath=cwdir+'/root-files/cal-tag-'+caltag+'/' 
if not os.path.isdir(rootpath):
   os.mkdir(rootpath)

# Copy root file
rootfile=cwdir+'/root-files/X0-'+name_alu+'-'+'reco.root'
shutil.copy(rootfile, rootpath)  

filename='root-files/cal-tag-default/X0-'+name_alu+'-'+'reco.root'
  
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


