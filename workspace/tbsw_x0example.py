"""
This is an example script to demonstrate how TBSW can be used to create an X0 image from a simulated
test beam experiment.

The script below simulates a test beam experiment where charged tracks from a monoenergetic beam 
cross a a misaligned pixel telescope containing six Mimosa 26 detector planes. Two data sets are 
simulated. A first 'air' run is done with no additional scatterer between the telescope arms. The 
'air' run is used to calibrate the telescope. In a second 'aluminium' run, a aluminium plate with 
a well known thickness profile is inserted in between the telescope arms. This second run is used 
to compute a X0 image from the reconstructed scattering angles. The known comparison between the 
reconstructed X0 image and the a priori known image is used to calibrate the beam energy and the 
angular resolution of the telescope. This second step completes the calibration of the telescope
for X0 imaging. 

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
rawfile_air = 'mc-air.slcio'
rawfile_alu = 'mc-alu.slcio'

# Defines the sequence of calibration steps. 
# XML steer files are taken from steerfiles. 
calpath = [ 
           'hotpixelkiller.xml' ,              
           'cluster-calibration-mc.xml',     
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

  




# Simulate a rawfile from a test beam experiment
# SimObj creates folder tmp-runs/name/ and populates it with 
# copies of steering files. After processing, the folder contains
# also logfiles. 

# Simulate Air data for telescope calibration 
SimObj_air = Simulation(steerfiles=steerfiles, name='mc-air-sim' )

# The following lines show how to change XML parameters before execution
xmlfile = SimObj_air.get_filename('simulation.xml')
override_xmlfileglobal(xmlfile=xmlfile, paramname='GearXMLFile', value='gear-air.xml') 
override_xmlfileglobal(xmlfile=xmlfile, paramname='MaxRecordNumber', value=200000) 

# Now start the simulation of air run
SimObj_air.simulate(path=['simulation.xml'], ofile=rawfile_air)  

# Al data for calibration of X0 image
SimObj_alu = Simulation(steerfiles=steerfiles, name='mc-alu-sim' )
xmlfile = SimObj_alu.get_filename('simulation.xml')
override_xmlfileglobal(xmlfile=xmlfile, paramname='GearXMLFile', value='gear.xml') 
override_xmlfileglobal(xmlfile=xmlfile, paramname='MaxRecordNumber', value=200000) 
SimObj_alu.simulate(path=['simulation.xml'], ofile=rawfile_alu)  

   
# Calibrate the telescope using the air rawfile. Creates a folder caltag 
# containing all calibrations 
CalObj = Calibration(steerfiles=steerfiles, name= 'mc-air-cal') 
CalObj.calibrate(path=calpath,ifile=rawfile_air,caltag=caltag)  
   
# Reconsruct the alu rawfile using caltag. Resulting root files are 
# written to folder root-files/
RecObj = Reconstruction(steerfiles=steerfiles, name='mc_alu-reco' )
RecObj.reconstruct(path=['reco.xml'],ifile=rawfile_alu,caltag=caltag) 

imagefilename='root-files/X0-mc-alu-default-reco-Uncalibrated-X0image.root'
imagecfgfilename='steering-files/x0-sim/image.cfg'
calibrationcfgfilename='steering-files/x0-sim/x0calibration.cfg'
deletetag=1

# Function which starts the imaging script
def x0imaging(filename,caltag,deletetag):

  flags='./tbsw_tools/x0imaging/GenerateImage.py -i '+filename+' -f '+imagecfgfilename+' -c '+caltag+' -d '+`deletetag`
  print('Starting X0 imaging')
  print(flags)
  subprocess.call(flags, shell=True)

  return None

# Function which starts the x0 calibration script
def x0calibration(filename,imagefilename,caltag):

  flags='./tbsw_tools/x0imaging/X0Calibration.py -i '+filename+' -f '+calibrationcfgfilename+' -m '+imagefilename+' -c '+caltag
  print('Starting X0 calibration')
  print(flags)
  subprocess.call(flags, shell=True)

  return None

# Some strings which will be needed in file operations later
# current directory
cwdir = os.getcwd()

# Base filename of the X0 root file
basefilename='X0-'+name_alu+'-'+'reco'

# Total path of X0 root file
filename='root-files/'+basefilename+'.root'

# Total path if the different kinds of X0 image files
imagefile=cwdir+'/root-files/'+basefilename+'-X0image.root'
uncalibratedimagefile=cwdir+'/root-files/'+basefilename+'-Uncalibrated-X0image.root'
calibratedimagefile=cwdir+'/root-files/'+basefilename+'-Calibrated-X0image.root'

# Total path to the different temporary work directories,
# in which the images are generated
tmpdir = cwdir+'/tmp-runs/'+basefilename+'-X0image'
uncaltmpdir = cwdir+'/tmp-runs/'+basefilename+'-UncalibratedX0image'
caltmpdir = cwdir+'/tmp-runs/'+basefilename+'-CalibratedX0image'

# Generate a uncalibrated X/X0 image
x0imaging(filename,'',deletetag)

# Rename the image file
if os.path.isfile(uncalibratedimagefile):
   os.remove(uncalibratedimagefile) 
if os.path.isfile(imagefile):
   shutil.move(imagefile,uncalibratedimagefile)

# Rename the directory in tmp-runs, in which the image was generated
if os.path.isdir(uncaltmpdir):
   shutil.rmtree(uncaltmpdir) 
if os.path.isdir(tmpdir):
   shutil.copytree(tmpdir, uncaltmpdir)

# Do a calibration of the angle resolution
x0calibration(filename,imagefilename,caltag)

# Generate a calibrated X/X0 image
x0imaging(filename,caltag,deletetag)

# Rename the image file
if os.path.isfile(calibratedimagefile):
   os.remove(calibratedimagefile) 
if os.path.isfile(imagefile):
   shutil.move(imagefile,calibratedimagefile)

# Rename the directory in tmp-runs, in which the image was generated
if os.path.isdir(caltmpdir):
   shutil.rmtree(caltmpdir) 
if os.path.isdir(tmpdir):
   shutil.copytree(tmpdir, caltmpdir)



