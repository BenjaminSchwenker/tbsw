"""
This is an example script to demonstrate how TBSW can be used to analyze test beam 
data using Python scripts.

The script below calibrates a single run. Afterwards it reconstructs all runs from a
runlist using the calibration data. In this example, one common set of steerfiles is
applied to all runs. The second part of the data processing uses Pythons 
multiprocessing.Pool to distribute runs to multiple cores and speed up processing time. 

Author: Ulf Stolzenberg <ulf.stolzenberg@phys.uni-goettingen.de>  
"""

from tbsw import *
import multiprocessing

# Path to steering files 
steerfiles = 'steering-files/x0-tb/'

# Create a list of parameters
MaxOccupancyList = (0.01,0.005,0.001)

# Create a run list for the calibration
RunList_cal = [
		    '/work1/rawdata/DESY_Oktober16/2GeV_air/run006972.raw',
          ]

# Create a run list for the reconstruction
RunList_reco = [
		    '/work1/rawdata/DESY_Oktober16/2GeV_air/run006972.raw',
		    '/work1/rawdata/DESY_Oktober16/2GeV_air/run006973.raw',
          ]

parameters_cal=[(x,y) for x in RunList_cal for y in MaxOccupancyList ]
print "The parameters for the reconstruction are: " 
print parameters_cal

parameters_reco=[(x,y) for x in RunList_reco for y in MaxOccupancyList ]
print "The parameters for the reconstruction are: " 
print parameters_reco


# Defines the sequence of calibration steps. 
# XML steer files are taken from steerfiles. 
calpath = [ 
           'hotpixelkiller.xml' ,              
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
           'correlator-iteration-2.xml' ,
           'kalmanalign-iteration-3.xml',
           'kalmanalign-iteration-4.xml',
           'kalmanalign-iteration-4.xml',
           'kalmanalign-iteration-4.xml',
           'telescope-dqm-iteration-2.xml',
         ]

def calibrate(pars):  

  # Tag for calibration data 
  caltag = 'air-MaxOccupancy'+str(int(10000*pars[1]))+'permille'
  
  # Base name for temporary folder created in tmp-runs/ 
  name = os.path.splitext(os.path.basename(pars[0]))[0] + '-' + caltag 

  # Calibrate the telescope using the rawfile. Creates a folder caltag 
  # containing all calibrations. 
  CalObj = Calibration(steerfiles=steerfiles, name=name + '-cal') 

  # The following lines show how to change parameters in copied 
  # XML steer files managed by CalObj 
  xmlfile = CalObj.get_filename('hotpixelkiller.xml')
  override_xml(xmlfile=xmlfile, xpath="./processor[@name='M26HotPixelKiller']/parameter[@name='MaxOccupancy']", value=pars[1])
  
  CalObj.calibrate(path=calpath,ifile=pars[0],caltag=caltag) 
  
  return None

if __name__ == '__main__':
  count = 2 #multiprocessing.cpu_count()
  pool = multiprocessing.Pool(processes=count)
  pool.map(calibrate, parameters_cal)

def reconstruct(pars):

  # Tag for calibration data 
  caltag = 'air-MaxOccupancy'+str(int(10000*pars[1]))+'permille'
  
  # Name for temporary folder created in tmp-runs/ 
  name = os.path.splitext(os.path.basename(pars[0]))[0] 

  # Reconsruct the rawfile using caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=name + '-reco-caltag-' + caltag )
  RecObj.reconstruct(path=['reco.xml'],ifile=pars[0],caltag=caltag) 

  return None

if __name__ == '__main__':
  count = 2       #multiprocessing.cpu_count()
  pool = multiprocessing.Pool(processes=count)
  pool.map(reconstruct, parameters_reco)
  

  


