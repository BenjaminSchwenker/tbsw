"""
This is an example script to demonstrate how TBSW can be used to analyze test beam 
data using Python scripts.

The script below calibrates and reconstructs all runs from a runlist. In this example, 
one common set of steerfiles is applied to all runs. The calibration of the telescope
is repeated for each run using the calibration steps given in the calpath. The data 
processing uses Pythons multiprocessing.Pool to distribute runs to multiple cores and 
speed up processing time. 


Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""


from tbsw import *
import multiprocessing
      
# Path to steering files 
steerfiles = 'steering-files/x0-tb/'  

# Sequence of calibration steps 
calpath = [ 
           'hotpixelkiller.xml' ,              
           'correlator.xml' ,
           'kalmanalign-iteration-1.xml',
           'kalmanalign-iteration-2.xml',
           'kalmanalign-iteration-2.xml',
           'kalmanalign-iteration-2.xml',
           'telescope-dqm.xml',
         ]


# Create a list of runs to be processed
runlist = [
            '/home/benjamin/Desktop/desy_nov_2015/complete/run000061.raw',
            '/home/benjamin/Desktop/desy_nov_2015/complete/run000062.raw',
            '/home/benjamin/Desktop/desy_nov_2015/complete/run000063.raw',
          ]

def calibrate_and_reconstruct(rawfile):
  
  # Name for temporary folder created in tmp-runs/ 
  name = os.path.splitext(os.path.basename(rawfile))[0]
  # Tag for calibration data 
  caltag = name + 'python'
  
  # Calibrate the telescope using the rawfile. Creates a folder caltag 
  # containing all calibrations. 
  CalObj = Calibration(steerfiles=steerfiles, name=name + '-cal') 
  CalObj.calibrate(path=calpath,ifile=rawfile,caltag=caltag)  
   
  # Reconsruct the rawfile using caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=name + '-reco' )
  RecObj.reconstruct(path=['reco.xml'],ifile=rawfile,caltag=caltag) 
  
  return None

if __name__ == '__main__':
  count = multiprocessing.cpu_count()
  pool = multiprocessing.Pool(processes=count)
  pool.map(calibrate_and_reconstruct, runlist)
