"""
This is an example script to demonstrate how TBSW can be used to analyze test beam 
data using Python scripts.

The data processing uses Pythons multiprocessing.Pool to distribute runs to multiple cores and 
speed up processing time. 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *
import tbsw_tools
import multiprocessing

# Path to steering files 
steerfiles_partone = 'steering-files/depfet-H6-tb/'
steerfiles_parttwo = 'steering-files/depfet-H6-tb-geoid16/'

# Defines the sequence of calibration steps. 
# XML steer files are taken from steerfiles
# folder. 
calpath = [ 
           'hotpixelkiller.xml' ,              
           'correlator.xml' ,
           'kalmanalign-iteration-1.xml',
           'kalmanalign-iteration-2.xml',
           'kalmanalign-iteration-2.xml',
           'kalmanalign-iteration-2.xml',
           'telescope-dqm.xml',
           #'cluster-calibration-dep.xml',
           #'kalmanalign-iteration-3.xml',
           #'kalmanalign-iteration-3.xml',
           #'kalmanalign-iteration-3.xml',
         ]


# List of runs from part one of H6 measurement. Runnumbers range from 5-113.
runlist_partone = [
                    #'/home/benjamin/Desktop/desy_nov_2015/complete/run000005.raw',
                    #'/home/benjamin/Desktop/desy_nov_2015/complete/run000113.raw',
                    '/home/benjamin/Desktop/desy_nov_2015/complete/run000065.raw',
                  ]

# List of runs from part two of H6 measurement. Runnumbers range from 437-537.
runlist_parttwo = [
                    #'/home/benjamin/Desktop/desy_nov_2015/complete/run000437.raw',
                    '/home/benjamin/Desktop/desy_nov_2015/complete/run000508.raw',
                  ]


def calibrate_and_reconstruct(params):
  
  rawfile, steerfiles = params

  # Tag for calibration data
  caltag = os.path.splitext(os.path.basename(rawfile))[0]
  
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
  CalObj = Calibration(steerfiles=steerfiles, name=caltag + '-cal') 
  CalObj.calibrate(path=calpath,ifile=rawfile,caltag=caltag)  
   
  # Reconsruct the rawfile using caltag. Resulting root files are 
  # written to folder root-files/
  RecObj = Reconstruction(steerfiles=steerfiles, name=caltag + '-reco' )
  RecObj.reconstruct(path=['reco.xml'],ifile=rawfile,caltag=caltag) 
  
  return None

if __name__ == '__main__':
  count = multiprocessing.cpu_count()
  pool = multiprocessing.Pool(processes=count)
  
  # Runs from first part must be processed with steerfiles_partone 
  params = [ (run, steerfiles_partone) for run in runlist_partone ]
  pool.map(calibrate_and_reconstruct, params)

  # Runs from second part must be processed with steerfiles_parttwo 
  params = [ (run, steerfiles_parttwo) for run in runlist_parttwo ]
  pool.map(calibrate_and_reconstruct, params)
  
  
  # Create some efficiency and residual plots from root files after 
  # reconstruction. 
  trackfiles =  glob.glob('root-files/Histos*.root')
   
  for trackfile in trackfiles: 
        
    ofile = 'Residuals-' + os.path.basename(trackfile)
    tbsw_tools.residuals.plot(inputfilename = trackfile, histofilename = ofile, basecut = "hasTrack==0")
      
    ofile = 'Efficiency-' + os.path.basename(trackfile)  
    tbsw_tools.efficiency.plot(inputfilename = trackfile, histofilename = ofile, basecut = "nTelTracks == 1 && nDutHits > 0 && cellU_fit >= 0 && cellU_fit < 64", matchcut="hasHit == 0", ucells=64, vcells=480)




