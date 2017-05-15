"""
Simple variant of testbeam.py written for reconstruction of ministrip data 
from Luises test beam. 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""


from tbsw import *
import multiprocessing
import residuals
import occupancy
import efficiency
import inpixel

      
# Path to steering files 
steerfiles = 'steering-files/luise-tb/'  

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
            '../workspace_desy_itk_luise_test/data_from_kristin/run000282-merger.slcio',
          ]

def calibrate_and_reconstruct(rawfile):
  
  # Name for temporary folder created in tmp-runs/ 
  name = os.path.splitext(os.path.basename(rawfile))[0]
  # Tag for calibration data 
  caltag = name + 'test'
  
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
  
  trackfiles =  glob.glob('root-files/*Histos*.root')
   
  for trackfile in trackfiles: 
        
    ofile = 'Residuals-' + os.path.basename(trackfile)
    residuals.plot(inputfilename = trackfile, histofilename = ofile, basecut = "hasTrack==0")
    
    ofile = 'Occupancy-' + os.path.basename(trackfile)
    occupancy.plot(inputfilename = trackfile, histofilename = ofile, ucells=128, vcells=1)
    
    ofile = 'Efficiency-' + os.path.basename(trackfile)  
    efficiency.plot(inputfilename = trackfile, histofilename = ofile, ucells=128, vcells=1)
    
    ofile = 'InpixEfficiency-' + os.path.basename(trackfile)  
    efficiency.plot_unit_inpix(inputfilename = trackfile, histofilename = ofile, basecut = "nTelTracks == 1", matchcut="hasHit == 0", upitch=0.0745, vpitch=0.5, ubins=20, vbins=10)
    
    ofile = 'UnitInpixEfficiency-' + os.path.basename(trackfile)
    efficiency.plot_inpix(inputfilename = trackfile, histofilename = ofile, basecut = "nTelTracks == 1", matchcut="hasHit == 0", usize=1.0, vsize=1.0, ubins=200, vbins=1)
    
    ofile = 'InpixSignal-' + os.path.basename(trackfile)
    inpixel.plot(inputfilename = trackfile, histofilename = ofile, usize=1.0, vsize=1.0, ubins=100, vbins=10)
    
    ofile = 'UnitInpixSignal-' + os.path.basename(trackfile)
    inpixel.plot_unit(inputfilename = trackfile, histofilename = ofile, upitch=0.0745, vpitch=0.5, ubins=20, vbins=10)
    
       
