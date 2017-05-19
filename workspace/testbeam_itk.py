"""
This is an example script to demonstrate how TBSW can be used to analyze test beam 
data using Python scripts.

The use case is a configuration the EUDET/AIDA telescope with two ministrip sensors 
installed as devices under test betweenn the EUDET arms. Input data are lcio files 
containing collections with Mimosa26 digits and digits from ministrip sensors.  

The script below calibrates and reconstructs all runs from a runlist. In this example, 
one common set of XML steerfiles is applied to all runs. The calibration of the data 
is done for each run using the calibration steps given in the calpath. The calibration 
results are stored in the folder cal-files. After the calibration, the main reconstruction
of the run is done using the calibration data prepared before. The reconstruction 
executes a single XML steerfile, by convention called reco.xml. The output are one or 
more root files containing root TTree objects for plotting by end users. As a last step, 
a couple of plotting scripts are started to create some pre-defined histograms from the 
trees using PyROOT. 

The data processing uses Pythons multiprocessing.Pool to distribute runs to multiple cores and 
speed up processing time. 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *
import tbsw_tools
import multiprocessing

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


# Create a list of runs to be processed using the same set of steerfiles, i.e.
# having same number of modules, beam energy etc. Small movements of sensors 
# can be tolerated since alignment is done run by run. 
runlist = [
            '../workspace_desy_itk_luise_test/data_from_kristin/run000282-merger.slcio',
          ]


def calibrate_and_reconstruct(rawfile):
  
  # Name for temporary folder created in tmp-runs/ 
  name = os.path.splitext(os.path.basename(rawfile))[0]
  # Tag for calibration data.
  caltag = name + 'python'
  
  # Calibrate of the run using beam data. Creates a folder cal-files/caltag 
  # containing all calibration data. 
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

  # Create some standard plots from root files after 
  # reconstruction. 
  trackfiles =  glob.glob('root-files/*Histos*.root')
   
  for trackfile in trackfiles: 
        
    ofile = 'Residuals-' + os.path.basename(trackfile)
    tbsw_tools.residuals.plot(inputfilename = trackfile, histofilename = ofile, basecut = "hasTrack==0")
    
    ofile = 'Occupancy-' + os.path.basename(trackfile)
    tbsw_tools.occupancy.plot(inputfilename = trackfile, histofilename = ofile, ucells=128, vcells=1)
    
    ofile = 'Efficiency-' + os.path.basename(trackfile)  
    tbsw_tools.efficiency.plot(inputfilename = trackfile, histofilename = ofile, ucells=128, vcells=1)
    
    ofile = 'InpixEfficiency-' + os.path.basename(trackfile)  
    tbsw_tools.efficiency.plot_unit_inpix(inputfilename = trackfile, histofilename = ofile, basecut = "nTelTracks == 1", matchcut="hasHit == 0", upitch=0.0745, vpitch=0.5, ubins=20, vbins=10)
    
    ofile = 'UnitInpixEfficiency-' + os.path.basename(trackfile)
    tbsw_tools.efficiency.plot_inpix(inputfilename = trackfile, histofilename = ofile, basecut = "nTelTracks == 1", matchcut="hasHit == 0", usize=1.0, vsize=1.0, ubins=200, vbins=1)
    
    ofile = 'InpixSignal-' + os.path.basename(trackfile)
    tbsw_tools.inpixel.plot(inputfilename = trackfile, histofilename = ofile, usize=1.0, vsize=1.0, ubins=100, vbins=10)
    
    ofile = 'UnitInpixSignal-' + os.path.basename(trackfile)
    tbsw_tools.inpixel.plot_unit(inputfilename = trackfile, histofilename = ofile, upitch=0.0745, vpitch=0.5, ubins=20, vbins=10)
