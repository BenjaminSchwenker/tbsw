"""
Simple script for plotting testbeam data using fully reconstructed root files 
from folder root-files. 

All DUT specific adjustements can be made in the DUTConfig dictionary. The 
plotter assumes pixel matrix with rectangular pixels.  

python3 histo-plotter-tj2.py --ifile=root-files/Histos-<something>.root

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

#from tbsw import *
import tbsw.residuals as residuals
import tbsw.efficiency as efficiency
import tbsw.inpixel as inpixel
import ROOT
import os
import glob

import argparse
parser = argparse.ArgumentParser(description="Perform plotting of test beam runs")
parser.add_argument('--ifile', dest='ifile', default='*', type=str, help='Input file pattern of runs to process')
args = parser.parse_args()
  

# Every plotting axis is given as a tuple (nbins,min,max) 
DUTConfig = { 'pitch_u' :          0.033,              # in mm 
              'pitch_v' :          0.033,              # in mm
              'residual_u_axis':   (151,-0.1,+0.1),    # in mm    
              'residual_v_axis':   (151,-0.1,+0.1),    # in mm 
              'charge_unit':       'ToT / LSB',   
              'seed_charge_axis':  (128,0,128),    
              'clus_charge_axis':  (256,0,256), 
              'ucell_axis':        (512,0,512),        
              'vcell_axis':        (512,0,512),       
              'sensor_u_axis':     (512,-0.5*512*0.033,0.5*512*0.033),
              'sensor_v_axis':     (512,-0.5*512*0.033,0.5*512*0.033), 
            }
       

    
for inputfilename in glob.glob(args.ifile): 
    
    # Open files with reconstructed run data 
    inputfile = ROOT.TFile( inputfilename, 'READ' )    
    
    # Create one histofile per run  
    histofile = ROOT.TFile( 'Plotter-' + os.path.basename(inputfilename), 'RECREATE', 'Histos created from file ' + inputfilename )
    
    # Add residual plots 
    residuals.plot(inputfile, histofile, basecut="hasTrack==0 && localChi2<20", Config=DUTConfig)
    
    # Add efficiency plots   
    efficiency.plot(inputfile, histofile, basecut="hasTestPixels==1", matchcut="hasHit==0 && localChi2<20", uaxis=(512,0,512), vaxis=(512,0,512))
    
    # Add superpixel in-pixel charge plots 
    inpixel.plot_superpixel(inputfile, histofile, pixeltype=0, upitch=DUTConfig['pitch_u'], vpitch=DUTConfig['pitch_v'], ubins=20, vbins=20, ufold=2, vfold=2)             
      
    # Add superpixel in-pixel efficiency plots 
    efficiency.plot_super_inpix(inputfile, histofile, basecut="maskedPixel==0", matchcut="hasHit==0 && localChi2<20", upitch=DUTConfig['pitch_u'], vpitch=DUTConfig['pitch_v'], ubins=20, vbins=20)

    # Compute efficiency (and error) in specified ROI
    efficiency.extract_roi(inputfile, basecut="hasRefHit==0 && maskedPixel==0", matchcut="hasHit==0 && localChi2<20")
    efficiency.extract_roi(inputfile, basecut="maskedPixel==0 && cellU_fit>1 && cellU_fit<220 && cellV_fit>190 && cellV_fit<200", matchcut="hasHit==0 && localChi2<20")


    # Make a pdf containing all plots 
    pdfName = os.path.splitext( os.path.basename( inputfilename ) )[0] + '.pdf' 
    residuals.make_pdf(histofile, pdfName)

    # Close all files 
    histofile.Write()
    histofile.Close()
    inputfile.Close()    
    
    
 
