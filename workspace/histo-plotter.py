"""
Simple script for plotting testbeam data using fully reconstructed root files 
from folder root-files. 

All DUT specific adjustements can be made in the DUTConfig dictionary. The 
plotter assumes pixel matrix with rectangular pixels.  

Usage: python histo-plotter.py

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *
import ROOT
import os

import argparse
parser = argparse.ArgumentParser(description="Perform plotting of test beam runs")
parser.add_argument('--pattern', dest='pattern', default='*', type=str, help='Pattern of run numbers to process')
args = parser.parse_args()
  

# Every plotting axis is given as a tuple (nbins,min,max) 
DUTConfig = { 'pitch_u' :          0.04,               # in mm 
              'pitch_v' :          0.036,              # in mm
              'residual_u_axis':   (151,-0.1,+0.1),    # in mm    
              'residual_v_axis':   (151,-0.1,+0.1),    # in mm 
              'charge_unit':       'ToT',   
              'seed_charge_axis':  (64,0,64),    
              'clus_charge_axis':  (300,0,300), 
              'ucell_axis':        (112,1,113),        
              'vcell_axis':        (224,1,225),       
              'sensor_u_axis':     (112,-0.5*112*0.04,0.5*112*0.04),
              'sensor_v_axis':     (224,-0.5*224*0.036,0.5*224*0.036), 
            }
       

    
for inputfilename in glob.glob('root-files/Histos-{}-reco.root'.format(args.pattern)): 
    
    # Open files with reconstructed run data 
    inputfile = ROOT.TFile( inputfilename, 'READ' )    
    
    # Create one histofile per run  
    histofile = ROOT.TFile( 'Plotter-' + os.path.basename(inputfilename), 'RECREATE', 'Histos created from file ' + inputfilename )
    
    # Add residual plots 
    residuals.plot(inputfile, histofile, basecut="hasTrack==0 && localChi2<20", Config=DUTConfig)
    
    # Add efficiency plots   
    efficiency.plot(inputfile, histofile, basecut="hasRefHit==0 && nDutDigits>0", matchcut="hasHit==0 && localChi2<20", uaxis=(112,1,113), vaxis=(224,1,225))
    
    # Add superpixel in-pixel charge plots 
    inpixel.plot_superpixel(inputfile, histofile, pixeltype=0, upitch=0.04, vpitch=0.036, ubins=20, vbins=20, ufold=2, vfold=2)             
      
    # Add superpixel in-pixel efficiency plots 
    efficiency.plot_super_inpix(inputfile, histofile, basecut="hasRefHit==0 && nDutDigits>0", matchcut="hasHit==0 && localChi2<20", upitch=0.04, vpitch=0.036, ubins=20, vbins=20)
    
    # Make a pdf containing all plots 
    pdfName = os.path.splitext( os.path.basename( inputfilename ) )[0] + '.pdf' 
    residuals.make_pdf(histofile, pdfName)

    # Close all files 
    histofile.Write()
    histofile.Close()
    inputfile.Close()    
    
    
 
