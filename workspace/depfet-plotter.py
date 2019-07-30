"""
Simple script for plotting depfet testbeam data using fully reconstructed root files 
from folder root-files. 

For example, this script can be run on the output root files from running the script 
depfet-reco.py. 

Usage: 

python depfet-plotter.py

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *

import argparse
parser = argparse.ArgumentParser(description="Perform plotting of a test beam runs")
parser.add_argument('--pattern', dest='pattern', default='*', type=str, help='Pattern of run numbers to process')
parser.add_argument('--mapping', dest='mapping', default='IF', type=str, help='IF, OF')
args = parser.parse_args()
  
    
for trackfile in glob.glob('root-files/Histos-PXD-{}-reco.root'.format(args.pattern)): 
        
    ofile = 'Residuals-LargePitch-' + os.path.basename(trackfile)
    residuals.plot(inputfilename=trackfile, histofilename=ofile, basecut="hasTrack==0 && cellV_hit<512", nbins=501, urange=400, vrange=400)

    ofile = 'Residuals-SmallPitch-' + os.path.basename(trackfile)
    residuals.plot(inputfilename=trackfile, histofilename=ofile, basecut="hasTrack==0 && cellV_hit>=512", nbins=501, urange=400, vrange=400)
    
    ofile = 'Efficiency-' + os.path.basename(trackfile)  
    efficiency.plot(inputfilename=trackfile, histofilename=ofile, basecut="hasRefHit==0 && hasTestPixels==1", matchcut="hasHit==0", uaxis=(250,0,250), vaxis=(768,0,768))

    # A single run will not give enough statistic for such a plot, you will need to chain tree from many runs for this
    # Anyway, here is the commond to obtain inpixel plots for the complete sensor area
    #ofile = 'InpixMaps-' + os.path.basename(trackfile)  
    #inpixel.plot(inputfilename=trackfile, histofilename=ofile, uaxis=(100,-4,-2), vaxis=(150,7,10))

    if args.mapping == 'OF':

      # Even a single run of 10min is enough to obtain a inpix plot for a generic cell of (ufold,vfold) superpixels. 
      ofile = 'InpixUnitMaps-OF-LargePitch-' + os.path.basename(trackfile)  
      inpixel.plot_unit(inputfilename=trackfile, histofilename=ofile, pixeltype=0, upitch=0.05, vpitch=0.085, ubins=30, vbins=30, ufold=2, vfold=2)             
      
      # Even a single run of 10min is enough to obtain a inpix plot for a generic cell of (ufold,vfold) superpixels. 
      ofile = 'InpixUnitMaps-OF-SmallPitch' + os.path.basename(trackfile)  
      inpixel.plot_unit(inputfilename=trackfile, histofilename=ofile, pixeltype=1, upitch=0.05, vpitch=0.070, ubins=30, vbins=30, ufold=2, vfold=2)         

    elif args.mapping == 'IF':

      # Even a single run of 10min is enough to obtain a inpix plot for a generic cell of (ufold,vfold) superpixels. 
      ofile = 'InpixUnitMaps-IF-LargePitch' + os.path.basename(trackfile)  
      inpixel.plot_unit(inputfilename=trackfile, histofilename=ofile, pixeltype=0, upitch=0.05, vpitch=0.060, ubins=30, vbins=30, ufold=2, vfold=2)             
      
      # Even a single run of 10min is enough to obtain a inpix plot for a generic cell of (ufold,vfold) superpixels. 
      ofile = 'InpixUnitMaps-IF-SmallPitch' + os.path.basename(trackfile)  
      inpixel.plot_unit(inputfilename=trackfile, histofilename=ofile, pixeltype=1, upitch=0.05, vpitch=0.055, ubins=30, vbins=30, ufold=2, vfold=2)      

    else: 
      
      print("No inpixel maps for PXD matrix produced since no valid mapping was specified")      
