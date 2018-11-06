"""
Script for plotting depfet testbeam data

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *

import argparse
parser = argparse.ArgumentParser(description="Perform plotting of a test beam runs")
parser.add_argument('--pattern', dest='pattern', default='*', type=str, help='Pattern of run numbers to process')
args = parser.parse_args()
  
  
for trackfile in glob.glob('root-files/Histos-FEI4-{}-reco.root'.format(args.pattern)): 
        
    ofile = 'Residuals-' + os.path.basename(trackfile)
    residuals.plot(inputfilename = trackfile, histofilename = ofile, basecut = "hasTrack==0", nbins=201, urange=400, vrange=400)
    
    ofile = 'Efficiency-' + os.path.basename(trackfile)  
    efficiency.plot(inputfilename=trackfile, histofilename=ofile, basecut="trackNHits==7 && iEvt>=0 && nDutDigits>0", matchcut="hasHit==0", uaxis=(80,0,80), vaxis=(336,0,336))
    

for trackfile in glob.glob('root-files/Histos-DEPH5-{}-reco.root'.format(args.pattern)): 
        
    ofile = 'Residuals-' + os.path.basename(trackfile)
    residuals.plot(inputfilename = trackfile, histofilename = ofile, basecut = "hasTrack==0", nbins=201, urange=400, vrange=400)
    
    ofile = 'Efficiency-' + os.path.basename(trackfile)  
    efficiency.plot(inputfilename=trackfile, histofilename=ofile, basecut="trackNHits==7 && iEvt>=0 && nDutDigits>4", matchcut="hasHit==0", uaxis=(10,0,64), vaxis=(10,0,64))
    
  
for trackfile in glob.glob('root-files/Histos-DEPBIG-{}-reco.root'.format(args.pattern)): 
        
    ofile = 'Residuals-' + os.path.basename(trackfile)
    residuals.plot(inputfilename=trackfile, histofilename=ofile, basecut="hasTrack==0", nbins=201, urange=400, vrange=400)
    
    ofile = 'Efficiency-' + os.path.basename(trackfile)  
    efficiency.plot(inputfilename=trackfile, histofilename=ofile, basecut="trackNHits==7 && iEvt>=0 && nDutDigits>4", matchcut="hasHit==0", uaxis=(100,0,250), vaxis=(100,0,768))

    # A single run will not give enough statistic for such a plot, you will need to chain tree from many runs for this
    # Anyway, here is the commond to obtain inpixel plots for the complete sensor area
    #ofile = 'InpixMaps-' + os.path.basename(trackfile)  
    #inpixel.plot(inputfilename=trackfile, histofilename=ofile, uaxis=(100,-4,-2), vaxis=(150,7,10))

    # Even a single run of 10min is enough to obtain a inpix plot for a generic cell of (ufold,vfold) superpixels. 
    ofile = 'InpixUnitMaps-' + os.path.basename(trackfile)  
    inpixel.plot_unit(inputfilename=trackfile, histofilename=ofile, upitch=0.05, vpitch=0.085, ubins=30, vbins=30, ufold=2, vfold=2)
