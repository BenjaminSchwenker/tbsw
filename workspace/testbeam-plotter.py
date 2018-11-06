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
    efficiency.plot(inputfilename=trackfile, histofilename=ofile, basecut="trackNHits==7 && iEvt>=0 && nDutDigits>0", matchcut="hasHit==0", uaxis=(10,0,64), vaxis=(10,0,64))
    
  
for trackfile in glob.glob('root-files/Histos-DEPBIG-{}-reco.root'.format(args.pattern)): 
        
    ofile = 'Residuals-' + os.path.basename(trackfile)
    residuals.plot(inputfilename=trackfile, histofilename=ofile, basecut="hasTrack==0", nbins=201, urange=400, vrange=400)
    
    ofile = 'Efficiency-' + os.path.basename(trackfile)  
    efficiency.plot(inputfilename=trackfile, histofilename=ofile, basecut="trackNHits==7 && iEvt>=0 && nDutDigits>0", matchcut="hasHit==0", uaxis=(100,0,250), vaxis=(100,0,768))
