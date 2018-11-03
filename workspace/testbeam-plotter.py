"""
Script for plotting depfet testbeam data

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *
  
  
for trackfile in glob.glob('root-files/Histos-FEI*.root'): 
        
    ofile = 'Residuals-' + os.path.basename(trackfile)
    residuals.plot(inputfilename = trackfile, histofilename = ofile, basecut = "hasTrack==0", nbins=201, urange=400, vrange=400)
    
    ofile = 'Efficiency-' + os.path.basename(trackfile)  
    #efficiency.plot(inputfilename=trackfile, histofilename=ofile, basecut="nTelTracks==1", matchcut="hasHit==0", uaxis=(64,0,64), vaxis=(64,0,64))
    

for trackfile in glob.glob('root-files/Histos-DEPH5*.root'): 
        
    ofile = 'Residuals-' + os.path.basename(trackfile)
    residuals.plot(inputfilename = trackfile, histofilename = ofile, basecut = "hasTrack==0", nbins=201, urange=400, vrange=400)
    
    ofile = 'Efficiency-' + os.path.basename(trackfile)  
    efficiency.plot(inputfilename=trackfile, histofilename=ofile, basecut="nTelTracks==1", matchcut="hasHit==0", uaxis=(64,0,64), vaxis=(64,0,64))
    
  
for trackfile in glob.glob('root-files/Histos-DEPBIG*.root'): 
        
    ofile = 'Residuals-' + os.path.basename(trackfile)
    residuals.plot(inputfilename=trackfile, histofilename=ofile, basecut="hasTrack==0", nbins=201, urange=400, vrange=400)
    
    ofile = 'Efficiency-' + os.path.basename(trackfile)  
    efficiency.plot(inputfilename=trackfile, histofilename=ofile, basecut="nTelTracks>1 && cellV_fit<450 && cellV_fit>300 && cellU_fit>50 && cellU_fit<220", matchcut="hasHit==0", uaxis=(50,0,250), vaxis=(50,0,768))


