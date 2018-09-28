"""
Script for plotting mini 3D diamond test beam data for Helge. 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *


trackfiles =  glob.glob('root-files/Histos-FEI*.root')  
  
for trackfile in trackfiles: 
        
    ofile = 'Residuals-' + os.path.basename(trackfile)
    residuals.plot(inputfilename = trackfile, histofilename = ofile, basecut = "hasTrack==0", urange=400, vrange=400)
    
    ofile = 'Efficiency-' + os.path.basename(trackfile)  
    efficiency.plot(inputfilename = trackfile, histofilename = ofile, ucells=400, vcells=400)
    
    


trackfiles =  glob.glob('root-files/Histos-DIA*.root')  
  
for trackfile in trackfiles: 
        
    ofile = 'Residuals-' + os.path.basename(trackfile)
    residuals.plot(inputfilename = trackfile, histofilename = ofile, basecut = "hasTrack==0", urange=400, vrange=400)
    
    ofile = 'Efficiency-' + os.path.basename(trackfile)  
    efficiency.plot(inputfilename = trackfile, histofilename = ofile, ucells=400, vcells=400)
    
    
