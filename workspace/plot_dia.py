"""
Script for plotting mini 3D diamond test beam data for Helge. 

Author: Benjamin Schwenker <benjamin.schwenker@phys.uni-goettingen.de>  
"""

from tbsw import *


trackfiles =  glob.glob('root-files/Histos-FEI*.root')  
  
for trackfile in trackfiles: 
        
    ofile = 'Residuals-' + os.path.basename(trackfile)
    residuals.plot(inputfilename = trackfile, histofilename = ofile, basecut = "hasTrack==0")
    
    ofile = 'Efficiency-' + os.path.basename(trackfile)  
    efficiency.plot(inputfilename = trackfile, histofilename = ofile, ucells=400, vcells=400)
    
    #ofile = 'InpixEfficiency-' + os.path.basename(trackfile)  
    #efficiency.plot_unit_inpix(inputfilename = trackfile, histofilename = ofile, basecut = "nTelTracks == 1", matchcut="hasHit == 0", upitch=0.0745, vpitch=0.5, ubins=20, vbins=10)
    
    #ofile = 'UnitInpixEfficiency-' + os.path.basename(trackfile)
    #efficiency.plot_inpix(inputfilename = trackfile, histofilename = ofile, basecut = "nTelTracks == 1", matchcut="hasHit == 0", usize=1.0, vsize=1.0, ubins=200, vbins=1)
    
    #ofile = 'InpixSignal-' + os.path.basename(trackfile)
    #inpixel.plot(inputfilename = trackfile, histofilename = ofile, usize=1.0, vsize=1.0, ubins=100, vbins=10)
    
    #ofile = 'UnitInpixSignal-' + os.path.basename(trackfile)
    #inpixel.plot_unit(inputfilename = trackfile, histofilename = ofile, upitch=0.0745, vpitch=0.5, ubins=20, vbins=10)


trackfiles =  glob.glob('root-files/Histos-DIA*.root')  
  
for trackfile in trackfiles: 
        
    ofile = 'Residuals-' + os.path.basename(trackfile)
    residuals.plot(inputfilename = trackfile, histofilename = ofile, basecut = "hasTrack==0")
    
    ofile = 'Efficiency-' + os.path.basename(trackfile)  
    efficiency.plot(inputfilename = trackfile, histofilename = ofile, ucells=400, vcells=400)
    
    #ofile = 'InpixEfficiency-' + os.path.basename(trackfile)  
    #efficiency.plot_unit_inpix(inputfilename = trackfile, histofilename = ofile, basecut = "nTelTracks == 1", matchcut="hasHit == 0", upitch=0.0745, vpitch=0.5, ubins=20, vbins=10)
    
    #ofile = 'UnitInpixEfficiency-' + os.path.basename(trackfile)
    #efficiency.plot_inpix(inputfilename = trackfile, histofilename = ofile, basecut = "nTelTracks == 1", matchcut="hasHit == 0", usize=1.0, vsize=1.0, ubins=200, vbins=1)
    
    #ofile = 'InpixSignal-' + os.path.basename(trackfile)
    #inpixel.plot(inputfilename = trackfile, histofilename = ofile, usize=1.0, vsize=1.0, ubins=100, vbins=10)
    
    #ofile = 'UnitInpixSignal-' + os.path.basename(trackfile)
    #inpixel.plot_unit(inputfilename = trackfile, histofilename = ofile, upitch=0.0745, vpitch=0.5, ubins=20, vbins=10)
