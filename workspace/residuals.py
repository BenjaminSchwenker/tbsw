from ROOT import TFile, TH1F, TH2F
from ROOT import gROOT, Double, TCut


def plot_residuals(inputfilename = None, histofilename = "MyHistos.root" , cutstring = "hasTrack==0"):
  
  # Analysis cuts 
  hitcut = TCut(cutstring)
  
  if inputfilename == None:
    return None
  
  
  rawfile = gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = TFile( inputfilename, 'READ' )
  
  histofile = gROOT.FindObject( histofilename )
  if histofile:
    histofile.Close()
  histofile = TFile( histofilename, 'RECREATE', 'Residual plots created from input file ' + inputfilename )
  
  # Get access to tracks 
  tracktree = rawfile.Get("Track")
  
  tracktree.Draw("cellV_fit:cellU_fit>>htrackspot", "", "goff")
  htrackspot = gROOT.FindObject("htrackspot")
  if htrackspot:
    htrackspot.SetTitle("") 
    htrackspot.GetXaxis().SetTitle("cell U [ID]")
    htrackspot.GetYaxis().SetTitle("cell V [cell ID]")
    htrackspot.GetZaxis().SetTitle("number of tracks")
    htrackspot.SetStats(0)

  tracktree.Draw("v_fit:u_fit>>hbeamspot", "","goff")
  hbeamspot = gROOT.FindObject("hbeamspot")
  if hbeamspot:
    hbeamspot.SetTitle("");
    hbeamspot.GetXaxis().SetTitle("u [mm]")
    hbeamspot.GetYaxis().SetTitle("v [mm]")
    hbeamspot.GetZaxis().SetTitle("number of tracks")
    hbeamspot.SetStats(0)

  
  #Get access to hits
  hittree = rawfile.Get("Hit")

  hittree.Draw("cellV_hit:cellU_hit>>hspot", hitcut,"goff")
  hspot = gROOT.FindObject("hspot")
  if hspot: 
    hspot.SetTitle("")
    hspot.GetXaxis().SetTitle("cell U [ID]")
    hspot.GetYaxis().SetTitle("cell V [ID]")
    hspot.GetZaxis().SetTitle("number of clusters")
    hspot.SetStats(0)


  hittree.Draw("(u_hit - u_fit)*1000 >>hu(101,-75.,+75.)", hitcut,"goff")
  hu = gROOT.FindObject("hu")
  if hu: 
    hu.SetTitle("")
    hu.GetXaxis().SetTitle("residuals u [#mum]")
    hu.GetYaxis().SetTitle("number of hits")
    hu.GetYaxis().SetTitleOffset(1.2)


  hittree.Draw("(v_hit - v_fit)*1000 >>hv(81,-75.,+75.)", hitcut,"goff")
  hv = gROOT.FindObject("hv")
  if hv:
    hv.SetTitle("")
    hv.GetXaxis().SetTitle("residuals v [#mum]")
    hv.GetYaxis().SetTitle("number of hits")
    hv.GetYaxis().SetTitleOffset(1.2)


  hittree.Draw("trackChi2 / trackNdof >> hchisqundof(100,0,10)", hitcut,"goff")
  hchisqundof = gROOT.FindObject("hchisqundof")
  if hchisqundof:
    hchisqundof.SetTitle("")
    hchisqundof.GetXaxis().SetTitle("#chi^2/ndof")
    hchisqundof.GetYaxis().SetTitle("number of tracks")
    hchisqundof.GetYaxis().SetTitleOffset(1.2)


  hittree.Draw("seedCharge >>hCharge(101,0,100)", hitcut,"goff");
  hCharge = gROOT.FindObject("hCharge")
  if hCharge:
    hCharge.SetTitle("")
    hCharge.GetXaxis().SetTitle("seed charge [ADU]")
    hCharge.GetYaxis().SetTitle("number of hits")
    hCharge.GetYaxis().SetTitleOffset(1.2)

  hittree.Draw("clusterCharge >>hCCharge(101,0,100)", hitcut,"goff")
  hCCharge = gROOT.FindObject("hCCharge")
  if hCCharge:
    hCCharge.SetTitle("")
    hCCharge.GetXaxis().SetTitle("cluster charge [ADU]")
    hCCharge.GetYaxis().SetTitle("number of hits")
    hCCharge.GetYaxis().SetTitleOffset(1.2)


  hittree.Draw("size >>hsize(20,0,20)", hitcut,"goff");
  hsize = gROOT.FindObject("hsize")
  if hsize:
    hsize.SetTitle("")
    hsize.GetXaxis().SetTitle("cluster size [pixels]")
    hsize.GetYaxis().SetTitle("number of hits")
    hsize.GetYaxis().SetTitleOffset(1.2)

  hittree.Draw("sizeU >>hsizeU(10,0,10)", hitcut,"goff")
  hsizeU = gROOT.FindObject("hsizeU")
  if hsizeU:
    hsizeU.SetTitle("")
    hsizeU.GetXaxis().SetTitle("size [ucells]")
    hsizeU.GetYaxis().SetTitle("number of hits")
    hsizeU.GetYaxis().SetTitleOffset(1.2)

  hittree.Draw("sizeV >>hsizeV(10,0,10)", hitcut,"goff")
  hsizeV = gROOT.FindObject("hsizeV")
  if hsizeV:
    hsizeV.SetTitle("")
    hsizeV.GetXaxis().SetTitle("size [vcells]")
    hsizeV.GetYaxis().SetTitle("number of hits")
    hsizeV.GetYaxis().SetTitleOffset(1.2)


  # write the tree into the output file and close the file
  histofile.Write()
  histofile.Close()
  rawfile.Close()



