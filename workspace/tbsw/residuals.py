from ROOT import TFile, TH1F, TH2F
from ROOT import gROOT, Double, TCut


def plot(inputfilename = None, histofilename = "ResidualHistos.root" , basecut = "hasTrack==0", nbins=101, urange=200, vrange=200):
  
  # Analysis cuts 
  hitcut = TCut(basecut)
  
  if inputfilename == None:
    return None
  
  
  rawfile = gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = TFile( inputfilename, 'READ' )
  
  histofile = gROOT.FindObject( histofilename )
  if histofile:
    histofile.Close()
  histofile = TFile( histofilename, 'RECREATE', 'Residual histos created from input file ' + inputfilename )
  
  # Get access to tracks 
  tracktree = rawfile.Get("Track")
  
  htrackspot = TH2F("htrackspot","",40,0,0,40,0,0)
  tracktree.Draw("cellV_fit:cellU_fit>>+htrackspot", "", "goff")
  #htrackspot.SetTitle("") 
  htrackspot.GetXaxis().SetTitle("cellU_fit [cell ID]")
  htrackspot.GetYaxis().SetTitle("cellV_fit [cell ID]")
  htrackspot.GetZaxis().SetTitle("number of tracks")
  htrackspot.SetStats(0)

  hbeamspot = TH2F("hbeamspot","",40,0,0,40,0,0)
  tracktree.Draw("v_fit:u_fit>>+hbeamspot", "","goff")
  hbeamspot.SetTitle("")
  hbeamspot.GetXaxis().SetTitle("u_fit [mm]")
  hbeamspot.GetYaxis().SetTitle("v_fit [mm]")
  hbeamspot.GetZaxis().SetTitle("number of tracks")
  hbeamspot.SetStats(0)
  
   
  #Get access to hits
  hittree = rawfile.Get("Hit")

  hspot = TH2F("hspot","",40,0,0,40,0,0)
  hittree.Draw("cellV_hit:cellU_hit>>hspot", hitcut,"goff")
  hspot.SetTitle("")
  hspot.GetXaxis().SetTitle("cellU_hit [cell ID]")
  hspot.GetYaxis().SetTitle("cellV_hit [cell ID]")
  hspot.GetZaxis().SetTitle("number of clusters")
  hspot.SetStats(0)

  hu = TH1F("hu","",nbins,-urange,+urange)
  hittree.Draw("(u_hit - u_fit)*1000 >>+hu", hitcut,"goff")
  hu.SetTitle("")
  hu.GetXaxis().SetTitle("residuals u [#mum]")
  hu.GetYaxis().SetTitle("number of hits")
  hu.GetYaxis().SetTitleOffset(1.2)

  hv = TH1F("hv","",nbins,-vrange,+vrange)
  hittree.Draw("(v_hit - v_fit)*1000 >>+hv", hitcut,"goff")
  hv.SetTitle("")
  hv.GetXaxis().SetTitle("residuals v [#mum]")
  hv.GetYaxis().SetTitle("number of hits")
  hv.GetYaxis().SetTitleOffset(1.2)
  
  huv = TH2F("huv","",nbins,-urange,+urange, nbins,-vrange,+vrange )
  hittree.Draw("(v_hit - v_fit)*1000:(u_hit - u_fit)*1000 >>+huv", hitcut,"goff")  
  huv.SetTitle("")
  huv.GetXaxis().SetTitle("residuals u [#mum]")
  huv.GetYaxis().SetTitle("residuals v [#mum]")
  huv.GetZaxis().SetTitle("number of hits")
  huv.SetStats(0)
  
  hchisqundof = TH1F("hchisqundof","",100,0,0)
  hittree.Draw("trackChi2 / trackNdof >> +hchisqundof", hitcut,"goff")
  hchisqundof.SetTitle("")
  hchisqundof.GetXaxis().SetTitle("#chi^2/ndof")
  hchisqundof.GetYaxis().SetTitle("number of tracks")
  hchisqundof.GetYaxis().SetTitleOffset(1.2)
  
  
  hCharge = TH1F("hCharge","",255,0,255)
  hittree.Draw("seedCharge >>+hCharge", hitcut,"goff")
  hCharge.SetTitle("")
  hCharge.GetXaxis().SetTitle("seed charge [ADU]")
  hCharge.GetYaxis().SetTitle("number of hits")
  hCharge.GetYaxis().SetTitleOffset(1.2)

  hCCharge = TH1F("hCCharge","",255,0,255)
  hittree.Draw("clusterCharge >>+hCCharge", hitcut,"goff")
  hCCharge.SetTitle("")
  hCCharge.GetXaxis().SetTitle("cluster charge [ADU]")
  hCCharge.GetYaxis().SetTitle("number of hits")
  hCCharge.GetYaxis().SetTitleOffset(1.2)

  hsize = TH1F("hsize","",10,0,10)
  hittree.Draw("size >>+hsize", hitcut,"goff")
  hsize.SetTitle("")
  hsize.GetXaxis().SetTitle("cluster size [pixels]")
  hsize.GetYaxis().SetTitle("number of hits")
  hsize.GetYaxis().SetTitleOffset(1.2)

  hsizeU = TH1F("hsizeU","",10,0,10)
  hittree.Draw("sizeU >>+hsizeU", hitcut,"goff")
  hsizeU.SetTitle("")
  hsizeU.GetXaxis().SetTitle("sizeU [cell ID]")
  hsizeU.GetYaxis().SetTitle("number of hits")
  hsizeU.GetYaxis().SetTitleOffset(1.2)

  hsizeV = TH1F("hsizeV","",10,0,10)
  hittree.Draw("sizeV >>+hsizeV", hitcut,"goff")
  hsizeV.SetTitle("")
  hsizeV.GetXaxis().SetTitle("sizeV [cell ID]")
  hsizeV.GetYaxis().SetTitle("number of hits")
  hsizeV.GetYaxis().SetTitleOffset(1.2)
   
  # write the tree into the output file and close the file
  histofile.Write()
  histofile.Close()
  rawfile.Close()



