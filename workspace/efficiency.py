from ROOT import TFile, TH1F, TH2F, TGraphAsymmErrors
from ROOT import gROOT, Double, TCut


def plot_unit_inpix(inputfilename = None, histofilename = "EfficiencyHistos.root", basecut = "nTelTracks == 1", matchcut="hasHit == 0", upitch=0.0, vpitch=0.0, ubins=10, vbins=10):
  
  if inputfilename == None:
    return None
  
  total_cut = TCut(basecut)
  pass_cut   = TCut(basecut) + TCut(matchcut)   
  
  rawfile = gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = TFile( inputfilename, 'READ' )
  
  histofile = gROOT.FindObject( histofilename )
  if histofile:
    histofile.Close()
  histofile = TFile( histofilename, 'RECREATE', 'Efficiency histos created from input file ' + inputfilename )
  
  # Get access to tracks 
  tree = rawfile.Get("Track")

  h_track_total = TH2F("h_track_total","",ubins,0,2*upitch,vbins,0,2*vpitch)
  tree.Draw("(v_fit-cellVCenter_fit)+%0.5f+(cellV_fit%%2)*%0.5f:(u_fit-cellUCenter_fit)+%0.5f+(cellU_fit%%2)*%0.5f >>+h_track_total" % (vpitch/2,vpitch,upitch/2,upitch),total_cut,"goff")
  h_track_total.GetXaxis().SetTitle("u_{m} [mm]")
  h_track_total.GetYaxis().SetTitle("v_{m} [mm]")
  h_track_total.GetZaxis().SetTitle("tracks")
  h_track_total.SetStats(0)

  h_track_pass = TH2F("h_track_pass","",ubins,0,2*upitch,vbins,0,2*vpitch)
  tree.Draw("(v_fit-cellVCenter_fit)+%0.5f+(cellV_fit%%2)*%0.5f:(u_fit-cellUCenter_fit)+%0.5f+(cellU_fit%%2)*%0.5f >>+h_track_pass" % (vpitch/2,vpitch,upitch/2,upitch),pass_cut,"goff")
  h_track_pass.GetXaxis().SetTitle("u_{m} [mm]")
  h_track_pass.GetYaxis().SetTitle("v_{m} [mm]")
  h_track_pass.GetZaxis().SetTitle("tracks")
  h_track_pass.SetStats(0) 

  h_efficiency =  h_track_pass
  h_efficiency.SetName("h_efficiency")
  h_efficiency.Divide(h_track_total)
  h_efficiency.SetTitle("")
  h_efficiency.SetXTitle("u_{m} [mm]")
  h_efficiency.SetYTitle("v_{m} [mm]")
  h_efficiency.SetZTitle("efficiency")
  h_efficiency.SetStats(0) 

  # write the tree into the output file and close the file
  histofile.Write()
  histofile.Close()
  rawfile.Close()


def plot_inpix(inputfilename = None, histofilename = "EfficiencyHistos.root", basecut = "nTelTracks == 1", matchcut="hasHit == 0", usize=0.0, vsize=0.0, ubins=10, vbins=10):
  
  if inputfilename == None:
    return None
  
  total_cut = TCut(basecut)
  pass_cut   = TCut(basecut) + TCut(matchcut)   
  
  rawfile = gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = TFile( inputfilename, 'READ' )
  
  histofile = gROOT.FindObject( histofilename )
  if histofile:
    histofile.Close()
  histofile = TFile( histofilename, 'RECREATE', 'Efficiency histos created from input file ' + inputfilename )
  
  # Get access to tracks 
  tree = rawfile.Get("Track")

  h_track_total = TH2F("h_track_total","",ubins,-usize,+usize,vbins,-vsize,vsize)
  tree.Draw(" v_fit:u_fit >> +h_track_total", total_cut , "goff" )
  h_track_total.GetXaxis().SetTitle("u_fit [mm]")
  h_track_total.GetYaxis().SetTitle("v_fit [mm]")
  h_track_total.GetZaxis().SetTitle("tracks")
  h_track_total.SetStats(0)
  
  h_track_pass = TH2F("h_track_pass","",ubins,-usize,usize,vbins,-vsize,+usize)
  tree.Draw(" v_fit:u_fit >> +h_track_pass" , pass_cut , "goff" )
  h_track_pass.GetXaxis().SetTitle("u_fit [mm]")
  h_track_pass.GetYaxis().SetTitle("v_fit [mm]")
  h_track_pass.GetZaxis().SetTitle("tracks")
  h_track_pass.SetStats(0) 
  
  h_efficiency =  h_track_pass
  h_efficiency.SetName("h_efficiency")
  h_efficiency.Divide(h_track_total)
  h_efficiency.SetTitle("")
  h_efficiency.SetXTitle("u_fit [mm]")
  h_efficiency.SetYTitle("v_fit [mm]")
  h_efficiency.SetZTitle("efficiency")
  h_efficiency.SetStats(0) 

  # write the tree into the output file and close the file
  histofile.Write()
  histofile.Close()
  rawfile.Close()




def plot(inputfilename = None, histofilename = "EfficiencyHistos.root", basecut = "nTelTracks == 1", matchcut="hasHit == 0", ucells=0, vcells=0 ):

  ubins = ucells
  vbins = vcells
   
  total_cut = TCut(basecut)
  pass_cut   = TCut(basecut) + TCut(matchcut)   
  
   
  if inputfilename == None:
    return None
  
  rawfile = gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = TFile( inputfilename, 'READ' )
  
  histofile = gROOT.FindObject( histofilename )
  if histofile:
    histofile.Close()
  histofile = TFile( histofilename, 'RECREATE', 'Efficiency histos created from input file ' + inputfilename )
  
  
  # Get access to tracks 
  tree = rawfile.Get("Track")

  # Compute efficiency for v cells 
  h_track_v_total = TH1F("h_track_v_total","",vbins,0,vcells)
  tree.Draw("cellV_fit >> +h_track_v_total", total_cut  ,  "goff")
  h_track_v_total.GetXaxis().SetTitle("cellV_fit [cellID]")
  h_track_v_total.GetYaxis().SetTitle("tracks") 
 
  h_track_v_pass = TH1F("h_track_v_pass","",vbins,0,vcells)
  tree.Draw("cellV_fit >> +h_track_v_pass"  ,  pass_cut , "goff" )
  h_track_v_pass.GetXaxis().SetTitle("cellV_fit [cellID]")
  h_track_v_pass.GetYaxis().SetTitle("tracks") 
   
  g_efficiency_v = TGraphAsymmErrors()
  g_efficiency_v.SetName("g_efficiency_v")
  g_efficiency_v.Divide(h_track_v_pass,h_track_v_total) 
  g_efficiency_v.SetTitle("")
  g_efficiency_v.GetXaxis().SetTitle("cellV_fit [cellID]")
  g_efficiency_v.GetYaxis().SetTitle("efficiency") 
  g_efficiency_v.Write()

  # Compute efficiency for u cells 
  h_track_u_total = TH1F("h_track_u_total","",ubins,0,ucells)
  tree.Draw("cellU_fit >> +h_track_u_total", total_cut ,  "goff")
  h_track_u_total.GetXaxis().SetTitle("cellU_fit [cellID]")
  h_track_u_total.GetYaxis().SetTitle("tracks")  
  

  h_track_u_pass = TH1F("h_track_u_pass","",ubins,0,ucells)
  tree.Draw("cellU_fit >> +h_track_u_pass"  , pass_cut , "goff" )
  h_track_u_pass.GetXaxis().SetTitle("cellU_fit [cellID]")
  h_track_u_pass.GetYaxis().SetTitle("tracks")  
   
  g_efficiency_u = TGraphAsymmErrors()
  g_efficiency_u.SetName("g_efficiency_u")
  g_efficiency_u.Divide(h_track_u_pass,h_track_u_total) 
  g_efficiency_u.SetTitle("")
  g_efficiency_u.GetXaxis().SetTitle("cellU_fit [cellID]")
  g_efficiency_u.GetYaxis().SetTitle("efficiency") 
  g_efficiency_u.Write()
  
  # Compute efficiency for u:v 
  h_track_total = TH2F("h_track_total","",ubins,0,ucells,vbins,0,vcells)
  tree.Draw("cellV_fit:cellU_fit >> +h_track_total", total_cut  ,  "goff")
  h_track_total.GetXaxis().SetTitle("cellU_fit [cellID]")
  h_track_total.GetYaxis().SetTitle("cellV_fit [cellID]")
  h_track_total.GetZaxis().SetTitle("tracks")
 
  h_track_pass = TH2F("h_track_pass","",ubins,0,ucells,vbins,0,vcells)
  tree.Draw("cellV_fit:cellU_fit >> +h_track_pass"  ,  pass_cut , "goff" )
  h_track_pass.GetXaxis().SetTitle("cellU_fit [cellID]")
  h_track_pass.GetYaxis().SetTitle("cellV_fit [cellID]")
  h_track_pass.GetZaxis().SetTitle("tracks")

  g_efficiency = h_track_pass
  g_efficiency.SetName("g_efficiency")
  g_efficiency.Divide(h_track_pass,h_track_total)
  g_efficiency.SetTitle("")
  g_efficiency.GetXaxis().SetTitle("cellU_fit [cellID]")
  g_efficiency.GetYaxis().SetTitle("cellV_fit [cellID]")
  g_efficiency.GetZaxis().SetTitle("efficiency") 
  g_efficiency.SetStats(0)
  g_efficiency.Write()

 
  # write the tree into the output file and close the file
  histofile.Write()
  histofile.Close()
  rawfile.Close()


