from ROOT import TFile, TH1F, TH2F, TGraphAsymmErrors
from ROOT import gROOT, Double, TCut


def plot_super_inpix(inputfile=None, histofile=None, basecut="", matchcut="hasHit==0", upitch=0.0, vpitch=0.0, ubins=10, vbins=10):
  
  if inputfile == None or histofile==None:
    return None
  
  total_cut = TCut(basecut)
  pass_cut   = TCut(basecut) + TCut(matchcut)   
  
  # Get access to tracks 
  tree = inputfile.Get("Track")

  h_super_track_total = TH2F("h_super_track_total","",ubins,0,2*upitch,vbins,0,2*vpitch)
  tree.Draw("(v_fit-cellVCenter_fit)+%0.5f+(cellV_fit%%2)*%0.5f:(u_fit-cellUCenter_fit)+%0.5f+(cellU_fit%%2)*%0.5f >>+h_super_track_total" % (vpitch/2,vpitch,upitch/2,upitch),total_cut,"goff")
  h_super_track_total.GetXaxis().SetTitle("u_{m} [mm]")
  h_super_track_total.GetYaxis().SetTitle("v_{m} [mm]")
  h_super_track_total.GetZaxis().SetTitle("tracks")
  h_super_track_total.SetStats(0)

  h_super_track_pass = TH2F("h_super_track_pass","",ubins,0,2*upitch,vbins,0,2*vpitch)
  tree.Draw("(v_fit-cellVCenter_fit)+%0.5f+(cellV_fit%%2)*%0.5f:(u_fit-cellUCenter_fit)+%0.5f+(cellU_fit%%2)*%0.5f >>+h_super_track_pass" % (vpitch/2,vpitch,upitch/2,upitch),pass_cut,"goff")
  h_super_track_pass.GetXaxis().SetTitle("u_{m} [mm]")
  h_super_track_pass.GetYaxis().SetTitle("v_{m} [mm]")
  h_super_track_pass.GetZaxis().SetTitle("tracks")
  h_super_track_pass.SetStats(0) 

  h2_super_efficiencymap =  h_super_track_pass
  h2_super_efficiencymap.SetName("h2_super_efficiencymap")
  h2_super_efficiencymap.Divide(h_super_track_total)
  h2_super_efficiencymap.SetTitle("SuperPixel inpixel efficiency")
  h2_super_efficiencymap.SetXTitle("u_{m} [mm]")
  h2_super_efficiencymap.SetYTitle("v_{m} [mm]")
  h2_super_efficiencymap.SetZTitle("efficiency")
  h2_super_efficiencymap.SetStats(0) 
  h2_super_efficiencymap.Write()
  

def plot_inpix(inputfile=None, histofile=None, basecut="", matchcut="hasHit==0", usize=0.0, vsize=0.0, ubins=10, vbins=10):
  
  if inputfile == None or histofile==None:
    return None
  
  total_cut = TCut(basecut)
  pass_cut   = TCut(basecut) + TCut(matchcut)   
  
  # Get access to tracks 
  tree = inputfile.Get("Track")

  h2_inpix_track_total = TH2F("h2_inpix_track_total","",ubins,-usize,+usize,vbins,-vsize,vsize)
  tree.Draw(" v_fit:u_fit >> +h2_inpix_track_total", total_cut , "goff" )
  h2_inpix_track_total.GetXaxis().SetTitle("u_fit [mm]")
  h2_inpix_track_total.GetYaxis().SetTitle("v_fit [mm]")
  h2_inpix_track_total.GetZaxis().SetTitle("tracks")
  h2_inpix_track_total.SetStats(0)
  
  h2_inpix_pass_total = TH2F("h2_inpix_pass_total","",ubins,-usize,usize,vbins,-vsize,+usize)
  tree.Draw(" v_fit:u_fit >> +h2_inpix_pass_total" , pass_cut , "goff" )
  h2_inpix_pass_total.GetXaxis().SetTitle("u_fit [mm]")
  h2_inpix_pass_total.GetYaxis().SetTitle("v_fit [mm]")
  h2_inpix_pass_total.GetZaxis().SetTitle("tracks")
  h2_inpix_pass_total.SetStats(0) 
  
  h2_inpix_efficiencymap =  h2_inpix_pass_total
  h2_inpix_efficiencymap.SetName("h2_inpix_efficiencymap")
  h2_inpix_efficiencymap.Divide(h2_inpix_track_total)
  h2_inpix_efficiencymap.SetTitle("")
  h2_inpix_efficiencymap.SetXTitle("u_fit [mm]")
  h2_inpix_efficiencymap.SetYTitle("v_fit [mm]")
  h2_inpix_efficiencymap.SetZTitle("efficiency")
  h2_inpix_efficiencymap.SetStats(0)
  h2_inpix_efficiencymap.Write()
   

def plot(inputfile=None, histofile=None, basecut="", matchcut="hasHit==0", uaxis=(250,0,250), vaxis=(768,0,768)):
     
  if inputfile == None or histofile==None:
    return None
  
  total_cut = TCut(basecut)
  pass_cut   = TCut(basecut) + TCut(matchcut)   
  
  # Get access to tracks 
  tree = inputfile.Get("Track")

  # Number of tracks 
  h_nTelTracks = TH1F("h_nTelTracks","",50,0,50)
  tree.Draw("nTelTracks >> +h_nTelTracks", ""  ,  "goff")
  h_nTelTracks.SetTitle("Number of tracks on DUT matrix") 
  h_nTelTracks.GetXaxis().SetTitle("number of telescope tracks")
  h_nTelTracks.GetYaxis().SetTitle("events") 

  # Number of tracks with reference hit
  h_nTelTracksWithRef = TH1F("h_nTelTracksWithRef","",50,0,50)
  tree.Draw("nTelTracks >> +h_nTelTracksWithRef", "hasRefHit==0"  ,  "goff")
  h_nTelTracksWithRef.SetTitle("Number of tracks with timing hit on DUT matrix") 
  h_nTelTracksWithRef.GetXaxis().SetTitle("number of telescope tracks with timing hit")
  h_nTelTracksWithRef.GetYaxis().SetTitle("events") 

  # Compute efficiency for v cells 
  h_track_v_total = TH1F("h_track_v_total","",vaxis[0],vaxis[1],vaxis[2])
  tree.Draw("cellV_fit >> +h_track_v_total", total_cut  ,  "goff")
  h_track_v_total.GetXaxis().SetTitle("cellV_fit [cellID]")
  h_track_v_total.GetYaxis().SetTitle("tracks") 
 
  h_track_v_pass = TH1F("h_track_v_pass","",vaxis[0],vaxis[1],vaxis[2])
  tree.Draw("cellV_fit >> +h_track_v_pass"  ,  pass_cut , "goff" )
  h_track_v_pass.GetXaxis().SetTitle("cellV_fit [cellID]")
  h_track_v_pass.GetYaxis().SetTitle("tracks") 
   
  g_efficiency_v = TGraphAsymmErrors()
  g_efficiency_v.SetName("g_efficiency_v")
  g_efficiency_v.Divide(h_track_v_pass,h_track_v_total) 
  g_efficiency_v.SetTitle("Hit efficiency profile")
  g_efficiency_v.GetXaxis().SetTitle("cellV_{fit} [cellID]")
  g_efficiency_v.GetYaxis().SetTitle("efficiency") 
  g_efficiency_v.Write()

  # Compute efficiency for u cells 
  h_track_u_total = TH1F("h_track_u_total","",uaxis[0],uaxis[1],uaxis[2])
  tree.Draw("cellU_fit >> +h_track_u_total", total_cut ,  "goff")
  h_track_u_total.GetXaxis().SetTitle("cellU_fit [cellID]")
  h_track_u_total.GetYaxis().SetTitle("tracks")  
  

  h_track_u_pass = TH1F("h_track_u_pass","",uaxis[0],uaxis[1],uaxis[2])
  tree.Draw("cellU_fit >> +h_track_u_pass"  , pass_cut , "goff" )
  h_track_u_pass.GetXaxis().SetTitle("cellU_fit [cellID]")
  h_track_u_pass.GetYaxis().SetTitle("tracks")  
   
  g_efficiency_u = TGraphAsymmErrors()
  g_efficiency_u.SetName("g_efficiency_u")
  g_efficiency_u.Divide(h_track_u_pass,h_track_u_total) 
  g_efficiency_u.SetTitle("Hit efficiency profile")
  g_efficiency_u.GetXaxis().SetTitle("cellU_{fit} [cellID]")
  g_efficiency_u.GetYaxis().SetTitle("efficiency") 
  g_efficiency_u.Write()
  
  # Compute efficiency for u:v 
  h_track_total = TH2F("h_track_total","",uaxis[0],uaxis[1],uaxis[2],vaxis[0],vaxis[1],vaxis[2])
  tree.Draw("cellV_fit:cellU_fit >> +h_track_total", total_cut  ,  "goff")
  h_track_total.GetXaxis().SetTitle("cellU_fit [cellID]")
  h_track_total.GetYaxis().SetTitle("cellV_fit [cellID]")
  h_track_total.GetZaxis().SetTitle("tracks")
 
  h_track_pass = TH2F("h_track_pass","",uaxis[0],uaxis[1],uaxis[2],vaxis[0],vaxis[1],vaxis[2])
  tree.Draw("cellV_fit:cellU_fit >> +h_track_pass"  ,  pass_cut , "goff" )
  h_track_pass.GetXaxis().SetTitle("cellU_fit [cellID]")
  h_track_pass.GetYaxis().SetTitle("cellV_fit [cellID]")
  h_track_pass.GetZaxis().SetTitle("tracks")

  h2_efficiencymap = h_track_pass      
  h2_efficiencymap.SetName("h2_efficiencymap")
  h2_efficiencymap.Divide(h_track_pass,h_track_total)
  h2_efficiencymap.SetTitle("Hit efficiency map")
  h2_efficiencymap.GetXaxis().SetTitle("cellU_{fit} [cellID]")
  h2_efficiencymap.GetYaxis().SetTitle("cellV_{fit} [cellID]")
  h2_efficiencymap.GetZaxis().SetTitle("efficiency") 
  h2_efficiencymap.SetStats(0)
  h2_efficiencymap.Write()

  
