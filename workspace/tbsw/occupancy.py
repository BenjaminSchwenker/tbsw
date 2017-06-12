from ROOT import TFile, TH1F, TH2F, TGraphAsymmErrors
from ROOT import gROOT, Double, TCut

def plot(inputfilename = None, histofilename = "OccupancyHistos.root", ucells=0,vcells=0):
  
  if inputfilename == None:
    return None
  
  rawfile = gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = TFile( inputfilename, 'READ' )
  
  histofile = gROOT.FindObject( histofilename )
  if histofile:
    histofile.Close()
  histofile = TFile( histofilename, 'RECREATE', 'Occupancy histos created from input file ' + inputfilename )
  
  
  # Get access to hits 
  tree = rawfile.Get("Hit")

  # Only count hits not matched to a track
  basecut = TCut("hasTrack == -1")

  # Number of events 
  events = rawfile.Get("Event").GetEntries()
   
  # Compute occupancy for v cells 
  h_occupancy_v = TH1F("h_occupancy_v","",vcells,0,vcells)
  tree.Draw("cellV_hit >> +h_occupancy_v", basecut ,  "goff")
  h_occupancy_v.Scale(1.0/(ucells*events ) ) 
  h_occupancy_v.SetTitle("")
  h_occupancy_v.GetXaxis().SetTitle("cellV [cellID]")
  h_occupancy_v.GetYaxis().SetTitle("occupancy") 
 
  # Compute occupancy for u cells 
  h_occupancy_u = TH1F("h_occupancy_u","",ucells,0,ucells)
  tree.Draw("cellU_hit >> +h_occupancy_u", basecut ,  "goff")
  h_occupancy_u.Scale(1.0/(vcells*events ) ) 
  h_occupancy_u.SetTitle("")
  h_occupancy_u.GetXaxis().SetTitle("cellU [cellID]")
  h_occupancy_u.GetYaxis().SetTitle("occupancy") 

  # Compute occupancy for u:v cells 
  h_occupancy = TH2F("h_occupancy","",ucells,0,ucells,vcells,0,vcells)
  tree.Draw("cellV_hit:cellU_hit >> +h_occupancy", basecut ,  "goff")
  h_occupancy.Scale(1.0/( 1.0*events ) ) 
  h_occupancy.SetTitle("")
  h_occupancy.GetXaxis().SetTitle("cellU [cellID]")
  h_occupancy.GetYaxis().SetTitle("cellV [cellID]")
  h_occupancy.GetYaxis().SetTitle("occupancy") 

 
  # write the tree into the output file and close the file
  histofile.Write()
  histofile.Close()
  rawfile.Close()



