from ROOT import TFile, TH1F, TH2F, TGraphAsymmErrors
from ROOT import gROOT, Double, TCut


NBINSU = 1152
NBINSV = 576
First = 25000 
Last = 100000

EVENTS = Last - First

basecut = TCut("hasTrack == -1 && iEvt > %d && iEvt < %d "% (First, Last))

 
def plot_occupancy(inputfilename = None, histofilename = "EfficiencyHistos.root"):
  
  if inputfilename == None:
    return None
  
  rawfile = gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = TFile( inputfilename, 'READ' )
  
  histofile = gROOT.FindObject( histofilename )
  if histofile:
    histofile.Close()
  histofile = TFile( histofilename, 'RECREATE', 'Efficiency plots created from input file ' + inputfilename )
  
  
  # Get access to hits 
  tree = rawfile.Get("Hit")
   
  # Compute occupancy for v cells 
   
  h_occupancy_v = TH1F("h_occupancy_v","",NBINSV,0,NBINSV)
  tree.Draw("cellV_hit >> +h_occupancy_v", basecut ,  "goff")
  h_occupancy_v.Scale(1.0/(NBINSU*EVENTS ) ) 
  h_occupancy_v.SetTitle("")
  h_occupancy_v.GetXaxis().SetTitle("cellV [cellID]")
  h_occupancy_v.GetYaxis().SetTitle("occupancy") 
 
  # Compute occupancy for u cells 
   
  h_occupancy_u = TH1F("h_occupancy_u","",NBINSU,0,NBINSU)
  tree.Draw("cellU_hit >> +h_occupancy_u", basecut ,  "goff")
  h_occupancy_u.Scale(1.0/(NBINSV*EVENTS ) ) 
  h_occupancy_u.SetTitle("")
  h_occupancy_u.GetXaxis().SetTitle("cellU [cellID]")
  h_occupancy_u.GetYaxis().SetTitle("occupancy") 

  # Compute occupancy for u:v cells 
   
  h_occupancy = TH2F("h_occupancy","",NBINSU,0,NBINSU,NBINSV,0,NBINSV)
  tree.Draw("cellV_hit:cellU_hit >> +h_occupancy", basecut ,  "goff")
  h_occupancy.Scale(1.0/( 1.0*EVENTS ) ) 
  h_occupancy.SetTitle("")
  h_occupancy.GetXaxis().SetTitle("cellU [cellID]")
  h_occupancy.GetYaxis().SetTitle("cellV [cellID]")
  h_occupancy.GetYaxis().SetTitle("occupancy") 

 
  # write the tree into the output file and close the file
  histofile.Write()
  histofile.Close()
  rawfile.Close()

if __name__ == '__main__':
  plot_occupancy(inputfilename = "root-files/Strip-ID6-Histos-run000282-merger-reco.root", histofilename = "OccupancyHistos.root")

