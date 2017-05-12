from ROOT import TFile, TH1F, TH2F, TGraphAsymmErrors
from ROOT import gROOT, Double, TCut

BINSU = 256
BINSV = 128

MAXU = 1152
MAXV = 576

 
# Analysis cuts 
basecut = TCut("trackNHits  == 6 && nTelTracks == 1")
matched = TCut("hasHit == 0 && seedCharge > 0") 
u_cut = TCut("cellU_fit >= 0 && cellU_fit < %d" % MAXU)
v_cut = TCut("cellV_fit >= 0 && cellV_fit < %d" % MAXV) 


def plot_efficiency(inputfilename = None, histofilename = "EfficiencyHistos.root"):
  
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
  
  
  # Get access to tracks 
  tree = rawfile.Get("Track")

  # Compute efficiency for rows 
  
  h_track_row_total = TH1F("h_track_row_total","",BINSV,0,MAXV)
  tree.Draw("cellV_fit >> +h_track_row_total", basecut + u_cut  ,  "goff")
 
  h_track_row_pass = TH1F("h_track_row_pass","",BINSV,0,MAXV)
  tree.Draw("cellV_fit >> +h_track_row_pass"  ,  basecut + u_cut + matched , "goff" )
   
  g_effi_row = TGraphAsymmErrors()
  g_effi_row.Divide(h_track_row_pass,h_track_row_total) 
  
  g_effi_row.SetTitle("")
  g_effi_row.GetXaxis().SetTitle("cellV [cellID]")
  g_effi_row.GetYaxis().SetTitle("efficiency") 
  g_effi_row.Write()

  # Compute efficiency for columns 

  h_track_col_total = TH1F("h_track_col_total","",BINSU,0,MAXU)
  tree.Draw("cellU_fit >> +h_track_col_total", basecut + v_cut  ,  "goff")
 
  h_track_col_pass = TH1F("h_track_col_pass","",BINSU,0,MAXU)
  tree.Draw("cellU_fit >> +h_track_col_pass"  ,  basecut + v_cut + matched , "goff" )
   
  g_effi_col = TGraphAsymmErrors()
  g_effi_col.Divide(h_track_col_pass,h_track_col_total) 
  
  g_effi_col.SetTitle("")
  g_effi_col.GetXaxis().SetTitle("cellU [cellID]")
  g_effi_col.GetYaxis().SetTitle("efficiency") 
  g_effi_col.Write()
  
  # Compute efficiency for cols:rows 

  h_track_total = TH2F("h_track_total","",BINSU,0,MAXU,BINSV,0,MAXV)
  tree.Draw("cellV_fit:cellU_fit >> +h_track_total", basecut  ,  "goff")
 
  h_track_pass = TH2F("h_track_pass","",BINSU,0,MAXU,BINSV,0,MAXV)
  tree.Draw("cellV_fit:cellU_fit >> +h_track_pass"  ,  basecut  + matched , "goff" )

  g_effi = h_track_pass
  g_effi.Divide(h_track_pass,h_track_total,1.0,1.0,"B")
   
  g_effi.SetTitle("")
  g_effi.GetXaxis().SetTitle("cellU [cellID]")
  g_effi.GetYaxis().SetTitle("cellV [cellID]")
  g_effi.GetZaxis().SetTitle("efficiency") 
  g_effi.Write()

 
  # write the tree into the output file and close the file
  histofile.Write()
  histofile.Close()
  rawfile.Close()

if __name__ == '__main__':
  plot_efficiency(inputfilename = "root-files/Strip-ID6-Histos-run000282-merger-reco.root", histofilename = "EfficiencyHistos.root")

