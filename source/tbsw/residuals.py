from ROOT import TFile, TH1F, TH2F
from ROOT import gROOT, Double, TCut



def make_pdf(histofile, pdfName='plots.pdf'):

    
  if histofile == None:
    return  
  
  import ROOT

  c1 = ROOT.TCanvas("c1","",10,10,1100,700)
  c1.SetRightMargin(0.2)
  c1.Print(pdfName+"(","pdf")
  
  for key in histofile.GetListOfKeys():
         
    cl = gROOT.GetClass(key.GetClassName())
    if cl.InheritsFrom("TH1"):
      
      print("Printing histo " + key.GetName())
       
      h1 = key.ReadObj().Clone()
      c1.Clear()
      c1.cd()
      c1.SetName(key.GetName())
      c1.SetTitle(key.GetTitle())
      
      if  cl.InheritsFrom("TH2"): 
        h1.Draw("colz")
      else: 
        h1.Draw()
      
      ROOT.gPad.Modified()
      ROOT.gPad.Update()   
      c1.Print(pdfName+"(","pdf")

    elif cl.InheritsFrom("TGraph"):
      
      print("Printing graph " + key.GetName())
       
      g1 = key.ReadObj().Clone()
      c1.Clear()
      c1.cd()
      c1.SetName(key.GetName())
      c1.SetTitle(key.GetTitle())
      g1.Draw()
      
      ROOT.gPad.Modified()
      ROOT.gPad.Update()   
      c1.Print(pdfName+"(","pdf")

  
  c1.Print(pdfName+")","pdf")




def plot(inputfile=None, histofile=None, basecut="hasTrack==0", Config=None):
  
  if inputfile == None or histofile == None:
    return

  if Config==None:    
    print('Missing Config objects for plotting')
    return 
  
  # Analysis cuts 
  hitcut = TCut(basecut)
  
  # Get access to tracks 
  tracktree = inputfile.Get("Track")
  
  nucells = Config['ucell_axis'][0]
  minucell = Config['ucell_axis'][1]
  maxucell = Config['ucell_axis'][2]
  nvcells = Config['vcell_axis'][0]
  minvcell = Config['vcell_axis'][1]
  maxvcell = Config['vcell_axis'][2]

  htrackspot = TH2F("htrackspot","htrackspot",nucells,minucell,maxucell, nvcells, minvcell, maxvcell)
  tracktree.Draw("cellV_fit:cellU_fit>>+htrackspot", "", "goff")
  htrackspot.SetTitle("Track intersections on DUT matrix") 
  htrackspot.GetXaxis().SetTitle("cellU_{fit} [cell ID]")
  htrackspot.GetYaxis().SetTitle("cellV_{fit} [cell ID]")
  htrackspot.GetZaxis().SetTitle("number of tracks")
  htrackspot.SetStats(0)
  htrackspot.Write()  

  nu = Config['sensor_u_axis'][0]
  minu = Config['sensor_u_axis'][1]
  maxu = Config['sensor_u_axis'][2]
  nv = Config['sensor_v_axis'][0]
  minv = Config['sensor_v_axis'][1]
  maxv = Config['sensor_v_axis'][2]

  hbeamspot = TH2F("hbeamspot","hbeamspot",nu,minu,maxu,nv,minv,maxv)
  tracktree.Draw("v_fit:u_fit>>+hbeamspot", "","goff")
  hbeamspot.SetTitle("Track intersections on DUT plane")
  hbeamspot.GetXaxis().SetTitle("u_{fit} [mm]")
  hbeamspot.GetYaxis().SetTitle("v_{fit} [mm]")
  hbeamspot.GetZaxis().SetTitle("number of tracks")
  hbeamspot.SetStats(0)
  hbeamspot.Write()
  
  #Get access to hits
  hittree = inputfile.Get("Hit")

  hspot = TH2F("hspot","",nucells,minucell,maxucell, nvcells, minvcell, maxvcell)
  hittree.Draw("cellV_hit:cellU_hit>>hspot", hitcut,"goff")
  hspot.SetTitle("DUT clusters on matrix")
  hspot.GetXaxis().SetTitle("cellU_{hit} [cell ID]")
  hspot.GetYaxis().SetTitle("cellV_{hit} [cell ID]")
  hspot.GetZaxis().SetTitle("number of clusters")
  hspot.SetStats(0)
  hspot.Write()


  nbinsu = Config['residual_u_axis'][0]
  minu = 1000*Config['residual_u_axis'][1]
  maxu = 1000*Config['residual_u_axis'][2]

  hres_u = TH1F("hres_u","",nbinsu,minu,maxu)
  hittree.Draw("(u_hit - u_fit)*1000 >>+hres_u", hitcut,"goff")
  hres_u.SetTitle("Unbiased DUT residuals u_{hit}-u_{fit}")
  hres_u.GetXaxis().SetTitle("u_{hit} - u_{fit} [#mum]")
  hres_u.GetYaxis().SetTitle("number of hits")
  hres_u.GetYaxis().SetTitleOffset(1.2)
  hres_u.Write()  

  nbinsv = Config['residual_u_axis'][0]
  minv = 1000*Config['residual_u_axis'][1]
  maxv = 1000*Config['residual_u_axis'][2]

  hres_v = TH1F("hres_v","",nbinsv,minv,maxv)
  hittree.Draw("(v_hit - v_fit)*1000 >>+hres_v", hitcut,"goff")
  hres_v.SetTitle("Unbiased DUT residuals v_{hit}-v_{fit}")
  hres_v.GetXaxis().SetTitle("v_{hit} - v_{fit} [#mum]")
  hres_v.GetYaxis().SetTitle("number of hits")
  hres_v.GetYaxis().SetTitleOffset(1.2)
  hres_v.Write()  

  hres_uv = TH2F("hres_uv","",nbinsu,minu,maxu, nbinsv,minv,maxv )
  hittree.Draw("(v_hit - v_fit)*1000:(u_hit - u_fit)*1000 >>+hres_uv", hitcut,"goff")  
  hres_uv.SetTitle("Unbiased 2D DUT residuals")
  hres_uv.GetXaxis().SetTitle("u_{hit} - u_{fit} [#mum]")
  hres_uv.GetYaxis().SetTitle("v_{hit} - v_{fit} [#mum]")
  hres_uv.GetZaxis().SetTitle("number of hits")
  hres_uv.SetStats(0)
  hres_uv.Write()  

  hchisqundof = TH1F("hchisqundof","",100,0,0)
  hittree.Draw("trackChi2 / trackNdof >> +hchisqundof", hitcut,"goff")
  hchisqundof.SetTitle("#chi^{2}/ndof")
  hchisqundof.GetXaxis().SetTitle("#chi^{2}/ndof")
  hchisqundof.GetYaxis().SetTitle("number of tracks")
  hchisqundof.GetYaxis().SetTitleOffset(1.2)
  hchisqundof.Write()  

  nseed = Config['seed_charge_axis'][0]
  seedmin = Config['seed_charge_axis'][1]
  seedmax = Config['seed_charge_axis'][2]  
  
  hCharge = TH1F("hCharge","",nseed,seedmin,seedmax)
  hittree.Draw("seedCharge >>+hCharge", hitcut,"goff")
  hCharge.SetTitle("Seed Charge")
  hCharge.GetXaxis().SetTitle("seed charge [{}]".format(Config['charge_unit']))
  hCharge.GetYaxis().SetTitle("number of hits")
  hCharge.GetYaxis().SetTitleOffset(1.2)
  hCharge.Write()

  nclu = Config['clus_charge_axis'][0]
  clumin = Config['clus_charge_axis'][1]
  clumax = Config['clus_charge_axis'][2]  

  hCCharge = TH1F("hCCharge","",nclu,clumin,clumax)
  hittree.Draw("clusterCharge >>+hCCharge", hitcut,"goff")
  hCCharge.SetTitle("Cluster Charge")
  hCCharge.GetXaxis().SetTitle("cluster charge [{}]".format(Config['charge_unit']))
  hCCharge.GetYaxis().SetTitle("number of hits")
  hCCharge.GetYaxis().SetTitleOffset(1.2)
  hCCharge.Write()

  hsize = TH1F("hsize","",10,0,10)
  hittree.Draw("size >>+hsize", hitcut,"goff")
  hsize.SetTitle("Cluster size")
  hsize.GetXaxis().SetTitle("cluster size [pixels]")
  hsize.GetYaxis().SetTitle("number of hits")
  hsize.GetYaxis().SetTitleOffset(1.2)
  hsize.Write()

  hsizeU = TH1F("hsizeU","",10,0,10)
  hittree.Draw("sizeU >>+hsizeU", hitcut,"goff")
  hsizeU.SetTitle("Cluster size projection")
  hsizeU.GetXaxis().SetTitle("sizeU [cell ID]")
  hsizeU.GetYaxis().SetTitle("number of hits")
  hsizeU.GetYaxis().SetTitleOffset(1.2)
  hsizeU.Write()

  hsizeV = TH1F("hsizeV","",10,0,10)
  hittree.Draw("sizeV >>+hsizeV", hitcut,"goff")
  hsizeV.SetTitle("Cluster size projection")
  hsizeV.GetXaxis().SetTitle("sizeV [cell ID]")
  hsizeV.GetYaxis().SetTitle("number of hits")
  hsizeV.GetYaxis().SetTitleOffset(1.2)
  hsizeV.Write() 



