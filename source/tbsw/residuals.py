from ROOT import TFile, TH1F, TH2F
from ROOT import gROOT, Double, TCut
import math



def make_pdf(histofile, pdfName='plots.pdf'):

    
  if histofile == None:
    return  
  
  import ROOT

  c1 = ROOT.TCanvas("c1","",10,10,1100,700)
  c1.SetRightMargin(0.2)
  
  
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
  
  def make_residual_histo(nbins, xmin, xmax, size_cut='sizeU>0', axis='u', label='all'):
    """ Helper function for residual making residual histo."""
    histo = TH1F("hres_{:s}_{:s}".format(axis, label),"",nbins,xmin,xmax)
    hittree.Draw("({:s}_hit - {:s}_fit)*1000 >>+hres_{:s}_{:s}".format(axis, axis, axis, label), hitcut+TCut(size_cut),"goff")
    histo.SetTitle("Unbiased DUT residuals {:s}_{{hit}}-{:s}_{{fit}} for {:s} clusters".format(axis, axis, label))
    histo.GetXaxis().SetTitle("{:s}_{{hit}} - {:s}_{{fit}} [#mum]".format(axis,axis))
    histo.GetYaxis().SetTitle("number of clusters")
    histo.GetYaxis().SetTitleOffset(1.2)
    histo.Write()
    return histo 
  
  def make_resolution_histo(residual_histos, tel_sigma=0, axis='u'):
    """ Helper function for making resolution histo."""
    histo = TH1F("hpoint_resolution_{:s}".format(axis),"",len(residual_histos),0 ,len(residual_histos))
    histo.SetTitle("DUT pointing resolution {:s} overview".format(axis))
    histo.GetXaxis().SetTitle("")
    histo.GetYaxis().SetTitle("DUT pointing resolution #sigma_{:s} [#mum]".format(axis))
    histo.GetYaxis().SetTitleOffset(1.2)    
    histo.SetStats(0)     

    bin=1
    for key, residual_histo in sorted(residual_histos.items(), key = lambda kv: kv[0][-1] )    :
      histo.GetXaxis().SetBinLabel(bin, key)
      if residual_histo.GetRMS()**2 - tel_sigma**2 > 0:
        dut_sigma = math.sqrt( residual_histo.GetRMS()**2 - tel_sigma**2 )
        dut_sigma_sigma = residual_histo.GetRMSError()
        histo.SetBinContent(bin, dut_sigma) 
        histo.SetBinError(bin, dut_sigma_sigma) 
      bin+=1
    
    histo.Write()
    return histo   
    
  nbinsu = Config['residual_u_axis'][0]
  minu = 1000*Config['residual_u_axis'][1]
  maxu = 1000*Config['residual_u_axis'][2]
  
  hfit_sigma_u = TH1F("hfit_sigma_u","hfit_sigma_u",200,0,maxu)
  hittree.Draw("u_fiterr*1000 >>+hfit_sigma_u", hitcut,"goff")
  hfit_sigma_u.SetTitle("Track intersection uncertainty #sigma_{u}")
  hfit_sigma_u.GetXaxis().SetTitle("#sigma_{u_{fit}} [#mum]")
  hfit_sigma_u.GetYaxis().SetTitle("number of tracks")
  hfit_sigma_u.GetYaxis().SetTitleOffset(1.2)
  hfit_sigma_u.Write()  

  hfit_cov_u = TH1F("hfit_cov_u","hfit_cov_u",200,0,maxu**2)
  hittree.Draw("u_fiterr*u_fiterr*1000*1000 >>+hfit_cov_u", hitcut,"goff")
  hfit_cov_u.SetTitle("Track intersection uncertainty #sigma_{u}^{2}")
  hfit_cov_u.GetXaxis().SetTitle("#sigma_{u_{fit}}^{2} [#mum^{2}]")
  hfit_cov_u.GetYaxis().SetTitle("number of tracks")
  hfit_cov_u.GetYaxis().SetTitleOffset(1.2)
  hfit_cov_u.Write()  
  
  residual_histos_u = {}
  residual_histos_u['sizeU>0'] = make_residual_histo(nbinsu, minu, maxu, size_cut='sizeU>0', axis='u', label='all') 
  residual_histos_u['sizeU==1'] = make_residual_histo(nbinsu, minu, maxu, size_cut='sizeU==1', axis='u', label='sizeU==1') 
  residual_histos_u['sizeU==2'] = make_residual_histo(nbinsu, minu, maxu, size_cut='sizeU==2', axis='u', label='sizeU==2') 
  #residual_histos_u['sizeU==3'] = make_residual_histo(nbinsu, minu, maxu, size_cut='sizeU==3', axis='u', label='sizeU==3') 
  
  make_resolution_histo(residual_histos_u, tel_sigma=math.sqrt(hfit_cov_u.GetMean()), axis='u')

  nbinsv = Config['residual_u_axis'][0]
  minv = 1000*Config['residual_u_axis'][1]
  maxv = 1000*Config['residual_u_axis'][2]

  hfit_sigma_v = TH1F("hfit_sigma_v","hfit_sigma_v",200,0,maxv)
  hittree.Draw("v_fiterr*1000 >>+hfit_sigma_v", hitcut,"goff")
  hfit_sigma_v.SetTitle("Track intersection uncertainty #sigma_{v}")
  hfit_sigma_v.GetXaxis().SetTitle("#sigma_{v_{fit}} [#mum]")
  hfit_sigma_v.GetYaxis().SetTitle("number of tracks")
  hfit_sigma_v.GetYaxis().SetTitleOffset(1.2)
  hfit_sigma_v.Write()  

  hfit_cov_v = TH1F("hfit_cov_v","hfit_cov_v",200,0,maxv**2)
  hittree.Draw("v_fiterr*v_fiterr*1000*1000 >>+hfit_cov_v", hitcut,"goff")
  hfit_cov_v.SetTitle("Track intersection uncertainty #sigma_{v}^{2}")
  hfit_cov_v.GetXaxis().SetTitle("#sigma_{v_{fit}}^{2} [#mum^{2}]")
  hfit_cov_v.GetYaxis().SetTitle("number of tracks")
  hfit_cov_v.GetYaxis().SetTitleOffset(1.2)
  hfit_cov_v.Write()  
  
  residual_histos_v = {}
  residual_histos_v['sizeV>0'] = make_residual_histo(nbinsv, minv, maxv, size_cut='sizeV>0', axis='v', label='all') 
  residual_histos_v['sizeV==1'] = make_residual_histo(nbinsv, minv, maxv, size_cut='sizeV==1', axis='v', label='sizeV==1') 
  residual_histos_v['sizeV==2'] = make_residual_histo(nbinsv, minv, maxv, size_cut='sizeV==2', axis='v', label='sizeV==2') 
  #residual_histos_v['sizeV==3'] = make_residual_histo(nbinsv, minv, maxv, size_cut='sizeV==3', axis='v', label='sizeV==3') 
   
  make_resolution_histo(residual_histos_v, tel_sigma=math.sqrt(hfit_cov_v.GetMean()), axis='v')
  
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



