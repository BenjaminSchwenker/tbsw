from ROOT import TFile, TH1F, TH2F
from ROOT import gROOT, Double, TCut
import numpy 
import math


def plot(inputfilename = None, histofilename = "InpixHistos.root", usize=0.0, vsize=0.0, ubins=10, vbins=10):
    
  if inputfilename == None:
    return None
  
  rawfile = gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = TFile( inputfilename, 'READ' )
  
  histofile = gROOT.FindObject( histofilename )
  if histofile:
    histofile.Close()
  histofile = TFile( histofilename, 'RECREATE', 'Inpixel histos created from input file ' + inputfilename )
  
  # Get access to hits 
  tree = rawfile.Get("Hit")
  
  # Temporary data objects   
  histofile.mkdir("tmp")
  histofile.cd("tmp")
  
  counter = numpy.zeros(shape=(ubins,vbins))  
  h2_cs = [[TH1F("hcs_%d_%d" % (iv,iu),"",400,0,400) for iv in range(vbins)] for iu in range(ubins)]
  h2_ss = [[TH1F("hss_%d_%d" % (iv,iu),"",400,0,400) for iv in range(vbins)] for iu in range(ubins)]
  h2_uu = [[TH1F("huu_%d_%d" % (iv,iu),"",10,0,10) for iv in range(vbins)] for iu in range(ubins)]
  h2_vv = [[TH1F("hvv_%d_%d" % (iv,iu),"",10,0,10) for iv in range(vbins)] for iu in range(ubins)]
  h2_px = [[TH1F("hpx_%d_%d" % (iv,iu),"",10,0,10) for iv in range(vbins)] for iu in range(ubins)]
    
  histofile.cd("")
  
  # Final histos 
  h2c = TH2F("h2c","h2c",ubins,-usize,usize,vbins,-vsize,+vsize)
  h2cs = TH2F("h2cs","h2cs",ubins,-usize,usize,vbins,-vsize,+vsize)
  h2ss = TH2F("h2ss","h2ss",ubins,-usize,usize,vbins,-vsize,+vsize)
  h2uu = TH2F("h2uu","h2uu",ubins,-usize,usize,vbins,-vsize,vsize)
  h2vv = TH2F("h2vv","h2vv",ubins,-usize,usize,vbins,-vsize,vsize)
  h2px = TH2F("h2px","h2px",ubins,-usize,usize,vbins,-vsize,vsize)
  
  for event in tree: 
    if event.hasTrack == 0: 
      m_u = event.u_fit  
      m_v = event.v_fit
      
      if m_u<-usize:
        continue
      if m_u>+usize:
        continue
      if m_v<-vsize:
        continue
      if m_v>+vsize:
        continue
        
      # Note: This must be c-style index 
      iu = h2cs.GetXaxis().FindBin(m_u)-1
      iv = h2cs.GetYaxis().FindBin(m_v)-1
      
      # Fill maps  
      counter[iu,iv] += 1
      h2_cs[iu][iv].Fill(event.clusterCharge)
      h2_ss[iu][iv].Fill(event.seedCharge)
      h2_uu[iu][iv].Fill(event.sizeU)
      h2_vv[iu][iv].Fill(event.sizeV)
      h2_px[iu][iv].Fill(event.size)
   
  # fill 2d histos 
  
  for iu in range(ubins):
    for iv in range(vbins):  
      mean_cs = h2_cs[iu][iv].GetMean()
      mean_ss= h2_ss[iu][iv].GetMean()
      mean_uu= h2_uu[iu][iv].GetMean()
      mean_vv= h2_vv[iu][iv].GetMean()
      mean_px= h2_px[iu][iv].GetMean()
        
      h2cs.SetBinContent(iu+1,iv+1,mean_cs)
      h2ss.SetBinContent(iu+1,iv+1,mean_ss)
      h2c.SetBinContent(iu+1,iv+1,(counter[iu,iv]))
      h2uu.SetBinContent(iu+1,iv+1,(mean_uu))
      h2vv.SetBinContent(iu+1,iv+1,(mean_vv))
      h2px.SetBinContent(iu+1,iv+1,(mean_px))


  # plot histograms    
       
  h2c.SetTitle("")
  h2c.GetXaxis().SetTitle("u_fit [mm]")
  h2c.GetYaxis().SetTitle("v_fit [mm]")
  h2c.GetZaxis().SetTitle("Number of tracks")
  h2c.SetStats(0)

  h2ss.SetTitle("")
  h2ss.GetXaxis().SetTitle("u_fit [mm]")
  h2ss.GetYaxis().SetTitle("v_fit [mm]")
  h2ss.GetZaxis().SetTitle("Mean Seed Signal [LSB]")
  h2ss.SetStats(0)

  h2cs.SetTitle("")
  h2cs.GetXaxis().SetTitle("u_fit [mm]")
  h2cs.GetYaxis().SetTitle("v_fit [mm]")
  h2cs.GetZaxis().SetTitle("Mean Cluster Signal [LSB]")
  h2cs.SetStats(0)

  h2uu.SetTitle("")
  h2uu.GetXaxis().SetTitle("u_fit [mm]")
  h2uu.GetYaxis().SetTitle("v_fit [mm]")
  h2uu.GetZaxis().SetTitle("Mean sizeU [cells]")
  h2uu.SetStats(0)

  h2vv.SetTitle("")
  h2vv.GetXaxis().SetTitle("u_fit [mm]")
  h2vv.GetYaxis().SetTitle("v_fit [mm]")
  h2vv.GetZaxis().SetTitle("Mean sizeV [cells]")
  h2vv.SetStats(0)

  h2px.SetTitle("")
  h2px.GetXaxis().SetTitle("u_fit [mm]")
  h2px.GetYaxis().SetTitle("v_fit [mm]")
  h2px.GetZaxis().SetTitle("Mean size [celss]")
  h2px.SetStats(0)

  # write the tree into the output file and close the file
  histofile.Write()
  histofile.Close()
  rawfile.Close()


def plot_unit(inputfilename = None, histofilename = "InpixHistos.root", upitch=0.0, vpitch=0.0, ubins=10, vbins=10):
  
  if inputfilename == None:
    return None
  
  rawfile = gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = TFile( inputfilename, 'READ' )
  
  histofile = gROOT.FindObject( histofilename )
  if histofile:
    histofile.Close()
  histofile = TFile( histofilename, 'RECREATE', 'Inpixel histos created from input file ' + inputfilename )
  
  # Get access to hits 
  tree = rawfile.Get("Hit")
 
  # Super pixel contains upix x vpix cells  
  upix = 2  
  vpix = 1 

  counter = numpy.zeros(shape=(ubins,vbins))
  
  histofile.mkdir("tmp")
  histofile.cd("tmp")

  h2_cs = [[TH1F("hcs_%d_%d" % (iv,iu),"",400,0,400) for iv in range(vbins)] for iu in range(ubins)]
  h2_ss = [[TH1F("hss_%d_%d" % (iv,iu),"",400,0,400) for iv in range(vbins)] for iu in range(ubins)]
  h2_uu = [[TH1F("huu_%d_%d" % (iv,iu),"",10,0,10) for iv in range(vbins)] for iu in range(ubins)]
  h2_vv = [[TH1F("hvv_%d_%d" % (iv,iu),"",10,0,10) for iv in range(vbins)] for iu in range(ubins)]
  h2_px = [[TH1F("hpx_%d_%d" % (iv,iu),"",10,0,10) for iv in range(vbins)] for iu in range(ubins)]
  
  
  histofile.cd("")

  for event in tree: 
    if event.hasTrack == 0: 
      m_u = (event.u_fit - event.cellUCenter_fit) 
      m_v = (event.v_fit - event.cellVCenter_fit)
      
      m_u += upitch/2 
      m_v += vpitch/2
      
      if m_u<0:
        continue
      if m_u>=upitch:
        continue
      if m_v<0:
        continue
      if m_v>=vpitch:
        continue
      
      # Pitch of inpixel bins 
      subpitchU = (upix*upitch)/ubins
      subpitchV = (vpix*vpitch)/vbins
      
      # Calculate inpixel bins 
      if event.cellU_fit < 0:
        continue 
      if event.cellV_fit < 0:
        continue 
       
      m_u+= (event.cellU_fit%upix) * upitch   
      m_v+= (event.cellV_fit%vpix) * vpitch
        
      iu = int( m_u / subpitchU )
      iv = int( m_v / subpitchV )
      
      # Fill maps 
      counter[iu,iv] += 1
      h2_cs[iu][iv].Fill(event.clusterCharge)
      h2_ss[iu][iv].Fill(event.seedCharge)
      h2_uu[iu][iv].Fill(event.sizeU)
      h2_vv[iu][iv].Fill(event.sizeV)
      h2_px[iu][iv].Fill(event.size)
     


  h2c =  TH2F("h2c","h2c",ubins,0,upix*upitch,vbins,0,vpix*vpitch)
  h2cs = TH2F("h2cs","h2cs",ubins,0,upix*upitch,vbins,0,vpix*vpitch)
  h2ss = TH2F("h2ss","h2ss",ubins,0,upix*upitch,vbins,0,vpix*vpitch)
  h2uu = TH2F("h2uu","h2uu",ubins,0,upix*upitch,vbins,0,vpix*vpitch)
  h2vv = TH2F("h2vv","h2vv",ubins,0,upix*upitch,vbins,0,vpix*vpitch)
  h2px = TH2F("h2px","h2px",ubins,0,upix*upitch,vbins,0,vpix*vpitch)
  
  # fill 2d histos 
  
  for iu in range(ubins):
    for iv in range(vbins):  
      mean_cs = h2_cs[iu][iv].GetMean()
      mean_ss= h2_ss[iu][iv].GetMean()
      mean_uu= h2_uu[iu][iv].GetMean()
      mean_vv= h2_vv[iu][iv].GetMean()      
      mean_px= h2_px[iu][iv].GetMean()
  
      h2cs.SetBinContent(iu+1,iv+1,mean_cs)
      h2ss.SetBinContent(iu+1,iv+1,mean_ss)
      h2c.SetBinContent(iu+1,iv+1,(counter[iu,iv]))
      h2uu.SetBinContent(iu+1,iv+1,mean_uu)
      h2vv.SetBinContent(iu+1,iv+1,mean_vv)
      h2px.SetBinContent(iu+1,iv+1,mean_px)

  # plot histograms    
       
  h2c.SetTitle("")
  h2c.GetXaxis().SetTitle("u_{m} [#mum]")
  h2c.GetYaxis().SetTitle("v_{m} [#mum]")
  h2c.GetZaxis().SetTitle("Number of tracks")
  h2c.SetStats(0)   

  h2ss.SetTitle("")
  h2ss.GetXaxis().SetTitle("u_{m} [#mum]")
  h2ss.GetYaxis().SetTitle("v_{m} [#mum]")
  h2ss.GetZaxis().SetTitle("Mean Seed Signal [LSB]")
  h2ss.SetStats(0)  

  h2cs.SetTitle("")
  h2cs.GetXaxis().SetTitle("u_{m} [#mum]")
  h2cs.GetYaxis().SetTitle("v_{m} [#mum]")
  h2cs.GetZaxis().SetTitle("Mean Cluster Signal [LSB]")
  h2cs.SetStats(0)

  h2ss_x = h2ss.ProjectionX("h2ss_x")
  h2ss_x.SetTitle("")
  h2ss_x.GetXaxis().SetTitle("u_{m} [#mum]")
  h2ss_x.GetYaxis().SetTitle("Mean Seed Signal [LSB]")

  h2ss_y = h2ss.ProjectionY("h2ss_y")
  h2ss_y.SetTitle("")
  h2ss_y.GetXaxis().SetTitle("v_{m} [#mum]")
  h2ss_y.GetYaxis().SetTitle("Mean Seed Signal [LSB]")
  

  h2uu.SetTitle("")
  h2uu.GetXaxis().SetTitle("u_{m} [#mum]")
  h2uu.GetYaxis().SetTitle("v_{m} [#mum]")
  h2uu.GetZaxis().SetTitle("Mean sizeU [uCells]")
  h2uu.SetStats(0)

  h2vv.SetTitle("")
  h2vv.GetXaxis().SetTitle("u_{m} [#mum]")
  h2vv.GetYaxis().SetTitle("v_{m} [#mum]")
  h2vv.GetZaxis().SetTitle("Mean sizeV [vCells]")
  h2vv.SetStats(0)

  h2px.SetTitle("")
  h2px.GetXaxis().SetTitle("u_{m} [#mum]")
  h2px.GetYaxis().SetTitle("v_{m} [#mum]")
  h2px.GetZaxis().SetTitle("Mean size [Pixel]")
  h2px.SetStats(0)

  # write the tree into the output file and close the file
  histofile.Write()
  histofile.Close()
  rawfile.Close()




