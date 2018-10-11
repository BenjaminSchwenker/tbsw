from ROOT import TFile, TCanvas, TF1, TH1F, TH2F, TGraphAsymmErrors, TLegend, TAxis
from ROOT import gROOT, Double, TCut
import math

def plot_alignment_parameters(inputfilename = None):

  if inputfilename == None:
    return None

  gROOT.Reset()

  # Cluster validation file name
  rawfile = gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = TFile( inputfilename, 'READ' )

  histoname="alignment/hxshift_diff"
  hxshift_tmp = rawfile.Get(histoname)

  hxshift = TH1F("hxshift", "hxshift", hxshift_tmp.GetNbinsX(), 0.5, hxshift_tmp.GetNbinsX()+0.5)
  for i in range(0,hxshift_tmp.GetNbinsX()):
    hxshift.SetBinContent(i,hxshift_tmp.GetBinContent(i))
  hxshift.SetTitle("")
  hxshift.SetYTitle("x_{nom.}-x_{align.} [mm]")
  hxshift.GetYaxis().SetTitleOffset(1.25)
  hxshift.SetXTitle("sensor position")
  hxshift.SetNdivisions(hxshift_tmp.GetNbinsX())

  c_xshifts = TCanvas( 'c_xshifts', 'c_xshifts', 200, 10, 700, 500 )
  hxshift.Draw()
  hxshift.SetStats(False)
  c_xshifts.Update()
  c_xshifts.SaveAs("alignment_shift_x.pdf")


  histoname="alignment/hyshift_diff"
  hyshift_tmp = rawfile.Get(histoname)

  hyshift = TH1F("hyshift", "hyshift", hyshift_tmp.GetNbinsX(), 0.5, hyshift_tmp.GetNbinsX()+0.5)
  for i in range(0,hyshift_tmp.GetNbinsX()):
    hyshift.SetBinContent(i,hyshift_tmp.GetBinContent(i))
  hxshift.SetTitle("")
  hyshift.SetYTitle("y_{nom.}-y_{align.} [mm]")
  hyshift.GetYaxis().SetTitleOffset(1.25)
  hyshift.SetXTitle("sensor position")
  hyshift.SetNdivisions(hyshift_tmp.GetNbinsX())

  c_yshifts = TCanvas( 'c_yshifts', 'c_yshifts', 200, 10, 700, 500 )
  hyshift.Draw()
  hyshift.SetStats(False)
  c_yshifts.Update()
  c_yshifts.SaveAs("alignment_shift_y.pdf")


  histoname="alignment/hzshift_diff"
  hzshift_tmp = rawfile.Get(histoname)

  hzshift = TH1F("hzshift", "hzshift", hzshift_tmp.GetNbinsX(), 0.5, hzshift_tmp.GetNbinsX()+0.5)
  for i in range(0,hzshift_tmp.GetNbinsX()):
    hzshift.SetBinContent(i,hzshift_tmp.GetBinContent(i))
  hzshift.SetTitle("")
  hzshift.SetYTitle("z_{nom.}-z_{align.} [mm]")
  hzshift.GetYaxis().SetTitleOffset(1.25)
  hzshift.SetXTitle("sensor position")
  hzshift.SetNdivisions(hzshift_tmp.GetNbinsX())

  c_zshifts = TCanvas( 'c_zshifts', 'c_zshifts', 200, 10, 700, 500 )
  hzshift.Draw()
  hzshift.SetStats(False)
  c_zshifts.Update()
  c_zshifts.SaveAs("alignment_shift_z.pdf")


  histoname="alignment/hzrot_diff"
  hzrot_tmp = rawfile.Get(histoname)

  hzrot = TH1F("hzrot", "hzrot", hzrot_tmp.GetNbinsX(), 0.5, hzrot_tmp.GetNbinsX()+0.5)
  for i in range(0,hzrot_tmp.GetNbinsX()):
    hzrot.SetBinContent(i,hzrot_tmp.GetBinContent(i))
  hzrot.SetTitle("")
  hzrot.SetYTitle("#gamma_{nom.}-#gamma_{align.} [mrad]")
  hzrot.GetYaxis().SetTitleOffset(1.0)
  hzrot.SetXTitle("sensor position")
  hzrot.Scale(1E3)
  hzrot.SetNdivisions(hzrot_tmp.GetNbinsX())


  c_zrots = TCanvas( 'c_zrots', 'c_zrots', 200, 10, 700, 500 )
  hzrot.Draw("hist")
  hzrot.SetStats(False)
  c_zrots.Update()
  c_zrots.SaveAs("alignment_rotation_z.pdf")



def plot_clusterDB_parameters(inputfilename = None):

  if inputfilename == None:
    return None

  oldlabels=["H1.0.0D0.0.0","H2.0.0D0.0.0D0.1.0","H2.0.0D0.0.0D1.0.0","H3.0.0D0.0.0D0.1.0D1.0.0","H3.0.0D0.0.0D0.1.0D1.1.0","H3.0.0D0.0.0D1.0.0D1.1.0","H3.0.0D0.1.0D1.0.0D1.1.0","H4.0.0D0.0.0D0.1.0D1.0.0D1.1.0"]
  newlabels=["1p","2pu","2pv","3p1","3p2","3p3","3p4","4p"]

  # Cluster validation file name
  rawfile = gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = TFile( inputfilename, 'READ' )

  histoname="hDB_Weight"
  hfractions_tmp = rawfile.Get(histoname)
  integral =  hfractions_tmp.Integral()
  hfractions_tmp.Scale(100.0/integral)
  hfractions = TH1F("hfractions", "cluster type fractions", len(oldlabels),0.5,len(oldlabels)+0.5)

  histoname="hDB_Sigma2_U"
  hsigma2_u_tmp = rawfile.Get(histoname)
  hsigma_u = TH1F("hsigma_u", "cluster u resolution ", len(oldlabels),0.5,len(oldlabels)+0.5)

  histoname="hDB_Sigma2_V"
  hsigma2_v_tmp = rawfile.Get(histoname)
  hsigma_v = TH1F("hsigma_v", "cluster v resolution ", len(oldlabels),0.5,len(oldlabels)+0.5)

  for ilabel,label in enumerate(oldlabels):

    hfractions.GetXaxis().SetBinLabel( ilabel+1, newlabels[ilabel] )
    hfractions.SetBinContent(ilabel+1,hfractions_tmp.GetBinContent(hfractions_tmp.GetXaxis().FindBin(label)))
    hfractions.SetBinError(ilabel+1,hfractions_tmp.GetBinError(hfractions_tmp.GetXaxis().FindBin(label)))

    hsigma_u.GetXaxis().SetBinLabel( ilabel+1, newlabels[ilabel] )
    hsigma_u.SetBinContent(ilabel+1,1000*math.sqrt(hsigma2_u_tmp.GetBinContent(hsigma2_u_tmp.GetXaxis().FindBin(label))))
    hsigma_u.SetBinError(ilabel+1,1000*hsigma2_u_tmp.GetBinError(hsigma2_u_tmp.GetXaxis().FindBin(label))/math.sqrt(hsigma2_u_tmp.GetBinContent(hsigma2_u_tmp.GetXaxis().FindBin(label))))

    hsigma_v.GetXaxis().SetBinLabel( ilabel+1, newlabels[ilabel] )
    hsigma_v.SetBinContent(ilabel+1,1000*math.sqrt(hsigma2_v_tmp.GetBinContent(hsigma2_v_tmp.GetXaxis().FindBin(label))))
    hsigma_v.SetBinError(ilabel+1,1000*hsigma2_v_tmp.GetBinError(hsigma2_v_tmp.GetXaxis().FindBin(label))/math.sqrt(hsigma2_v_tmp.GetBinContent(hsigma2_v_tmp.GetXaxis().FindBin(label))))

  hfractions.SetStats(False)
  hfractions.GetXaxis().SetTitle("cluster type")
  hfractions.GetXaxis().SetTitleOffset(0.88)
  hfractions.GetXaxis().SetTitleSize(0.055)
  hfractions.GetXaxis().SetLabelSize(0.07)
  hfractions.GetYaxis().SetTitle("fraction")
  hfractions.GetYaxis().SetTitleOffset(0.80)
  hfractions.GetYaxis().SetTitleSize(0.055)
  hfractions.GetYaxis().SetLabelSize(0.05)

  hsigma_u.SetStats(False)
  hsigma_u.GetXaxis().SetTitle("cluster type")
  hsigma_u.GetXaxis().SetTitleOffset(0.88)
  hsigma_u.GetXaxis().SetTitleSize(0.055)
  hsigma_u.GetXaxis().SetLabelSize(0.07)
  hsigma_u.GetYaxis().SetTitle("#sigma_{u}[#mum]")
  hsigma_u.GetYaxis().SetTitleOffset(0.80)
  hsigma_u.GetYaxis().SetTitleSize(0.055)
  hsigma_u.GetYaxis().SetLabelSize(0.05)

  hsigma_v.SetStats(False)
  hsigma_v.GetXaxis().SetTitle("cluster type")
  hsigma_v.GetXaxis().SetTitleOffset(0.88)
  hsigma_v.GetXaxis().SetTitleSize(0.055)
  hsigma_v.GetXaxis().SetLabelSize(0.07)
  hsigma_v.GetYaxis().SetTitle("#sigma_{v}[#mum]")
  hsigma_v.GetYaxis().SetTitleOffset(0.80)
  hsigma_v.GetYaxis().SetTitleSize(0.055)
  hsigma_v.GetYaxis().SetLabelSize(0.05)


  c1 = TCanvas( 'c1', 'c1', 200, 10, 700, 500 )
  hfractions.Draw("histE")
  c1.Update()
  c1.SaveAs("cluster_fractions.pdf")

  c2 = TCanvas( 'c2', 'c2', 200, 10, 700, 500 )
  hsigma_u.Draw("histE")
  c2.Update()
  c2.SaveAs("cluster_sigma_u.pdf")

  c3 = TCanvas( 'c3', 'c3', 200, 10, 700, 500 )
  hsigma_v.Draw("histE")
  c3.Update()
  c3.SaveAs("cluster_sigma_v.pdf")



def plot_pulls(inputfilename = None):

  if inputfilename == None:
    return None

  gROOT.Reset()

  # Cluster validation file name
  rawfile = gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = TFile( inputfilename, 'READ' )

  nsensors=6
  sensorIDs=[0,1,2,4,5,6]
  colors=[1,2,4,6,7,8]

  leg_u = TLegend(.73,.73,.97,.97)
  c_pullu = TCanvas( 'c_pullu', 'c_pullu', 200, 10, 700, 500 )

  for isens in range(0,nsensors):

    sensorname='Sensor'+str(sensorIDs[isens])
    histoname=sensorname+"/hpull_resU_sensor"+str(sensorIDs[isens])
    hpullu = rawfile.Get(histoname)

    hpullu.GetXaxis().SetTitle("track pull")
    hpullu.GetYaxis().SetTitle("number of tracks[10^{3}]")
    hpullu.SetTitle("track pull (u direction)")

#    integral =  hpullu.Integral()
    hpullu.Scale(1E-3)

    hpullu.SetLineColor(colors[isens])

    if isens > 0:
      hpullu.Draw("Hsame")
    else:
      hpullu.Draw()

    hpullu.SetStats(False)
    leg_u.AddEntry(hpullu,sensorname,"LE")
    c_pullu.Update()

  leg_u.Draw()
  c_pullu.SaveAs("trackpulls_u.pdf")


  leg_v = TLegend(.73,.73,.97,.97)
  c_pullv = TCanvas( 'c_pullv', 'c_pullv', 200, 10, 700, 500 )

  for isens in range(0,nsensors):

    sensorname='Sensor'+str(sensorIDs[isens])

    histoname=sensorname+"/hpull_resV_sensor"+str(sensorIDs[isens])
    hpullv = rawfile.Get(histoname)

    hpullv.GetXaxis().SetTitle("track pull")
    hpullv.GetYaxis().SetTitle("number of tracks[10^{3}]")
    hpullv.SetTitle("track pull (v direction)")

#   integral =  hpullv.Integral()
    hpullv.Scale(1E-3)

    hpullv.SetLineColor(colors[isens])

    if isens > 0:
      hpullv.Draw("Hsame")
    else:
      hpullv.Draw()

    hpullv.SetStats(False)
    leg_v.AddEntry(hpullv,sensorname,"LE")
    c_pullv.Update()

  leg_v.Draw()
  c_pullv.SaveAs("trackpulls_v.pdf")



def plot_trackDQM(inputfilename = None):

  if inputfilename == None:
    return None

  gROOT.Reset()

  # Cluster validation file name
  rawfile = gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = TFile( inputfilename, 'READ' )

  histoname="hchi2prob"
  hpvalues = rawfile.Get(histoname)

  histoname="hntracks"
  hntracks = rawfile.Get(histoname)

  histoname="Sensor2/BeamProfile/hhitmap_sensor2"
  hhitmap = rawfile.Get(histoname)

  hpvalues.Scale(1E-3)
  hpvalues.GetXaxis().SetTitle("track p value")
  hpvalues.GetYaxis().SetTitle("number of tracks[10^{3}]")

  hntracks.Scale(1E-3)
  hntracks.GetYaxis().SetTitle("number of events[10^{3}]")

  hhitmap.GetZaxis().SetTitleOffset(0.95)


  c_pvalues = TCanvas( 'c_pvalues', 'c_pvalues', 200, 10, 700, 500 )
  hpvalues.Draw("hist")
  hpvalues.SetStats(False)
  c_pvalues.Update()
  c_pvalues.SaveAs("track_pvalues.pdf")

  c_ntracks = TCanvas( 'c_ntracks', 'c_ntracks', 200, 10, 700, 500 )
  hntracks.Draw("hist")
  hntracks.SetStats(False)
  c_ntracks.Update()
  c_ntracks.SaveAs("ntracks.pdf")

  c_hitmap = TCanvas( 'c_hitmap', 'c_hitmap', 200, 10, 700, 500 )
  c_hitmap.SetRightMargin(0.13)
  hhitmap.Draw("colz")
  hhitmap.SetStats(False)
  c_hitmap.Update()
  c_hitmap.SaveAs("hitmap.pdf")



def plot_anglereco_DQM(inputfilename = None):

  if inputfilename == None:
    return None

  gROOT.Reset()

  rawfile = gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = TFile( inputfilename, 'READ' )

  tree = rawfile.Get("MSCTree")

  h_theta_u_reso = TH1F()
  tree.Draw("1E6*sqrt(theta1_var) >>h_theta_u_reso")
  h_theta_u_reso= tree.GetHistogram()
  h_theta_u_reso.Scale(1E-3)
  h_theta_u_reso.GetXaxis().SetTitle("#theta_{u} reconstruction error[#murad]")
  h_theta_u_reso.GetYaxis().SetTitle("reconstructed angles[10^{3}]")
  h_theta_u_reso.SetTitle("")
  h_theta_u_reso.SetStats(False)

  h_theta_v_reso = TH1F()
  tree.Draw("1E6*sqrt(theta2_var) >>h_theta_v_reso")
  h_theta_v_reso= tree.GetHistogram()
  h_theta_v_reso.Scale(1E-3)
  h_theta_v_reso.GetXaxis().SetTitle("#theta_{v} reconstruction error[#murad]")
  h_theta_v_reso.GetYaxis().SetTitle("reconstructed angles[10^{3}]")
  h_theta_v_reso.SetTitle("")
  h_theta_v_reso.SetStats(False)

  c_theta_u_reso = TCanvas( 'c_theta_u_reso', 'c_theta_u_reso', 200, 10, 700, 500 )
  c_theta_u_reso.SetRightMargin(0.13)
  h_theta_u_reso.Draw("hist")
  h_theta_u_reso.SetStats(False)
  c_theta_u_reso.Update()
  c_theta_u_reso.SaveAs("theta_u_reso.pdf")

  c_theta_v_reso = TCanvas( 'c_theta_v_reso', 'c_theta_v_reso', 200, 10, 700, 500 )
  c_theta_v_reso.SetRightMargin(0.13)
  h_theta_v_reso.Draw("hist")
  h_theta_v_reso.SetStats(False)
  c_theta_v_reso.Update()
  c_theta_v_reso.SaveAs("theta_v_reso.pdf")




