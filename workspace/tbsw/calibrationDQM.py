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
  c_xshifts.SaveAs("hxshift.pdf")


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
  c_yshifts.SaveAs("hyshift.pdf")


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
  c_zshifts.SaveAs("hzshift.pdf")


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
  c_zrots.SaveAs("hzrots.pdf")



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
  colors=[0,1,2,4,6,7]

  c_resu = TCanvas( 'c_resu', 'c_resu', 200, 10, 700, 500 )

  for isens in range(0,nsensors):

    sensorname='Sensor'+str(sensorIDs[isens])

    histoname=sensorname+"/hresU_sensor"+str(sensorIDs[isens])
    hresu = rawfile.Get(histoname)

    integral =  hresu.Integral()
    hresu.Scale(1.0/integral)

    hresu.SetMaximum(0.045)
    hresu.SetMinimum(0.0)
    hresu.SetLineColor(colors[isens])

    if isens > 0:
      hresu.Draw("Hsame")
    else:
      hresu.Draw()
    hresu.SetStats(False)

    c_resu.Update()
  c_resu.SaveAs("hresu.pdf")


  c_resv = TCanvas( 'c_resv', 'c_resv', 200, 10, 700, 500 )
  leg = TLegend(.73,.32,.97,.53)

  for isens in range(0,nsensors):

    sensorname='Sensor'+str(sensorIDs[isens])

    histoname=sensorname+"/hresV_sensor"+str(sensorIDs[isens])
    hresv = rawfile.Get(histoname)

    integral =  hresv.Integral()
    hresv.Scale(1.0/integral)

    hresv.SetMaximum(0.045)
    hresv.SetMinimum(0.0)
    hresv.SetLineColor(colors[isens])

    if isens > 0:
      hresv.Draw("Hsame")
    else:
      hresv.Draw()

    hresv.SetStats(False)


    leg.AddEntry(hresv,sensorname,"H")
    c_resv.Update()

  leg.Draw()
  c_resv.SaveAs("hresv.pdf")


  c_pullu = TCanvas( 'c_pullu', 'c_pullu', 200, 10, 700, 500 )

  for isens in range(0,nsensors):

    sensorname='Sensor'+str(sensorIDs[isens])

    histoname=sensorname+"/hpull_resU_sensor"+str(sensorIDs[isens])
    hpullu = rawfile.Get(histoname)

    integral =  hpullu.Integral()
    hpullu.Scale(1.0/integral)

    #hpullu.SetMaximum(0.045)
    #hpullu.SetMinimum(0.0)
    hpullu.SetLineColor(colors[isens])

    if isens > 0:
      hpullu.Draw("Hsame")
    else:
      hpullu.Draw()
    hpullu.SetStats(False)

    c_pullu.Update()
  c_pullu.SaveAs("hpullu.pdf")


  c_pullv = TCanvas( 'c_pullv', 'c_pullv', 200, 10, 700, 500 )

  for isens in range(0,nsensors):

    sensorname='Sensor'+str(sensorIDs[isens])

    histoname=sensorname+"/hpull_resV_sensor"+str(sensorIDs[isens])
    hpullv = rawfile.Get(histoname)

    integral =  hpullv.Integral()
    hpullv.Scale(1.0/integral)

    #hpullv.SetMaximum(0.045)
    #hpullv.SetMinimum(0.0)
    hpullv.SetLineColor(colors[isens])

    if isens > 0:
      hpullv.Draw("Hsame")
    else:
      hpullv.Draw()

    hpullv.SetStats(False)
    c_pullv.Update()
  c_pullv.SaveAs("hpullv.pdf")



def plot_pvalues(inputfilename = None):

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

  c_pvalues = TCanvas( 'c_pvalues', 'c_pvalues', 200, 10, 700, 500 )
  hpvalues.Draw()
  hpvalues.SetStats(False)
  c_pvalues.Update()
  c_pvalues.SaveAs("hpvalues.pdf")



