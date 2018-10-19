import os
import shutil
import glob
import math

import ROOT
import math


canvases = {}
def plot_align_parameter(histoname, paramname, unit, scale,uselatex, inputfilename = None):

  if inputfilename == None:
    return None

  ROOT.gROOT.Reset()

  rawfile = ROOT.gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = ROOT.TFile( inputfilename, 'READ' )

  diffhisto_tmp = rawfile.Get(histoname)

  diffhisto = ROOT.TH1F("h{:s}shift", "h{:s}shift", diffhisto_tmp.GetNbinsX(), 0.5, diffhisto_tmp.GetNbinsX()+0.5)
  for i in range(0,diffhisto_tmp.GetNbinsX()):
    diffhisto.SetBinContent(i,diffhisto_tmp.GetBinContent(i))
  diffhisto.SetTitle("")

  ytitle='{:s}'.format(paramname)+'_{nom.}-'+'{:s}'.format(paramname)+'_{align.} '+'[{:s}]'.format(unit)
  if uselatex:
    ytitle='#{:s}'.format(paramname)+'_{nom.}-'+'#{:s}'.format(paramname)+'_{align.} '+'[{:s}]'.format(unit)
  diffhisto.SetYTitle(ytitle)
  diffhisto.GetYaxis().SetTitleOffset(1.25)
  diffhisto.SetXTitle("sensor position")
  diffhisto.Scale(scale)
  diffhisto.SetNdivisions(diffhisto.GetNbinsX())

  canvases[histoname] = ROOT.TCanvas( '{:s}'.format(histoname), '{:s}'.format(histoname), 200, 10, 700, 500 )
  diffhisto.Draw("hist")
  diffhisto.SetStats(False)
  canvases[histoname].Update()
  canvases[histoname].SaveAs("alignment_shift_{:s}.pdf".format(paramname))


def plot_alignment_parameters(inputfilename = None):

  if inputfilename == None:
    return None

  ROOT.gROOT.Reset()

  histoname="alignment/hxshift_diff"
  plot_align_parameter(histoname,'x','mm',1,False,inputfilename)

  histoname="alignment/hyshift_diff"
  plot_align_parameter(histoname,'y','mm',1,False,inputfilename)

  histoname="alignment/hzshift_diff"
  plot_align_parameter(histoname,'z','mm',1,False,inputfilename)

  histoname="alignment/hzrot_diff"
  plot_align_parameter(histoname,'gamma','mrad',1E3,True,inputfilename)


def select_labels(histo, tmphisto, calculateSigma, oldlabels, newlabels):
  for ilabel,label in enumerate(oldlabels):
    histo.GetXaxis().SetBinLabel( ilabel+1, newlabels[ilabel] )

    if calculateSigma:
      histo.SetBinContent(ilabel+1,1000*math.sqrt(tmphisto.GetBinContent(tmphisto.GetXaxis().FindBin(label))))
      histo.SetBinError(ilabel+1,1000*tmphisto.GetBinError(tmphisto.GetXaxis().FindBin(label))/math.sqrt(tmphisto.GetBinContent(tmphisto.GetXaxis().FindBin(label))))

    else:
      histo.SetBinContent(ilabel+1,tmphisto.GetBinContent(tmphisto.GetXaxis().FindBin(label)))
      histo.SetBinError(ilabel+1,tmphisto.GetBinError(tmphisto.GetXaxis().FindBin(label)))


def plot_clusterDB_parameter(histo, paramname, ytitle):

  ROOT.gROOT.Reset()

  histo.SetStats(False)
  histo.GetXaxis().SetTitle("cluster type")
  histo.GetXaxis().SetTitleOffset(0.88)

  histo.GetXaxis().SetTitleSize(0.055)
  histo.GetXaxis().SetLabelSize(0.07)
  histo.GetYaxis().SetTitle(ytitle)
  histo.GetYaxis().SetTitleOffset(0.80)
  histo.GetYaxis().SetTitleSize(0.055)
  histo.GetYaxis().SetLabelSize(0.05)

  canvases[paramname] = ROOT.TCanvas( '{:s}'.format(paramname), '{:s}'.format(paramname), 200, 10, 700, 500 )
  histo.Draw("HE")
  histo.SetStats(False)
  canvases[paramname].Update()
  canvases[paramname].SaveAs("cluster_{:s}.pdf".format(paramname))



def plot_clusterDB_parameters(inputfilename = None):

  if inputfilename == None:
    return None

  oldlabels=["H1.0.0D0.0.0","H2.0.0D0.0.0D0.1.0","H2.0.0D0.0.0D1.0.0","H3.0.0D0.0.0D0.1.0D1.0.0","H3.0.0D0.0.0D0.1.0D1.1.0","H3.0.0D0.0.0D1.0.0D1.1.0","H3.0.0D0.1.0D1.0.0D1.1.0","H4.0.0D0.0.0D0.1.0D1.0.0D1.1.0"]
  newlabels=["1p","2pu","2pv","3p1","3p2","3p3","3p4","4p"]

  rawfile = ROOT.gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = ROOT.TFile( inputfilename, 'READ' )

  histoname="hDB_Weight"
  hfractions_tmp = rawfile.Get(histoname)
  integral =  hfractions_tmp.Integral()
  hfractions_tmp.Scale(100.0/integral)
  hfractions = ROOT.TH1F("hfractions", "cluster type fractions", len(oldlabels),0.5,len(oldlabels)+0.5)

  histoname="hDB_Sigma2_U"
  hsigma2_u_tmp = rawfile.Get(histoname)
  hsigma_u = ROOT.TH1F("hsigma_u", "cluster u resolution ", len(oldlabels),0.5,len(oldlabels)+0.5)

  histoname="hDB_Sigma2_V"
  hsigma2_v_tmp = rawfile.Get(histoname)
  hsigma_v = ROOT.TH1F("hsigma_v", "cluster v resolution ", len(oldlabels),0.5,len(oldlabels)+0.5)

  select_labels(hfractions, hfractions_tmp, False, oldlabels, newlabels)
  select_labels(hsigma_u, hsigma2_u_tmp, True, oldlabels, newlabels)
  select_labels(hsigma_v, hsigma2_v_tmp, True, oldlabels, newlabels)

  plot_clusterDB_parameter(hfractions, 'fractions', 'fraction')
  plot_clusterDB_parameter(hsigma_u, 'sigma_u', '#sigma_{u}[#mum]')
  plot_clusterDB_parameter(hsigma_v, 'sigma_v', '#sigma_{v}[#mum]')



def plot_pulls(inputfilename = None):

  if inputfilename == None:
    return None

  ROOT.gROOT.Reset()

  rawfile = ROOT.gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = ROOT.TFile( inputfilename, 'READ' )

  nsensors=6
  sensorIDs=[0,1,2,4,5,6]
  colors=[1,2,4,6,7,8]

  leg_u = ROOT.TLegend(.82,.35,.97,.97)
  leg_u.SetTextSize(.03)
  c_pullu = ROOT.TCanvas( 'c_pullu', 'c_pullu', 200, 10, 700, 500 )

  for isens in range(0,nsensors):

    sensorname='Sensor'+str(sensorIDs[isens])
    histoname=sensorname+"/hpull_resU_sensor"+str(sensorIDs[isens])
    hpullu = rawfile.Get(histoname)

    hpullu.GetXaxis().SetTitle("track pull")
    hpullu.GetYaxis().SetTitle("number of tracks[10^{3}]")
    hpullu.SetTitle("track pull (u direction)")

    hpullu.Scale(1E-3)

    hpullu.SetLineColor(colors[isens])

    if isens > 0:
      hpullu.Draw("Hsame")
    else:
      hpullu.Draw()

    hpullu.SetStats(False)
    entry='#splitline{#splitline{'+sensorname+'}{mean: '+str(round(hpullu.GetMean(),2))+'}}{RMS: '+str(round(hpullu.GetRMS(),2))+'}'
    leg_u.AddEntry(hpullu,entry,"L")
    c_pullu.Update()

  leg_u.Draw()
  c_pullu.SaveAs("trackpulls_u.pdf")


  leg_v = ROOT.TLegend(.82,.35,.97,.97)
  leg_v.SetTextSize(.03)
  c_pullv = ROOT.TCanvas( 'c_pullv', 'c_pullv', 200, 10, 700, 500 )

  for isens in range(0,nsensors):

    sensorname='Sensor'+str(sensorIDs[isens])

    histoname=sensorname+"/hpull_resV_sensor"+str(sensorIDs[isens])
    hpullv = rawfile.Get(histoname)

    hpullv.GetXaxis().SetTitle("track pull")
    hpullv.GetYaxis().SetTitle("number of tracks[10^{3}]")
    hpullv.SetTitle("track pull (v direction)")

    hpullv.Scale(1E-3)

    hpullv.SetLineColor(colors[isens])

    if isens > 0:
      hpullv.Draw("Hsame")
    else:
      hpullv.Draw()

    hpullv.SetStats(False)
    entry='#splitline{#splitline{'+sensorname+'}{mean: '+str(round(hpullv.GetMean(),2))+'}}{RMS: '+str(round(hpullv.GetRMS(),2))+'}'
    leg_v.AddEntry(hpullv,entry,"L")
    c_pullv.Update()

  leg_v.Draw()
  c_pullv.SaveAs("trackpulls_v.pdf")



def plot_trackDQM(inputfilename = None):

  if inputfilename == None:
    return None

  ROOT.gROOT.Reset()

  rawfile = ROOT.gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = ROOT.TFile( inputfilename, 'READ' )

  histoname="hchi2prob"
  hpvalues = rawfile.Get(histoname)

  histoname="hntracks"
  hntracks = rawfile.Get(histoname)

  hpvalues.Scale(1E-3)
  hpvalues.GetXaxis().SetTitle("track p value")
  hpvalues.GetYaxis().SetTitle("number of tracks[10^{3}]")

  hntracks.Scale(1E-3)
  hntracks.GetYaxis().SetTitle("number of events[10^{3}]")


  c_pvalues = ROOT.TCanvas( 'c_pvalues', 'c_pvalues', 200, 10, 700, 500 )
  hpvalues.Draw("hist")
  hpvalues.SetStats(False)
  c_pvalues.Update()
  c_pvalues.SaveAs("track_pvalues.pdf")

  c_ntracks = ROOT.TCanvas( 'c_ntracks', 'c_ntracks', 200, 10, 700, 500 )
  hntracks.Draw("hist")
  hntracks.SetStats(False)
  c_ntracks.Update()
  c_ntracks.SaveAs("ntracks.pdf")


def plot_anglereco_DQM(inputfilename = None):

  if inputfilename == None:
    return None

  ROOT.gROOT.Reset()

  rawfile = ROOT.gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = ROOT.TFile( inputfilename, 'READ' )

  tree = rawfile.Get("MSCTree")

  h_theta_u_reso = ROOT.TH1F()
  tree.Draw("1E6*sqrt(theta1_var) >>h_theta_u_reso")
  h_theta_u_reso= tree.GetHistogram()
  h_theta_u_reso.Scale(1E-3)
  h_theta_u_reso.GetXaxis().SetTitle("#theta_{u} reconstruction error[#murad]")
  h_theta_u_reso.GetYaxis().SetTitle("reconstructed angles[10^{3}]")
  h_theta_u_reso.SetTitle("")
  h_theta_u_reso.SetStats(False)

  h_theta_v_reso = ROOT.TH1F()
  tree.Draw("1E6*sqrt(theta2_var) >>h_theta_v_reso")
  h_theta_v_reso= tree.GetHistogram()
  h_theta_v_reso.Scale(1E-3)
  h_theta_v_reso.GetXaxis().SetTitle("Angle reconstruction error[#murad]")
  h_theta_v_reso.GetYaxis().SetTitle("reconstructed angles[10^{3}]")
  h_theta_v_reso.SetTitle("")
  h_theta_v_reso.SetLineColor(2)
  h_theta_v_reso.SetStats(False)

  h_trackmap = ROOT.TH2F()
  tree.Draw("v:u >>h_trackmap")
  h_trackmap= tree.GetHistogram()
  h_trackmap.GetXaxis().SetTitle("u[mm]") 
  h_trackmap.GetYaxis().SetTitle("v[mm]")
  h_trackmap.GetZaxis().SetTitle("number of tracks")
  h_trackmap.GetZaxis().SetTitleOffset(1.1)
  NoTracks=str(round(1E-6*tree.GetEntries(),1))+' million tracks'
  h_trackmap.SetTitle(NoTracks)
  h_trackmap.SetStats(False)

  c_theta_reso = ROOT.TCanvas( 'c_theta_u_reso', 'c_theta_u_reso', 200, 10, 700, 500 )
  leg = ROOT.TLegend(.13,.57,.43,.93)
  leg.SetTextSize(.038)

  c_theta_reso.SetRightMargin(0.13)
  h_theta_u_reso.Draw("hist")
  h_theta_v_reso.Draw("histsame")
  h_theta_u_reso.SetStats(False)

  description='#splitline{#theta_{u} resolution}{#splitline{mean: '+str(round(h_theta_u_reso.GetMean(),1))+' #murad}{RMS: '+str(round(h_theta_u_reso.GetRMS(),1))+' #murad}}'
  leg.AddEntry(h_theta_u_reso,description,"L")

  description='#splitline{#theta_{v} resolution}{#splitline{mean: '+str(round(h_theta_v_reso.GetMean(),1))+' #murad}{RMS: '+str(round(h_theta_v_reso.GetRMS(),1))+' #murad}}'
  leg.AddEntry(h_theta_v_reso,description,"L")
  c_theta_reso.Update()

  leg.Draw()
  c_theta_reso.SaveAs("theta_reso.pdf")

  c_trackmap = ROOT.TCanvas( 'c_trackmap', 'c_trackmap', 200, 10, 700, 500 )
  c_trackmap.SetRightMargin(0.13)
  h_trackmap.Draw("colz")
  c_trackmap.SaveAs("beamspot.pdf")
	

def plot_x0image_DQM(inputfilename = None):

  if inputfilename == None:
    return None

  ROOT.gROOT.Reset()

  rawfile = ROOT.gROOT.FindObject( inputfilename )
  if rawfile:
    rawfile.Close()
  rawfile = ROOT.TFile( inputfilename, 'READ' )

  histoname="x0_image"
  hx0 = rawfile.Get(histoname)

  c_x0 = ROOT.TCanvas( 'c_x0', 'c_x0', 200, 10, 700, 500 )
  maxbin=hx0.GetMaximumBin()
  maxvalue=hx0.GetBinContent(maxbin)
  hx0.SetMaximum(1.2*maxvalue)
  hx0.Draw("colz")
  hx0.SetStats(False)

  c_x0.Update()
  c_x0.SaveAs("x0image.pdf")


def calibration_DQMPlots(caltag,UseclusterDB):

  caldir='tmp-runs/'+caltag+'-cal'

  DQMfilename=caldir+'/TelescopeDQM2.root'
  clusterDBfilename=caldir+'/localDB/clusterDB-M26.root'

  # remember current working dir 
  fullpath = os.getcwd() 
    
  # create results dir if not exist
  if not os.path.isdir(fullpath+'/results'):
      os.mkdir(fullpath+'/results')

  # create telescopeDQM dir if not exist
  if not os.path.isdir(fullpath+'/results/telescopeDQM'):
      os.mkdir(fullpath+'/results/telescopeDQM')

  # Name of directory
  workdir = 'results/telescopeDQM/'+caltag+'-cal'

  # remove olddir if exists 
  if os.path.isdir(workdir):
    shutil.rmtree(workdir)
    
  # create dir and change directory
  os.mkdir(workdir) 
  os.chdir(workdir)

  if os.path.isfile(fullpath+'/'+DQMfilename):

    # Create file links in the current work dir
    os.symlink(fullpath+'/'+DQMfilename,'TelescopeDQM')
    os.symlink(fullpath+'/'+clusterDBfilename,'clusterDB')

  # Generate DQM plots
  plot_alignment_parameters(inputfilename = 'TelescopeDQM') 
  if UseclusterDB:
    plot_clusterDB_parameters(inputfilename = 'clusterDB') 
  plot_pulls(inputfilename = 'TelescopeDQM') 
  plot_trackDQM(inputfilename = 'TelescopeDQM') 

  os.chdir(fullpath)

def anglereco_DQMPlots(runspec,caltag):

  recofile='root-files/X0-'+runspec+caltag+'-reco.root'
  paramsDQM = (caltag,recofile)

  # remember current working dir 
  fullpath = os.getcwd() 
    
  # create results dir if not exist
  if not os.path.isdir(fullpath+'/results'):
      os.mkdir(fullpath+'/results')

  # create telescopeDQM dir if not exist
  if not os.path.isdir(fullpath+'/results/anglerecoDQM'):
      os.mkdir(fullpath+'/results/anglerecoDQM')

  # Name of directory
  workdir = 'results/anglerecoDQM/'+runspec+caltag

  # remove olddir if exists 
  if os.path.isdir(workdir):
    shutil.rmtree(workdir)
    
  # create dir and change directory
  os.mkdir(workdir) 
  os.chdir(workdir)

  if os.path.isfile(fullpath+'/'+recofile):

    # Create file link in the current work dir
    os.symlink(fullpath+'/'+recofile,'anglereco')

  plot_anglereco_DQM(inputfilename = 'anglereco') 

  os.chdir(fullpath)


def x0image_Plots(namespec,caltag,calibrated,x0tag):

  calibrationflag='Uncalibrated'
  if calibrated:
    calibrationflag='Calibrated'

  imagefile='tmp-runs/X0-'+namespec+caltag+'-reco-'+calibrationflag+'-X0image/X0-completeimage.root'

  # remember current working dir 
  fullpath = os.getcwd() 
    
  # create results dir if not exist
  if not os.path.isdir(fullpath+'/results'):
      os.mkdir(fullpath+'/results')

  # create telescopeDQM dir if not exist
  if not os.path.isdir(fullpath+'/results/images'):
      os.mkdir(fullpath+'/results/images')

  # Name of directory
  workdir = 'results/images/'+namespec+caltag+'-'+calibrationflag
  if calibrated:
    workdir = 'results/images/'+namespec+caltag+'-'+calibrationflag+'-'+x0tag

  # remove olddir if exists 
  if os.path.isdir(workdir):
    shutil.rmtree(workdir)
    
  # create dir and change directory
  os.mkdir(workdir) 
  os.chdir(workdir)

  if os.path.isfile(fullpath+'/'+imagefile):

    # Create file link in the current work dir
    os.symlink(fullpath+'/'+imagefile,'x0image')

  plot_x0image_DQM(inputfilename = 'x0image') 

  os.chdir(fullpath)


def x0calibration_DQMPlots(namespec,caltag,x0tag):

  # remember current working dir 
  fullpath = os.getcwd() 
    
  # create results dir if not exist
  if not os.path.isdir(fullpath+'/results'):
      os.mkdir(fullpath+'/results')

  # create telescopeDQM dir if not exist
  if not os.path.isdir(fullpath+'/results/x0calibrationDQM'):
      os.mkdir(fullpath+'/results/x0calibrationDQM')

  # Name of directory
  workdir = 'results/x0calibrationDQM/'+namespec+caltag+'-reco-'+x0tag

  # remove olddir if exists 
  if os.path.isdir(workdir):
    shutil.rmtree(workdir)
    
  # create dir and change directory
  os.mkdir(workdir)

  dirname='tmp-runs/X0-'+namespec+caltag+'-reco-'+x0tag+'-X0Calibration/'

  for pdffile in glob.glob(dirname+'*.pdf'): 
    shutil.copy(pdffile, os.path.join(workdir+'/',os.path.basename(pdffile))) 


