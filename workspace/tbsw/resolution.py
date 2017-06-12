from ROOT import TFile, TH1F, TH2F, TString, TMath, TVectorD
from ROOT import gROOT, Double
import re

# List of types which should be contained in a complete cluster DB
standard_types = ["1D0.0","2D0.0D0.1","2D0.0D1.0","3D0.0D0.1D1.0","3D0.0D0.1D1.1",
                  "3D0.0D1.0D1.1","3D0.1D1.0D1.1","4D0.0D0.1D1.0D1.1"
                 ]


def get_labeltype(label=None): 
  return_type = ""
  for tok in re.split('D',label):
    addr = re.split('\.',tok)
    if len(addr) == 1:
     return_type += tok    
    else: 
      return_type += "D" + addr[0] + '.' + addr[1]   
    
  return return_type
      
    
 
def compare_clusterdb(dblist = [], histofilename = "MyHistos.root" , compare_type="2D0.0D0.1"):
  
  if dblist == []:
    return None  
  
  dict_sigu     = {}
  dict_sigv     = {}
  dict_nlabels  = {}
  dict_ntypes   = {}
  dict_coverage = {}
  dict_trafo    = {}
  
  for dbfilename in dblist:
    
    # open cluster db 
    dbfile = gROOT.FindObject( dbfilename )
    if dbfile:
      dbfile.Close()
    dbfile = TFile( dbfilename, 'READ' )
  
    # read data for comparison plots 
    swADCSteps = dbfile.Get("DB_swADCSteps")  
    adclabel = '' 
    if swADCSteps:  
      for index in range(0, swADCSteps.GetNrows()):
        adclabel =  adclabel + 'S' + str(int(swADCSteps[index])) 
    dict_trafo[dbfilename] = adclabel    

    print ("db sw adc ", dict_trafo[dbfilename])
      
    coverageobj = dbfile.Get("hDB_Coverage")
    if coverageobj: 
      dict_coverage[dbfilename] =  coverageobj.GetBinContent(1)
    else: 
      dict_coverage[dbfilename] = -1    
    
    print ("db coverage ", dict_coverage[dbfilename])

    weights = dbfile.Get("hDB_Weight")
    dict_nlabels[dbfilename] = weights.GetNbinsX()
    
    print ("db labels ", dict_nlabels[dbfilename])  
 
    typeset = []    
    for bin in range(1,weights.GetNbinsX()+1):
      current_type = get_labeltype(label=weights.GetXaxis().GetBinLabel(bin))    
      if not current_type in typeset:   
        typeset.append(current_type) 
    
    dict_ntypes[dbfilename] = len(typeset)
    
    print ("db types ", dict_ntypes[dbfilename])
    
    
    histo_sigma2U = dbfile.Get("hDB_Sigma2_U")
    histo_sigma2V = dbfile.Get("hDB_Sigma2_V")

    weightedSigma2U = 0.0
    weightedSigma2V = 0.0
    labelNorm = 0.0
     
    for bin in range(1,weights.GetNbinsX()+1):  
      current_type = get_labeltype(label=weights.GetXaxis().GetBinLabel(bin))    
      if current_type == compare_type:  
        w = weights.GetBinContent(bin)
        labelNorm += w
        weightedSigma2U += w*histo_sigma2U.GetBinContent(bin)
        weightedSigma2V += w*histo_sigma2V.GetBinContent(bin)
    
    if labelNorm >0:
      weightedSigma2U/=labelNorm
      weightedSigma2V/=labelNorm
      dict_sigu[dbfilename] = TMath.Sqrt(weightedSigma2U)
      dict_sigv[dbfilename] = TMath.Sqrt(weightedSigma2V)
    else:
      dict_sigu[dbfilename] = 0
      dict_sigv[dbfilename] = 0
        
    print ("db sigma2 u ", dict_sigu[dbfilename])  
    print ("db sigma2 v ", dict_sigv[dbfilename])  

    dbfile.Close()
  
  histofile = gROOT.FindObject( histofilename )
  if histofile:
    histofile.Close()
  histofile = TFile( histofilename, 'RECREATE', 'Resolution plots created from ' + dbfilename )

  # summary histograms on type resolution 
  histofile.cd("")
 
  NDB = len(dblist)

  hcoverage = TH1F("hcoverage","",NDB,0,NDB)
  hcoverage.SetStats(0)
  hcoverage.SetFillColor(38)
  hcoverage.SetYTitle("cluster coverage [%]")  
  hcoverage.SetXTitle("cluster db")
  
  hntypes = TH1F("hntypes","",NDB,0,NDB)
  hntypes.SetStats(0)
  hntypes.SetFillColor(38)
  hntypes.SetYTitle("number of cluster types")  
  hntypes.SetXTitle("cluster db")

  hnlabels = TH1F("hnlabels","",NDB,0,NDB)
  hnlabels.SetStats(0)
  hnlabels.SetFillColor(38)
  hnlabels.SetYTitle("number of cluster labels")  
  hnlabels.SetXTitle("cluster db")  

  hsigmau = TH1F("hsigmau","",NDB,0,NDB)
  hsigmau.SetStats(0)
  hsigmau.SetFillColor(38)
  hsigmau.SetYTitle("cluster sigma u [mm]")  
  hsigmau.SetXTitle("cluster db")

  hsigmav = TH1F("hsigmav","",NDB,0,NDB)
  hsigmav.SetStats(0)
  hsigmav.SetFillColor(38)
  hsigmav.SetYTitle("cluster sigma v [mm]") 
  hsigmav.SetXTitle("cluster db")  
 
  for j, dbfilename in enumerate(dblist):
       
    hcoverage.GetXaxis().SetBinLabel(j+1, str(dict_trafo[dbfilename]))
    hntypes.GetXaxis().SetBinLabel(j+1, str(dict_trafo[dbfilename])) 
    hnlabels.GetXaxis().SetBinLabel(j+1, str(dict_trafo[dbfilename])) 
    hsigmau.GetXaxis().SetBinLabel(j+1, str(dict_trafo[dbfilename])) 
    hsigmav.GetXaxis().SetBinLabel(j+1, str(dict_trafo[dbfilename])) 
    
    hcoverage.SetBinContent(j+1, dict_coverage[dbfilename])  
    hntypes.SetBinContent(j+1, dict_ntypes[dbfilename])  
    hnlabels.SetBinContent(j+1, dict_nlabels[dbfilename])  
    hsigmau.SetBinContent(j+1, dict_sigu[dbfilename])  
    hsigmav.SetBinContent(j+1, dict_sigv[dbfilename])   

  histofile.Write()
  histofile.Close()   
       
  
    
    
def plot_clusterdb(dbfilename = None, histofilename = "MyHistos.root" ):

  if dbfilename == None:
    return None
    
  dbfile = gROOT.FindObject( dbfilename )
  if dbfile:
    dbfile.Close()
  dbfile = TFile( dbfilename, 'READ' )
  
  h_Weight  = gROOT.FindObject("hDB_Weight")
  h_Weight.SetDirectory(0)
  
  h_U = gROOT.FindObject("hDB_U")
  h_U.SetDirectory(0)
  
  h_V = gROOT.FindObject("hDB_V")
  h_V.SetDirectory(0)
  
  h_Var_U = gROOT.FindObject("hDB_Sigma2_U")
  h_Var_U.SetDirectory(0)
 
  h_Var_V = gROOT.FindObject("hDB_Sigma2_V")
  h_Var_V.SetDirectory(0)
  
  dbfile.Close()
  
  # First, we want to compute a list of all different cluster 
  # types found in the db
  
  typeset = []
  
  for bin in range(1,h_Weight.GetNbinsX()+1):
      
    # The bin label decodes the cluster shape. The  
    # label contains tokens seperated by "D". The 
    # first token is the cluster size, all other 
    # are digits. Each digits contains three tokens
    # seperated by "."; nameley iu, iv and signal.
    # These can be converted to integers. 
        
    label = h_Weight.GetXaxis().GetBinLabel(bin)
    
    current_type = ""
    for tok in re.split('D',label):
      addr = re.split('\.',tok)
      if len(addr) == 1:
        current_type += tok    
      else: 
        current_type += "D" + addr[0] + '.' + addr[1]   
          
    if not current_type in typeset:   
      typeset.append(current_type) 
    
  print("Number of labels in clusterDB is ", h_Weight.GetNbinsX() )
  print("Number of types in clusterDB is ", len(typeset) )
  
  
  histofile = gROOT.FindObject( histofilename )
  if histofile:
    histofile.Close()
  histofile = TFile( histofilename, 'RECREATE', 'Resolution plots created from ' + dbfilename )
  
  histomap_weight = {}
  histomap_u = {}
  histomap_v = {}
  histomap_sigu = {}
  histomap_sigv = {}
  
  for currenttype in  typeset:
  
    map_weight = {}
    map_u = {}
    map_v = {}
    map_sigu = {}
    map_sigv = {}
    
    for bin in range(1,h_Weight.GetNbinsX()+1):
    
      # The bin label decodes the cluster shape. The  
      # label contains tokens seperated by "D". The 
      # first token is the cluster size, all other 
      # are digits. Each digits contains three tokens
      # seperated by "."; nameley iu, iv and signal.
      # These can be converted to integers. 
       
      label = h_Weight.GetXaxis().GetBinLabel(bin)
      
      # We compute the type string by stripping all 
      # signal information from the label 
       
      current_type = ""
      for tok in re.split('D',label):
        addr = re.split('\.',tok)
        if len(addr) == 1:
          current_type += tok    
        else: 
          current_type += "D" + addr[0] + '.' + addr[1]   
      
      
      if current_type == currenttype:
        map_weight[label] = h_Weight.GetBinContent(bin)
        map_u[label] = h_U.GetBinContent(bin)
        map_v[label] = h_V.GetBinContent(bin)
        map_sigu[label] = TMath.Sqrt(h_Var_U.GetBinContent(bin))
        map_sigv[label] = TMath.Sqrt(h_Var_V.GetBinContent(bin)) 
      
  
    
    histofile.cd("")
    histofile.mkdir(currenttype)
    histofile.cd(currenttype)
    
    # these are the reprocessed histos for viewing
    LABELS = len(map_weight)
    
    histomap_weight[currenttype] = TH1F("hweight_"+currenttype,"",LABELS,0,LABELS)
    histomap_weight[currenttype].SetStats(0)
    histomap_weight[currenttype].SetFillColor(38)
    histomap_weight[currenttype].SetYTitle("weight")

    histomap_u[currenttype] = TH1F("hu_"+currenttype,"",LABELS,0,LABELS)
    histomap_u[currenttype].SetStats(0)
    histomap_u[currenttype].SetFillColor(38)
    histomap_u[currenttype].SetYTitle("offset u [mm]")
    
    histomap_v[currenttype] = TH1F("hv_"+currenttype,"",LABELS,0,LABELS)
    histomap_v[currenttype].SetStats(0)
    histomap_v[currenttype].SetFillColor(38)
    histomap_v[currenttype].SetYTitle("offset v [mm]")
    
    histomap_sigu[currenttype] = TH1F("hsigu_"+currenttype,"",LABELS,0,LABELS)
    histomap_sigu[currenttype].SetStats(0)
    histomap_sigu[currenttype].SetFillColor(38)
    histomap_sigu[currenttype].SetYTitle("cluster sigma u [mm]")
  
    histomap_sigv[currenttype] = TH1F("hsigv_"+currenttype,"",LABELS,0,LABELS)
    histomap_sigv[currenttype].SetStats(0)
    histomap_sigv[currenttype].SetFillColor(38)
    histomap_sigv[currenttype].SetYTitle("cluster sigma v [mm]")  

    
    for i, label in enumerate( map_weight.keys() ):
      
      histomap_weight[currenttype].SetBinContent(i+1, map_weight[label])
      histomap_u[currenttype].SetBinContent(i+1, map_u[label])
      histomap_v[currenttype].SetBinContent(i+1, map_v[label]) 
      histomap_sigu[currenttype].SetBinContent(i+1, map_sigu[label])  
      histomap_sigv[currenttype].SetBinContent(i+1, map_sigv[label]) 
      
      histomap_weight[currenttype].GetXaxis().SetBinLabel(i+1, label)
      histomap_u[currenttype].GetXaxis().SetBinLabel(i+1, label)
      histomap_v[currenttype].GetXaxis().SetBinLabel(i+1, label) 
      histomap_sigu[currenttype].GetXaxis().SetBinLabel(i+1, label)  
      histomap_sigv[currenttype].GetXaxis().SetBinLabel(i+1, label) 
    


  
  # summary histograms on type resolution 
  histofile.cd("")
  
  TYPES = len(typeset)

  htypes_sigu = TH1F("htypes_sigu","",TYPES,0,TYPES)
  htypes_sigu.SetStats(0)
  htypes_sigu.SetFillColor(38)
  htypes_sigu.SetYTitle("weighted cluster sigma u [mm]")  
  
  htypes_sigv = TH1F("htypes_sigv","",TYPES,0,TYPES)
  htypes_sigv.SetStats(0)
  htypes_sigv.SetFillColor(38)
  htypes_sigv.SetYTitle("weighted cluster sigma v [mm]")  
  
  htypes_weight = TH1F("htypes_weight","",TYPES,0,TYPES)
  htypes_weight.SetStats(0)
  htypes_weight.SetFillColor(38)
  htypes_weight.SetYTitle("weight")  

  for j, currenttype in enumerate( typeset ):
  
    htypes_sigu.GetXaxis().SetBinLabel(j+1, currenttype)  
    htypes_sigv.GetXaxis().SetBinLabel(j+1, currenttype)
    htypes_weight.GetXaxis().SetBinLabel(j+1, currenttype) 
  
    weightedTypeVarU = 0.0
    weightedTypeVarV = 0.0
    typeNorm = 0.0
     
    for bin in range(1,histomap_weight[currenttype].GetNbinsX()+1):  
      w = histomap_weight[currenttype].GetBinContent(bin)
      typeNorm += w
     
      weightedTypeVarU += w*TMath.Power(histomap_sigu[currenttype].GetBinContent(bin),2)
      weightedTypeVarV += w*TMath.Power(histomap_sigv[currenttype].GetBinContent(bin),2)
    
    htypes_weight.SetBinContent(j+1, typeNorm)

    if typeNorm >0: 
      weightedTypeVarU/=typeNorm
      weightedTypeVarV/=typeNorm
      htypes_sigu.SetBinContent(j+1, TMath.Sqrt(weightedTypeVarU))  
      htypes_sigv.SetBinContent(j+1, TMath.Sqrt(weightedTypeVarV))  
    else:
      htypes_sigu.SetBinContent(j+1, 0)  # invalid
      htypes_sigv.SetBinContent(j+1, 0)  # invalid
    
  
  # Standardized plots
  
  SDTYPES = len(standard_types)
  
  hsdtypes_sigu = TH1F("hsdtypes_sigu","",SDTYPES,0,SDTYPES)
  hsdtypes_sigu.SetStats(0)
  hsdtypes_sigu.SetFillColor(38)
  hsdtypes_sigu.SetYTitle("weighted cluster sigma u [mm]")
  
  hsdtypes_sigv = TH1F("hsdtypes_sigv","",SDTYPES,0,SDTYPES)
  hsdtypes_sigv.SetStats(0)
  hsdtypes_sigv.SetFillColor(38)
  hsdtypes_sigv.SetYTitle("weighted cluster sigma v [mm]")  

  hsdtypes_weight = TH1F("hsdtypes_weight","",SDTYPES,0,SDTYPES)
  hsdtypes_weight.SetStats(0)
  hsdtypes_weight.SetFillColor(38)
  hsdtypes_weight.SetYTitle("weight")  
  
  for j, currenttype in enumerate( standard_types ):
   
    hsdtypes_sigu.GetXaxis().SetBinLabel(j+1, currenttype)
    hsdtypes_sigv.GetXaxis().SetBinLabel(j+1, currenttype) 
    hsdtypes_weight.GetXaxis().SetBinLabel(j+1, currenttype) 
    
    weightedTypeVarU = 0.0
    weightedTypeVarV = 0.0
    typeNorm = 0.0
   
    if currenttype in histomap_weight:
      for bin in range(1,histomap_weight[currenttype].GetNbinsX()+1):  
      
        w = histomap_weight[currenttype].GetBinContent(bin)
        typeNorm += w
        
        weightedTypeVarU += w*TMath.Power(histomap_sigu[currenttype].GetBinContent(bin),2)
        weightedTypeVarV += w*TMath.Power(histomap_sigv[currenttype].GetBinContent(bin),2)
     
    
    hsdtypes_weight.SetBinContent(j+1, typeNorm)

    if typeNorm >0:
      weightedTypeVarU/=typeNorm
      weightedTypeVarV/=typeNorm
      hsdtypes_sigu.SetBinContent(j+1, TMath.Sqrt(weightedTypeVarU))  
      hsdtypes_sigv.SetBinContent(j+1, TMath.Sqrt(weightedTypeVarV))  
    else: 
      hsdtypes_sigu.SetBinContent(j+1, 0)  # invalid
      hsdtypes_sigv.SetBinContent(j+1, 0)  # invalid
    
  # summary histograms on overall resolution

  hweighted_sigma_sensor =  TH1F("hweighted_sigma_sensor","",2,0,2)
  hweighted_sigma_sensor.SetStats(0)
  hweighted_sigma_sensor.SetFillColor(38)
  hweighted_sigma_sensor.SetYTitle("cluster sigma [mm]") 
  
  hweighted_sigma_sensor.GetXaxis().SetBinLabel(1, "sigma u")  
  hweighted_sigma_sensor.GetXaxis().SetBinLabel(2, "sigma v") 

  weightedVarU = 0.0
  weightedVarV = 0.0
  norm = 0.0

  
  for bin in range(1,h_Weight.GetNbinsX()+1):  
    w = h_Weight.GetBinContent(bin)
    norm += w
  
    weightedVarU += w*h_Var_U.GetBinContent(bin)
    weightedVarV += w*h_Var_V.GetBinContent(bin)
  

  print( "Number of tracks used for calibration is " , norm )
  
  if norm >0:
    weightedVarU/=norm
    weightedVarV/=norm

    hweighted_sigma_sensor.SetBinContent(1, TMath.Sqrt(weightedVarU))  
    hweighted_sigma_sensor.SetBinContent(2, TMath.Sqrt(weightedVarV)) 
    
    print("Weighted clusterDB sigmaU [mm]: ", TMath.Sqrt(weightedVarU) )
    print("Weighted clusterDB sigmaV [mm]: ", TMath.Sqrt(weightedVarV) )
  else:
    hweighted_sigma_sensor.SetBinContent(1, 0)  # invalid
    hweighted_sigma_sensor.SetBinContent(2, 0)  # invalid
  

  histofile.Write()
  histofile.Close()

                      


