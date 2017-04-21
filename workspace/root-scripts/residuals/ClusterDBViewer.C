{ 

// Input file name
TString inputfile = "cal-files/run65-s8/clusterDB-DEP2.root";

// Output file name
TString histofile = "ViewClusterDB-s8.root";

//------------------------------------------------------------------------


//
// get histos from db
TFile *ftb = new TFile(inputfile);

TH1F * h_Weight = (TH1F*) ftb->Get("hDB_Weight");
h_Weight->SetDirectory(0);

TH1F * h_U = (TH1F*) ftb->Get("hDB_U");
h_U->SetDirectory(0);

TH1F * h_V = (TH1F*) ftb->Get("hDB_V");
h_V->SetDirectory(0);

TH1F * h_Var_U = (TH1F*) ftb->Get("hDB_Sigma2_U");
h_Var_U->SetDirectory(0);

TH1F * h_Var_V = (TH1F*) ftb->Get("hDB_Sigma2_V");
h_Var_V->SetDirectory(0);

// Close root  file
ftb->Close();

//
// First, we want to compute a list of all different cluster 
// types found in the db

std::set<TString> typeset;

for (int bin = 1; bin<= h_Weight->GetNbinsX(); bin++) {
  
  // The bin label decodes the cluster shape. The  
  // label contains tokens seperated by "D". The 
  // first token is the cluster size, all other 
  // are digits. Each digits contains three tokens
  // seperated by "."; nameley iu, iv and signal.
  // These can be converted to integers. 
  
  TString label(h_Weight->GetXaxis()->GetBinLabel(bin));
  
  // We compute the type string by stripping all 
  // signal information from the label 
  TString type("");
  
  TString tok;
  Ssiz_t from = 0;
  while (label.Tokenize(tok, from, "D")) {
    // Analyse tok
    if (TPRegexp("^\\b(\\d+).(\\d+).(\\d+)\\b$").MatchB(tok)) {
      TPRegexp("^\\b(\\d+).(\\d+).(\\d+)\\b$").Substitute(tok,"$1.$2");
      type.Append("D"+tok);
    }
    if (TPRegexp("^\\b(\\d+)\\b$").MatchB(tok)) {
      type.Append(tok);
    }
  }
  typeset.insert(type); 
}

cout << "Number of labels in clusterDB is " << h_Weight->GetNbinsX() << endl;
cout << "Number of types in clusterDB is " << typeset.size() << endl;

//
// output root file with final cluster db histos 
TFile *fDB = new TFile(histofile,"RECREATE");

std::map<TString, TH1F *> histomap_weight;
std::map<TString, TH1F *> histomap_u;
std::map<TString, TH1F *> histomap_v;
std::map<TString, TH1F *> histomap_sigu;
std::map<TString, TH1F *> histomap_sigv;

for (auto currenttype : typeset ) {
  
 
  std::map<TString, float> map_weight;
  std::map<TString, float> map_u;
  std::map<TString, float> map_v;
  std::map<TString, float> map_sigu;
  std::map<TString, float> map_sigv;
    

  for (int bin = 1; bin<= h_Weight->GetNbinsX(); bin++) {
  
    // The bin label decodes the cluster shape. The  
    // label contains tokens seperated by "D". The 
    // first token is the cluster size, all other 
    // are digits. Each digits contains three tokens
    // seperated by "."; nameley iu, iv and signal.
    // These can be converted to integers. 
     
    TString label(h_Weight->GetXaxis()->GetBinLabel(bin));
    
    // We compute the type string by stripping all 
    // signal information from the label 
    TString type("");
    
    TString tok;
    Ssiz_t from = 0;
    while (label.Tokenize(tok, from, "D")) {
      // Analyse tok
      if (TPRegexp("^\\b(\\d+).(\\d+).(\\d+)\\b$").MatchB(tok)) {
        TPRegexp("^\\b(\\d+).(\\d+).(\\d+)\\b$").Substitute(tok,"$1.$2");
        type.Append("D"+tok);
      }
      if (TPRegexp("^\\b(\\d+)\\b$").MatchB(tok)) {
        type.Append(tok);
      }
    }
    
    if ( type == currenttype ) {
      map_weight[label] = h_Weight->GetBinContent(bin);
      map_u[label] = h_U->GetBinContent(bin);
      map_v[label] = h_V->GetBinContent(bin);
      map_sigu[label] = TMath::Sqrt(h_Var_U->GetBinContent(bin));
      map_sigv[label] = TMath::Sqrt(h_Var_V->GetBinContent(bin)); 
    }   
  }
  
  fDB->cd("");
  fDB->mkdir(currenttype);
  fDB->cd(currenttype);

  //
  // these are the reprocessed histos for viewing
  int LABELS = map_weight.size();
  
  histomap_weight[currenttype] = new TH1F(TString("hweight_")+currenttype,"",LABELS,0,LABELS);
  histomap_weight[currenttype]->SetStats(0);
  histomap_weight[currenttype]->SetFillColor(38);
  histomap_weight[currenttype]->SetYTitle("weight");  

  histomap_u[currenttype] = new TH1F(TString("hu_")+currenttype,"",LABELS,0,LABELS);
  histomap_u[currenttype]->SetStats(0);
  histomap_u[currenttype]->SetFillColor(38);
  histomap_u[currenttype]->SetYTitle("offset u [mm]");  

  histomap_v[currenttype] =  new TH1F(TString("hv_")+currenttype,"",LABELS,0,LABELS);
  histomap_v[currenttype]->SetStats(0);
  histomap_v[currenttype]->SetFillColor(38);
  histomap_v[currenttype]->SetYTitle("offset v [mm]");  
  
  histomap_sigu[currenttype] = new TH1F(TString("hsigu_")+currenttype,"",LABELS,0,LABELS);
  histomap_sigu[currenttype]->SetStats(0);
  histomap_sigu[currenttype]->SetFillColor(38);
  histomap_sigu[currenttype]->SetYTitle("cluster sigma u [mm]");  
  
  histomap_sigv[currenttype] = new TH1F(TString("hsigv_")+currenttype,"",LABELS,0,LABELS);
  histomap_sigv[currenttype]->SetStats(0);
  histomap_sigv[currenttype]->SetFillColor(38);
  histomap_sigv[currenttype]->SetYTitle("cluster sigma v [mm]");  

  int i = 0; 
  for (auto labelAndWeight : map_weight) {
    i++;
    TString label = labelAndWeight.first;
    histomap_weight[currenttype]->SetBinContent(i, map_weight[label]);
    histomap_u[currenttype]->SetBinContent(i, map_u[label]);  
    histomap_v[currenttype]->SetBinContent(i, map_v[label]); 
    histomap_sigu[currenttype]->SetBinContent(i, map_sigu[label]);  
    histomap_sigv[currenttype]->SetBinContent(i, map_sigv[label]); 
  
    histomap_weight[currenttype]->GetXaxis()->SetBinLabel(i, label);
    histomap_u[currenttype]->GetXaxis()->SetBinLabel(i, label);  
    histomap_v[currenttype]->GetXaxis()->SetBinLabel(i, label); 
    histomap_sigu[currenttype]->GetXaxis()->SetBinLabel(i, label);  
    histomap_sigv[currenttype]->GetXaxis()->SetBinLabel(i, label); 
  }

}

//
// summary histograms on type resolution 
fDB->cd("");

int TYPES = typeset.size();

htypes_sigu = new TH1F("htypes_sigu","",TYPES,0,TYPES);
htypes_sigu->SetStats(0);
htypes_sigu->SetFillColor(38);
htypes_sigu->SetYTitle("weighted cluster sigma u [mm]");  
  
htypes_sigv = new TH1F("htypes_sigv","",TYPES,0,TYPES);
htypes_sigv->SetStats(0);
htypes_sigv->SetFillColor(38);
htypes_sigv->SetYTitle("weighted cluster sigma v [mm]");  

int j = 0;
for (auto currenttype : typeset ) {
  j++;

  htypes_sigu->GetXaxis()->SetBinLabel(j, currenttype);  
  htypes_sigv->GetXaxis()->SetBinLabel(j, currenttype); 
  
  double weightedTypeVarU = 0;
  double weightedTypeVarV = 0;
  double typeNorm = 0;
  
  for (int bin = 1; bin<= histomap_weight[currenttype]->GetNbinsX(); bin++) {
    
    double w = histomap_weight[currenttype]->GetBinContent(bin);
    typeNorm += w;
     
    weightedTypeVarU += w*TMath::Power(histomap_sigu[currenttype]->GetBinContent(bin),2);
    weightedTypeVarV += w*TMath::Power(histomap_sigv[currenttype]->GetBinContent(bin),2);
  }

  if (typeNorm >0) {
    weightedTypeVarU/=typeNorm;
    weightedTypeVarV/=typeNorm;
    htypes_sigu->SetBinContent(j, TMath::Sqrt(weightedTypeVarU));  
    htypes_sigv->SetBinContent(j, TMath::Sqrt(weightedTypeVarV)); 
  } else {
    htypes_sigu->SetBinContent(j, 0);  // invalid
    htypes_sigv->SetBinContent(j, 0);  // invalid
  }
}

//
// summary histograms on overall resolution

hweighted_sigma_sensor = new TH1F("hweighted_sigma_sensor","",2,0,2);
hweighted_sigma_sensor->SetStats(0);
hweighted_sigma_sensor->SetFillColor(38);
hweighted_sigma_sensor->SetYTitle("cluster sigma [mm]");  

hweighted_sigma_sensor->GetXaxis()->SetBinLabel(1, "sigma u");  
hweighted_sigma_sensor->GetXaxis()->SetBinLabel(2, "sigma v"); 

double weightedVarU = 0;
double weightedVarV = 0;
double norm = 0;

for (int bin = 1; bin<= h_Weight->GetNbinsX(); bin++) {

  double w = h_Weight->GetBinContent(bin);
  norm += w;
  
  weightedVarU += w*h_Var_U->GetBinContent(bin);
  weightedVarV += w*h_Var_V->GetBinContent(bin);
}

cout << "Number of tracks used for calibration is " << norm << endl;

if (norm >0) {
  weightedVarU/=norm;
  weightedVarV/=norm;

  hweighted_sigma_sensor->SetBinContent(1, TMath::Sqrt(weightedVarU));  
  hweighted_sigma_sensor->SetBinContent(2, TMath::Sqrt(weightedVarV)); 

  cout << "Weighted clusterDB sigmaU [mm]: " << TMath::Sqrt(weightedVarU) << endl;
  cout << "Weighted clusterDB sigmaV [mm]: " << TMath::Sqrt(weightedVarV) << endl;
} else {
  hweighted_sigma_sensor->SetBinContent(1, 0);  // invalid
  hweighted_sigma_sensor->SetBinContent(2, 0);  // invalid
}

fDB->Write();
fDB->Close();

}
