{ 

//
// get histos from db
TFile *ftb = new TFile("cal-files/vxd/clusterDB-MC.root");

TH1F * h_ID = (TH1F*) ftb->Get("hDB_ID");
h_ID->SetDirectory(0);

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

cout << "total number of labels in clusterDB is " << h_ID->GetNbinsX() << endl;

//
// First, we want to compute a list of all different cluster 
// types found in the db

std::set<TString> typeset;

for (int bin = 1; bin<= h_ID->GetNbinsX(); bin++) {
  
  // The bin label decodes the cluster shape. The  
  // label contains tokens seperated by "D". The 
  // first token is the cluster size, all other 
  // are digits. Each digits contains three tokens
  // seperated by "."; nameley iu, iv and signal.
  // These can be converted to integers. 
  
  TString label(h_ID->GetXaxis()->GetBinLabel(bin));
  
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

cout << "total number of types in clusterDB is " << typeset.size() << endl;

//
// output root file with final cluster db histos 
TFile *fDB = new TFile("ClusterDBViewer.root","RECREATE");

std::map<TString, TH1F *> histomap_id;
std::map<TString, TH1F *> histomap_u;
std::map<TString, TH1F *> histomap_v;
std::map<TString, TH1F *> histomap_sigu;
std::map<TString, TH1F *> histomap_sigv;

for (auto currenttype : typeset ) {
  
 
  std::map<TString, float> map_id;
  std::map<TString, float> map_u;
  std::map<TString, float> map_v;
  std::map<TString, float> map_sigu;
  std::map<TString, float> map_sigv;
    

  for (int bin = 1; bin<= h_ID->GetNbinsX(); bin++) {
  
    // The bin label decodes the cluster shape. The  
    // label contains tokens seperated by "D". The 
    // first token is the cluster size, all other 
    // are digits. Each digits contains three tokens
    // seperated by "."; nameley iu, iv and signal.
    // These can be converted to integers. 
     
    TString label(h_ID->GetXaxis()->GetBinLabel(bin));
    
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
    
    if (type == currenttype ) {
      map_id[label] = h_ID->GetBinContent(bin);
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
  int LABELS = map_id.size();
  
  histomap_id[currenttype] = new TH1F(TString("hid_")+currenttype,"",LABELS,0,LABELS);
  histomap_id[currenttype]->SetStats(0);
  histomap_id[currenttype]->SetFillColor(38);
  histomap_id[currenttype]->SetYTitle("weight");  

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
  for (auto idAndWeight : map_id) {
    i++;
    TString ID = idAndWeight.first;
    histomap_id[currenttype]->SetBinContent(i, map_id[ID]);
    histomap_u[currenttype]->SetBinContent(i, map_u[ID]);  
    histomap_v[currenttype]->SetBinContent(i, map_v[ID]); 
    histomap_sigu[currenttype]->SetBinContent(i, map_sigu[ID]);  
    histomap_sigv[currenttype]->SetBinContent(i, map_sigv[ID]); 
  
    histomap_id[currenttype]->GetXaxis()->SetBinLabel(i, ID);
    histomap_u[currenttype]->GetXaxis()->SetBinLabel(i, ID);  
    histomap_v[currenttype]->GetXaxis()->SetBinLabel(i, ID); 
    histomap_sigu[currenttype]->GetXaxis()->SetBinLabel(i, ID);  
    histomap_sigv[currenttype]->GetXaxis()->SetBinLabel(i, ID); 
  }

}


fDB->Write();
fDB->Close();

}
