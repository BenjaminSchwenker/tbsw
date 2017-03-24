{ 



// Define a cluster shape 
//TRegexp e("^1D0\\.0\\.[0-9]+$");
TRegexp e("^2D0\\.0\\.[0-9]+D0\\.1\\.[0-9]+$");

//
// get histos from db
TFile *ftb = new TFile("clusterDB-DEP2-run000065-run65.root");

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


//
// Start reprocess clusters (filter + sort)

auto mycomp = [](const TString&a, const TString& b) { 
  
  TObjArray *ta = a.Tokenize("D"); 
  TObjArray *tb = b.Tokenize("D");

  TString a1 = ((TObjString *)(ta->At(1)))->String();
  TString b1 = ((TObjString *)(tb->At(1)))->String();
   
  TObjArray *ta1 = a1.Tokenize("."); 
  TObjArray *tb1 = b1.Tokenize("."); 
  
  int a_charge =  ((TObjString *)(ta1->At(2)))->String().Atoi();  
  int b_charge =  ((TObjString *)(tb1->At(2)))->String().Atoi();

  return a_charge < b_charge;   
};


//
// temporary maps to sort cluster db histos
std::map<TString, float, decltype(mycomp)> map_id(mycomp);
std::map<TString, float, decltype(mycomp)> map_u(mycomp);
std::map<TString, float, decltype(mycomp)> map_v(mycomp);
std::map<TString, float, decltype(mycomp)> map_sigu(mycomp);
std::map<TString, float, decltype(mycomp)> map_sigv(mycomp);


cout << "total number of entries in clusterDB is " << h_ID->GetNbinsX() << endl;

for (int bin = 1; bin<= h_ID->GetNbinsX(); bin++) {
  
  TString x(h_ID->GetXaxis()->GetBinLabel(bin));
   
  TObjArray *tx = x.Tokenize("D");
  int nDigits =  ((TObjString *)(tx->At(0)))->String().Atoi();
   
  if (x.Contains(e)) {
  //if (nDigits == 1) {
    map_id[x] = h_ID->GetBinContent(bin);
    map_u[x] = h_U->GetBinContent(bin);
    map_v[x] = h_V->GetBinContent(bin);
    map_sigu[x] = TMath::Sqrt(h_Var_U->GetBinContent(bin));
    map_sigv[x] = TMath::Sqrt(h_Var_V->GetBinContent(bin));    
  }
}

//
// cluster db file  
TFile *fDB = new TFile("ClusterDBViewer.root","RECREATE");

//
// these are the reprocessed histos for viewing
int NCLUSTERS = map_id.size();

TH1F *hid = new TH1F("hid","",NCLUSTERS,0,NCLUSTERS);
hid->SetStats(0);
hid->SetFillColor(38);
hid->SetYTitle("weight");  

TH1F *hu = new TH1F("hu","",NCLUSTERS,0,NCLUSTERS);
hu->SetStats(0);
hu->SetFillColor(38);
hu->SetYTitle("offset u [mm]");  

TH1F *hv = new TH1F("hv","",NCLUSTERS,0,NCLUSTERS);
hv->SetStats(0);
hv->SetFillColor(38);
hv->SetYTitle("offset v [mm]");  

TH1F *hsigu = new TH1F("hsigu","",NCLUSTERS,0,NCLUSTERS);
hsigu->SetStats(0);
hsigu->SetFillColor(38);
hsigu->SetYTitle("cluster sigma u [mm]");  

TH1F *hsigv = new TH1F("hsigv","",NCLUSTERS,0,NCLUSTERS);
hsigv->SetStats(0);
hsigv->SetFillColor(38);
hsigv->SetYTitle("cluster sigma v [mm]");  


int i = 0; 
for (auto idAndWeight : map_id) {
  i++;
  TString ID = idAndWeight.first;
  hid->SetBinContent(i, map_id[ID]);
  hu->SetBinContent(i, map_u[ID]);  
  hv->SetBinContent(i, map_v[ID]); 
  hsigu->SetBinContent(i, map_sigu[ID]);  
  hsigv->SetBinContent(i, map_sigv[ID]); 
  
  hid->GetXaxis()->SetBinLabel(i, ID);
  hu->GetXaxis()->SetBinLabel(i, ID);  
  hv->GetXaxis()->SetBinLabel(i, ID); 
  hsigu->GetXaxis()->SetBinLabel(i, ID);  
  hsigv->GetXaxis()->SetBinLabel(i, ID); 
}

fDB->Write();
fDB->Close();


}
