{ 

gROOT->Reset(); 
gROOT->SetStyle("Plain"); 

gStyle->SetOptFit();
gStyle->SetPalette(1);
gStyle->SetCanvasColor(0);
gStyle->SetTitleFillColor(0);
gStyle->SetPaintTextFormat("1f");
gStyle->SetOptStat(0);


//
// event data 


TFile *fmc50 = new TFile("Histos-50mu.root");
TTree *tmc50track = (TTree*) fmc50->Get("Track");

TFile *fmc75 = new TFile("Histos-75mu.root");
TTree *tmc75track = (TTree*) fmc75->Get("Track");
 
TFile *ftb = new TFile("Histos-71074.root");
TTree *ttbtrack = (TTree*) ftb->Get("Track");

// ----------------------------------------------------------------------

int n = 50; 
double LSB = 180; 

TH1D hpass("hpass","hpass", 1,0,32); 
TH1D htotal("htotal","htotal",1,0,32); 

// loop over tb data 

double seedcut[n]= {0};
double eseedcut[n]= {0};
double effi[n] =  {0};
double error[n]=  {0};

for (int i = 0; i<n; ++i) {

  hpass.Reset();
  htotal.Reset();
  
  seedcut[i]= i*LSB;
  eseedcut[i] =  i*10;
  cout << "seed cut " << seedcut[i] << " electrons "<< endl;  
  
  std::stringstream TotalSelector;
  TotalSelector << "               1x1Quality == 0 && chi2 < 100  && row_fit > 0 && row_fit < 63 && col_fit > 0 && col_fit < 31 && startgate != 0 && startgate != 15 && iEvt > 25000 && nTelTracks == 1"; 
  
  std::stringstream PassSelector;
  PassSelector << " hasHit == 0 && 1x1Quality == 0 && chi2 < 100  && row_fit > 0 && row_fit < 63 && col_fit > 0 && col_fit < 31 && startgate != 0 && startgate != 15 && iEvt > 25000 && nTelTracks == 1 && 180*seedCharge > " <<  seedcut[i];
   
   
  ttbtrack->Draw( "col_fit >> htotal" , TotalSelector.str().c_str() );
  ttbtrack->Draw( "col_fit >> hpass", PassSelector.str().c_str()  );
 
  TGraphAsymmErrors heffi( &hpass, &htotal, "w" ); 
 
  error[i]= heffi.GetErrorY(0);
  effi[i] = heffi.GetY()[0]; 
   
  cout << "effi " << effi[i] << " +/- " << error[i] << endl; 
  
}

// loop over 50mu Si data 

double seedcut_mc50[n]= {0};
double eseedcut_mc50[n]= {0};
double effi_mc50[n] =  {0};
double error_mc50[n]=  {0};

for (int i = 0; i<n; ++i) {

  hpass.Reset();
  htotal.Reset();
 
  seedcut_mc50[i]= i*LSB; 
  cout << "seed cut " << seedcut_mc50[i] << " electrons "<< endl; 
  

  std::stringstream TotalSelector;
  TotalSelector << "               row_fit > 0 && row_fit < 63 && col_fit > 0 && col_fit < 31 && nTelTracks == 1"; 
  
  std::stringstream PassSelector;
  PassSelector << " hasHit == 0 && row_fit > 0 && row_fit < 63 && col_fit > 0 && col_fit < 31 && nTelTracks == 1 && seedCharge > " << seedcut_mc50[i];
   
  
 
  tmc50track->Draw( "col_fit >> htotal" , TotalSelector.str().c_str() );
  tmc50track->Draw( "col_fit >> hpass", PassSelector.str().c_str()  );
 
  TGraphAsymmErrors heffi( &hpass, &htotal, "w" ); 
 
  
  error_mc50[i]= heffi.GetErrorY(0);
  effi_mc50[i] = heffi.GetY()[0]; 
    
  cout << "effi " << effi_mc50[i] << " +/- " << error_mc50[i] << endl; 

}


double seedcut_mc75[n]= {0};
double eseedcut_mc75[n]= {0};
double effi_mc75[n] =  {0};
double error_mc75[n]=  {0};


for (int i = 0; i<n; ++i) {

  hpass.Reset();
  htotal.Reset();
 
  seedcut_mc75[i]= i*LSB; 
  cout << "seed cut " << seedcut_mc75[i] << " electrons "<< endl; 
  

  std::stringstream TotalSelector;
  TotalSelector << "               row_fit > 0 && row_fit < 63 && col_fit > 0 && col_fit < 31 && nTelTracks == 1"; 
  
  std::stringstream PassSelector;
  PassSelector << " hasHit == 0 && row_fit > 0 && row_fit < 63 && col_fit > 0 && col_fit < 31 && nTelTracks == 1 && seedCharge > " << seedcut_mc75[i];
   
  
 
  tmc75track->Draw( "col_fit >> htotal" , TotalSelector.str().c_str() );
  tmc75track->Draw( "col_fit >> hpass", PassSelector.str().c_str()  );
 
  TGraphAsymmErrors heffi( &hpass, &htotal, "w" ); 
 
  
  error_mc75[i]= heffi.GetErrorY(0);
  effi_mc75[i] = heffi.GetY()[0]; 
    
  cout << "effi " << effi_mc75[i] << " +/- " << error_mc75[i] << endl; 

}





TCanvas * c01  = new TCanvas("c01","c01",800,400);
c01->SetLeftMargin(0.2);
c01->SetRightMargin(0.1);
c01->SetTopMargin(0.1);
c01->SetBottomMargin(0.16);


TH1F *vSizeFrame = c01.DrawFrame(0,0,9000,1);
vSizeFrame->SetXTitle("seed threshold [e]");
vSizeFrame->SetYTitle("efficiency");
vSizeFrame->GetXaxis()->SetTitleSize(0.06);
vSizeFrame->GetYaxis()->SetTitleSize(0.06);
vSizeFrame->GetXaxis()->SetLabelSize(0.06);
vSizeFrame->GetYaxis()->SetLabelSize(0.06);
//vSizeFrame->GetXaxis()->SetTitleOffset(1.3);
vSizeFrame->GetYaxis()->SetTitleOffset(1.7);



gr = new TGraphErrors(n,seedcut,effi,eseedcut,error);
gr->SetName("gr");
gr->SetMarkerColor(kBlue);
gr->SetMarkerStyle(21);
gr->Draw("LP");

grmc50 = new TGraphErrors(n,seedcut_mc50,effi_mc50,eseedcut_mc50,error_mc50);
grmc50->SetName("grmc50");
grmc50->SetMarkerColor(kRed);
grmc50->SetMarkerStyle(20);
grmc50->Draw("LP");

grmc75 = new TGraphErrors(n,seedcut_mc75,effi_mc75,eseedcut_mc75,error_mc75);
grmc75->SetName("grmc75");
grmc75->SetMarkerColor(kGreen);
grmc75->SetMarkerStyle(23);
grmc75->Draw("LP");

TLegend* l3 = new TLegend(0.4,0.32,0.6,0.6);
l3->SetFillColor(kWhite); 
l3->SetBorderSize(0);
l3->AddEntry("gr",   "TB 2012 H4.1.15","lep");
l3->AddEntry("grmc50", "Digitizer [3GeV e^{-}perp, 50#mum]","lep");
l3->AddEntry("grmc75", "Digitizer [3GeV e^{-}perp, 75#mum]","lep");
l3->Draw();


c01->Modified();
c01->Update();




}
