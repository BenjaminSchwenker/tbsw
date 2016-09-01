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
// Parameters
int n = 50; 
double LSB = 180; 
double GAIN = 180; 

// 
// Analysis cuts 
TCut basecut = "1x1Quality == 0 && chi2 < 100 && iEvt > 25000 && nTelTracks == 1 && startgate != 0 && startgate != 15";
TCut roicut  = "col_fit > 0 && col_fit < 31 && row_fit > 0 && row_fit < 63";

//
// event data 
TFile *ftb = new TFile("Histos-71074.root");
TTree *ttbtrack = (TTree*) ftb->Get("Track");

// ----------------------------------------------------------------------

double seedcut[n]= {0};
double eseedcut[n]= {0};
double effi[n] =  {0};
double error[n]=  {0};

TH1D hpass("hpass","hpass", 1,0,32); 
TH1D htotal("htotal","htotal",1,0,32); 


for (int i = 0; i<n; ++i) {

  hpass.Reset();
  htotal.Reset();
  
  seedcut[i]= i*LSB;
  eseedcut[i] =  i*10;
  cout << "seed cut " << seedcut[i] << " electrons "<< endl;  

  TCut matched = Form(" hasHit == 0 && %f*seedCharge > %f",GAIN,seedcut[i]);
  
  ttbtrack->Draw( "col_fit >> htotal" , basecut && roicut );
  ttbtrack->Draw( "col_fit >> hpass" ,  basecut && roicut && matched);
  
  TGraphAsymmErrors heffi( &hpass, &htotal, "w" ); 
 
  error[i]= heffi.GetErrorY(0);
  effi[i] = heffi.GetY()[0]; 
   
  cout << "effi " << effi[i] << " +/- " << error[i] << endl; 
  
}



TCanvas * c01  = new TCanvas("c01","c01",600,400);
c01->SetLeftMargin(0.2);
c01->SetRightMargin(0.1);
c01->SetTopMargin(0.1);
c01->SetBottomMargin(0.16);


TH1F *vSizeFrame = c01->DrawFrame(0,0,12000,1);
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

c01->Modified();
c01->Update();

TCanvas * c02  = new TCanvas("c02","c02",600,400);
c02->SetLeftMargin(0.2);
c02->SetRightMargin(0.1);
c02->SetTopMargin(0.1);
c02->SetBottomMargin(0.16);


TH1F *vSizeFrame2 = c02->DrawFrame(0,0.80,2400,1);
vSizeFrame2->SetXTitle("seed threshold [e]");
vSizeFrame2->SetYTitle("efficiency");
vSizeFrame2->GetXaxis()->SetTitleSize(0.06);
vSizeFrame2->GetYaxis()->SetTitleSize(0.06);
vSizeFrame2->GetXaxis()->SetLabelSize(0.06);
vSizeFrame2->GetYaxis()->SetLabelSize(0.06);
//vSizeFrame2->GetXaxis()->SetTitleOffset(1.3);
vSizeFrame2->GetYaxis()->SetTitleOffset(1.7);



gr2 = new TGraphErrors(n,seedcut,effi,eseedcut,error);
gr2->SetName("gr2");
gr2->SetMarkerColor(kBlue);
gr2->SetMarkerStyle(21);
gr2->Draw("LP");

c02->Modified();
c02->Update();



}
