{ 

gROOT->Reset(); 
gROOT->SetStyle("Plain"); 

gStyle->SetOptFit();
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetTitleFillColor(0);
gStyle->SetCanvasColor(0);
gStyle->SetPaintTextFormat("1f");

//
// MC data 
TFile *fmc = new TFile("root-final/tuneA.root");
TTree *tmc = (TTree*) fmc->Get("Hit");

//
// TB data 
TFile *ftb = new TFile("root-final/tbtiltscan.root");
TTree *ttb = (TTree*) ftb->Get("Hit");

TCut TBcut = "hasTrack == 0 && chi2 <100 && ndof == 8 && iRun == 409";
TCut MCcut = "hasTrack == 0 && iRun == 4";

Double_t telerru =  4.04;  // 2.26; 2.79; 3.24; 4.05; 5.5; 6.04; 6.62; 9.14 
Double_t telerrv =  3.5;    // 2.26; 2.75; 3.04; 3.5; 4.15; 4.2; 4.3; 4.7

//Int_t iRun_TB = 435;  // 0deg 
//Int_t iRun_MC = 0; 

//Int_t iRun_TB = 420;  // 10deg
//Int_t iRun_MC = 1;  

//Int_t iRun_TB = 416;  // 20deg
//Int_t iRun_MC = 2;  

//Int_t iRun_TB = 409;  // 30deg
//Int_t iRun_MC = 4;  

//Int_t iRun_TB = 393;  // 40deg
//Int_t iRun_MC = 6;  

//Int_t iRun_TB = 390;  // 45deg
//Int_t iRun_MC = 7;  

//Int_t iRun_TB = 384;  // 50deg
//Int_t iRun_MC = 8;  

//Int_t iRun_TB = 378;  // 60deg
//Int_t iRun_MC = 9;  



// ----------------------------------------------------------------------

float noiseU;

TBranch *bnoiseU = tmc->Branch("noiseU",&noiseU,"noiseU/F");
Long64_t nentries = tmc->GetEntries();
for (Long64_t i=0;i<nentries;i++) {
  noiseU = telerru*gRandom->Gaus(0,1); 
  bnoiseU->Fill();
}

float noiseV;

TBranch *bnoiseV = tmc->Branch("noiseV",&noiseV,"noiseV/F");
Long64_t nentries = tmc->GetEntries();
for (Long64_t i=0;i<nentries;i++) {
  noiseV = telerrv*gRandom->Gaus(0,1);
  bnoiseV->Fill();
}
  
// ----------------------------------------------------------------------


Double_t Norm; 

TCanvas * c9 = new TCanvas("c9","c9",500,500);   
c9->SetLeftMargin(0.22);
//c9->SetRightMargin(0.3);
c9->SetTopMargin(0.1);
c9->SetBottomMargin(0.16);


ttb->Draw("(u_hit - u_fit)*1000 >>htb(251,-105.,+105.)", TBcut );
htb->SetName("htb"); 
htb->SetXTitle("u residual [#mum]");
htb->SetYTitle("probability [1/#mum]"); 
htb->SetTitle(""); 
htb->SetLineColor(kBlack);
htb->SetLineWidth(2);
htb->GetYaxis()->SetTitleOffset(1.3);

htb->Sumw2();
Norm = htb->GetEntries()*htb->GetBinWidth(1);
if ( Norm != 0 )   htb->Scale( 1.0/Norm ); 


htb->GetXaxis()->SetTitleSize(0.06);
htb->GetYaxis()->SetTitleSize(0.06);
htb->GetXaxis()->SetLabelSize(0.06);
htb->GetYaxis()->SetLabelSize(0.06);
htb->GetYaxis()->SetTitleOffset(1.8);



tmc->Draw("(v_hit - v_fit)*1000 + noiseV >>hmc(251,-105.,+105.)", MCcut);   // u-v swapped
hmc->SetName("hmc"); 
hmc->SetXTitle("residual u [#mum]");
hmc->SetYTitle("probability [1/#mum]"); 
hmc->SetTitle(""); 
hmc->SetOption("hist");
hmc->SetFillColor(kGray);
hmc->GetYaxis()->SetTitleOffset(1.3);

hmc->Sumw2();
Norm = hmc->GetEntries()*hmc->GetBinWidth(1);
if ( Norm != 0 )   hmc->Scale( 1.0/Norm ); 


hmc->GetXaxis()->SetTitleSize(0.06);
hmc->GetYaxis()->SetTitleSize(0.06);
hmc->GetXaxis()->SetLabelSize(0.06);
hmc->GetYaxis()->SetLabelSize(0.06);
hmc->GetYaxis()->SetTitleOffset(1.8);


hmc->Draw();
htb->Draw("same");

//TLegend* l3 = new TLegend(0.6,0.5,0.85,0.7);
TLegend* l3 = new TLegend(0.45,0.7,0.85,0.88);
l3->SetFillColor(kWhite);
l3->SetBorderSize(0);
l3->AddEntry("hmc","Simulation: RMS=XY#mum","f");
l3->AddEntry("htb","Module A: RMS=XY#mum","lep");
l3->Draw();


TCanvas * c10 = new TCanvas("c10","c10",500,500);   
c10->SetLeftMargin(0.22);
//c10->SetRightMargin(0.3);
c10->SetTopMargin(0.1);
c10->SetBottomMargin(0.16);

ttb->Draw("(v_hit - v_fit)*1000 >>htbv(251,-105.,+105.)", TBcut);
htbv->SetName("htbv"); 
htbv->SetXTitle("residual v [#mum]");
htbv->SetYTitle("probability [1/#mum]"); 
htbv->SetTitle(""); 
htbv->SetLineColor(kBlack);
htbv->SetLineWidth(2);
htbv->GetYaxis()->SetTitleOffset(1.3);

htbv->Sumw2();
Norm = htbv->GetEntries()*htbv->GetBinWidth(1);
if ( Norm != 0 )   htbv->Scale( 1.0/Norm ); 


htbv->GetXaxis()->SetTitleSize(0.06);
htbv->GetYaxis()->SetTitleSize(0.06);
htbv->GetXaxis()->SetLabelSize(0.06);
htbv->GetYaxis()->SetLabelSize(0.06);
htbv->GetYaxis()->SetTitleOffset(1.8);


tmc->Draw("(u_hit - u_fit)*1000 + noiseU >>hmcv(251,-105.,+105.)", MCcut);   // u-v swapped
hmcv->SetName("hmcv"); 
hmcv->SetXTitle("residual v [#mum]");
hmcv->SetYTitle("probability [1/#mum]"); 
hmcv->SetTitle(""); 
hmcv->SetOption("hist");
hmcv->SetFillColor(kGray);
hmcv->GetYaxis()->SetTitleOffset(1.3);

hmcv->Sumw2();
Norm = hmcv->GetEntries()*hmcv->GetBinWidth(1);
if ( Norm != 0 )   hmcv->Scale( 1.0/Norm ); 


hmcv->GetXaxis()->SetTitleSize(0.06);
hmcv->GetYaxis()->SetTitleSize(0.06);
hmcv->GetXaxis()->SetLabelSize(0.06);
hmcv->GetYaxis()->SetLabelSize(0.06);
hmcv->GetYaxis()->SetTitleOffset(1.8);

hmcv->Draw();
htbv->Draw("same");

//TLegend* l4 = new TLegend(0.6,0.5,0.85,0.7);
TLegend* l4 = new TLegend(0.45,0.7,0.85,0.88);
l4->SetFillColor(kWhite); 
l4->SetBorderSize(0);
l4->AddEntry("hmcv","Simulation: RMS=XY#mum ","f");
l4->AddEntry("htbv","Module A: RMS=XY#mum","lep");
l4->Draw();





TCanvas c12;

c12->SetLeftMargin(0.1);
c12->SetRightMargin(0.15);


tmc->Draw("((u_hit - u_fit)*1000 + noiseU):((v_hit - v_fit)*1000 + noiseV) >>hmc2d(71,-65.,+65.,71,-65.,+65.)", MCcut);  // u-v swapped

hmc2d->SetName("hmc2d"); 
hmc2d->SetXTitle("MC residual u [#mum]");
hmc2d->SetYTitle("MC residual v [#mum]"); 
//hmc2d->SetZTitle("tracks"); 
hmc2d->SetZTitle("probability [1/#mum^{2}]"); 
hmc2d->SetTitle(""); 
//htb2d->GetYaxis()->SetTitleOffset(1.3);

Norm = hmc2d->GetEntries()*hmc2d->GetXaxis()->GetBinWidth(1)*hmc2d->GetYaxis()->GetBinWidth(1);
if ( Norm != 0 )   hmc2d->Scale( 1.0/Norm ); 


hmc2d->Draw("colz");

TCanvas c11;

c11->SetLeftMargin(0.1);
c11->SetRightMargin(0.15);


ttb->Draw("(v_fit - v_hit)*1000:(u_fit - u_hit)*1000  >>htb2d(71,-65.,+65.,71,-65.,+65.)", TBcut);

htb2d->SetName("htb2d"); 
htb2d->SetXTitle("TB residual u [#mum]");
htb2d->SetYTitle("TB residual v [#mum]"); 
//htb2d->SetZTitle("tracks"); 
htb2d->SetZTitle("probability [1/#mum^{2}]"); 
htb2d->SetTitle(""); 
//htb2d->GetYaxis()->SetTitleOffset(1.3);

Norm = htb2d->GetEntries()*htb2d->GetXaxis()->GetBinWidth(1)*htb2d->GetYaxis()->GetBinWidth(1);
if ( Norm != 0 )   htb2d->Scale( 1.0/Norm ); 



htb2d->Draw("colz");


}
