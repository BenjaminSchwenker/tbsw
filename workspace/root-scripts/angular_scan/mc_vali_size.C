{ 

gROOT->Reset(); 
gROOT->SetStyle("Plain"); 
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetTitleFillColor(0);


//
// MC data 
TFile *fmc = new TFile("root-final/tuneA.root");
TTree *tmc = (TTree*) fmc->Get("Hit");

//
// TB data 
TFile *ftb = new TFile("root-final/tbtiltscan.root");
TTree *ttb = (TTree*) ftb->Get("Hit");

TCut TBcut = "hasTrack == 0 && chi2 <100 && ndof == 8 && iRun == 378";
TCut MCcut = "hasTrack == 0 && iRun == 9";

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

Double_t Norm; 

TCanvas * c9 = new TCanvas("c9","c9",600,400);   
c9->SetLeftMargin(0.13);
//c9->SetRightMargin(0.3);
c9->SetTopMargin(0.1);
c9->SetBottomMargin(0.16);


ttb->Draw("sizeCol >>htb(7,0,7)", TBcut );
htb->SetName("htb"); 
htb->SetXTitle("column size [cols]");
htb->SetYTitle("probability [1/cols]"); 
htb->SetTitle(""); 
htb->SetLineColor(kBlack);
htb->SetLineWidth(2);
htb->GetYaxis()->SetTitleOffset(1.2);
Norm = htb->GetEntries()*htb->GetBinWidth(1);
if ( Norm != 0 )   htb->Scale( 1.0/Norm ); 





tmc->Draw("sizeRow >>hmc(7,0,7)", MCcut);
hmc->SetName("hmc"); 
hmc->SetXTitle("column size [cols]");
hmc->SetYTitle("probability [1/cols]"); 
hmc->SetTitle(""); 
hmc->SetLineColor(14);
hmc->SetFillColor(14);
hmc->SetFillStyle(3001); 
hmc->GetYaxis()->SetTitleOffset(1.2);
Norm = hmc->GetEntries()*hmc->GetBinWidth(1);
if ( Norm != 0 )   hmc->Scale( 1.0/Norm ); 


hmc->Draw();
htb->Draw("same");


TLegend* l3 = new TLegend(0.6,0.5,0.85,0.7);
l3->SetFillColor(kWhite); 
l3->SetBorderSize(0);
l3->AddEntry("hmc","Simulation","f");
l3->AddEntry("htb","Test Beam","l");
l3->Draw();


TCanvas * c10 = new TCanvas("c10","c10",600,400);   
c10->SetLeftMargin(0.13);
//c10->SetRightMargin(0.3);
c10->SetTopMargin(0.1);
c10->SetBottomMargin(0.16);


ttb->Draw("sizeRow >>htbr(7,0,7)", TBcut );
htbr->SetName("htbr"); 
htbr->SetXTitle("row size [cols]");
htbr->SetYTitle("probability [1/rows]"); 
htbr->SetTitle(""); 
htbr->SetLineColor(kBlack);
htbr->SetLineWidth(2);
htbr->GetYaxis()->SetTitleOffset(1.2);
Norm = htbr->GetEntries()*htbr->GetBinWidth(1);
if ( Norm != 0 )   htbr->Scale( 1.0/Norm ); 



tmc->Draw("sizeCol >>hmcr(7,0,7)", MCcut);
hmcr->SetName("hmcr"); 
hmcr->SetXTitle("row size [rows]");
hmcr->SetYTitle("probability [1/rows]"); 
hmcr->SetTitle(""); 
hmcr->SetLineColor(14);
hmcr->SetFillColor(14);
hmcr->SetFillStyle(3001); 
hmcr->GetYaxis()->SetTitleOffset(1.2);
Norm = hmcr->GetEntries()*hmcr->GetBinWidth(1);
if ( Norm != 0 )   hmcr->Scale( 1.0/Norm ); 


hmcr->Draw();
htbr->Draw("same");


TLegend* l4 = new TLegend(0.6,0.5,0.85,0.7);
l4->SetFillColor(kWhite); 
l4->SetBorderSize(0);
l4->AddEntry("hmcr","Simulation","f");
l4->AddEntry("htbr","Test Beam","l");
l4->Draw();





}
