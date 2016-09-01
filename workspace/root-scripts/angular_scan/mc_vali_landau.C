{ 

gROOT->Reset(); 
gROOT->SetStyle("Plain"); 
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetTitleFillColor(0);
gStyle->SetPaintTextFormat("1f");

//
// MC data 
TFile *fmc = new TFile("root-final/tune_B.root");
TTree *tmc = (TTree*) fmc->Get("Hit");

//
// TB data 
TFile *ftb = new TFile("root-final/tbtiltscan.root");
TTree *ttb = (TTree*) ftb->Get("Hit");

TCut TBcut = "hasTrack == 0 && chi2 <100 && ndof == 8 && iRun == 378 && clusterCharge > 30";
TCut MCcut = "hasTrack == 0 && iRun == 9 && clusterCharge > 30";

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


ttb->Draw("clusterCharge >>htb(200,0,200)", TBcut);
htb->SetName("htb"); 
htb->SetXTitle("cluster signal [ADU]");
htb->SetYTitle("probability [1/ADU]"); 
htb->SetTitle(""); 
htb->SetLineColor(kBlack);
htb->SetLineWidth(2);

htb->Sumw2();
Norm = htb->GetEntries()*htb->GetBinWidth(1);
if ( Norm != 0 )   htb->Scale( 1.0/Norm ); 


htb->GetXaxis()->SetTitleSize(0.06);
htb->GetYaxis()->SetTitleSize(0.06);
htb->GetXaxis()->SetLabelSize(0.06);
htb->GetYaxis()->SetLabelSize(0.06);
//htb->GetYaxis()->SetTitleOffset(1.2);


tmc->Draw("clusterCharge >>hmc(200,0,200)", MCcut);
hmc->SetName("hmc"); 
hmc->SetXTitle("cluster signal [ADU]");
hmc->SetYTitle("probability [1/ADU]"); 
hmc->SetTitle(""); 
hmc->SetOption("hist");
hmc->SetFillColor(kGray);

hmc->Sumw2();
Norm = hmc->GetEntries()*hmc->GetBinWidth(1);
if ( Norm != 0 )   hmc->Scale( 1.0/Norm ); 


hmc->GetXaxis()->SetTitleSize(0.06);
hmc->GetYaxis()->SetTitleSize(0.06);
hmc->GetXaxis()->SetLabelSize(0.06);
hmc->GetYaxis()->SetLabelSize(0.06);
//hmc->GetYaxis()->SetTitleOffset(1.2);


hmc->Draw();
htb->Draw("e same");



int bin1, bin2 , bin; 
double fwhm, max; 

max= hmc->GetMaximum();
bin1 = hmc->FindFirstBinAbove(max/2);
bin2 = hmc->FindLastBinAbove(max/2);
fwhm = hmc->GetBinCenter(bin2) - hmc->GetBinCenter(bin1);

bin = hmc->GetMaximumBin(); 
cout << "width " << hmc->GetBinWidth(1) << endl; 
cout << "max bin " << hmc->GetBinCenter(bin)   << endl; 
cout << "MC FWHM is " << fwhm << endl; 

max= htb->GetMaximum();
bin1 = htb->FindFirstBinAbove(max/2);
bin2 = htb->FindLastBinAbove(max/2);
fwhm = htb->GetBinCenter(bin2) - htb->GetBinCenter(bin1);

cout << "TB FWHM is " << fwhm << endl; 



TLegend* l3 = new TLegend(0.6,0.6,0.88,0.8);
l3->SetFillColor(kWhite); 
l3->SetBorderSize(0);
l3->AddEntry("hmc","MC 4GeV e^{-}","f");
l3->AddEntry("htb","TB 4GeV e^{-}","lep");
l3->Draw();


//--------- seed charge 

TCanvas * c10 = new TCanvas("c10","c10",600,400);   
c10->SetLeftMargin(0.13);
//c10->SetRightMargin(0.3);
c10->SetTopMargin(0.1);
c10->SetBottomMargin(0.16);


ttb->Draw("seedCharge >>htbs(200,0,200)", TBcut);
htbs->SetName("htbs"); 
htbs->SetXTitle("seed signal [ADU]");
htbs->SetYTitle("probability [1/ADU]"); 
htbs->SetTitle(""); 
htbs->SetLineColor(kBlack);
htbs->SetLineWidth(2);

htbs->Sumw2();
Norm = htbs->GetEntries()*htbs->GetBinWidth(1);
if ( Norm != 0 )   htbs->Scale( 1.0/Norm ); 


htbs->GetXaxis()->SetTitleSize(0.06);
htbs->GetYaxis()->SetTitleSize(0.06);
htbs->GetXaxis()->SetLabelSize(0.06);
htbs->GetYaxis()->SetLabelSize(0.06);
//htbs->GetYaxis()->SetTitleOffset(1.2);


tmc->Draw("seedCharge >>hmcs(200,0,200)", MCcut);
hmcs->SetName("hmcs"); 
hmcs->SetXTitle("seed signal [ADU]");
hmcs->SetYTitle("probability [1/ADU]"); 
hmcs->SetTitle(""); 
hmcs->SetOption("hist");
hmcs->SetFillColor(kGray);

hmcs->Sumw2();
Norm = hmcs->GetEntries()*hmcs->GetBinWidth(1);
if ( Norm != 0 )   hmcs->Scale( 1.0/Norm ); 


hmcs->GetXaxis()->SetTitleSize(0.06);
hmcs->GetYaxis()->SetTitleSize(0.06);
hmcs->GetXaxis()->SetLabelSize(0.06);
hmcs->GetYaxis()->SetLabelSize(0.06);
//hmcs->GetYaxis()->SetTitleOffset(1.2);


hmcs->Draw();
htbs->Draw("e same");




TLegend* l4 = new TLegend(0.6,0.6,0.88,0.8);
l4->SetFillColor(kWhite); 
l4->SetBorderSize(0);
l4->AddEntry("hmcs","MC 4GeV e^{-}","f");
l4->AddEntry("htbs","TB 4GeV e^{-}","lep");
l4->Draw();




}
