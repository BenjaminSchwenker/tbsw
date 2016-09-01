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
int BINS = 50; 
double LSB = 50;
double GAIN = 175;  
int nEvt = 560000 - 25000;

// 
// Analysis cuts 
TCut basecut = "hasTrack == -1 && nTelTracks == 1 && iEvt > 25000 && iEvt < 560000";
TCut roicut  = "col_hit > 1 && col_hit < 30";

//
// Event data 
TFile *ftb = new TFile("Histos-71074.root");
TTree *ttbhit = (TTree*) ftb->Get("Hit");

// ----------------------------------------------------------------------



cout << "Using data from " << nEvt << " triggers" << endl; 


TH1D * hocc = new TH1D("hocc","hocc",BINS,0,BINS*LSB);
TH1D htmp("htmp","htmp", 1,0,64); 

for (int bin = 1; bin<=BINS; ++bin) {
  
  htmp.Reset();

  double seedcut =  hocc->GetXaxis()->GetBinCenter(bin); 
  cout << "seed cut " << seedcut << endl; 

  TCut matched = Form(" %f*seedCharge > %f",GAIN,seedcut); 

  ttbhit->Draw( "row_hit >> htmp" , basecut && roicut && matched );
  htmp.Scale(1.0/(28*64*nEvt) ); 

  double occ = htmp.GetBinContent(1); 
  hocc->SetBinContent(bin,occ);  
  cout << "occupancy " << occ << endl; 
  
}


TCanvas * c01  = new TCanvas("c01","c01",800,400);
c01->SetLeftMargin(0.13);
c01->SetRightMargin(0.1);
c01->SetTopMargin(0.1);
c01->SetBottomMargin(0.16);

hocc->SetTitle("");
hocc->SetXTitle("seed threshold [e]");
hocc->SetYTitle("occupancy");
hocc->GetXaxis()->SetTitleSize(0.06);
hocc->GetYaxis()->SetTitleSize(0.06);
hocc->GetXaxis()->SetLabelSize(0.06);
hocc->GetYaxis()->SetLabelSize(0.06);
//hocc->GetXaxis()->SetTitleOffset(1.3);
//hocc->GetYaxis()->SetTitleOffset(1.3);

hocc->Draw();

TLegend* l1 = new TLegend(0.45,0.2,0.65,0.35);
l1->SetFillColor(kWhite); 
l1->SetBorderSize(0);
l1->AddEntry("hocc","digitizer","p");
l1->Draw();





}
