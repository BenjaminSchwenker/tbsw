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
// Analysis cuts 
TCut basecut = "hasTrack == -1";
TCut seedcut = "seedCharge > 0";
TCut col_cut = "";
TCut row_cut = "";

double thr = 5; 
double tocc = 0.5 + 0.5*TMath::Erf(-thr/TMath::Sqrt(2));

cout << "Cut at " << thr << "xNoise gives noise occupancy " << tocc << endl;  
  

//
// event data 
//TFile *ftb = new TFile("root-final/tbtiltscan.root");
TFile *ftb = new TFile("root-final/tuneA.root");
TTree *ttb = (TTree*) ftb->Get("Hit");



TCanvas * c1  = new TCanvas("c1","c1",600,400);
c1->SetLeftMargin(0.2);
c1->SetRightMargin(0.1);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.16);

ttb->Draw(" row_hit >> hocc_row(64,0,64)" , basecut && col_cut );
hocc_row->Scale(1.0/(32) ); 

hocc_row->SetTitle("");
hocc_row->GetXaxis()->SetTitle("rows");
hocc_row->GetYaxis()->SetTitle("occupancy");
hocc_row->GetXaxis()->SetRangeUser(0,64);
//hocc_row->GetYaxis()->SetRangeUser(0.95,1.0);
hocc_row->GetXaxis()->SetTitleSize(0.06);
hocc_row->GetYaxis()->SetTitleSize(0.06);
hocc_row->GetXaxis()->SetLabelSize(0.06);
hocc_row->GetYaxis()->SetLabelSize(0.06);
//hocc_row->GetXaxis()->SetTitleOffset(1.3);
hocc_row->GetYaxis()->SetTitleOffset(1.7);
hocc_row->Draw();
 


TCanvas * c2  = new TCanvas("c2","c2",600,400);
c2->SetLeftMargin(0.2);
c2->SetRightMargin(0.1);
c2->SetTopMargin(0.1);
c2->SetBottomMargin(0.16);

ttb->Draw(" col_hit >> hocc_col(64,0,64)" , basecut && row_cut );
hocc_col->Scale(1.0/64 ); 

hocc_col->SetTitle("");
hocc_col->GetXaxis()->SetTitle("columns");
hocc_col->GetYaxis()->SetTitle("occupancy");
hocc_col->GetXaxis()->SetRangeUser(0,64);
//hocc_col->GetYaxis()->SetRangeUser(0.95,1.0);
hocc_col->GetXaxis()->SetTitleSize(0.06);
hocc_col->GetYaxis()->SetTitleSize(0.06);
hocc_col->GetXaxis()->SetLabelSize(0.06);
hocc_col->GetYaxis()->SetLabelSize(0.06);
//hocc_col->GetXaxis()->SetTitleOffset(1.3);
hocc_col->GetYaxis()->SetTitleOffset(1.7);
hocc_col->Draw();
 




TCanvas * c3  = new TCanvas("c3","c3",600,400);
c3->SetLeftMargin(0.14);
c3->SetRightMargin(0.26);
c3->SetTopMargin(0.1);
c3->SetBottomMargin(0.16);

ttb->Draw(" row_hit:col_hit >> hocc_2d(64,0,64,64,0,64)" , basecut );
//hocc_2d->Scale(1.0/(1.0*nEvt) ); 


hocc_2d->SetTitle("");
hocc_2d->GetXaxis()->SetTitle("columns");
hocc_2d->GetYaxis()->SetTitle("rows");
hocc_2d->GetZaxis()->SetTitle("occupancy");
//hocc_2d->GetZaxis()->SetRangeUser(0.9,1.0);
hocc_2d->GetXaxis()->SetTitleSize(0.06);
hocc_2d->GetYaxis()->SetTitleSize(0.06);
hocc_2d->GetZaxis()->SetTitleSize(0.06);
hocc_2d->GetXaxis()->SetLabelSize(0.06);
hocc_2d->GetYaxis()->SetLabelSize(0.06);
hocc_2d->GetZaxis()->SetLabelSize(0.06);
hocc_2d->GetYaxis()->SetTitleOffset(1.15);
hocc_2d->GetZaxis()->SetTitleOffset(1.7);
hocc_2d->Draw("colz");




}
