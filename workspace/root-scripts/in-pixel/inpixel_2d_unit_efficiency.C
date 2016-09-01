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
TFile *ftb = new TFile("Histos.root");
TTree *ttbtrack = (TTree*) ftb->Get("Track");


// 
// Analysis cuts 
TCut basecut = "1x1Quality == 0 && chi2 < 100 && iEvt > 25000 && nTelTracks == 1";
TCut seedcut = "seedCharge > 10";
TCut matched = "hasHit == 0"; 
TCut col_cut = "col_fit >= 0 && col_fit < 32";
TCut row_cut = "row_fit >= 0 && row_fit < 63"; 



// ----------------------------------------------------------------------


TCanvas * c3 = new TCanvas("c3","c3",400,400); 

c3->SetLeftMargin(0.22);
c3->SetRightMargin(0.3);
c3->SetTopMargin(0.1);
c3->SetBottomMargin(0.16);

ttbtrack->Draw(" (v_fit-v_pixel)*1000+37.5+(row_fit%2)*75:(u_fit-u_pixel)*1000+25+(col_fit%2)*50 >> htotal(10,0,100,20,0,150)" , basecut && col_cut && row_cut );
ttbtrack->Draw(" (v_fit-v_pixel)*1000+37.5+(row_fit%2)*75:(u_fit-u_pixel)*1000+25+(col_fit%2)*50 >> hpass(10, 0,100,20,0,150)" , basecut && col_cut && row_cut && seedcut && matched );



TH2D * effi = hpass->Clone(); 
effi->Divide(htotal); 

effi->GetXaxis()->SetNdivisions(2,kFALSE);
effi->GetYaxis()->SetNdivisions(2,kFALSE);
effi->SetTitle("");
effi->SetXTitle("u [#mum]");
effi->SetYTitle("v [#mum]");
effi->SetZTitle("efficiency");
effi->GetXaxis()->SetTitleSize(0.06);
effi->GetYaxis()->SetTitleSize(0.06);
effi->GetZaxis()->SetTitleSize(0.06);
effi->GetXaxis()->SetLabelSize(0.06);
effi->GetYaxis()->SetLabelSize(0.06);
effi->GetZaxis()->SetLabelSize(0.06);
effi->GetZaxis()->SetTitleOffset(1.7);
effi->GetYaxis()->SetTitleOffset(1.3);
effi->GetZaxis()->SetRangeUser(0.97,1.0);
effi->Draw("colz"); 






}
