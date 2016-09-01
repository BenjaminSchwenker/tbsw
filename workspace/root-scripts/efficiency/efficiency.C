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
TCut basecut = "1x1Quality == 0 && chi2 < 100 && iEvt > 25000 && nTelTracks == 1 && nhits>0";
TCut seedcut = "seedCharge > 0";
TCut matched = "hasHit == 0"; 
TCut col_cut = "col_fit >= 0 && col_fit < 32";
TCut row_cut = "row_fit >= 0 && row_fit < 63"; 

//
// event data 
TFile *ftb = new TFile("root-final/Histos-run000331.root");
TTree *ttbtrack = (TTree*) ftb->Get("Track");


TCanvas * c1  = new TCanvas("c1","c1",600,400);
c1->SetLeftMargin(0.2);
c1->SetRightMargin(0.1);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.16);


ttbtrack->Draw(" row_fit >> htotal_row(64,0,64)" , basecut && col_cut );
ttbtrack->Draw(" row_fit >> hpass_row(64,0,64)" ,  basecut && col_cut && seedcut && matched );

TGraphAsymmErrors * effiRow = new TGraphAsymmErrors ( hpass_row, htotal_row, "w" ); 

effiRow->SetTitle("");
effiRow->GetXaxis()->SetTitle("rows");
effiRow->GetYaxis()->SetTitle("efficiency");
effiRow->GetXaxis()->SetRangeUser(0,64);
effiRow->GetYaxis()->SetRangeUser(0.98,1.0);
effiRow->GetXaxis()->SetTitleSize(0.06);
effiRow->GetYaxis()->SetTitleSize(0.06);
effiRow->GetXaxis()->SetLabelSize(0.06);
effiRow->GetYaxis()->SetLabelSize(0.06);
//effiRow->GetXaxis()->SetTitleOffset(1.3);
effiRow->GetYaxis()->SetTitleOffset(1.7);

hpass_row->Draw();
effiRow->Draw("a*"); 



TLegend* l1 = new TLegend(0.45,0.2,0.85,0.45);
l1->SetFillColor(kWhite); 
l1->SetBorderSize(0);
l1->AddEntry("effiRow","H4.1.15 CERN 2012","pe");
l1->Draw();


TCanvas * c2  = new TCanvas("c2","c2",600,400);
c2->SetLeftMargin(0.2);
c2->SetRightMargin(0.1);
c2->SetTopMargin(0.1);
c2->SetBottomMargin(0.16);

ttbtrack->Draw(" col_fit >> htotal_col(32,0,32)" , basecut && row_cut );
ttbtrack->Draw(" col_fit >> hpass_col(32,0,32)" ,  basecut && row_cut && seedcut && matched );

TGraphAsymmErrors * effiCol = new TGraphAsymmErrors ( hpass_col, htotal_col, "w" ); 

effiCol->SetTitle("");
effiCol->GetXaxis()->SetTitle("columns");
effiCol->GetYaxis()->SetTitle("efficiency");
effiCol->GetYaxis()->SetRangeUser(0.98,1.0);
effiCol->GetXaxis()->SetRangeUser(0,32);
effiCol->GetXaxis()->SetTitleSize(0.06);
effiCol->GetYaxis()->SetTitleSize(0.06);
effiCol->GetXaxis()->SetLabelSize(0.06);
effiCol->GetYaxis()->SetLabelSize(0.06);
//effiCol->GetXaxis()->SetTitleOffset(1.3);
effiCol->GetYaxis()->SetTitleOffset(1.7);


hpass_col->Draw();
effiCol->Draw("a*"); 


TLegend* l2 = new TLegend(0.45,0.2,0.85,0.45);
l2->SetFillColor(kWhite); 
l2->SetBorderSize(0);
l2->AddEntry("effiCol","H4.1.15 CERN 2012","pe");
l2->Draw();


TCanvas * c3  = new TCanvas("c3","c3",600,600);
c3->SetLeftMargin(0.14);
c3->SetRightMargin(0.26);
c3->SetTopMargin(0.1);
c3->SetBottomMargin(0.16);

ttbtrack->Draw(" row_fit:col_fit >> htotal(64,0,64,64,0,64)" , basecut );
ttbtrack->Draw(" row_fit:col_fit >> hpass(64,0,64,64,0,64)" ,  basecut && seedcut && matched );

TH2D * effi = hpass->Clone(); 
effi->Divide(htotal); 


effi->SetTitle("");
effi->GetXaxis()->SetTitle("columns");
effi->GetYaxis()->SetTitle("rows");
effi->GetZaxis()->SetTitle("efficiency");
effi->GetZaxis()->SetRangeUser(0.0,1.0);
effi->GetXaxis()->SetTitleSize(0.06);
effi->GetYaxis()->SetTitleSize(0.06);
effi->GetZaxis()->SetTitleSize(0.06);
effi->GetXaxis()->SetLabelSize(0.06);
effi->GetYaxis()->SetLabelSize(0.06);
effi->GetZaxis()->SetLabelSize(0.06);
effi->GetYaxis()->SetTitleOffset(1.15);
effi->GetZaxis()->SetTitleOffset(1.7);
effi->Draw("colz");




}
