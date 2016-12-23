{ 

gROOT->Reset(); 
gROOT->SetStyle("Plain"); 

gStyle->SetOptFit();
gStyle->SetPalette(1);
gStyle->SetCanvasColor(0);
gStyle->SetTitleFillColor(0);
gStyle->SetPaintTextFormat("1f");
gStyle->SetOptStat(0);


Int_t BINSU = 256;
Int_t BINSV = 128;

Int_t MAXU = 1152;
Int_t MAXV = 576;



// 
// Analysis cuts 
TCut basecut = "trackChi2 < 100 && trackNHits  == 6 &&  iEvt > 25000 && nTelTracks == 1";
TCut seedcut = "seedCharge > 0";
TCut matched = "hasHit == 0"; 
TCut u_cut = Form("cellU_fit >= 0 && cellU_fit < %d",MAXU);
TCut v_cut = Form("cellV_fit >= 0 && cellV_fit < %d",MAXV); 




//
// event data 
TFile *ftb = new TFile("root-files/DUT-Histos-dutmc.root");
TTree *ttbtrack = (TTree*) ftb->Get("Track");


TCanvas * c1  = new TCanvas("c1","c1",600,400);
c1->SetLeftMargin(0.2);
c1->SetRightMargin(0.1);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.16);


ttbtrack->Draw(Form("cellV_fit >> htotal_row(%d,0,%d)",BINSV,MAXV) , basecut && u_cut );
ttbtrack->Draw(Form("cellV_fit >> hpass_row(%d,0,%d)",BINSV,MAXV) ,  basecut && u_cut && seedcut && matched );



TH1F* pass_row = dynamic_cast<TH1F*> (hpass_row); 
TH1F* total_row = dynamic_cast<TH1F*> (htotal_row); 



TGraphAsymmErrors * effiRow = new TGraphAsymmErrors ( pass_row, total_row, "w" ); 



effiRow->SetTitle("");
effiRow->GetXaxis()->SetTitle("cellV [cellID]");
effiRow->GetYaxis()->SetTitle("efficiency");
effiRow->GetXaxis()->SetRangeUser(0,MAXV);
effiRow->GetYaxis()->SetRangeUser(0.0,1.0);
effiRow->GetXaxis()->SetTitleSize(0.06);
effiRow->GetYaxis()->SetTitleSize(0.06);
effiRow->GetXaxis()->SetLabelSize(0.06);
effiRow->GetYaxis()->SetLabelSize(0.06);
//effiRow->GetXaxis()->SetTitleOffset(1.3);
effiRow->GetYaxis()->SetTitleOffset(1.7);

pass_row->Draw();
effiRow->Draw("a*"); 



TCanvas * c2  = new TCanvas("c2","c2",600,400);
c2->SetLeftMargin(0.2);
c2->SetRightMargin(0.1);
c2->SetTopMargin(0.1);
c2->SetBottomMargin(0.16);

ttbtrack->Draw(Form("cellU_fit  >> htotal_col(%d,0,%d)",BINSU,MAXU) , basecut && v_cut );
ttbtrack->Draw(Form("cellU_fit  >> hpass_col(%d,0,%d)",BINSU,MAXU) ,  basecut && v_cut && seedcut && matched );




TH1F* pass_col = dynamic_cast<TH1F*> (hpass_col); 
TH1F* total_col = dynamic_cast<TH1F*> (htotal_col); 




TGraphAsymmErrors * effiCol = new TGraphAsymmErrors ( pass_col, total_col, "w" ); 

effiCol->SetTitle("");
effiCol->GetXaxis()->SetTitle("cellU [cellID]");
effiCol->GetYaxis()->SetTitle("efficiency");
effiCol->GetYaxis()->SetRangeUser(0.98,1.0);
effiCol->GetXaxis()->SetRangeUser(0,MAXU);
effiCol->GetXaxis()->SetTitleSize(0.06);
effiCol->GetYaxis()->SetTitleSize(0.06);
effiCol->GetXaxis()->SetLabelSize(0.06);
effiCol->GetYaxis()->SetLabelSize(0.06);
//effiCol->GetXaxis()->SetTitleOffset(1.3);
effiCol->GetYaxis()->SetTitleOffset(1.7);


hpass_col->Draw();
effiCol->Draw("a*"); 


TCanvas * c3  = new TCanvas("c3","c3",600,600);
c3->SetLeftMargin(0.14);
c3->SetRightMargin(0.26);
c3->SetTopMargin(0.1);
c3->SetBottomMargin(0.16);

ttbtrack->Draw(Form(" cellV_fit:cellU_fit >> htotal(%d,0,%d,%d,0,%d)", BINSU,MAXU,BINSV,MAXV) , basecut );
ttbtrack->Draw(Form(" cellV_fit:cellU_fit >> hpass(%d,0,%d,%d,0,%d)",BINSU,MAXU,BINSV,MAXV) ,  basecut && seedcut && matched );

TH2D * effi = hpass->Clone(); 
effi->Divide(htotal); 


effi->SetTitle("");
effi->GetXaxis()->SetTitle("cellU [cellID]");
effi->GetYaxis()->SetTitle("cellV [cellID]");
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
