{ 

gROOT->Reset(); 
gROOT->SetStyle("Plain"); 

gStyle->SetOptFit();
gStyle->SetPalette(1);
gStyle->SetCanvasColor(0);
gStyle->SetTitleFillColor(0);
gStyle->SetPaintTextFormat("1f");
gStyle->SetOptStat(0);




Int_t NBINSU = 1152;
Int_t NBINSV = 576;
Int_t First = 25000; 
Int_t Last = 100000;



// 
// Analysis cuts 
TCut basecut = Form("hasTrack == -1 && iEvt > %d && iEvt < %d ",First, Last);

Int_t EVENTS = Last - First; 

//
// event data 
TFile *ftb = new TFile("root-files/DUT-Histos-dutmc.root");
TTree *ttb = (TTree*) ftb->Get("Hit");


TCanvas * c1  = new TCanvas("c1","c1",600,400);
c1->SetLeftMargin(0.2);
c1->SetRightMargin(0.1);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.16);

ttb->Draw(Form("cellV_hit >> hocc_row(%d,0,%d)",NBINSV,NBINSV) , basecut );
hocc_row->Scale(1.0/(NBINSU*EVENTS ) ); 

hocc_row->SetTitle("");
hocc_row->GetXaxis()->SetTitle("cellV [cellID]");
hocc_row->GetYaxis()->SetTitle("occupancy");
hocc_row->GetXaxis()->SetRangeUser(0,NBINSV);
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

ttb->Draw(Form(" cellU_hit >> hocc_col(%d,0,%d)",NBINSU,NBINSU) , basecut  );
hocc_col->Scale(1.0/(NBINSV*EVENTS)); 

hocc_col->SetTitle("");
hocc_col->GetXaxis()->SetTitle("cellU [cellID]");
hocc_col->GetYaxis()->SetTitle("occupancy");
hocc_col->GetXaxis()->SetRangeUser(0,NBINSU);
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

ttb->Draw(Form(" cellV_hit:cellU_hit >> hocc_2d(%d,0,%d,%d,0,%d)",NBINSU,NBINSU,NBINSV,NBINSV) , basecut );
hocc_2d->Scale(1.0/(1.0*EVENTS) ); 


hocc_2d->SetTitle("");
hocc_2d->GetXaxis()->SetTitle("cellU [cellID]");
hocc_2d->GetYaxis()->SetTitle("cellV [cellID]");
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
