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
TCut hitcut = "chi2 < 100  && hasTrack == 0 && clusterQuality == 0";
TCut trkcut = "chi2 < 100 && clusterQuality == 0";

//
// event data 
TFile *ftb = new TFile("Histos.root");
TTree *ttb = (TTree*) ftb->Get("Hit");


// ----------------------------------------------------------------------

Double_t Norm; 

TCanvas * c0  = new TCanvas("c0","c0",600,400);
c0->SetLeftMargin(0.2);
c0->SetRightMargin(0.2);
c0->SetTopMargin(0.1);
c0->SetBottomMargin(0.16);

ttb->Draw("row_hit:col_hit>>hspot(256,0,256,768,0,768)", hitcut);
hspot->SetName("hspot"); 
hspot->SetTitle("");
hspot->SetStats( false );  
hspot->GetXaxis()->SetTitle("columns");
hspot->GetYaxis()->SetTitle("rows");
hspot->GetZaxis()->SetTitle("number of matched hits");
hspot->GetXaxis()->SetTitleSize(0.08);
hspot->GetYaxis()->SetTitleSize(0.08);
hspot->GetZaxis()->SetTitleSize(0.08);
hspot->GetXaxis()->SetLabelSize(0.08);
hspot->GetYaxis()->SetLabelSize(0.08);
hspot->GetZaxis()->SetLabelSize(0.08);
hspot->Draw("colz");

TCanvas * c1  = new TCanvas("c1","c1",600,400);
c1->SetLeftMargin(0.2);
c1->SetRightMargin(0.2);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.16);

ttb->Draw("row_hit:col_hit>>hspot2(256,0,256,768,0,768)", trkcut);
hspot2->SetName("hspot2"); 
hspot2->SetTitle("");
hspot2->SetStats( false );  
hspot2->GetXaxis()->SetTitle("columns");
hspot2->GetYaxis()->SetTitle("rows");
hspot2->GetZaxis()->SetTitle("number of hits");
hspot2->GetXaxis()->SetTitleSize(0.08);
hspot2->GetYaxis()->SetTitleSize(0.08);
hspot2->GetZaxis()->SetTitleSize(0.08);
hspot2->GetXaxis()->SetLabelSize(0.08);
hspot2->GetYaxis()->SetLabelSize(0.08);
hspot2->GetZaxis()->SetLabelSize(0.08);
hspot2->Draw("colz");


TCanvas * c3  = new TCanvas("c3","c3",600,400);
c3->SetLeftMargin(0.2);
c3->SetRightMargin(0.1);
c3->SetTopMargin(0.1);
c3->SetBottomMargin(0.16);

ttb->Draw("(u_hit - u_fit)*1000 >>hu(201,-475.,+475.)", hitcut);
hu->SetName("hu"); 
hu->SetTitle("");
hu->GetXaxis()->SetTitle("residuals u [#mum]");
hu->GetYaxis()->SetTitle("number of hits");
hu->GetXaxis()->SetTitleSize(0.06);
hu->GetYaxis()->SetTitleSize(0.06);
hu->GetYaxis()->SetTitleOffset(1.5);
hu->GetXaxis()->SetLabelSize(0.06);
hu->GetYaxis()->SetLabelSize(0.06);
hu->Draw("");

TCanvas * c31  = new TCanvas("c31","c31",600,400);
c31->SetLeftMargin(0.2);
c31->SetRightMargin(0.1);
c31->SetTopMargin(0.1);
c31->SetBottomMargin(0.16);

ttb->Draw("(u_hit - u_fit)*1000:col_hit >>hucol(250,0,250,40,-550.,+550.)", hitcut);
hucol->SetName("hucol"); 
hucol->SetTitle("");
hucol->GetYaxis()->SetTitle("residuals u [#mum]");
hucol->GetXaxis()->SetTitle("columns");
hucol->GetZaxis()->SetTitle("number of hits");
hucol->GetXaxis()->SetTitleSize(0.06);
hucol->GetYaxis()->SetTitleSize(0.06);
hucol->GetYaxis()->SetTitleOffset(1.5);
hucol->GetXaxis()->SetLabelSize(0.06);
hucol->GetYaxis()->SetLabelSize(0.06);
hucol->SetStats( false );  
hucol->Draw("colz");

TCanvas * c32  = new TCanvas("c32","c32",600,400);
c32->SetLeftMargin(0.2);
c32->SetRightMargin(0.1);
c32->SetTopMargin(0.1);
c32->SetBottomMargin(0.16);

ttb->Draw("(u_hit - u_fit)*1000:row_hit >>hurow(768,0,768,40,-550.,+550.)", hitcut);
hurow->SetName("hurow"); 
hurow->SetTitle("");
hurow->GetYaxis()->SetTitle("residuals u [#mum]");
hurow->GetXaxis()->SetTitle("rows");
hurow->GetZaxis()->SetTitle("number of hits");
hurow->GetXaxis()->SetTitleSize(0.06);
hurow->GetYaxis()->SetTitleSize(0.06);
hurow->GetYaxis()->SetTitleOffset(1.5);
hurow->GetXaxis()->SetLabelSize(0.06);
hurow->GetYaxis()->SetLabelSize(0.06);
hurow->SetStats( false );  
hurow->Draw("colz");



TCanvas * c4  = new TCanvas("c4","c4",600,400);
c4->SetLeftMargin(0.2);
c4->SetRightMargin(0.1);
c4->SetTopMargin(0.1);
c4->SetBottomMargin(0.16);

ttb->Draw("(v_hit - v_fit)*1000 >>hv(201,-475.,+475.)", hitcut);
hv->SetName("hv"); 
hv->SetTitle("");
hv->GetXaxis()->SetTitle("residuals v [#mum]");
hv->GetYaxis()->SetTitle("number of hits");
hv->GetXaxis()->SetTitleSize(0.06);
hv->GetYaxis()->SetTitleSize(0.06);
hv->GetYaxis()->SetTitleOffset(1.5);
hv->GetXaxis()->SetLabelSize(0.06);
hv->GetYaxis()->SetLabelSize(0.06);
hv->Draw("");

TCanvas * c41  = new TCanvas("c41","c41",600,400);
c41->SetLeftMargin(0.2);
c41->SetRightMargin(0.1);
c41->SetTopMargin(0.1);
c41->SetBottomMargin(0.16);

ttb->Draw("(v_hit - v_fit)*1000:col_hit >>hvcol(250,0,250,40,-550.,+550.)", hitcut);
hvcol->SetName("hvcol"); 
hvcol->SetTitle("");
hvcol->GetYaxis()->SetTitle("residuals v [#mum]");
hvcol->GetXaxis()->SetTitle("columns");
hvcol->GetZaxis()->SetTitle("number of hits");
hvcol->GetXaxis()->SetTitleSize(0.06);
hvcol->GetYaxis()->SetTitleSize(0.06);
hvcol->GetYaxis()->SetTitleOffset(1.5);
hvcol->GetXaxis()->SetLabelSize(0.06);
hvcol->GetYaxis()->SetLabelSize(0.06);
hvcol->SetStats( false );  
hvcol->Draw("colz");

TCanvas * c42  = new TCanvas("c42","c42",600,400);
c42->SetLeftMargin(0.2);
c42->SetRightMargin(0.1);
c42->SetTopMargin(0.1);
c42->SetBottomMargin(0.16);

ttb->Draw("(v_hit - v_fit)*1000:row_hit >>hvrow(768,0,768,40,-550.,+550.)", hitcut);
hvrow->SetName("hvrow"); 
hvrow->SetTitle("");
hvrow->GetYaxis()->SetTitle("residuals v [#mum]");
hvrow->GetXaxis()->SetTitle("rows");
hvrow->GetZaxis()->SetTitle("number of hits");
hvrow->GetXaxis()->SetTitleSize(0.06);
hvrow->GetYaxis()->SetTitleSize(0.06);
hvrow->GetYaxis()->SetTitleOffset(1.5);
hvrow->GetXaxis()->SetLabelSize(0.06);
hvrow->GetYaxis()->SetLabelSize(0.06);
hvrow->SetStats( false );  
hvrow->Draw("colz");


}
