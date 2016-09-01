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
TCut roicut = "";
TCut trkcut = "";

//
// event data 
TFile *ftb = new TFile("Histos.root");

TTree *ttb = (TTree*) ftb->Get("Hit");
TTree *ttbtrack = (TTree*) ftb->Get("Track");


// ----------------------------------------------------------------------

Double_t Norm; 

TCanvas * c0  = new TCanvas("c0","c0",600,400);
c0->SetLeftMargin(0.2);
c0->SetRightMargin(0.1);
c0->SetTopMargin(0.1);
c0->SetBottomMargin(0.16);

ttb->Draw("row_hit:col_hit>>hspot", hitcut);
hspot->SetName("hspot"); 
hspot->GetXaxis()->SetTitle("columns");
hspot->GetYaxis()->SetTitle("rows");
hspot->GetZaxis()->SetTitle("number of hits");
hspot->GetXaxis()->SetTitleSize(0.08);
hspot->GetYaxis()->SetTitleSize(0.08);
hspot->GetZaxis()->SetTitleSize(0.08);
hspot->GetXaxis()->SetLabelSize(0.08);
hspot->GetYaxis()->SetLabelSize(0.08);
hspot->GetZaxis()->SetLabelSize(0.08);
hspot->Draw("colz");

TCanvas * c1  = new TCanvas("c1","c1",600,400);
c1->SetLeftMargin(0.2);
c1->SetRightMargin(0.1);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.16);

ttbtrack->Draw("row_fit:col_fit>>htrackspot", trkcut && roicut);
htrackspot->SetName("htrackspot"); 
htrackspot->GetXaxis()->SetTitle("columns");
htrackspot->GetYaxis()->SetTitle("rows");
htrackspot->GetZaxis()->SetTitle("number of tracks");
htrackspot->GetXaxis()->SetTitleSize(0.08);
htrackspot->GetYaxis()->SetTitleSize(0.08);
htrackspot->GetZaxis()->SetTitleSize(0.08);
htrackspot->GetXaxis()->SetLabelSize(0.08);
htrackspot->GetYaxis()->SetLabelSize(0.08);
htrackspot->GetZaxis()->SetLabelSize(0.08);
htrackspot->Draw("colz");



TCanvas * c2  = new TCanvas("c2","c2",600,400);
c2->SetLeftMargin(0.2);
c2->SetRightMargin(0.1);
c2->SetTopMargin(0.1);
c2->SetBottomMargin(0.16);

ttbtrack->Draw("v_fit:u_fit>>hbeamspot", trkcut);
hbeamspot->SetName("hbeamspot"); 
hbeamspot->GetXaxis()->SetTitle("u [mm]");
hbeamspot->GetYaxis()->SetTitle("v [mm]");
hbeamspot->GetZaxis()->SetTitle("number of tracks");
hbeamspot->GetXaxis()->SetTitleSize(0.08);
hbeamspot->GetYaxis()->SetTitleSize(0.08);
hbeamspot->GetZaxis()->SetTitleSize(0.08);
hbeamspot->GetXaxis()->SetLabelSize(0.08);
hbeamspot->GetYaxis()->SetLabelSize(0.08);
hbeamspot->GetZaxis()->SetLabelSize(0.08);
 
hbeamspot->Draw("colz");



TCanvas * c3  = new TCanvas("c3","c3",600,400);
c3->SetLeftMargin(0.2);
c3->SetRightMargin(0.1);
c3->SetTopMargin(0.1);
c3->SetBottomMargin(0.16);

ttb->Draw("(u_hit - u_fit)*1000 >>hu(101,-75.,+75.)", hitcut);
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



TCanvas * c4  = new TCanvas("c4","c4",600,400);
c4->SetLeftMargin(0.2);
c4->SetRightMargin(0.1);
c4->SetTopMargin(0.1);
c4->SetBottomMargin(0.16);

ttb->Draw("(v_hit - v_fit)*1000 >>hv(81,-75.,+75.)", hitcut);
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

TCanvas * c5  = new TCanvas("c5","c5",600,400);
c5->SetLeftMargin(0.2);
c5->SetRightMargin(0.1);
c5->SetTopMargin(0.1);
c5->SetBottomMargin(0.16);

ttb->Draw("chi2 / ndof >> hchisqundof(100,0,10)", hitcut);
hchisqundof->SetName("hchisqundof"); 
hchisqundof->SetTitle("");
hchisqundof->GetXaxis()->SetTitle("#chi^2/ndof");
hchisqundof->GetYaxis()->SetTitle("number of tracks");
hchisqundof->GetXaxis()->SetTitleSize(0.06);
hchisqundof->GetYaxis()->SetTitleSize(0.06);
hchisqundof->GetYaxis()->SetTitleOffset(1.5);
hchisqundof->GetXaxis()->SetLabelSize(0.06);
hchisqundof->GetYaxis()->SetLabelSize(0.06);
hchisqundof->Draw();


TCanvas * c6  = new TCanvas("c6","c6",600,400);
c6->SetLeftMargin(0.2);
c6->SetRightMargin(0.1);
c6->SetTopMargin(0.1);
c6->SetBottomMargin(0.16);

ttb->Draw("chi2pred >>hchi2pred(100,0,20)", hitcut);
hchi2pred->SetName("hchi2pred"); 
hchi2pred->SetTitle("");
hchi2pred->GetXaxis()->SetTitle("local #chi^2");
hchi2pred->GetYaxis()->SetTitle("number of tracks");
hchi2pred->GetXaxis()->SetTitleSize(0.06);
hchi2pred->GetYaxis()->SetTitleSize(0.06);
hchi2pred->GetYaxis()->SetTitleOffset(1.5);
hchi2pred->GetXaxis()->SetLabelSize(0.06);
hchi2pred->GetYaxis()->SetLabelSize(0.06);
hchi2pred->Draw();

TCanvas * c7  = new TCanvas("c7","c7",600,400);
c7->SetLeftMargin(0.2);
c7->SetRightMargin(0.1);
c7->SetTopMargin(0.1);
c7->SetBottomMargin(0.16);


ttb->Draw("seedCharge >>hCharge(101,0,100)", hitcut);
hCharge->SetName("hCharge"); 
hCharge->SetTitle("");
hCharge->GetXaxis()->SetTitle("seed charge [ADU]");
hCharge->GetYaxis()->SetTitle("number of hits");
hCharge->GetXaxis()->SetTitleSize(0.06);
hCharge->GetYaxis()->SetTitleSize(0.06);
hCharge->GetYaxis()->SetTitleOffset(1.5);
hCharge->GetXaxis()->SetLabelSize(0.06);
hCharge->GetYaxis()->SetLabelSize(0.06);
hCharge->Draw();


TCanvas * c8  = new TCanvas("c8","c8",600,400);
c8->SetLeftMargin(0.2);
c8->SetRightMargin(0.1);
c8->SetTopMargin(0.1);
c8->SetBottomMargin(0.16);


ttb->Draw("clusterCharge >>hCCharge(101,0,100)", hitcut);
hCCharge->SetName("hCCharge"); 
hCCharge->SetTitle("");
hCCharge->GetXaxis()->SetTitle("cluster charge [ADU]");
hCCharge->GetYaxis()->SetTitle("number of hits");
hCCharge->GetXaxis()->SetTitleSize(0.06);
hCCharge->GetYaxis()->SetTitleSize(0.06);
hCCharge->GetYaxis()->SetTitleOffset(1.5);
hCCharge->GetXaxis()->SetLabelSize(0.06);
hCCharge->GetYaxis()->SetLabelSize(0.06);
hCCharge->Draw();


TCanvas * c10  = new TCanvas("c10","c10",600,400);
c10->SetLeftMargin(0.2);
c10->SetRightMargin(0.1);
c10->SetTopMargin(0.1);
c10->SetBottomMargin(0.16);

ttb->Draw("size >>hsize(20,0,20)", hitcut);
hsize->SetName("hsize"); 
hsize->SetTitle("");
hsize->GetXaxis()->SetTitle("cluster size [pixels]");
hsize->GetYaxis()->SetTitle("number of hits");
hsize->GetXaxis()->SetTitleSize(0.06);
hsize->GetYaxis()->SetTitleSize(0.06);
hsize->GetYaxis()->SetTitleOffset(1.5);
hsize->GetXaxis()->SetLabelSize(0.06);
hsize->GetYaxis()->SetLabelSize(0.06);
hsize->Draw();


TCanvas * c11  = new TCanvas("c11","c11",600,400);
c11->SetLeftMargin(0.2);
c11->SetRightMargin(0.1);
c11->SetTopMargin(0.1);
c11->SetBottomMargin(0.16);

ttb->Draw("sizeCol >>hcol(10,0,10)", hitcut);
hcol->SetName("hcol"); 
hcol->SetTitle("");
hcol->GetXaxis()->SetTitle("proj. cluster size [pixels]");
hcol->GetYaxis()->SetTitle("number of hits");
hcol->GetXaxis()->SetTitleSize(0.06);
hcol->GetYaxis()->SetTitleSize(0.06);
hcol->GetYaxis()->SetTitleOffset(1.5);
hcol->GetXaxis()->SetLabelSize(0.06);
hcol->GetYaxis()->SetLabelSize(0.06);


ttb->Draw("sizeRow >>hrow(10,0,10)", hitcut);
hrow->SetName("hrow"); 
hrow->SetTitle("");
hrow->GetXaxis()->SetTitle("proj. cluster size [pixels]");
hrow->GetYaxis()->SetTitle("number of hits");
hrow->GetXaxis()->SetTitleSize(0.06);
hrow->GetYaxis()->SetTitleSize(0.06);
hrow->GetYaxis()->SetTitleOffset(1.5);
hrow->GetXaxis()->SetLabelSize(0.06);
hrow->GetYaxis()->SetLabelSize(0.06);
hrow->SetLineColor(4);
hrow->SetLineStyle(2);

hcol->Draw();
hrow->Draw("same");

TLegend* lc = new TLegend(0.6,0.5,0.85,0.7);
lc->SetFillColor(kWhite); 
lc->SetBorderSize(0);
lc->AddEntry("hcol","Columns","l");
lc->AddEntry("hrow","Rows","l");
lc->Draw();

}
