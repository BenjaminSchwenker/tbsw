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
TCut hitcut = "";


//
// event data 
TFile *ftb = new TFile("root-files/MC-Histos-dutmc.root");
TTree *ttb = (TTree*) ftb->Get("Cluster");



// ----------------------------------------------------------------------

Double_t Norm; 

TCanvas * c0  = new TCanvas("c0","c0",600,400);
c0->SetLeftMargin(0.2);
c0->SetRightMargin(0.1);
c0->SetTopMargin(0.1);
c0->SetBottomMargin(0.16);

ttb->Draw("trackPosV:trackPosU>>hspot", hitcut);
hspot->SetName(""); 
hspot->GetXaxis()->SetTitle("track position U [mm]");
hspot->GetYaxis()->SetTitle("track position V [mm]");
hspot->GetZaxis()->SetTitle("number of clusters");
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

ttb->Draw("clusterStartCellIdV:clusterStartCellIdU>>htrackspot", hitcut );
htrackspot->SetName("htrackspot"); 
htrackspot->GetXaxis()->SetTitle("start cell U [cell ID]");
htrackspot->GetYaxis()->SetTitle("start cell V [cell ID]");
htrackspot->GetZaxis()->SetTitle("number of clusters");
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

ttb->Draw("clusterStartCellPosV:clusterStartCellPosU>>hbeamspot", hitcut);
hbeamspot->SetName("hbeamspot"); 
hbeamspot->GetXaxis()->SetTitle("start position U [mm]");
hbeamspot->GetYaxis()->SetTitle("start position V [mm]");
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

ttb->Draw("(clusterPosU - trackPosU)*1000 >>hu(201,-55.,+55.)", hitcut);
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

ttb->Draw("(clusterPosV - trackPosV)*1000 >>hv(201,-55.,+55.)", hitcut);
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

ttb->Draw("clusterSize >>hsize(20,0,20)", hitcut);
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

ttb->Draw("clusterSizeU >>hcol(10,0,10)", hitcut);
hcol->SetName("hcol"); 
hcol->SetTitle("");
hcol->GetXaxis()->SetTitle("size [cells]");
hcol->GetYaxis()->SetTitle("number of hits");
hcol->GetXaxis()->SetTitleSize(0.06);
hcol->GetYaxis()->SetTitleSize(0.06);
hcol->GetYaxis()->SetTitleOffset(1.5);
hcol->GetXaxis()->SetLabelSize(0.06);
hcol->GetYaxis()->SetLabelSize(0.06);


ttb->Draw("clusterSizeV >>hrow(10,0,10)", hitcut);
hrow->SetName("hrow"); 
hrow->SetTitle("");
hrow->GetXaxis()->SetTitle("size [cells]");
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
lc->AddEntry("hcol","sizeU","l");
lc->AddEntry("hrow","sizeV","l");
lc->Draw();



TCanvas * c31  = new TCanvas("c31","c31",600,400);
c31->SetLeftMargin(0.2);
c31->SetRightMargin(0.1);
c31->SetTopMargin(0.1);
c31->SetBottomMargin(0.16);


ttb->Draw("(clusterPosU - trackPosU)*1000 >>hu1(201,-25.,+25.)", hitcut && "clusterSizeU==1");
hu1->SetName("hu1"); 
hu1->SetTitle("");
hu1->GetXaxis()->SetTitle("residuals u [#mum]");
hu1->GetYaxis()->SetTitle("number of clusters");
hu1->GetXaxis()->SetTitleSize(0.06);
hu1->GetYaxis()->SetTitleSize(0.06);
hu1->GetYaxis()->SetTitleOffset(1.5);
hu1->GetXaxis()->SetLabelSize(0.06);
hu1->GetYaxis()->SetLabelSize(0.06);


ttb->Draw("(clusterPosU - trackPosU)*1000 >>hu2(201,-25.,+25.)", hitcut && "clusterSizeU==2");
hu2->SetName("hu2"); 
hu2->SetTitle("");
hu2->GetXaxis()->SetTitle("residuals u [#mum]");
hu2->GetYaxis()->SetTitle("number of clusters");
hu2->GetXaxis()->SetTitleSize(0.06);
hu2->GetYaxis()->SetTitleSize(0.06);
hu2->GetYaxis()->SetTitleOffset(1.5);
hu2->GetXaxis()->SetLabelSize(0.06);
hu2->GetYaxis()->SetLabelSize(0.06);
hu2->SetLineColor(4);
hu2->SetLineStyle(2);


hu1->Draw();
hu2->Draw("same");

double rmsu1 = hu1->GetRMS();
double rmsu2 = hu2->GetRMS();

TLegend* lcu = new TLegend(0.6,0.5,0.85,0.7);
lcu->SetFillColor(kWhite); 
lcu->SetBorderSize(0);
lcu->AddEntry("hu1",Form("size=1: RMS=%.1f", rmsu1),"l");
lcu->AddEntry("hu2",Form("size=2: RMS=%.1f", rmsu2),"l");
lcu->Draw();

TCanvas * c41  = new TCanvas("c41","c41",600,400);
c41->SetLeftMargin(0.2);
c41->SetRightMargin(0.1);
c41->SetTopMargin(0.1);
c41->SetBottomMargin(0.16);


ttb->Draw("(clusterPosV - trackPosV)*1000 >>hv1(201,-25.,+25.)", hitcut && "clusterSizeV==1");
hv1->SetName("hu1"); 
hv1->SetTitle("");
hv1->GetXaxis()->SetTitle("residuals v [#mum]");
hv1->GetYaxis()->SetTitle("number of clusters");
hv1->GetXaxis()->SetTitleSize(0.06);
hv1->GetYaxis()->SetTitleSize(0.06);
hv1->GetYaxis()->SetTitleOffset(1.5);
hv1->GetXaxis()->SetLabelSize(0.06);
hv1->GetYaxis()->SetLabelSize(0.06);


ttb->Draw("(clusterPosV - trackPosV)*1000 >>hv2(201,-25.,+25.)", hitcut && "clusterSizeV==2");
hv2->SetName("hv2"); 
hv2->SetTitle("");
hv2->GetXaxis()->SetTitle("residuals v [#mum]");
hv2->GetYaxis()->SetTitle("number of clusters");
hv2->GetXaxis()->SetTitleSize(0.06);
hv2->GetYaxis()->SetTitleSize(0.06);
hv2->GetYaxis()->SetTitleOffset(1.5);
hv2->GetXaxis()->SetLabelSize(0.06);
hv2->GetYaxis()->SetLabelSize(0.06);
hv2->SetLineColor(4);
hv2->SetLineStyle(2);


hv1->Draw();
hv2->Draw("same");

double rmsv1 = hv1->GetRMS();
double rmsv2 = hv2->GetRMS();

TLegend* lcv = new TLegend(0.6,0.5,0.85,0.7);
lcv->SetFillColor(kWhite); 
lcv->SetBorderSize(0);
lcv->AddEntry("hv1",Form("size=1: RMS=%.1f", rmsv1),"l");
lcv->AddEntry("hv2",Form("size=2: RMS=%.1f", rmsv2),"l");
lcv->Draw();



}
