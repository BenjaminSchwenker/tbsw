{ 

gROOT->Reset(); 
gROOT->SetStyle("Plain"); 

gStyle->SetOptFit();
gStyle->SetPalette(1);
gStyle->SetCanvasColor(0);
gStyle->SetTitleFillColor(0);
gStyle->SetPaintTextFormat("1f");
//gStyle->SetOptStat(0);

// 
// Analysis cuts 
TCut hitcut = "chi2 < 50  && hasTrack == 0 && hitQuality == 0";
TCut roicut = "";
TCut trkcut = "";

//
// event data 
TFile *ftb = new TFile("root-files/L3-Histos-run000096_svdpxdteldigits_complete.root");

TTree *ttb = (TTree*) ftb->Get("Hit");
TTree *ttbtrack = (TTree*) ftb->Get("Track");


// ----------------------------------------------------------------------

Double_t Norm; 

TCanvas * c0  = new TCanvas("c0","c0",600,400);
c0->SetLeftMargin(0.2);
c0->SetRightMargin(0.1);
c0->SetTopMargin(0.1);
c0->SetBottomMargin(0.16);

ttb->Draw("cellV_hit:cellU_hit>>hspot", hitcut);
hspot->SetName("hspot"); 
hspot->GetXaxis()->SetTitle("uCells");
hspot->GetYaxis()->SetTitle("vCells");
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

ttbtrack->Draw("cellV_fit:cellU_fit>>htrackspot", trkcut && roicut);
htrackspot->SetName("htrackspot"); 
htrackspot->GetXaxis()->SetTitle("uCells");
htrackspot->GetYaxis()->SetTitle("vCells");
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


ttb->Draw("seedChargeU >>hChargeU(201,0,201)", hitcut);
hChargeU->SetName("hChargeU"); 
hChargeU->SetTitle("");
hChargeU->GetXaxis()->SetTitle("seed charge - u strips [ADU]");
hChargeU->GetYaxis()->SetTitle("number of hits");
hChargeU->GetXaxis()->SetTitleSize(0.06);
hChargeU->GetYaxis()->SetTitleSize(0.06);
hChargeU->GetYaxis()->SetTitleOffset(1.5);
hChargeU->GetXaxis()->SetLabelSize(0.06);
hChargeU->GetYaxis()->SetLabelSize(0.06);
hChargeU->Draw();

TCanvas * c71  = new TCanvas("c71","c71",600,400);
c71->SetLeftMargin(0.2);
c71->SetRightMargin(0.1);
c71->SetTopMargin(0.1);
c71->SetBottomMargin(0.16);


ttb->Draw("seedChargeV >>hChargeV(201,0,201)", hitcut);
hChargeV->SetName("hChargeV"); 
hChargeV->SetTitle("");
hChargeV->GetXaxis()->SetTitle("seed charge - v strips [ADU]");
hChargeV->GetYaxis()->SetTitle("number of hits");
hChargeV->GetXaxis()->SetTitleSize(0.06);
hChargeV->GetYaxis()->SetTitleSize(0.06);
hChargeV->GetYaxis()->SetTitleOffset(1.5);
hChargeV->GetXaxis()->SetLabelSize(0.06);
hChargeV->GetYaxis()->SetLabelSize(0.06);
hChargeV->Draw();


TCanvas * c8  = new TCanvas("c8","c8",600,400);
c8->SetLeftMargin(0.2);
c8->SetRightMargin(0.1);
c8->SetTopMargin(0.1);
c8->SetBottomMargin(0.16);


ttb->Draw("clusterChargeU >>hCChargeU(201,0,201)", hitcut);
hCChargeU->SetName("hCChargeU"); 
hCChargeU->SetTitle("");
hCChargeU->GetXaxis()->SetTitle("cluster charge - u strips [ADU]");
hCChargeU->GetYaxis()->SetTitle("number of hits");
hCChargeU->GetXaxis()->SetTitleSize(0.06);
hCChargeU->GetYaxis()->SetTitleSize(0.06);
hCChargeU->GetYaxis()->SetTitleOffset(1.5);
hCChargeU->GetXaxis()->SetLabelSize(0.06);
hCChargeU->GetYaxis()->SetLabelSize(0.06);
hCChargeU->Draw();

TCanvas * c81  = new TCanvas("c81","c81",600,400);
c81->SetLeftMargin(0.2);
c81->SetRightMargin(0.1);
c81->SetTopMargin(0.1);
c81->SetBottomMargin(0.16);


ttb->Draw("clusterChargeV >>hCChargeV(201,0,201)", hitcut);
hCChargeV->SetName("hCChargeV"); 
hCChargeV->SetTitle("");
hCChargeV->GetXaxis()->SetTitle("cluster charge - v strips [ADU]");
hCChargeV->GetYaxis()->SetTitle("number of hits");
hCChargeV->GetXaxis()->SetTitleSize(0.06);
hCChargeV->GetYaxis()->SetTitleSize(0.06);
hCChargeV->GetYaxis()->SetTitleOffset(1.5);
hCChargeV->GetXaxis()->SetLabelSize(0.06);
hCChargeV->GetYaxis()->SetLabelSize(0.06);
hCChargeV->Draw();


TCanvas * c10  = new TCanvas("c10","c10",600,400);
c10->SetLeftMargin(0.2);
c10->SetRightMargin(0.1);
c10->SetTopMargin(0.1);
c10->SetBottomMargin(0.16);

ttb->Draw("size >>hsize(20,0,20)", hitcut);
hsize->SetName("hsize"); 
hsize->SetTitle("");
hsize->GetXaxis()->SetTitle("uSize + vSize [cells]");
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

ttb->Draw("sizeU >>hcol(10,0,10)", hitcut);
hcol->SetName("hcol"); 
hcol->SetTitle("");
hcol->GetXaxis()->SetTitle("cluster size [cells]");
hcol->GetYaxis()->SetTitle("number of hits");
hcol->GetXaxis()->SetTitleSize(0.06);
hcol->GetYaxis()->SetTitleSize(0.06);
hcol->GetYaxis()->SetTitleOffset(1.5);
hcol->GetXaxis()->SetLabelSize(0.06);
hcol->GetYaxis()->SetLabelSize(0.06);


ttb->Draw("sizeV >>hrow(10,0,10)", hitcut);
hrow->SetName("hrow"); 
hrow->SetTitle("");
hrow->GetXaxis()->SetTitle("cluster size [cells]");
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
lc->AddEntry("hcol","uSize","l");
lc->AddEntry("hrow","vSize","l");
lc->Draw();

}
