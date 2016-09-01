{ 

 
gROOT->Reset(); 
gROOT->SetStyle("Plain"); 
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetTitleFillColor(0);
gStyle->SetPaintTextFormat("1f");

// TB data
TFile *ftb = new TFile("NoiseDB-M26.root"); 
TTree *eventtree = (TTree*) ftb->Get("Event");  

TCanvas * c1 = new TCanvas("c1","c1",1000,600);   
c1->SetLeftMargin(0.13);
//c1->SetRightMargin(0.3);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.16);

c1->Divide(2,3); 

c1->cd(1); 

eventtree->Draw("iEvt:nhits", "det == 0"); 
TGraph *adcVsTime0 = new TGraph( eventtree->GetSelectedRows(), eventtree->GetV1(), eventtree->GetV2() );
adcVsTime0->SetTitle("Firing Pixels Vs. iEvt Mod0");
adcVsTime0->GetXaxis()->SetTitle("iEvt");
adcVsTime0->GetYaxis()->SetTitle("Pixels");
adcVsTime0->Draw("a*");

c1->cd(2); 

eventtree->Draw("iEvt:nhits", "det == 1"); 
TGraph *adcVsTime1 = new TGraph( eventtree->GetSelectedRows(), eventtree->GetV1(), eventtree->GetV2() );
adcVsTime1->SetTitle("Firing Pixels Vs. iEvt Mod1");
adcVsTime1->GetXaxis()->SetTitle("iEvt");
adcVsTime1->GetYaxis()->SetTitle("Pixels");
adcVsTime1->Draw("a*");

c1->cd(3); 

eventtree->Draw("iEvt:nhits", "det == 2"); 
TGraph *adcVsTime2 = new TGraph( eventtree->GetSelectedRows(), eventtree->GetV1(), eventtree->GetV2() );
adcVsTime2->SetTitle("Firing Pixels Vs. iEvt Mod2");
adcVsTime2->GetXaxis()->SetTitle("iEvt");
adcVsTime2->GetYaxis()->SetTitle("Pixels");
adcVsTime2->Draw("a*");

c1->cd(4); 

eventtree->Draw("iEvt:nhits", "det == 3"); 
TGraph *adcVsTime3 = new TGraph( eventtree->GetSelectedRows(), eventtree->GetV1(), eventtree->GetV2() );
adcVsTime3->SetTitle("Firing Pixels Vs. iEvt Mod3");
adcVsTime3->GetXaxis()->SetTitle("iEvt");
adcVsTime3->GetYaxis()->SetTitle("Pixels");
adcVsTime3->Draw("a*");

c1->cd(5); 

eventtree->Draw("iEvt:nhits", "det == 4"); 
TGraph *adcVsTime4 = new TGraph( eventtree->GetSelectedRows(), eventtree->GetV1(), eventtree->GetV2() );
adcVsTime4->SetTitle("Firing Pixels Vs. iEvt Mod4");
adcVsTime4->GetXaxis()->SetTitle("iEvt");
adcVsTime4->GetYaxis()->SetTitle("Pixels");
adcVsTime4->Draw("a*");

c1->cd(6); 

eventtree->Draw("iEvt:nhits", "det == 5"); 
TGraph *adcVsTime5 = new TGraph( eventtree->GetSelectedRows(), eventtree->GetV1(), eventtree->GetV2() );
adcVsTime5->SetTitle("Firing Pixels Vs. iEvt Mod5");
adcVsTime5->GetXaxis()->SetTitle("iEvt");
adcVsTime5->GetYaxis()->SetTitle("Pixels");
adcVsTime5->Draw("a*");



TCanvas * c2  = new TCanvas("c2","c2",600,400);
c2->SetLeftMargin(0.2);
//c2->SetRightMargin(0.26);
c2->SetTopMargin(0.1);
c2->SetBottomMargin(0.16);

eventtree->Draw("nhits >>hhit0(200,0,200)", "det == 0" );
hhit0->SetName("hhit0");
hhit0->SetXTitle("hits per event");
hhit0->SetYTitle("events"); 
hhit0->SetTitle(""); 
//hhit0->GetYaxis()->SetRangeUser(0,3000);
hhit0->GetXaxis()->SetTitleSize(0.06);
hhit0->GetYaxis()->SetTitleSize(0.06);
hhit0->GetXaxis()->SetLabelSize(0.06);
hhit0->GetYaxis()->SetLabelSize(0.06);
hhit0->GetYaxis()->SetTitleOffset(1.2);
hhit0->SetLineColor(1);


eventtree->Draw("nhits >>hhit1(200,0,200)", "det == 1" );
hhit1->SetName("hhit1");
hhit1->SetLineColor(2);


eventtree->Draw("nhits >>hhit2(200,0,200)", "det == 2" );
hhit2->SetName("hhit2");
hhit2->SetLineColor(3);

eventtree->Draw("nhits >>hhit3(200,0,200)", "det == 3" );
hhit3->SetName("hhit3");
hhit3->SetLineColor(4);


eventtree->Draw("nhits >>hhit4(200,0,200)", "det == 4" );
hhit4->SetName("hhit4");
hhit4->SetLineColor(5);


eventtree->Draw("nhits >>hhit5(200,0,200)", "det == 5" );
hhit5->SetName("hhit5");
hhit5->SetLineColor(6);


hhit0->Draw();
hhit1->Draw("same");
hhit2->Draw("same");
hhit3->Draw("same");
hhit4->Draw("same");
hhit5->Draw("same");

TLegend* lc = new TLegend(0.6,0.5,0.85,0.7);
lc->SetFillColor(0); 
lc->SetBorderSize(0);
lc->AddEntry("hhit0","M26 plane 0","l");
lc->AddEntry("hhit1","M26 plane 1","l"); 
lc->AddEntry("hhit2","M26 plane 2","l");
lc->AddEntry("hhit3","M26 plane 3","l");
lc->AddEntry("hhit4","M26 plane 4","l");
lc->AddEntry("hhit5","M26 plane 5","l");
lc->Draw();
 

TCanvas * c3  = new TCanvas("c3","c3",600,400);
c3->SetLeftMargin(0.2);
//c3->SetRightMargin(0.26);
c3->SetTopMargin(0.1);
c3->SetBottomMargin(0.16);




eventtree->Draw("ngoodhits >>hhit0g(200,0,200)", "det == 0 && iEvt > 20000" );
hhit0g->SetName("hhit0g");
hhit0g->SetXTitle("good hits per event");
hhit0g->SetYTitle("events"); 
hhit0g->SetTitle(""); 
//hhit0g->GetYaxis()->SetRangeUser(0,3000);
hhit0g->GetXaxis()->SetTitleSize(0.06);
hhit0g->GetYaxis()->SetTitleSize(0.06);
hhit0g->GetXaxis()->SetLabelSize(0.06);
hhit0g->GetYaxis()->SetLabelSize(0.06);
hhit0g->GetYaxis()->SetTitleOffset(1.15);
hhit0g->SetLineColor(1);


eventtree->Draw("ngoodhits >>hhit1g(200,0,200)", "det == 1 && iEvt > 20000" );
hhit1g->SetName("hhit1g");
hhit1g->SetLineColor(2);


eventtree->Draw("ngoodhits >>hhit2g(200,0,200)", "det == 2 && iEvt > 20000" );
hhit2g->SetName("hhit2g");
hhit2g->SetLineColor(3);


eventtree->Draw("ngoodhits >>hhit3g(200,0,200)", "det == 3 && iEvt > 20000" );
hhit3g->SetName("hhit3g");
hhit3g->SetLineColor(4);


eventtree->Draw("ngoodhits >>hhit4g(200,0,200)", "det == 4 && iEvt > 20000" );
hhit4g->SetName("hhit4g");
hhit4g->SetLineColor(5);


eventtree->Draw("ngoodhits >>hhit5g(200,0,200)", "det == 5 && iEvt > 20000" );
hhit5g->SetName("hhit5g");
hhit5g->SetLineColor(6);


hhit0g->Draw();
hhit1g->Draw("same");
hhit2g->Draw("same");
hhit3g->Draw("same");
hhit4g->Draw("same");
hhit5g->Draw("same");

TLegend* lcg = new TLegend(0.6,0.5,0.85,0.7);
lcg->SetFillColor(0); 
lcg->SetBorderSize(0);
lcg->Draw();
lcg->AddEntry("hhit0g","M26 plane 0","l");
lcg->AddEntry("hhit1g","M26 plane 1","l"); 
lcg->AddEntry("hhit2g","M26 plane 2","l");
lcg->AddEntry("hhit3g","M26 plane 3","l");
lcg->AddEntry("hhit4g","M26 plane 4","l");
lcg->AddEntry("hhit5g","M26 plane 5","l");


}


