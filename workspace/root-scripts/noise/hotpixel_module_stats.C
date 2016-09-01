{ 

 
gROOT->Reset(); 
gROOT->SetStyle("Plain"); 
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetTitleFillColor(0);
gStyle->SetPaintTextFormat("1f"); 

// TB data
TFile *ftb = new TFile("NoiseDB-DEP.root");

TTree *occupancytree = (TTree*) ftb->Get("Occupancy"); 
TTree *eventtree = (TTree*) ftb->Get("Event");  

// Parameters
int modID =11;
int npixelX = 128; 
int npixelY = 128; 

// --------------------------------------------------------------------------------------------------------

TCut modcut = Form("det == %d ", modID);
TCut roicut = ""; 

double xmin = 0 - 0.5;
double xmax = npixelX -1 + 0.5;

double ymin = 0 - 0.5;
double ymax = npixelY -1 + 0.5;



TCanvas * c1 = new TCanvas("c1","c1",600,400);   
c1->SetLeftMargin(0.13);
//c1->SetRightMargin(0.3);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.16);

eventtree->Draw("iEvt:nhits", modcut); 
TGraph *adcVsTime = new TGraph( eventtree->GetSelectedRows(), eventtree->GetV1(), eventtree->GetV2() );
adcVsTime->SetTitle("");
adcVsTime->GetXaxis()->SetTitle("iEvt");
adcVsTime->GetYaxis()->SetTitle("Pixels");
adcVsTime->Draw("a*");

adcVsTime->GetXaxis()->SetTitleSize(0.06);
adcVsTime->GetYaxis()->SetTitleSize(0.06);
adcVsTime->GetXaxis()->SetLabelSize(0.06);
adcVsTime->GetYaxis()->SetLabelSize(0.06);
adcVsTime->GetYaxis()->SetTitleOffset(1.16);

TCanvas * c2 = new TCanvas("c2","c2",600,400);   
c2->SetLeftMargin(0.13);
c2->SetRightMargin(0.2);
c2->SetTopMargin(0.1);
c2->SetBottomMargin(0.16);


// here status maps 
TH2D * statusMap = new TH2D("statusMap","statusMap",npixelX,xmin,xmax,npixelY,ymin,ymax);


occupancytree->Draw("px_x:px_y:status", modcut ); 
for (int i = 0; i<occupancytree->GetSelectedRows(); i++) {
  statusMap->Fill(occupancytree->GetV1()[i],occupancytree->GetV2()[i],occupancytree->GetV3()[i]);
} 

statusMap->SetTitle("");
statusMap->GetXaxis()->SetTitle("columns");
statusMap->GetYaxis()->SetTitle("rows");
statusMap->GetZaxis()->SetTitle("status code (>0: masked)");
statusMap->GetXaxis()->SetTitleSize(0.06);
statusMap->GetYaxis()->SetTitleSize(0.06);
statusMap->GetZaxis()->SetTitleSize(0.06);
statusMap->GetXaxis()->SetLabelSize(0.06);
statusMap->GetYaxis()->SetLabelSize(0.06);
statusMap->GetZaxis()->SetLabelSize(0.06);
statusMap->GetYaxis()->SetTitleOffset(1.16);
statusMap->GetZaxis()->SetTitleOffset(1.7);
statusMap->Draw("colz");



TCanvas * c3 = new TCanvas("c3","c3",600,400);   
c3->SetLeftMargin(0.13);
c3->SetRightMargin(0.23);
c3->SetTopMargin(0.1);
c3->SetBottomMargin(0.16);


TH2D * freqMap = new TH2D("freqMap","Firing Frequency",npixelX,xmin,xmax,npixelY,ymin,ymax);

occupancytree->Draw("px_x:px_y:hitFreq", modcut   ); 
for (int i = 0; i<occupancytree->GetSelectedRows(); i++) {
  freqMap->Fill(occupancytree->GetV1()[i],occupancytree->GetV2()[i],occupancytree->GetV3()[i]);
} 

freqMap->SetTitle("");
freqMap->GetXaxis()->SetTitle("columns");
freqMap->GetYaxis()->SetTitle("rows");
freqMap->GetZaxis()->SetTitle("occupancy");
freqMap->GetXaxis()->SetTitleSize(0.06);
freqMap->GetYaxis()->SetTitleSize(0.06);
freqMap->GetZaxis()->SetTitleSize(0.06);
freqMap->GetXaxis()->SetLabelSize(0.06);
freqMap->GetYaxis()->SetLabelSize(0.06);
freqMap->GetZaxis()->SetLabelSize(0.06);
freqMap->GetYaxis()->SetTitleOffset(1.15);
freqMap->GetZaxis()->SetTitleOffset(1.6);
freqMap->Draw("colz");


TCanvas * c4 = new TCanvas("c4","c4",600,400);   
c4->SetLeftMargin(0.13);
//c4->SetRightMargin(0.3);
c4->SetTopMargin(0.1);
c4->SetBottomMargin(0.16);
c4->SetLogy();


occupancytree->Draw(" status >>hstat", modcut  );
hstat->SetName("hstat");
hstat->SetTitle("");
hstat->GetXaxis()->SetTitle("status code (>0: masked)");
hstat->GetYaxis()->SetTitle("pixels");
hstat->GetXaxis()->SetTitleSize(0.06);
hstat->GetYaxis()->SetTitleSize(0.06);
hstat->GetXaxis()->SetLabelSize(0.06);
hstat->GetYaxis()->SetLabelSize(0.06);
hstat->GetYaxis()->SetTitleOffset(1.16);
hstat->Draw();


TCanvas * c5 = new TCanvas("c5","c5",750,400);   
c5->SetLeftMargin(0.13);
//c5->SetRightMargin(0.3);
c5->SetTopMargin(0.1);
c5->SetBottomMargin(0.16);
c5->SetLogy();


occupancytree->Draw(" hitFreq >>hfreq(1000000,0,100)",  modcut && roicut );
hfreq->SetName("hfreq");
hfreq->SetTitle("");  
hfreq->SetXTitle("occupancy");
hfreq->SetYTitle("pixels"); 
hfreq->GetXaxis()->SetRangeUser(0,0.01);
hfreq->GetXaxis()->SetTitleSize(0.06);
hfreq->GetYaxis()->SetTitleSize(0.06);
hfreq->GetXaxis()->SetLabelSize(0.06);
hfreq->GetYaxis()->SetLabelSize(0.06);
//hfreq->GetYaxis()->SetTitleOffset(1.2);
hfreq->Draw(); 




}


