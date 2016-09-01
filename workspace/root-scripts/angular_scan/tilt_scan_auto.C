{ 

gROOT->Reset(); 
gROOT->SetStyle("Plain"); 


// display mode
gStyle->SetOptFit();
gStyle->SetPalette(1);
gStyle->SetCanvasColor(0);
gStyle->SetTitleFillColor(0);
gStyle->SetPaintTextFormat("1f");


//
// MC data 
TFile *fmc = new TFile("root-final/tune_B.root");
TTree *tmc = (TTree*) fmc->Get("Hit");

//
// TB data 
TFile *ftb = new TFile("root-final/tbtiltscan.root");
TTree *ttb = (TTree*) ftb->Get("Hit");

Double_t gain_MC = 1; // 170; 

Int_t nMC= 10;

Double_t tilt_MC[nMC] =  {0,10,20,25,30,35,40,45,50,60};
Double_t etilt_MC[nMC] = {0};

Double_t size_MC[nMC] = { 0 };
Double_t esize_MC[nMC]= { 0 };

Double_t sizeCol_MC[nMC] = { 0 };
Double_t esizeCol_MC[nMC]= { 0 };

Double_t sizeRow_MC[nMC] = { 0 };
Double_t esizeRow_MC[nMC]= { 0 };

Double_t charge_MC[nMC] = { 0 };
Double_t echarge_MC[nMC]= { 0 };

Double_t resoU_MC[nMC]   = { 0 };
Double_t eresoU_MC[nMC]  = { 0 };

Double_t resoV_MC[nMC]   = { 0 };
Double_t eresoV_MC[nMC]  = { 0 };

Int_t nTB= 8;

Double_t tilt_TB[nTB] =  {0,10,20,30,40,45,50,60};
Double_t etilt_TB[nTB] = {0};

Double_t size_TB[nTB] = { 0 };
Double_t esize_TB[nTB]= { 0 };

Double_t sizeCol_TB[nTB] = { 0 };
Double_t esizeCol_TB[nTB]= { 0 };

Double_t sizeRow_TB[nTB] = { 0 };
Double_t esizeRow_TB[nTB]= { 0 };

Double_t charge_TB[nTB] = { 0 };
Double_t echarge_TB[nTB]= { 0 };

Double_t resoU_TB[nTB]   = { 0 };
Double_t eresoU_TB[nTB]  = { 0 };

Double_t resoV_TB[nTB]   = { 0 };
Double_t eresoV_TB[nTB]  = { 0 };


std::vector<int> run_list(500,-1);
run_list[435] = 0;
run_list[420] = 1;
run_list[416] = 2;
run_list[409] = 3;
run_list[393] = 4;
run_list[391] = 5;
run_list[390] = 5;
run_list[384] = 6;
run_list[378] = 7;


// ----------------------------------------------------------------------
// Do not modify below here

// declare/init histos

TH1F *hmc_size[nMC];
TH1F *hmc_sizecol[nMC];
TH1F *hmc_sizerow[nMC];
TH1F *hmc_charge[nMC];
TH1F *hmc_resou[nMC];
TH1F *hmc_resov[nMC];


for (int i=0; i<nMC; i++){
  hmc_size[i] = new TH1F("","",10,0,10);
  hmc_sizecol[i] = new TH1F("","",10,0,10);
  hmc_sizerow[i] = new TH1F("","",10,0,10);
  hmc_charge[i] = new TH1F("","",200,0,200);
  hmc_resou[i] = new TH1F("","",101,-75.,+75.);
  hmc_resov[i] = new TH1F("","",101,-75.,+75.);
}


TH1F *htb_size[nTB];
TH1F *htb_sizecol[nTB];
TH1F *htb_sizerow[nTB];
TH1F *htb_charge[nTB];
TH1F *htb_residu[nTB];
TH1F *htb_residv[nTB];
TH1F *htb_telu[nTB];
TH1F *htb_telv[nTB];


for (int i=0; i<nTB; i++){
  htb_size[i] = new TH1F("","",10,0,10);
  htb_sizecol[i] = new TH1F("","",10,0,10);
  htb_sizerow[i] = new TH1F("","",10,0,10);
  htb_charge[i] = new TH1F("","",200,0,200);
  htb_residu[i] = new TH1F("","",251,-105.,+105.);
  htb_residv[i] = new TH1F("","",251,-105.,+105.);
  htb_telu[i] = new TH1F("","",551,0,50);
  htb_telv[i] = new TH1F("","",551,0,50);
}

// fill mc histos

Int_t _hasTrack;
Int_t _iRun;
Int_t _ndof; 
Int_t _col_hit;
Int_t _row_hit;
Int_t _size; 
Int_t _sizeCol;
Int_t _sizeRow;
Double_t _chi2;
Double_t  _u_fit; 
Double_t  _v_fit; 
Double_t  _u_hit; 
Double_t  _v_hit; 
Double_t  _clusterCharge; 
Double_t _u_fiterr; 
Double_t _v_fiterr;

tmc->SetBranchAddress("hasTrack",&_hasTrack);
tmc->SetBranchAddress("iRun",&_iRun);
tmc->SetBranchAddress("col_hit",&_col_hit);
tmc->SetBranchAddress("row_hit",&_row_hit);
tmc->SetBranchAddress("size",&_size);
tmc->SetBranchAddress("sizeCol",&_sizeCol);
tmc->SetBranchAddress("sizeRow",&_sizeRow);
tmc->SetBranchAddress("chi2",&_chi2);
tmc->SetBranchAddress("u_fit",&_u_fit);
tmc->SetBranchAddress("v_fit",&_v_fit);
tmc->SetBranchAddress("u_hit",&_u_hit);
tmc->SetBranchAddress("v_hit",&_v_hit);
tmc->SetBranchAddress("clusterCharge",&_clusterCharge);


for(int i=0; i< tmc->GetEntries(); i++){
    tmc->GetEntry(i);
    if (  _hasTrack == 0 ){
      
      // Fill histos
      hmc_size[_iRun]->Fill(_size);
      hmc_sizerow[_iRun]->Fill(_sizeCol);  //  mc col - row swapped
      hmc_sizecol[_iRun]->Fill(_sizeRow);  //  mc col - row swapped
      hmc_charge[_iRun]->Fill(_clusterCharge/gain_MC);
      hmc_resov[_iRun]->Fill(  1000*(_u_hit-_u_fit) );  //  mc col - row swapped
      hmc_resou[_iRun]->Fill(  1000*(_v_hit-_v_fit) );  //  mc col - row swapped
      
    }
}


cout << "MC summary: " << endl; 

for (int i=0; i<nMC ; ++i) {
  
  size_MC[i] = hmc_size[i]->GetMean();
  esize_MC[i]= 0;

  sizeCol_MC[i] = hmc_sizecol[i]->GetMean();
  esizeCol_MC[i]= 0;

  sizeRow_MC[i] = hmc_sizerow[i]->GetMean();
  esizeRow_MC[i]= 0;

  charge_MC[i] = hmc_charge[i]->GetMaximumBin();
  echarge_MC[i]= 0;

  resoU_MC[i]   = hmc_resou[i]->GetRMS();
  eresoU_MC[i]  = 0;

  resoV_MC[i]   = hmc_resov[i]->GetRMS();
  eresoV_MC[i]  = 0;

  cout << "   theta: "   << tilt_MC[i] << endl;
  cout << "   charge: "  << charge_MC[i] << endl;  
  cout << "   size: "    << size_MC[i] << endl; 
  cout << "   sizeCol: " << sizeCol_MC[i] << endl; 
  cout << "   sizeRow: " << sizeRow_MC[i] << endl; 
  cout << "   resoU: "   << resoU_MC[i] << endl;  
  cout << "   resoV: "   << resoV_MC[i]  << endl; 

  
}

// fill tb histos 

ttb->SetBranchAddress("hasTrack",&_hasTrack);
ttb->SetBranchAddress("iRun",&_iRun);
ttb->SetBranchAddress("ndof",&_ndof);
ttb->SetBranchAddress("col_hit",&_col_hit);
ttb->SetBranchAddress("row_hit",&_row_hit);
ttb->SetBranchAddress("size",&_size);
ttb->SetBranchAddress("sizeCol",&_sizeCol);
ttb->SetBranchAddress("sizeRow",&_sizeRow);
ttb->SetBranchAddress("chi2",&_chi2);
ttb->SetBranchAddress("u_fit",&_u_fit);
ttb->SetBranchAddress("v_fit",&_v_fit);
ttb->SetBranchAddress("u_hit",&_u_hit);
ttb->SetBranchAddress("v_hit",&_v_hit);
ttb->SetBranchAddress("clusterCharge",&_clusterCharge);
ttb->SetBranchAddress("u_fiterr",&_u_fiterr);
ttb->SetBranchAddress("v_fiterr",&_v_fiterr);


for(int i=0; i< ttb->GetEntries(); i++){
    ttb->GetEntry(i);
    if (  _hasTrack == 0 && _chi2 <100 && _ndof == 8){

      int value = run_list[_iRun];
      
      // Fill histos
      htb_size[value]->Fill(_size);
      htb_sizecol[value]->Fill(_sizeCol);
      htb_sizerow[value]->Fill(_sizeRow);
      htb_charge[value]->Fill(_clusterCharge);
      htb_residu[value]->Fill(1000*(_u_hit-_u_fit));
      htb_residv[value]->Fill(1000*(_v_hit-_v_fit));
      htb_telu[value]->Fill(1000*_u_fiterr);
      htb_telv[value]->Fill(1000*_v_fiterr);     
      
    }
}


cout << "TB summary: " << endl; 

for (int i=0; i<nTB ; ++i) {
  
  size_TB[i] = htb_size[i]->GetMean();
  esize_TB[i]= 0;

  sizeCol_TB[i] = htb_sizecol[i]->GetMean();
  esizeCol_TB[i]= 0;

  sizeRow_TB[i] = htb_sizerow[i]->GetMean();
  esizeRow_TB[i]= 0;

  charge_TB[i] = htb_charge[i]->GetMaximumBin();
  echarge_TB[i]= 1;

  double tu = htb_telu[i]->GetMean();
  double ru = htb_residu[i]->GetRMS();
  resoU_TB[i] = TMath::Sqrt( ru*ru - tu*tu ); 
  eresoU_TB[i]  = 0;

  double tv = htb_telv[i]->GetMean();
  double rv = htb_residv[i]->GetRMS();
  resoV_TB[i] = TMath::Sqrt( rv*rv - tv*tv ); 
  eresoV_TB[i]  = 0;

  cout << "   theta: "   << tilt_TB[i] << endl;
  cout << "   charge: "  << charge_TB[i] << endl;  
  cout << "   size: "    << size_TB[i] << endl; 
  cout << "   sizeCol: " << sizeCol_TB[i] << endl; 
  cout << "   sizeRow: " << sizeRow_TB[i] << endl; 
  cout << "   residU: "  << ru << endl;  
  cout << "   residV: "  << rv  << endl; 
  cout << "   telU: "    << tu << endl;  
  cout << "   telV: "    << tv  << endl; 
  cout << "   resoU: "   << resoU_TB[i] << endl;  
  cout << "   resoV: "   << resoV_TB[i]  << endl; 

  
}


// start plotting


TCanvas * c1 = new TCanvas("c1","c1",600,400);   
c1->SetLeftMargin(0.15);
//c1->SetRightMargin(0.3);
c1->SetTopMargin(0.1);
c1->SetBottomMargin(0.16);


TH1F *vFrame1 = c1.DrawFrame(0,1,70,4);
vFrame1->SetXTitle("#theta [degree]");
vFrame1->SetYTitle("proj. cluster size [px]");
vFrame1->GetXaxis()->SetTitleSize(0.06);
vFrame1->GetYaxis()->SetTitleSize(0.06);
vFrame1->GetXaxis()->SetLabelSize(0.06);
vFrame1->GetYaxis()->SetLabelSize(0.06);
//vFrame1->GetYaxis()->SetTitleOffset(1.8);


gr_sizeCol_MC = new TGraphErrors(nMC,tilt_MC,sizeCol_MC,etilt_MC,esizeCol_MC);
gr_sizeCol_MC->SetName("gr_sizeCol_MC");
gr_sizeCol_MC->SetLineColor(kBlack);
gr_sizeCol_MC->SetLineStyle(2);
gr_sizeCol_MC->SetLineWidth(2);

gr_sizeRow_MC = new TGraphErrors(nMC,tilt_MC,sizeRow_MC,etilt_MC,esizeRow_MC);
gr_sizeRow_MC->SetName("gr_sizeRow_MC");
gr_sizeRow_MC->SetLineColor(kRed);
gr_sizeRow_MC->SetLineStyle(5);
gr_sizeRow_MC->SetLineWidth(2);

gr_sizeCol_TB = new TGraphErrors(nTB,tilt_TB,sizeCol_TB,etilt_TB,esizeCol_TB);
gr_sizeCol_TB->SetName("gr_sizeCol_TB");
gr_sizeCol_TB->SetMarkerColor(kBlack);
gr_sizeCol_TB->SetMarkerStyle(22);

gr_sizeRow_TB = new TGraphErrors(nTB,tilt_TB,sizeRow_TB,etilt_TB,esizeRow_TB);
gr_sizeRow_TB->SetName("gr_sizeRow_TB");
gr_sizeRow_TB->SetMarkerColor(kRed);
gr_sizeRow_TB->SetMarkerStyle(21);

gr_sizeCol_MC->Draw("L");
gr_sizeRow_MC->Draw("L");
gr_sizeCol_TB->Draw("P");
gr_sizeRow_TB->Draw("P");

TLegend* l_size = new TLegend(0.18,0.6,0.50,0.85);
l_size->SetFillColor(kWhite); 
l_size->SetBorderSize(0);
l_size->AddEntry("gr_sizeCol_TB","TB column size ","p");
l_size->AddEntry("gr_sizeCol_MC","MC column size ","l");
l_size->AddEntry("gr_sizeRow_TB","TB row size","p");
l_size->AddEntry("gr_sizeRow_MC","MC row size","l");
l_size->Draw();

c1.Modified();
c1.Update();


TCanvas * c2 = new TCanvas("c2","c2",600,400);   
c2->SetLeftMargin(0.15);
//c2->SetRightMargin(0.3);
c2->SetTopMargin(0.1);
c2->SetBottomMargin(0.16);


TH1F *vFrame2 = c2.DrawFrame(0,1,70,100);
vFrame2->SetXTitle("#theta [degree]");
vFrame2->SetYTitle("MPV cluster charge [ADU]");
vFrame2->GetXaxis()->SetTitleSize(0.06);
vFrame2->GetYaxis()->SetTitleSize(0.06);
vFrame2->GetXaxis()->SetLabelSize(0.06);
vFrame2->GetYaxis()->SetLabelSize(0.06);
//vFrame2->GetYaxis()->SetTitleOffset(1.8);


gr_charge_MC = new TGraphErrors(nMC,tilt_MC,charge_MC,etilt_MC,echarge_MC);
gr_charge_MC->SetName("gr_charge_MC");
gr_charge_MC->SetLineColor(kRed);
gr_charge_MC->SetLineStyle(2);
gr_charge_MC->SetLineWidth(2);
gr_charge_MC->Draw("LP");

gr_charge_TB = new TGraphErrors(nTB,tilt_TB,charge_TB,etilt_TB,echarge_TB);
gr_charge_TB->SetName("gr_charge_TB");
gr_charge_TB->SetMarkerColor(kRed);
gr_charge_TB->SetMarkerStyle(22);
gr_charge_TB->Draw("P");

TLegend* l_charge = new TLegend(0.18,0.6,0.50,0.85);
l_charge->SetFillColor(kWhite); 
l_charge->SetBorderSize(0);
l_charge->AddEntry("gr_charge_TB","TB cluster charge","p");
l_charge->AddEntry("gr_charge_MC","MC cluster charge","lp");
l_charge->Draw();

c2.Modified();
c2.Update();


TCanvas * c3 = new TCanvas("c3","c3",600,400);   
c3->SetLeftMargin(0.15);
//c3->SetRightMargin(0.3);
c3->SetTopMargin(0.1);
c3->SetBottomMargin(0.16);


TH1F *vFrame3 = c3.DrawFrame(0,1,70,35);
vFrame3->SetXTitle("#theta [degree]");
vFrame3->SetYTitle("spatial resolution [#mum]");
vFrame3->GetXaxis()->SetTitleSize(0.06);
vFrame3->GetYaxis()->SetTitleSize(0.06);
vFrame3->GetXaxis()->SetLabelSize(0.06);
vFrame3->GetYaxis()->SetLabelSize(0.06);
//vFrame3->GetYaxis()->SetTitleOffset(1.8);


gr_resoU_MC = new TGraphErrors(nMC,tilt_MC,resoU_MC,etilt_MC,eresoU_MC);
gr_resoU_MC->SetName("gr_resoU_MC");
gr_resoU_MC->SetLineColor(kBlack);
gr_resoU_MC->SetLineStyle(2);
gr_resoU_MC->SetLineWidth(2);
gr_resoU_MC->Draw("L");


gr_resoV_MC = new TGraphErrors(nMC,tilt_MC,resoV_MC,etilt_MC,eresoV_MC);
gr_resoV_MC->SetName("gr_resoV_MC");
gr_resoV_MC->SetLineColor(kRed);
gr_resoV_MC->SetLineStyle(5);
gr_resoV_MC->SetLineWidth(2);
gr_resoV_MC->Draw("L");

gr_resoU_TB = new TGraphErrors(nTB,tilt_TB,resoU_TB,etilt_TB,eresoU_TB);
gr_resoU_TB->SetName("gr_resoU_TB");
gr_resoU_TB->SetMarkerColor(kBlack);
gr_resoU_TB->SetMarkerStyle(22);
gr_resoU_TB->Draw("P");


gr_resoV_TB = new TGraphErrors(nTB,tilt_TB,resoV_TB,etilt_TB,eresoV_TB);
gr_resoV_TB->SetName("gr_resoV_TB");
gr_resoV_TB->SetMarkerColor(kRed);
gr_resoV_TB->SetMarkerStyle(21);
gr_resoV_TB->Draw("P");


TLegend* l_reso = new TLegend(0.18,0.6,0.50,0.85);
l_reso->SetFillColor(kWhite); 
l_reso->SetBorderSize(0);
l_reso->AddEntry("gr_resoU_TB","TB resolution u ","p");
l_reso->AddEntry("gr_resoU_MC","MC resolution u ","l");
l_reso->AddEntry("gr_resoV_TB","TB resoultion v ","p");
l_reso->AddEntry("gr_resoV_MC","MC resolution v ","l");
l_reso->Draw();

c3.Modified();
c3.Update();


TCanvas * c9 = new TCanvas("c9","c9",600,400);   
c9->SetLeftMargin(0.13);
//c9->SetRightMargin(0.3);
c9->SetTopMargin(0.1);
c9->SetBottomMargin(0.16);

hmc_charge[9]->Draw();


}
