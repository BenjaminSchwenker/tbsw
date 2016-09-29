{ 

gROOT->Reset(); 
gROOT->SetStyle("Plain"); 

gStyle->SetOptFit();
gStyle->SetPalette(1);
gStyle->SetCanvasColor(0);
gStyle->SetTitleFillColor(0);
gStyle->SetPaintTextFormat("1f");
gStyle->SetOptStat(0);

// Event data 
TFile *ftb = new TFile("root-files/L3-Histos-run000096_svdpxdteldigits_complete.root");
       
//
// Sensor paramaters

double PITCHX = 50; // in um 
double PITCHY = 160; // in um 

// 
// Super pixel contains XPIX x YPIX pixel cells  

int XPIX = 1;  
int YPIX = 1;  

// 
// Super pixel is divided into bins along x and y

int BINSX = 10; // nbins in X 30
int BINSY = 10; // nbins in Y 40

// 
// Signal cuts 

double CHI2MAX = 50; 


// ----------------------------------------------------------------------

Double_t _chi2;
Int_t _hasTrack;
Double_t  _u_fit; 
Double_t  _v_fit; 
Double_t _cellCenterU;
Double_t _cellCenterV;
Int_t _cellU_fit;
Int_t _cellV_fit;
Int_t _size; 
Double_t _clusterChargeU;
Double_t _clusterChargeV;


TTree * mytuple = (TTree*) ftb->Get("Hit");  

mytuple->SetBranchAddress("chi2",&_chi2);
mytuple->SetBranchAddress("hasTrack",&_hasTrack);
mytuple->SetBranchAddress("u_fit",&_u_fit);
mytuple->SetBranchAddress("v_fit",&_v_fit);
mytuple->SetBranchAddress("cellU_fit",&_cellU_fit);
mytuple->SetBranchAddress("cellV_fit",&_cellV_fit);
mytuple->SetBranchAddress("cellCenterU",&_cellCenterU);
mytuple->SetBranchAddress("cellCenterV",&_cellCenterV);
mytuple->SetBranchAddress("size",&_size);
mytuple->SetBranchAddress("clusterChargeU",&_clusterChargeU);
mytuple->SetBranchAddress("clusterChargeV",&_clusterChargeV);



// 2d histos


vector<vector<int>> counter; 
counter.resize(BINSX,vector<int>(BINSY, 0));

vector<vector<TH1F*>> histo_cc; 
histo_cc.resize(BINSX,vector<TH1F*>(BINSY, NULL));

vector<vector<TH1F*>> histo_rr; 
histo_rr.resize(BINSX,vector<TH1F*>(BINSY, NULL));

vector<vector<TH1F*>> histo_px; 
histo_px.resize(BINSX,vector<TH1F*>(BINSY, NULL));


for (int u=0; u<BINSX; u++){
  for (int w=0; w<BINSY; w++){
      histo_cc[u][w] = new TH1F("","",200,0,200);
      histo_rr[u][w] = new TH1F("","",200,0,200);
      histo_px[u][w] = new TH1F("","",10,0,10);
  }
}


// 1d histos

vector<TH1F*>histo_c (BINSX, NULL);



for (int w=0; w<BINSX; w++){
  histo_c[w] = new TH1F("","",200,0,200);
}

vector<TH1F*>histo_r (BINSY, NULL);




for (int w=0; w<BINSY; w++){
  histo_r[w] = new TH1F("","",200,0,200);
}




int binx, biny;
float m_x, m_y;
  
for(int i=0; i< mytuple->GetEntries(); i++){
    mytuple->GetEntry(i);
    if (  _hasTrack == 0 && _chi2 < CHI2MAX ){
      
      // m_x and m_y are in pixel hit positions from track fit 
      // origin (0,0) is pixel centre, unit is mm 
      
      // Rescale unit from mm to um. Move origin to lower left corner
      // of pixel.
      
      m_x = (_u_fit - _cellCenterU)*1000; // rescale lenght unit: mm -> um 
      m_y = (_v_fit - _cellCenterV)*1000;
     
      m_x += PITCHX/2; // move origin to corner of pixel 
      m_y += PITCHY/2;
      
      if (m_x<0) continue;
      if (m_x>=PITCHX) continue;
      if (m_y<0) continue;
      if (m_y>=PITCHY) continue;
      
      // Pitch of inpixel bins 
      double subpitchX = (XPIX*PITCHX)/BINSX; 
      double subpitchY = (YPIX*PITCHY)/BINSY; 
      
      // Calculate inpixel bins  
      if (_cellU_fit < 0) continue; 
      if (_cellV_fit < 0) continue; 
       
      m_x+= (_cellU_fit%XPIX) * PITCHX;   
      m_y+= (_cellV_fit%YPIX) * PITCHY;
      
      binx = floor( m_x / subpitchX );
      biny = floor( m_y / subpitchY );
      
      

      // Fill maps 
      
      counter[binx][biny]++;
      
      histo_px[binx][biny]->Fill(_size);
      histo_cc[binx][biny]->Fill(_clusterChargeU);
      histo_rr[binx][biny]->Fill(_clusterChargeV);

      

      // Fill 1d histos
      histo_c[binx]->Fill(_clusterChargeU);
      histo_r[biny]->Fill(_clusterChargeV);
      
    }
}






  
TH2F * h2c = new TH2F("h2c","h2c",BINSX,0,XPIX*PITCHX,BINSY,0,YPIX*PITCHY);
TH2F * h2px = new TH2F("h2px","h2px",BINSX,0,XPIX*PITCHX,BINSY,0,YPIX*PITCHY);
TH2F * h2ccm = new TH2F("h2ccm","h2ccm",BINSX,0,XPIX*PITCHX,BINSY,0,YPIX*PITCHY);
TH2F * h2rrm = new TH2F("h2rrm","h2rrm",BINSX,0,XPIX*PITCHX,BINSY,0,YPIX*PITCHY);



for (int i=0; i<BINSX; i++){
  for (int j=0; j<BINSY; j++){
      if (counter[i][j]!=0){
        
        std::cout<< "column " << i << "   row " << j << endl;
       
        double mean_px= histo_px[i][j]->GetMean();
        double mean_cc= histo_cc[i][j]->GetMean();
        double mean_rr= histo_rr[i][j]->GetMean();

         
        h2px->SetBinContent(i+1,j+1,mean_px);
        h2ccm->SetBinContent(i+1,j+1,mean_cc);
        h2rrm->SetBinContent(i+1,j+1,mean_rr);
	h2c->SetBinContent(i+1,j+1,(counter[i][j]));
      }
  }
}


// create 1d histos

TH1F * hmean_c = new TH1F("hmean_c","hmean_c",BINSX,0,XPIX*PITCHX);
TH1F * hmean_r = new TH1F("hmean_r","hmean_r",BINSY,0,YPIX*PITCHY);



for (int i=0; i<BINSX; i++){
  double m = histo_c[i]->GetMean();
  double m_err = histo_c[i]->GetMeanError(); 
  
  hmean_c->SetBinContent(i+1,m);
  hmean_c->SetBinError(i+1,m_err);

 
  
}

for (int i=0; i<BINSY; i++){
  double m = histo_r[i]->GetMean();
  double m_err = histo_r[i]->GetMeanError(); 
  
  hmean_r->SetBinContent(i+1,m);
  hmean_r->SetBinError(i+1,m_err);

  
}

 
TCanvas * c3 = new TCanvas("c3","c3",400,400);
  
c3->SetLeftMargin(0.22);
c3->SetRightMargin(0.3);
c3->SetTopMargin(0.1);
c3->SetBottomMargin(0.16);
h2c->GetXaxis()->SetNdivisions(2,kFALSE);
h2c->GetYaxis()->SetNdivisions(2,kFALSE);
h2c->GetZaxis()->SetNdivisions(5,kTRUE);
h2c->SetTitle("");
h2c->GetXaxis()->SetTitle("u [#mum]");
h2c->GetYaxis()->SetTitle("v [#mum]");
h2c->GetZaxis()->SetTitle("Number Tracks");
h2c->GetXaxis()->SetTitleSize(0.06);
h2c->GetYaxis()->SetTitleSize(0.06);
h2c->GetZaxis()->SetTitleSize(0.06);
h2c->GetXaxis()->SetLabelSize(0.06);
h2c->GetYaxis()->SetLabelSize(0.06);
h2c->GetZaxis()->SetLabelSize(0.06);
h2c->GetYaxis()->SetTitleOffset(1.3);
h2c->GetZaxis()->SetTitleOffset(1.7);
h2c->Draw("colz");


TCanvas * c4 = new TCanvas("c4","c4",400,400);
c4->SetLeftMargin(0.2);
c4->SetRightMargin(0.28);
c4->SetTopMargin(0.1);
c4->SetBottomMargin(0.16);
h2px->GetXaxis()->SetNdivisions(2,kFALSE);
h2px->GetYaxis()->SetNdivisions(2,kFALSE);
h2px->GetZaxis()->SetNdivisions(5,kTRUE);
h2px->SetTitle("");
h2px->GetXaxis()->SetTitle("u [#mum]");
h2px->GetYaxis()->SetTitle("v [#mum]");
h2px->GetZaxis()->SetTitle("Mean cluster size [pixels]");
h2px->GetZaxis()->SetRangeUser(1,5);
h2px->GetXaxis()->SetTitleSize(0.06);
h2px->GetYaxis()->SetTitleSize(0.06);
h2px->GetZaxis()->SetTitleSize(0.06);
h2px->GetXaxis()->SetLabelSize(0.06);
h2px->GetYaxis()->SetLabelSize(0.06);
h2px->GetZaxis()->SetLabelSize(0.06);
h2px->GetYaxis()->SetTitleOffset(1.3);
h2px->GetZaxis()->SetTitleOffset(1.4);
h2px->Draw("colz");

TCanvas * c5 = new TCanvas("c5","c5",400,400);
c5->SetLeftMargin(0.2);
c5->SetRightMargin(0.28);
c5->SetTopMargin(0.1);
c5->SetBottomMargin(0.16);
h2ccm->GetXaxis()->SetNdivisions(2,kFALSE);
h2ccm->GetYaxis()->SetNdivisions(2,kFALSE);
h2ccm->GetZaxis()->SetNdivisions(5,kTRUE);
h2ccm->SetTitle("");
h2ccm->GetXaxis()->SetTitle("u [#mum]");
h2ccm->GetYaxis()->SetTitle("v [#mum]");
h2ccm->GetZaxis()->SetTitle("Mean u charge");
//h2ccm->GetZaxis()->SetRangeUser(1,3);
h2ccm->GetXaxis()->SetTitleSize(0.06);
h2ccm->GetYaxis()->SetTitleSize(0.06);
h2ccm->GetZaxis()->SetTitleSize(0.06);
h2ccm->GetXaxis()->SetLabelSize(0.06);
h2ccm->GetYaxis()->SetLabelSize(0.06);
h2ccm->GetZaxis()->SetLabelSize(0.06);
h2ccm->GetYaxis()->SetTitleOffset(1.3);
h2ccm->GetZaxis()->SetTitleOffset(1.4);
h2ccm->Draw("colz");

TCanvas * c6 = new TCanvas("c6","c6",400,400);
c6->SetLeftMargin(0.2);
c6->SetRightMargin(0.28);
c6->SetTopMargin(0.1);
c6->SetBottomMargin(0.16);
h2rrm->GetXaxis()->SetNdivisions(2,kFALSE);
h2rrm->GetYaxis()->SetNdivisions(2,kFALSE);
h2rrm->GetZaxis()->SetNdivisions(5,kTRUE);
h2rrm->SetTitle("");
h2rrm->GetXaxis()->SetTitle("u [#mum]");
h2rrm->GetYaxis()->SetTitle("v [#mum]");
h2rrm->GetZaxis()->SetTitle("Mean v charge");
//h2rrm->GetZaxis()->SetRangeUser(1,3);
h2rrm->GetXaxis()->SetTitleSize(0.06);
h2rrm->GetYaxis()->SetTitleSize(0.06);
h2rrm->GetZaxis()->SetTitleSize(0.06);
h2rrm->GetXaxis()->SetLabelSize(0.06);
h2rrm->GetYaxis()->SetLabelSize(0.06);
h2rrm->GetZaxis()->SetLabelSize(0.06);
h2rrm->GetYaxis()->SetTitleOffset(1.3);
h2rrm->GetZaxis()->SetTitleOffset(1.4);
h2rrm->Draw("colz");



TCanvas * c9 = new TCanvas("c9","c9",400,400);
c9->SetLeftMargin(0.2);
//c9->SetRightMargin(0.28);
c9->SetTopMargin(0.1);
c9->SetBottomMargin(0.16);
//hmean_c->GetXaxis()->SetNdivisions(2,kFALSE);
hmean_c->SetTitle("");
hmean_c->GetXaxis()->SetTitle("u [#mum]");
hmean_c->GetYaxis()->SetTitle("Mean u charge");
hmean_c->GetXaxis()->SetTitleSize(0.06);
hmean_c->GetYaxis()->SetTitleSize(0.06);
hmean_c->GetXaxis()->SetLabelSize(0.06);
hmean_c->GetYaxis()->SetLabelSize(0.06);
hmean_c->GetYaxis()->SetTitleOffset(1.35);
hmean_c->Draw();



TCanvas * c11 = new TCanvas("c11","c11",400,400);
c11->SetLeftMargin(0.2);
//c11->SetRightMargin(0.28);
c11->SetTopMargin(0.1);
c11->SetBottomMargin(0.16);
//hmean_r->GetXaxis()->SetNdivisions(2,kFALSE);
hmean_r->SetTitle("");
hmean_r->GetXaxis()->SetTitle("v^{x} [#mum]");
hmean_r->GetYaxis()->SetTitle("Mean v charge");
hmean_r->GetXaxis()->SetTitleSize(0.06);
hmean_r->GetYaxis()->SetTitleSize(0.06);
hmean_r->GetXaxis()->SetLabelSize(0.06);
hmean_r->GetYaxis()->SetLabelSize(0.06);
hmean_r->GetYaxis()->SetTitleOffset(1.35);
hmean_r->Draw();




}
