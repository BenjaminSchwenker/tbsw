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
TFile *ftb = new TFile("Histos.root");

//
// Sensor paramaters

double PITCHX = 50; // in um 
double PITCHY = 75; // in um 

// 
// Super pixel contains XPIX x YPIX pixel cells  

int XPIX = 2;  
int YPIX = 2;  

// 
// Super pixel is divided into bins along x and y

int BINSX = 10; // nbins in X 30
int BINSY = 20; // nbins in Y 40

// 
// Make projections into superpixel at for this bins 

int PROJ_X = 5; 
int PROJ_Y = 10;

// 
// Signal cuts 

double CHI2MAX = 50; 
double CLUSTERLOW = 5; 
double CLUSTERHIGH = 60; 
double SEEDLOW = 5; 
double SEEDHIGH = 60; 


// ----------------------------------------------------------------------

Double_t _chi2;
Int_t _hasTrack;
Double_t  _u_fit; 
Double_t  _v_fit; 
Double_t _u_pixel;
Double_t _v_pixel;
Int_t _col_fit;
Int_t _row_fit;
Double_t _clusterCharge;
Double_t _seedCharge;


TTree * mytuple = (TTree*) ftb->Get("Hit");  

mytuple->SetBranchAddress("chi2",&_chi2);
mytuple->SetBranchAddress("hasTrack",&_hasTrack);
mytuple->SetBranchAddress("u_fit",&_u_fit);
mytuple->SetBranchAddress("v_fit",&_v_fit);
mytuple->SetBranchAddress("col_fit",&_col_fit);
mytuple->SetBranchAddress("row_fit",&_row_fit);
mytuple->SetBranchAddress("u_pixel",&_u_pixel);
mytuple->SetBranchAddress("v_pixel",&_v_pixel);
mytuple->SetBranchAddress("clusterCharge",&_clusterCharge);
mytuple->SetBranchAddress("seedCharge",&_seedCharge);

// 2d histos

vector<vector<int>> counter; 
counter.resize(BINSX,vector<int>(BINSY, 0));


vector<vector<TH1F*>> histo_cs; 
histo_cs.resize(BINSX,vector<TH1F*>(BINSY, NULL));

vector<vector<TH1F*>> histo_ss; 
histo_ss.resize(BINSX,vector<TH1F*>(BINSY, NULL));

for (int u=0; u<BINSX; u++){
  for (int w=0; w<BINSY; w++){
      histo_cs[u][w] = new TH1F("","",CLUSTERHIGH-CLUSTERLOW,CLUSTERLOW,CLUSTERHIGH);
      histo_ss[u][w] = new TH1F("","",SEEDHIGH-SEEDLOW,SEEDLOW,SEEDHIGH);
  }
}




// 1d histos 

vector<TH1F*>histo_cs_y (BINSY, NULL);
vector<TH1F*>histo_ss_y (BINSY, NULL);



for (int w=0; w<BINSY; w++){
  histo_cs_y[w] = new TH1F("","",CLUSTERHIGH-CLUSTERLOW,CLUSTERLOW,CLUSTERHIGH);
  histo_ss_y[w] = new TH1F("","",SEEDHIGH-SEEDLOW,SEEDLOW,SEEDHIGH);
}


int binx, biny;
float m_x, m_y;
  
for(int i=0; i< mytuple->GetEntries(); i++){
    mytuple->GetEntry(i);
    if ( _hasTrack == 0 && _chi2 < CHI2MAX ){
      
      // m_x and m_y are in pixel hit positions from track fit 
      // origin (0,0) is pixel centre, unit is mm 
      
      // Rescale unit from mm to um. Move origin to lower left corner
      // of pixel.
      
      m_x = (_u_fit - _u_pixel)*1000; // rescale lenght unit: mm -> um 
      m_y = (_v_fit - _v_pixel)*1000;

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
      if (_col_fit < 0) continue; 
      if (_row_fit < 0) continue; 
       
      m_x+= (_col_fit%XPIX) * PITCHX;   
      m_y+= (_row_fit%YPIX) * PITCHY;
        
      binx = floor( m_x / subpitchX );
      biny = floor( m_y / subpitchY );
      
      // Fill maps 
      
      counter[binx][biny]++;
      histo_cs[binx][biny]->Fill(_clusterCharge);
      histo_ss[binx][biny]->Fill(_seedCharge);
      histo_cs_y[biny]->Fill(_clusterCharge);
      histo_ss_y[biny]->Fill(_seedCharge); 
    }
}

// create 2d histos

TH2F * h2c = new TH2F("h2c","h2c",BINSX,0,XPIX*PITCHX,BINSY,0,YPIX*PITCHY);
TH2F * h2cs = new TH2F("h2cs","h2cs",BINSX,0,XPIX*PITCHX,BINSY,0,YPIX*PITCHY);
TH2F * h2ss = new TH2F("h2ss","h2ss",BINSX,0,XPIX*PITCHX,BINSY,0,YPIX*PITCHY);


// fill 2d histos 

for (int i=0; i<BINSX; i++){
  for (int j=0; j<BINSY; j++){
      if (counter[i][j]!=0){
        
        std::cout<< "column " << i << "   row " << j << endl;
        double mean_cs= histo_cs[i][j]->GetMean();
        double mean_ss= histo_ss[i][j]->GetMean();
        
        h2cs->SetBinContent(i+1,j+1,mean_cs);
	h2ss->SetBinContent(i+1,j+1,mean_ss);
	h2c->SetBinContent(i+1,j+1,(counter[i][j]));
      }
  }
}
  

// create 1d histos


TH1F * h1cs_y = new TH1F("h1cs_y","h1cs_y",BINSY,0,YPIX*PITCHY);
TH1F * h1ss_y = new TH1F("h1ss_y","h1ss_y",BINSY,0,YPIX*PITCHY); 


// fill 1d histos

for (int w=0; w<BINSY; w++){
  std::cout<< "row " << w << endl;
  
  double mean_ss= histo_ss_y[w]->GetMean();
  double mean_cs= histo_cs_y[w]->GetMean();
  double mean_ss_error=histo_ss_y[w]->GetMeanError();
  double mean_cs_error=histo_cs_y[w]->GetMeanError();

  h1cs_y->SetBinContent(w+1,mean_cs);
  h1cs_y->SetBinError(w+1,mean_cs_error);
  
  h1ss_y->SetBinContent(w+1,mean_ss);
  h1ss_y->SetBinError(w+1,mean_ss_error);  
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
h2c->GetXaxis()->SetTitle("u_{p} [#mum]");
h2c->GetYaxis()->SetTitle("v_{p} [#mum]");
h2c->GetZaxis()->SetTitle("Number of tracks");
//h2c->GetZaxis()->SetRangeUser(1600,4500);
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
c4->SetLeftMargin(0.22);
c4->SetRightMargin(0.3);
c4->SetTopMargin(0.1);
c4->SetBottomMargin(0.16);
h2ss->GetXaxis()->SetNdivisions(2,kFALSE);
h2ss->GetYaxis()->SetNdivisions(2,kFALSE);
h2ss->GetZaxis()->SetNdivisions(5,kTRUE);
h2ss->SetTitle("");
h2ss->GetXaxis()->SetTitle("u_{p} [#mum]");
h2ss->GetYaxis()->SetTitle("v_{p} [#mum]");
h2ss->GetZaxis()->SetTitle("Mean Seed Signal [LSB]");
//h2ss->GetZaxis()->SetRangeUser(1600,4500);
h2ss->GetXaxis()->SetTitleSize(0.06);
h2ss->GetYaxis()->SetTitleSize(0.06);
h2ss->GetZaxis()->SetTitleSize(0.06);
h2ss->GetXaxis()->SetLabelSize(0.06);
h2ss->GetYaxis()->SetLabelSize(0.06);
h2ss->GetZaxis()->SetLabelSize(0.06);
h2ss->GetYaxis()->SetTitleOffset(1.3);
h2ss->GetZaxis()->SetTitleOffset(1.7);
h2ss->Draw("colz");

TCanvas * c5 = new TCanvas("c5","c5",400,400);
c5->SetLeftMargin(0.2);
c5->SetRightMargin(0.28);
c5->SetTopMargin(0.1);
c5->SetBottomMargin(0.16);
h2cs->GetXaxis()->SetNdivisions(2,kFALSE);
h2cs->GetYaxis()->SetNdivisions(2,kFALSE);
h2cs->GetZaxis()->SetNdivisions(5,kTRUE);
h2cs->SetTitle("");
h2cs->GetXaxis()->SetTitle("u_{p} [#mum]");
h2cs->GetYaxis()->SetTitle("v_{p} [#mum]");
h2cs->GetZaxis()->SetTitle("Mean Cluster Signal [LSB]");
//h2cs->GetZaxis()->SetRangeUser(3500,4500);
h2cs->GetXaxis()->SetTitleSize(0.06);
h2cs->GetYaxis()->SetTitleSize(0.06);
h2cs->GetZaxis()->SetTitleSize(0.06);
h2cs->GetXaxis()->SetLabelSize(0.06);
h2cs->GetYaxis()->SetLabelSize(0.06);
h2cs->GetZaxis()->SetLabelSize(0.06);
h2cs->GetYaxis()->SetTitleOffset(1.3);
h2cs->GetZaxis()->SetTitleOffset(1.7);
h2cs->Draw("colz");




TCanvas * c100 = new TCanvas("c100","c100",400,400);
c100->SetLeftMargin(0.2);
c100->SetRightMargin(0.1);
c100->SetTopMargin(0.1);
c100->SetBottomMargin(0.16);
c100->SetGrid(0,1);

TH1D *h2ss_x = h2ss->ProjectionX("h2ss_x", PROJ_X, PROJ_X);


h2ss_x->SetFillColor(kGray);

h2ss_x->GetXaxis()->SetNdivisions(4,kFALSE);
h2ss_x->SetTitle("");
h2ss_x->GetXaxis()->SetTitle("u_{p} [#mum]");
h2ss_x->GetYaxis()->SetTitle("Mean Seed Signal [LSB]");
h2ss_x->GetXaxis()->SetTitleSize(0.06);
h2ss_x->GetYaxis()->SetTitleSize(0.06);
h2ss_x->GetYaxis()->SetTitleOffset(1.5);
h2ss_x->GetXaxis()->SetLabelSize(0.06);
h2ss_x->GetYaxis()->SetLabelSize(0.06);
h2ss_x->Draw();

TCanvas * c110 = new TCanvas("c110","c110",400,400);
c110->SetLeftMargin(0.2);
c110->SetRightMargin(0.1);
c110->SetTopMargin(0.1);
c110->SetBottomMargin(0.16);
c110->SetGrid(0,1);

TH1D *h2ss_y = h2ss->ProjectionY("h2ss_y", PROJ_Y, PROJ_Y);


h2ss_y->SetFillColor(kGray);

h2ss_y->GetXaxis()->SetNdivisions(4,kFALSE);
h2ss_y->SetTitle("");
h2ss_y->GetXaxis()->SetTitle("v_{p} [#mum]");
h2ss_y->GetYaxis()->SetTitle("Mean Seed Signal [LSB]");
h2ss_y->GetXaxis()->SetTitleSize(0.06);
h2ss_y->GetYaxis()->SetTitleSize(0.06);
h2ss_y->GetYaxis()->SetTitleOffset(1.5);
h2ss_y->GetXaxis()->SetLabelSize(0.06);
h2ss_y->GetYaxis()->SetLabelSize(0.06);
h2ss_y->Draw();




}
