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
// Number of sensor columns and rows 

int COLS = 64; 
int ROWS = 64; 

//
// Pixel pitch 

double PITCHX = 50; // in um 
double PITCHY = 55; // in um 

//
// Map area 

double SIZEX = 1600; // 1600  // in um 
double SIZEY = 3520; // 3520  // in um 

// 
// Number of bins 

int BINSX = 512;   // 800
int BINSY = 1024;  // 1600

// 
// Signal cuts 

double CHI2MAX = 50; 
double CLUSTERLOW = 5; 
double CLUSTERHIGH = 60; 
double SEEDLOW = 5; 
double SEEDHIGH = 60; 


// ----------------------------------------------------------------------


// create 2d histos

TH2F * h2matrix = new TH2F("h2matrix","h2matrix",COLS,0,COLS,ROWS,0,ROWS);
TH2F * h2c = new TH2F("h2c","h2c",BINSX,-SIZEX/2,+SIZEX/2,BINSY,-SIZEY/2,+SIZEY/2);
TH2F * h2cs = new TH2F("h2cs","h2cs",BINSX,-SIZEX/2,+SIZEX/2,BINSY,-SIZEY/2,+SIZEY/2);
TH2F * h2ss = new TH2F("h2ss","h2ss",BINSX,-SIZEX/2,+SIZEX/2,BINSY,-SIZEY/2,+SIZEY/2);


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

vector<vector<int>> matrix; 
matrix.resize(BINSX,vector<int>(BINSY, 0));


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





int binx, biny;
float m_x, m_y;
  
for(int i=0; i< mytuple->GetEntries(); i++){
    mytuple->GetEntry(i);
    if ( _hasTrack == 0 && _chi2 < CHI2MAX  ){
      
      // m_x and m_y are in pixel hit positions from track fit 
      // origin (0,0) is pixel centre, unit is mm 
      
      // Rescale unit from mm to um
      
      m_x = _u_fit *1000; 
      m_y = _v_fit *1000;
 
      //if (_col_fit<0 )    continue;
      //if (_col_fit>=COLS) continue;
      //if (_row_fit<0)     continue;
      //if (_row_fit>=ROWS) continue;

      if (m_x<-SIZEX/2) continue;
      if (m_x>+SIZEX/2) continue;
      if (m_y<-SIZEY/2) continue;
      if (m_y>+SIZEY/2) continue;
        
      // note: this must be c-style index 
      binx = h2cs->GetXaxis()->FindBin(m_x)-1;
      biny = h2cs->GetYaxis()->FindBin(m_y)-1;

      // Fill maps 

      matrix[_col_fit][_row_fit]++;    
      counter[binx][biny]++;
      histo_cs[binx][biny]->Fill(_clusterCharge);
      histo_ss[binx][biny]->Fill(_seedCharge);
      
    }
}


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
  
for (int i=0; i<COLS; i++){
  for (int j=0; j<ROWS; j++){
       
	h2matrix->SetBinContent(i+1,j+1,(matrix[i][j]));
     
  }
}



 
TCanvas * c2 = new TCanvas("c2","c2",400,400);
c2->SetLeftMargin(0.22);
c2->SetRightMargin(0.3);
c2->SetTopMargin(0.1);
c2->SetBottomMargin(0.16);

h2matrix->GetXaxis()->SetNdivisions(2,kFALSE);
h2matrix->GetYaxis()->SetNdivisions(2,kFALSE);
h2matrix->GetZaxis()->SetNdivisions(5,kTRUE);
h2matrix->SetTitle("");
h2matrix->GetXaxis()->SetTitle("column");
h2matrix->GetYaxis()->SetTitle("row");
h2matrix->GetZaxis()->SetTitle("Number of hits");
//h2c->GetZaxis()->SetRangeUser(1600,4500);
h2matrix->GetXaxis()->SetTitleSize(0.06);
h2matrix->GetYaxis()->SetTitleSize(0.06);
h2matrix->GetZaxis()->SetTitleSize(0.06);
h2matrix->GetXaxis()->SetLabelSize(0.06);
h2matrix->GetYaxis()->SetLabelSize(0.06);
h2matrix->GetZaxis()->SetLabelSize(0.06);
h2matrix->GetYaxis()->SetTitleOffset(1.3);
h2matrix->GetZaxis()->SetTitleOffset(1.7);
h2matrix->Draw("colz");


 
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





}
