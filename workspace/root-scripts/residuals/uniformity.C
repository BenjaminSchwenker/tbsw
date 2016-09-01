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


int UCELLS = 250; 
int VCELLS = 768; 

int BINSX = 12; 
int BINSY = 12; 


// 
// Signal cuts 


double CLUSTERLOW = 5; 
double CLUSTERHIGH = 60; 
 


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


TTree * hits = (TTree*) ftb->Get("Hit");  

hits->SetBranchAddress("chi2",&_chi2);
hits->SetBranchAddress("hasTrack",&_hasTrack);
hits->SetBranchAddress("u_fit",&_u_fit);
hits->SetBranchAddress("v_fit",&_v_fit);
hits->SetBranchAddress("col_fit",&_col_fit);
hits->SetBranchAddress("row_fit",&_row_fit);
hits->SetBranchAddress("u_pixel",&_u_pixel);
hits->SetBranchAddress("v_pixel",&_v_pixel);
hits->SetBranchAddress("clusterCharge",&_clusterCharge);
hits->SetBranchAddress("seedCharge",&_seedCharge);

// 2d histos

vector<vector<int>> counter; 
counter.resize(BINSX,vector<int>(BINSY, 0));

vector<vector<TH1F*>> histo_cs; 
histo_cs.resize(BINSX,vector<TH1F*>(BINSY, NULL));

for (int u=0; u<BINSX; u++){
  for (int w=0; w<BINSY; w++){
      histo_cs[u][w] = new TH1F("","",100,CLUSTERLOW,CLUSTERHIGH);
  }
}



int binx, biny;

for(int i=0; i< hits->GetEntries(); i++){
    hits->GetEntry(i);
    if ( _hasTrack == 0  ){

      binx = floor( _col_fit / (UCELLS / BINSX) );
      biny = floor( _row_fit / (VCELLS / BINSY) );
      
      if (binx<0) continue;
      if (binx>=BINSX) continue;
      if (biny<0) continue;
      if (biny>=BINSY) continue;

      // Fill maps 
      
      counter[binx][biny]++;
      histo_cs[binx][biny]->Fill(_clusterCharge);
     
    }
}



// create 2d histos

TH2F * h2c = new TH2F("h2c","h2c",BINSX,0,UCELLS,BINSY,0,VCELLS);
TH2F * h2cs = new TH2F("h2cs","h2cs",BINSX,0,UCELLS,BINSY,0,VCELLS);


// fill 2d histos 

for (int i=0; i<BINSX; i++){
  for (int j=0; j<BINSY; j++){
      if (counter[i][j]!=0){
        
        std::cout<< "column " << i << "   row " << j << endl;
        double mean_cs= histo_cs[i][j]->GetMean();
        int total = counter[i][j];
        
        h2cs->SetBinContent(i+1,j+1,mean_cs);
	h2c->SetBinContent(i+1,j+1,total);
      }
  }
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
h2c->GetXaxis()->SetTitle("column");
h2c->GetYaxis()->SetTitle("row");
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
  


TCanvas * c5 = new TCanvas("c5","c5",400,400);
c5->SetLeftMargin(0.2);
c5->SetRightMargin(0.28);
c5->SetTopMargin(0.1);
c5->SetBottomMargin(0.16);
h2cs->GetXaxis()->SetNdivisions(2,kFALSE);
h2cs->GetYaxis()->SetNdivisions(2,kFALSE);
h2cs->GetZaxis()->SetNdivisions(5,kTRUE);
h2cs->SetTitle("");
h2cs->GetXaxis()->SetTitle("column");
h2cs->GetYaxis()->SetTitle("row");
h2cs->GetZaxis()->SetTitle("Mean Cluster Signal [ADU]");
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
