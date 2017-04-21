{ 

gROOT->Reset(); 
gROOT->SetStyle("Plain"); 

gStyle->SetOptFit();
gStyle->SetPalette(1);
gStyle->SetCanvasColor(0);
gStyle->SetTitleFillColor(0);
gStyle->SetPaintTextFormat("1f");
gStyle->SetOptStat(0);


// Simulate energy loss by bremsstrahlung (Bethe Heitler theory)

double Ei = 4; 

double t = 0.1;  // thickness in units X0 
  
t = 30000.0/305000;
      
double c = t/TMath::Log(2);

TH1D * h1e = new TH1D("h1e","h1e",200,0,1.4*Ei);


for (int i = 0; i< 100000; i++) {

  double x = gRandom->Rndm(1);

  double u = ROOT::Math::gamma_quantile(x,c,1.0);
            
  double z = TMath::Exp(-u);    

  double Eism = gRandom->Gaus(Ei,0.2);

  double Ef = z*Eism; 


  h1e->Fill(Ef);
}
 
TCanvas * c3 = new TCanvas("c3","c3",400,400);
c3->SetLeftMargin(0.22);
c3->SetRightMargin(0.3);
c3->SetTopMargin(0.1);
c3->SetBottomMargin(0.16);


h1e->Draw();


cout << "mean is " << h1e->GetMean() << endl;
cout << "rms is " << h1e->GetRMS() << endl;



}
