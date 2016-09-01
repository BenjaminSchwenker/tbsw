{

TFile* f = new TFile("./vxdCheckPlot.root"); // open the file

TCanvas* c1 = new TCanvas();

c1->Divide(3,3); // divide up canvas into 3x3
c1->SetTitle("VXD Checkplots") ;

//f->cd("vxdMap"); // navigate to the vxdMap directory

TTree* vxdMap          = (TTree*)gDirectory->Get("vxdMap") ;            // get the trees
TTree* vxdDetail       = (TTree*)gDirectory->Get("vxdDetailMap") ;
TTree* vxdPoints       = (TTree*)gDirectory->Get("vxdPointsOnSurface") ;
TTree* vxdIntersection = (TTree* )gDirectory->Get("vxdIntersectionPoints") ;
TH2F htemp ;

TCut cutSens = "sensitive==1" ;
TCut cutLadd = "sensitive==0" ;

c1->cd(1);                 // upper left
vxdMap->Draw("x:y" , cutLadd) ; // plot ladder
htemp = (TH2F*)gPad->GetPrimitive("htemp") ;
htemp->SetTitle("LoRes Overview Map of Ladders") ;


c1->cd(2);                 // upper middle
vxdMap->Draw("x:y" , cutSens) ; // plot sensitive
htemp = (TH2F*)gPad->GetPrimitive("htemp") ;
htemp->SetTitle("LoRes Overview Map of Sensitive") ;


c1->cd(4) ;                // middle left
vxdDetail->Draw("x:y", cutLadd) ;
htemp = (TH2F*)gPad->GetPrimitive("htemp") ;
htemp->SetTitle("HiRes Overview Map of Ladders") ;


c1->cd(5) ;               //middle right
vxdDetail->Draw("x:y", cutSens) ;
htemp = (TH2F*)gPad->GetPrimitive("htemp") ;
htemp->SetTitle("HiRes Overview Map of Sensitive") ;


c1->cd(7) ;               //lower left
vxdPoints->Draw("x:y") ;
htemp = (TH2F*)gPad->GetPrimitive("htemp") ;
htemp->SetTitle("Points on Surface") ;


c1->cd(8) ;               //lower right
distance->Draw() ;

c1->cd(3) ;               //upper right
vxdIntersection->Draw("x:y") ;
htemp = (TH2F*)gPad->GetPrimitive("htemp") ;
htemp->SetTitle("Intersection Points") ;

c1->Update(); // update canvas (draw all the pads)

}
