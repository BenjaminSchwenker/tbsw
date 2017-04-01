{ 

// 
// Analysis cuts 
TCut hitcut = "trackChi2 < 100  && hasTrack == 0 && clusterQuality == 0 && iEvt < 100000 ";
TCut trkcut = "";

// Input file name
TString inputfile = "root-files/H6-S2-Histos-DB-run000065-run65-s8.root";

// Output file name
TString histofile = "ResidualHistos.root";


// ----------------------------------------------------------------------

TFile *ftb = new TFile(inputfile);
TFile *fnew = new TFile(histofile,"RECREATE");

//
// event data 
TTree *ttb = (TTree*) ftb->Get("Hit");
TTree *ttbtrack = (TTree*) ftb->Get("Track");

ttb->Draw("cellV_hit:cellU_hit>>hspot", hitcut,"goff");
hspot->SetName("hhitmap_clusters_cells"); 
hspot->SetTitle(""); 
hspot->GetXaxis()->SetTitle("cell U [cell ID]");
hspot->GetYaxis()->SetTitle("cell V [cell ID]");
hspot->GetZaxis()->SetTitle("number of clusters");
hspot->SetStats(0);

ttbtrack->Draw("cellV_fit:cellU_fit>>htrackspot", trkcut, "goff");
htrackspot->SetName("hhitmap_tracks_uv"); 
htrackspot->SetTitle(""); 
htrackspot->GetXaxis()->SetTitle("cell U [cell ID]");
htrackspot->GetYaxis()->SetTitle("cell V [cell ID]");
htrackspot->GetZaxis()->SetTitle("number of tracks");
htrackspot->SetStats(0);

ttbtrack->Draw("v_fit:u_fit>>hbeamspot", trkcut,"goff");
hbeamspot->SetName("hhitmap_tracks_uv"); 
hbeamspot->SetTitle("");
hbeamspot->GetXaxis()->SetTitle("u [mm]");
hbeamspot->GetYaxis()->SetTitle("v [mm]");
hbeamspot->GetZaxis()->SetTitle("number of tracks");
hbeamspot->SetStats(0);

ttb->Draw("(u_hit - u_fit)*1000 >>hu(101,-75.,+75.)", hitcut,"goff");
hu->SetName("hresidual_u"); 
hu->SetTitle("");
hu->GetXaxis()->SetTitle("residuals u [#mum]");
hu->GetYaxis()->SetTitle("number of hits");
hu->GetYaxis()->SetTitleOffset(1.2);

ttb->Draw("(v_hit - v_fit)*1000 >>hv(81,-75.,+75.)", hitcut,"goff");
hv->SetName("hresidual_v"); 
hv->SetTitle("");
hv->GetXaxis()->SetTitle("residuals v [#mum]");
hv->GetYaxis()->SetTitle("number of hits");
hv->GetYaxis()->SetTitleOffset(1.2);

ttb->Draw("trackChi2 / trackNdof >> hchisqundof(100,0,10)", hitcut,"goff");
hchisqundof->SetName("hchisqundof"); 
hchisqundof->SetTitle("");
hchisqundof->GetXaxis()->SetTitle("#chi^2/ndof");
hchisqundof->GetYaxis()->SetTitle("number of tracks");
hchisqundof->GetYaxis()->SetTitleOffset(1.2);

ttb->Draw("seedCharge >>hCharge(101,0,100)", hitcut,"goff");
hCharge->SetName("hSeedCharge"); 
hCharge->SetTitle("");
hCharge->GetXaxis()->SetTitle("seed charge [ADU]");
hCharge->GetYaxis()->SetTitle("number of hits");
hCharge->GetYaxis()->SetTitleOffset(1.2);

ttb->Draw("clusterCharge >>hCCharge(101,0,100)", hitcut,"goff");
hCCharge->SetName("hClusterCharge"); 
hCCharge->SetTitle("");
hCCharge->GetXaxis()->SetTitle("cluster charge [ADU]");
hCCharge->GetYaxis()->SetTitle("number of hits");
hCCharge->GetYaxis()->SetTitleOffset(1.2);

ttb->Draw("size >>hsize(20,0,20)", hitcut,"goff");
hsize->SetName("hsize"); 
hsize->SetTitle("");
hsize->GetXaxis()->SetTitle("cluster size [pixels]");
hsize->GetYaxis()->SetTitle("number of hits");
hsize->GetYaxis()->SetTitleOffset(1.2);

ttb->Draw("sizeU >>hsizeU(10,0,10)", hitcut,"goff");
hsizeU->SetName("hsizeU"); 
hsizeU->SetTitle("");
hsizeU->GetXaxis()->SetTitle("size [ucells]");
hsizeU->GetYaxis()->SetTitle("number of hits");
hsizeU->GetYaxis()->SetTitleOffset(1.2);

ttb->Draw("sizeV >>hsizeV(10,0,10)", hitcut,"goff");
hsizeV->SetName("hsizeV"); 
hsizeV->SetTitle("");
hsizeV->GetXaxis()->SetTitle("size [vcells]");
hsizeV->GetYaxis()->SetTitle("number of hits");
hsizeV->GetYaxis()->SetTitleOffset(1.2);
hsizeV->SetLineColor(4);
hsizeV->SetLineStyle(2);



fnew->Write();
fnew->Close();
ftb->Close();


}
