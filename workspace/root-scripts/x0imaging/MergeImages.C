#include <fstream>
#include <algorithm>
using namespace std ;


// This script is used to create a map of a plane in a test beam telescope. The input is a TTree including 
// MSC projected scattering angle distributions and reconstruction errors.
int MergeImages()
{
	gSystem->Load("libProof.so");
	gSystem->Load("libTreePlayer.so");

	gROOT->Reset(); 
	//gROOT->SetStyle("Plain"); 

	// display mode
	gStyle->SetPalette(1);
	gStyle->SetOptStat(11111111);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);

    gROOT->ForceStyle();

    // Read config file
    //------------------
    TEnv mEnv("x0merge.cfg");

	// Number of Input files (image parts)
	const int num_usplits=mEnv.GetValue("usplits", 1);
	const int num_vsplits=mEnv.GetValue("vsplits", 1);

	// TString for the input root file name
	TString filenameB[num_usplits][num_vsplits];

	// Input root files
	TFile *X0file[num_usplits][num_vsplits]; 

	// TString for the input root file name
	TString filenameA;
	filenameA=mEnv.GetValue("x0filename", "X0-merge");

	for(int i=0;i<num_usplits;i++) 
	{
		for(int j=0;j<num_vsplits;j++) 
		{
			filenameB[i][j].Form("-part-%i-%i.root",i+1,j+1); 
			TString filename=filenameA+filenameB[i][j];
			X0file[i][j] = new TFile(filename, "READ");
		}
	}

	// Histogras from the X0 image parts
	TH2F * hX0map_aid[num_usplits][num_vsplits];		   	// X0 images
	TH2F * hX0errmap_aid[num_usplits][num_vsplits];			// statistical error images
	TH2F * hX0relerrmap_aid[num_usplits][num_vsplits];		// relative X0 error images
	TH2F * hchi2map_aid[num_usplits][num_vsplits];			// fit chi2 images
	TH2F * hprobmap1_aid[num_usplits][num_vsplits];			// fit prob images of first scattering angle
	TH2F * hprobmap2_aid[num_usplits][num_vsplits];			// fit prob images of second scattering angle
	TH2F * hprobmapsum_aid[num_usplits][num_vsplits];		// fit prob images of merged scattering angle distribution
	TH2F * hmeanmap1_aid[num_usplits][num_vsplits];			// images of mean value of first scattering angle distribution
	TH2F * hmeanmap2_aid[num_usplits][num_vsplits];			// images of mean value of second scattering angle distribution
	TH2F * huresidualmeanmap_aid[num_usplits][num_vsplits];		// u residual images
	TH2F * hvresidualmeanmap_aid[num_usplits][num_vsplits];		// v residual images
	TH2F * htrackchi2map_aid[num_usplits][num_vsplits];		// track chi2 images
	TH2F * hnummap_aid[num_usplits][num_vsplits];			// hit map
	TH2F * hmommap_aid[num_usplits][num_vsplits];			// Momentum images

	for(int i=0;i<num_usplits;i++) 
	{
		for(int j=0;j<num_vsplits;j++) 
		{
			hX0map_aid[i][j]=(TH2F*)X0file[i][j]->Get("mapping/result/hX0map");
			hX0errmap_aid[i][j]=(TH2F*)X0file[i][j]->Get("mapping/result/hX0errmap");
			hX0relerrmap_aid[i][j]=(TH2F*)X0file[i][j]->Get("mapping/result/hX0relerrmap");
			hchi2map_aid[i][j]=(TH2F*)X0file[i][j]->Get("mapping/result/hchi2map");
			hprobmap1_aid[i][j]=(TH2F*)X0file[i][j]->Get("mapping/result/hprobmap1");
			hprobmap2_aid[i][j]=(TH2F*)X0file[i][j]->Get("mapping/result/hprobmap2");
			hprobmapsum_aid[i][j]=(TH2F*)X0file[i][j]->Get("mapping/result/hprobmapsum");
			hmeanmap1_aid[i][j]=(TH2F*)X0file[i][j]->Get("mapping/result/hmeanmap1");
			hmeanmap2_aid[i][j]=(TH2F*)X0file[i][j]->Get("mapping/result/hmeanmap2");
			huresidualmeanmap_aid[i][j]=(TH2F*)X0file[i][j]->Get("mapping/result/huresidualmeanmap");
			hvresidualmeanmap_aid[i][j]=(TH2F*)X0file[i][j]->Get("mapping/result/hvresidualmeanmap");
			hnummap_aid[i][j]=(TH2F*)X0file[i][j]->Get("mapping/result/hnummap");
			hmommap_aid[i][j]=(TH2F*)X0file[i][j]->Get("mapping/result/hmommap");
		}
	}

	// Results file
	TString resultsfilename=mEnv.GetValue("resultsfilename", "X0-completeimage");
	resultsfilename+=".root";
	TFile *Resultsfile = new TFile(resultsfilename, "RECREATE");	

	// PIXEL Size of the image
	double upixelsize = mEnv.GetValue("upixelsize", 100); // in µm
	double vpixelsize = mEnv.GetValue("vpixelsize", 100); // in µm

	// u minimum and v maximum value (in mm)
	double umin=mEnv.GetValue("umin", -10.0);
	double vmax=mEnv.GetValue("vmax", 5.0);

	// maximal number of pixels per partial image
	int max_u_pixels=mEnv.GetValue("maxupixels", 100);
	int max_v_pixels=mEnv.GetValue("maxvpixels", 50);

	// u and v length of the map (in mm)
	double ulength=num_usplits*upixelsize*max_u_pixels/1000.0;
	double vlength=num_vsplits*vpixelsize*max_v_pixels/1000.0;

	// Print map parameters
	cout<<"Length of total image in u direction: "<<ulength<<" mm"<<endl;
	cout<<"Length of total image in v direction: "<<vlength<<" mm"<<endl;

	// umax and vmin can be computed from the parameters already implemented 
	double umax=umin+ulength;
	double vmin=vmax-vlength;

	// Number of Rows and Columns of the sensor map
	int numcol = num_usplits*max_u_pixels;//max(max_u_pixels,ulength/upixelsize*1000); //1000 because pixelsize is in µm
	int numrow = num_vsplits*max_v_pixels;//max(max_v_pixels,vlength/vpixelsize*1000); //1000 because pixelsize is in µm

	// Print map parameters
	cout<<"Column and Row values of the whole area:"<<endl;
	cout<<"Lowest Column value: "<<0<<endl;
	cout<<"Highest Column value: "<<numcol-1<<endl;
	cout<<"Lowest Row value: "<<0<<endl;
	cout<<"Highest Row value: "<<numrow-1<<endl;

	cout<<endl;

	cout<<"Minimal u value:"<<umin<<" mm"<<endl;
	cout<<"Max. u value:"<<umax<<" mm"<<endl;
	cout<<"Minimal v value:"<<vmin<<" mm"<<endl;
	cout<<"Max. v value:"<<vmax<<" mm"<<endl;


	cout<<"Maximal number of u pixels per partial X0file:"<<max_u_pixels<<endl;
	cout<<"Maximal number of v pixels per partial X0file:"<<max_v_pixels<<endl;

	Resultsfile->cd("");

	// X0 map
	TH2F * hX0map = new TH2F("hX0map","hX0map",numcol,umin,umax,numrow,vmin,vmax);
    hX0map->SetStats(kFALSE);
        hX0map->SetMaximum(8);
        hX0map->SetMinimum(0);
        hX0map->GetXaxis()->SetTitle("u [mm]");
        hX0map->GetYaxis()->SetTitle("v [mm]");
        hX0map->GetZaxis()->SetTitle("X/X0 [%]");
        hX0map->GetZaxis()->SetTitleSize(0.02);
        hX0map->GetZaxis()->SetLabelSize(0.02);

	// X0 statistical error map (absolute value)
	TH2F * hX0errmap = new TH2F("hX0errmap","hX0errmap",numcol,umin,umax,numrow,vmin,vmax);
        hX0errmap->SetStats(kFALSE);
        hX0errmap->SetMaximum(5);
        hX0errmap->SetMinimum(0);
        hX0errmap->GetXaxis()->SetTitle("u [mm]");
        hX0errmap->GetYaxis()->SetTitle("v [mm]");
        hX0errmap->GetZaxis()->SetTitle("X/X0 [%]");
        hX0errmap->GetZaxis()->SetTitleSize(0.02);
        hX0errmap->GetZaxis()->SetLabelSize(0.02);

	// X0 statistical error map (relative value)
	TH2F * hX0relerrmap = new TH2F("hX0relerrmap","hX0relerrmap",numcol,umin,umax,numrow,vmin,vmax);
        hX0relerrmap->SetStats(kFALSE);
        hX0relerrmap->SetMaximum(100);
        hX0relerrmap->SetMinimum(0);
        hX0relerrmap->GetXaxis()->SetTitle("u [mm]");
        hX0relerrmap->GetYaxis()->SetTitle("v [mm]");
        hX0relerrmap->GetZaxis()->SetTitle("rel. error [%]");
        hX0relerrmap->GetZaxis()->SetTitleOffset(1.4);
        hX0relerrmap->GetZaxis()->SetTitleSize(0.02);
        hX0relerrmap->GetZaxis()->SetLabelSize(0.02);

	// Fit Chi2 map
	TH2F * hchi2map = new TH2F("hchi2map","hchi2map",numcol,umin,umax,numrow,vmin,vmax);
        hchi2map->SetStats(kFALSE);
        hchi2map->GetXaxis()->SetTitle("u [mm]");
        hchi2map->GetYaxis()->SetTitle("v [mm]");
        hchi2map->GetZaxis()->SetTitle("chi2");
        hchi2map->GetZaxis()->SetTitleSize(0.02);
        hchi2map->GetZaxis()->SetLabelSize(0.02);

	// Fit Probability distribution for first angle dist
	TH1F * hprobhisto1 = new TH1F("hprobhisto1","hprobhisto1",50,0.0,1.0);
	hprobhisto1->SetStats(kFALSE);
        hprobhisto1->GetXaxis()->SetTitle("u [mm]");
        hprobhisto1->GetYaxis()->SetTitle("fit p1 value");

	// Fit probability map for first angle dist
	TH2F * hprobmap1 = new TH2F("hprobmap1","hprobmap1",numcol,umin,umax,numrow,vmin,vmax);
	hprobmap1->SetStats(kFALSE);
        hprobmap1->GetXaxis()->SetTitle("u [mm]");
        hprobmap1->GetYaxis()->SetTitle("v [mm]");
        hprobmap1->GetZaxis()->SetTitle("fit p1 value");
        hprobmap1->GetZaxis()->SetTitleSize(0.02);
        hprobmap1->GetZaxis()->SetLabelSize(0.02);


	// Fit Probability distribution for second angle dist
	TH1F * hprobhisto2 = new TH1F("hprobhisto2","hprobhisto2",50,0.0,1.0);
	hprobhisto2->SetStats(kFALSE);
        hprobhisto2->GetXaxis()->SetTitle("u [mm]");
        hprobhisto2->GetYaxis()->SetTitle("fit p2 value");

	// Fit probability map for second angle dist
	TH2F * hprobmap2 = new TH2F("hprobmap2","hprobmap2",numcol,umin,umax,numrow,vmin,vmax);
	hprobmap2->SetStats(kFALSE);
        hprobmap2->GetXaxis()->SetTitle("u [mm]");
        hprobmap2->GetYaxis()->SetTitle("v [mm]");
        hprobmap2->GetZaxis()->SetTitle("fit p2 value");
        hprobmap2->GetZaxis()->SetTitleSize(0.02);
        hprobmap2->GetZaxis()->SetLabelSize(0.02);

	// Fit Probability distribution for merged angle dist
	TH1F * hprobhistosum = new TH1F("hprobhistosum","hprobhistosum",50,0.0,1.0);
	hprobhistosum->SetStats(kFALSE);
        hprobhistosum->GetXaxis()->SetTitle("u [mm]");
        hprobhistosum->GetYaxis()->SetTitle("fit psum value");

	// Fit probability map for merged angle dist
	TH2F * hprobmapsum = new TH2F("hprobmapsum","hprobmapsum",numcol,umin,umax,numrow,vmin,vmax);
	hprobmapsum->SetStats(kFALSE);
        hprobmapsum->GetXaxis()->SetTitle("u [mm]");
        hprobmapsum->GetYaxis()->SetTitle("v [mm]");
        hprobmapsum->GetZaxis()->SetTitle("fit psum value");
        hprobmapsum->GetZaxis()->SetTitleSize(0.02);
        hprobmapsum->GetZaxis()->SetLabelSize(0.02);

	// Fit mean value of first distribution
	TH2F * hmeanmap1 = new TH2F("hmeanmap1","hmeanmap1",numcol,umin,umax,numrow,vmin,vmax);
	hmeanmap1->SetStats(kFALSE);
        hmeanmap1->GetXaxis()->SetTitle("u [mm]");
        hmeanmap1->GetYaxis()->SetTitle("v [mm]");
        hmeanmap1->GetZaxis()->SetTitle("theta1 mean value[rad]");
        hmeanmap1->GetZaxis()->SetTitleSize(0.02);
        hmeanmap1->GetZaxis()->SetLabelSize(0.02);

	// Fit mean value of second distribution
	TH2F * hmeanmap2 = new TH2F("hmeanmap2","hmeanmap2",numcol,umin,umax,numrow,vmin,vmax);
	hmeanmap2->SetStats(kFALSE);
        hmeanmap2->GetXaxis()->SetTitle("u [mm]");
        hmeanmap2->GetYaxis()->SetTitle("v [mm]");
        hmeanmap2->GetZaxis()->SetTitle("theta2 mean value[rad]");
        hmeanmap2->GetZaxis()->SetTitleSize(0.02);
        hmeanmap2->GetZaxis()->SetLabelSize(0.02);


	// Fit mean value of u residuals
	TH2F * huresidualmeanmap = new TH2F("huresidualmeanmap","huresidualmeanmap",numcol,umin,umax,numrow,vmin,vmax);
	huresidualmeanmap->SetStats(kFALSE);
        huresidualmeanmap->GetXaxis()->SetTitle("u [mm]");
        huresidualmeanmap->GetYaxis()->SetTitle("v [mm]");
        huresidualmeanmap->GetZaxis()->SetTitle("u residual[µm]");
        huresidualmeanmap->GetZaxis()->SetTitleSize(0.02);
        huresidualmeanmap->GetZaxis()->SetLabelSize(0.02);

	// Fit mean value of v residuals
	TH2F * hvresidualmeanmap = new TH2F("hvresidualmeanmap","hvresidualmeanmap",numcol,umin,umax,numrow,vmin,vmax);
	hvresidualmeanmap->SetStats(kFALSE);
        hvresidualmeanmap->GetXaxis()->SetTitle("u [mm]");
        hvresidualmeanmap->GetYaxis()->SetTitle("v [mm]");
        hvresidualmeanmap->GetZaxis()->SetTitle("v residual[µm]");
        hvresidualmeanmap->GetZaxis()->SetTitleSize(0.02);
        hvresidualmeanmap->GetZaxis()->SetLabelSize(0.02);
	
	// #Tracks map
	TH2F * hnummap = new TH2F("hnummap","hnummap",numcol,umin,umax,numrow,vmin,vmax);
        hnummap->SetStats(kFALSE);
        hnummap->GetXaxis()->SetTitle("u [mm]");
        hnummap->GetYaxis()->SetTitle("v [mm]");
        hnummap->GetZaxis()->SetTitle("number of tracks");
        hnummap->GetZaxis()->SetTitleOffset(1.4);
        hnummap->GetZaxis()->SetTitleSize(0.02);
        hnummap->GetZaxis()->SetLabelSize(0.02);

	// Momentum map
	TH2F * hmommap = new TH2F("hmommap","hmommap",numcol,umin,umax,numrow,vmin,vmax);
        hmommap->SetStats(kFALSE);
        hmommap->GetXaxis()->SetTitle("u [mm]");
        hmommap->GetYaxis()->SetTitle("v [mm]");
        hmommap->GetZaxis()->SetTitle("momentum [GeV/c]");
        hmommap->GetZaxis()->SetTitleOffset(1.4);
        hmommap->GetZaxis()->SetTitleSize(0.02);
        hmommap->GetZaxis()->SetLabelSize(0.02);

	for(int col=0; col<numcol; col++)
	{
		for(int row=0; row<numrow; row++)
		{
			// Copy entries of first images
			hX0map->SetBinContent(col+1,row+1,hX0map_aid[int(col/max_u_pixels)][num_vsplits-(int(row/max_v_pixels)+1)]->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			hX0errmap->SetBinContent(col+1,row+1,hX0errmap_aid[int(col/max_u_pixels)][num_vsplits-(int(row/max_v_pixels)+1)]->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			hX0relerrmap->SetBinContent(col+1,row+1,hX0relerrmap_aid[int(col/max_u_pixels)][num_vsplits-(int(row/max_v_pixels)+1)]->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			hchi2map->SetBinContent(col+1,row+1,hchi2map_aid[int(col/max_u_pixels)][num_vsplits-(int(row/max_v_pixels)+1)]->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			hprobmap1->SetBinContent(col+1,row+1,hprobmap1_aid[int(col/max_u_pixels)][num_vsplits-(int(row/max_v_pixels)+1)]->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			hprobmap2->SetBinContent(col+1,row+1,hprobmap2_aid[int(col/max_u_pixels)][num_vsplits-(int(row/max_v_pixels)+1)]->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			hprobmapsum->SetBinContent(col+1,row+1,hprobmapsum_aid[int(col/max_u_pixels)][num_vsplits-(int(row/max_v_pixels)+1)]->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			hmeanmap1->SetBinContent(col+1,row+1,hmeanmap1_aid[int(col/max_u_pixels)][num_vsplits-(int(row/max_v_pixels)+1)]->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			hmeanmap2->SetBinContent(col+1,row+1,hmeanmap2_aid[int(col/max_u_pixels)][num_vsplits-(int(row/max_v_pixels)+1)]->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			huresidualmeanmap->SetBinContent(col+1,row+1,huresidualmeanmap_aid[int(col/max_u_pixels)][num_vsplits-(int(row/max_v_pixels)+1)]->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			hvresidualmeanmap->SetBinContent(col+1,row+1,hvresidualmeanmap_aid[int(col/max_u_pixels)][num_vsplits-(int(row/max_v_pixels)+1)]->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			hnummap->SetBinContent(col+1,row+1,hnummap_aid[int(col/max_u_pixels)][num_vsplits-(int(row/max_v_pixels)+1)]->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			hmommap->SetBinContent(col+1,row+1,hmommap_aid[int(col/max_u_pixels)][num_vsplits-(int(row/max_v_pixels)+1)]->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
		}
	}
	
	hX0map->Write();
	hX0errmap->Write();
	hX0relerrmap->Write();
	hchi2map->Write();
	hprobmap1->Write();
	hprobmap2->Write();
	hprobmapsum->Write();
	hmeanmap1->Write();
	hmeanmap2->Write();
	huresidualmeanmap->Write();
	hvresidualmeanmap->Write();

	hnummap->Write();
	hmommap->Write();

	Resultsfile->Close();

        return 0;

	
}
