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
	TString filenameB;

	// Input root files
	TFile* X0file; 

	// TString for the input root file name
	TString filenameA;
	filenameA=mEnv.GetValue("x0filename", "X0-merge");

	// Aid Histograms from the X0 image parts
	TH2F * x0_image_aid;		   		            // X0 images
	TH2F * x0err_image_aid;				            // statistical error images
	TH2F * x0relerr_image_aid;			            // relative X0 error images
	TH2F * fit1chi2ndof_image_aid;				    // fit chi2 images
	TH2F * fit2chi2ndof_image_aid;				    // fit chi2 images
	TH2F * fitsumchi2ndof_image_aid;			    // fit chi2 images
	TH2F * fit1prob_image_aid;				        // fit prob images of first scattering angle
	TH2F * fit2prob_image_aid;				        // fit prob images of second scattering angle
	TH2F * fitsumprob_image_aid;			        // fit prob images of merged scattering angle distribution
	TH2F * theta1mean_image_aid;				    // images of mean value of first scattering angle distribution
	TH2F * theta2mean_image_aid;				    // images of mean value of second scattering angle distribution
	TH2F * correctedtheta1mean_image_aid;	        // images of mean value of first scattering angle distribution
	TH2F * correctedtheta2mean_image_aid;	        // images of mean value of second scattering angle distribution
	TH2F * uresidualmean_image_aid;		            // u residual images
	TH2F * vresidualmean_image_aid;		            // v residual images
	TH2F * htrackchi2map_aid;			            // track chi2 images
	TH2F * beamspot_aid;				            // hit map
	TH2F * BE_image_aid;				            // Momentum images

	TH2F * vertex_w_mean_image_aid;				    // Vertex w mean images
	TH2F * vertex_w_rms_image_aid;				    // Vertex w RMS images
	TH2F * vertex_chi2_image_aid;				    // Vertex chi2 images
	TH2F * vertex_multiplicity_image_aid;			// Vertex multiplicity images

	TH2F * res_u_vtx_trk_image_aid;					// Vertex trk u mean residual images
	TH2F * res_v_vtx_trk_image_aid;					// Vertex trk v mean residual images
	TH2F * res_u_rms_vtx_trk_image_aid;					// Vertex trk u residual rms images
	TH2F * res_v_rms_vtx_trk_image_aid;					// Vertex trk v residual rms images

	TH1F * fit1prob_histo_aid;	    	            // fit prob histo of first scattering angle
	TH1F * fit2prob_histo_aid;			            // fit prob histo of second scattering angle
	TH1F * fitsumprob_histo_aid;			        // fit prob histo of merged scattering angle distribution

	TH2F * hscatt_theta1_vs_resu_aid;	            // theta vs residual scatter plot
	TH2F * hscatt_theta1_vs_resv_aid;	            // theta vs residual scatter plot
	TH2F * hscatt_theta2_vs_resu_aid;	            // theta vs residual scatter plot
	TH2F * hscatt_theta2_vs_resv_aid;	            // theta vs residual scatter plot


    // Fit Probability distribution for first angle dist
	TH1F * fit1prob_histo = new TH1F("fit1prob_histo","fit1prob_histo",50,0.0,1.0);
	fit1prob_histo->SetStats(kFALSE);
    fit1prob_histo->GetXaxis()->SetTitle("fit1 p value");
    fit1prob_histo->GetYaxis()->SetTitle("number of fits");

	// Fit Probability distribution for second angle dist
	TH1F * fit2prob_histo = new TH1F("fit2prob_histo","fit2prob_histo",50,0.0,1.0);
	fit2prob_histo->SetStats(kFALSE);
    fit2prob_histo->GetXaxis()->SetTitle("fit2 p value");
    fit2prob_histo->GetYaxis()->SetTitle("number of fits");

	// Fit Probability distribution for merged angle dist
	TH1F * fitsumprob_histo = new TH1F("fitsumprob_histo","fitsumprob_histo",50,0.0,1.0);
	fitsumprob_histo->SetStats(kFALSE);
    fitsumprob_histo->GetXaxis()->SetTitle("fitsum p value");
    fitsumprob_histo->GetYaxis()->SetTitle("number of fits");

	// Scatter theta1 vs residual u
	TH2F * hscatt_theta1_vs_resu = new TH2F("hscatt_theta1_vs_resu","hscatt_theta1_vs_resu",100,-0.1,0.1,150,-0.005,0.005);
    hscatt_theta1_vs_resu->SetStats(kFALSE);
    hscatt_theta1_vs_resu->GetXaxis()->SetTitle("u residual[mm]");
    hscatt_theta1_vs_resu->GetYaxis()->SetTitle("theta1[rad]");
    hscatt_theta1_vs_resu->GetZaxis()->SetTitle("Number of tracks");
    hscatt_theta1_vs_resu->GetZaxis()->SetTitleOffset(1.4);
    hscatt_theta1_vs_resu->GetZaxis()->SetTitleSize(0.02);
    hscatt_theta1_vs_resu->GetZaxis()->SetLabelSize(0.02);

	// Scatter theta2 vs residual u
	TH2F * hscatt_theta2_vs_resu = new TH2F("hscatt_theta2_vs_resu","hscatt_theta2_vs_resu",100,-0.1,0.1,150,-0.005,0.005);
    hscatt_theta2_vs_resu->SetStats(kFALSE);
    hscatt_theta2_vs_resu->GetXaxis()->SetTitle("u residual[mm]");
    hscatt_theta2_vs_resu->GetYaxis()->SetTitle("theta2[rad]");
    hscatt_theta2_vs_resu->GetZaxis()->SetTitle("Number of tracks");
    hscatt_theta2_vs_resu->GetZaxis()->SetTitleOffset(1.4);
    hscatt_theta2_vs_resu->GetZaxis()->SetTitleSize(0.02);
    hscatt_theta2_vs_resu->GetZaxis()->SetLabelSize(0.02);

	// Scatter theta1 vs residual v
	TH2F * hscatt_theta1_vs_resv = new TH2F("hscatt_theta1_vs_resv","hscatt_theta1_vs_resv",100,-0.1,0.1,150,-0.005,0.005);
    hscatt_theta1_vs_resv->SetStats(kFALSE);
    hscatt_theta1_vs_resv->GetXaxis()->SetTitle("u residual[mm]");
    hscatt_theta1_vs_resv->GetYaxis()->SetTitle("theta1[rad]");
    hscatt_theta1_vs_resv->GetZaxis()->SetTitle("Number of tracks");
    hscatt_theta1_vs_resv->GetZaxis()->SetTitleOffset(1.4);
    hscatt_theta1_vs_resv->GetZaxis()->SetTitleSize(0.02);
    hscatt_theta1_vs_resv->GetZaxis()->SetLabelSize(0.02);

	// Scatter theta2 vs residual v
	TH2F * hscatt_theta2_vs_resv = new TH2F("hscatt_theta2_vs_resv","hscatt_theta2_vs_resv",100,-0.1,0.1,150,-0.005,0.005);
    hscatt_theta2_vs_resv->SetStats(kFALSE);
    hscatt_theta2_vs_resv->GetXaxis()->SetTitle("u residual[mm]");
    hscatt_theta2_vs_resv->GetYaxis()->SetTitle("theta2[rad]");
    hscatt_theta2_vs_resv->GetZaxis()->SetTitle("Number of tracks");
    hscatt_theta2_vs_resv->GetZaxis()->SetTitleOffset(1.4);
    hscatt_theta2_vs_resv->GetZaxis()->SetTitleSize(0.02);
    hscatt_theta2_vs_resv->GetZaxis()->SetLabelSize(0.02);



	// Get the root file of the first partial image
	// and copy the scatter plots
	filenameB.Form("-part-1-1.root"); 

	TString filename=filenameA+filenameB;
    X0file = new TFile(filename, "READ");

	hscatt_theta1_vs_resu_aid=(TH2F*)X0file->Get("mapping/result/hscatt_theta1_vs_resu");
	hscatt_theta1_vs_resv_aid=(TH2F*)X0file->Get("mapping/result/hscatt_theta1_vs_resv");
	hscatt_theta2_vs_resu_aid=(TH2F*)X0file->Get("mapping/result/hscatt_theta2_vs_resu");
	hscatt_theta2_vs_resv_aid=(TH2F*)X0file->Get("mapping/result/hscatt_theta2_vs_resv");

	hscatt_theta1_vs_resu=hscatt_theta1_vs_resu_aid;
	hscatt_theta1_vs_resv=hscatt_theta1_vs_resv_aid;
	hscatt_theta2_vs_resu=hscatt_theta2_vs_resu_aid;
	hscatt_theta2_vs_resv=hscatt_theta2_vs_resv_aid;

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
	TH2F * x0_image = new TH2F("x0_image","x0_image",numcol,umin,umax,numrow,vmin,vmax);
    x0_image->SetStats(kFALSE);
    x0_image->SetMaximum(10);
    x0_image->SetMinimum(0);
    x0_image->GetXaxis()->SetTitle("u [mm]");
    x0_image->GetYaxis()->SetTitle("v [mm]");
    x0_image->GetZaxis()->SetTitle("X/X0 [%]");
    x0_image->GetZaxis()->SetTitleSize(0.02);
    x0_image->GetZaxis()->SetLabelSize(0.02);

	// X0 statistical error map (absolute value)
	TH2F * x0err_image = new TH2F("x0err_image","x0err_image",numcol,umin,umax,numrow,vmin,vmax);
    x0err_image->SetStats(kFALSE);
    x0err_image->SetMaximum(10);
    x0err_image->SetMinimum(0);
    x0err_image->GetXaxis()->SetTitle("u [mm]");
    x0err_image->GetYaxis()->SetTitle("v [mm]");
    x0err_image->GetZaxis()->SetTitle("X/X0 [%]");
    x0err_image->GetZaxis()->SetTitleSize(0.02);
    x0err_image->GetZaxis()->SetLabelSize(0.02);

	// X0 statistical error map (relative value)
	TH2F * x0relerr_image = new TH2F("x0relerr_image","x0relerr_image",numcol,umin,umax,numrow,vmin,vmax);
    x0relerr_image->SetStats(kFALSE);
    x0relerr_image->SetMaximum(100);
    x0relerr_image->SetMinimum(0);
    x0relerr_image->GetXaxis()->SetTitle("u [mm]");
    x0relerr_image->GetYaxis()->SetTitle("v [mm]");
    x0relerr_image->GetZaxis()->SetTitle("rel. error [%]");
    x0relerr_image->GetZaxis()->SetTitleOffset(1.4);
    x0relerr_image->GetZaxis()->SetTitleSize(0.02);
    x0relerr_image->GetZaxis()->SetLabelSize(0.02);

	// Fit Chi2 image of the merged angle distribution
	TH2F * fitsumchi2ndof_image = new TH2F("fitsumchi2ndof_image","fitsumchi2ndof_image",numcol,umin,umax,numrow,vmin,vmax);
    fitsumchi2ndof_image->SetStats(kFALSE);
    fitsumchi2ndof_image->GetXaxis()->SetTitle("u [mm]");
    fitsumchi2ndof_image->GetYaxis()->SetTitle("v [mm]");
    fitsumchi2ndof_image->GetZaxis()->SetTitle("fitsum chi2");
    fitsumchi2ndof_image->GetZaxis()->SetTitleSize(0.02);
    fitsumchi2ndof_image->GetZaxis()->SetLabelSize(0.02);

	// Fit Chi2 image of the merged angle distribution
	TH2F * fit1chi2ndof_image = new TH2F("fit1chi2ndof_image","fit1chi2ndof_image",numcol,umin,umax,numrow,vmin,vmax);
    fit1chi2ndof_image->SetStats(kFALSE);
    fit1chi2ndof_image->GetXaxis()->SetTitle("u [mm]");
    fit1chi2ndof_image->GetYaxis()->SetTitle("v [mm]");
    fit1chi2ndof_image->GetZaxis()->SetTitle("fit1 chi2");
    fit1chi2ndof_image->GetZaxis()->SetTitleSize(0.02);
    fit1chi2ndof_image->GetZaxis()->SetLabelSize(0.02);

	// Fit Chi2 image of the merged angle distribution
	TH2F * fit2chi2ndof_image = new TH2F("fit2chi2ndof_image","fit2chi2ndof_image",numcol,umin,umax,numrow,vmin,vmax);
    fit2chi2ndof_image->SetStats(kFALSE);
    fit2chi2ndof_image->GetXaxis()->SetTitle("u [mm]");
    fit2chi2ndof_image->GetYaxis()->SetTitle("v [mm]");
    fit2chi2ndof_image->GetZaxis()->SetTitle("fit2 chi2");
    fit2chi2ndof_image->GetZaxis()->SetTitleSize(0.02);
    fit2chi2ndof_image->GetZaxis()->SetLabelSize(0.02);

	// Fit probability map for first angle dist
	TH2F * fit1prob_image = new TH2F("fit1prob_image","fit1prob_image",numcol,umin,umax,numrow,vmin,vmax);
	fit1prob_image->SetStats(kFALSE);
    fit1prob_image->GetXaxis()->SetTitle("u [mm]");
    fit1prob_image->GetYaxis()->SetTitle("v [mm]");
    fit1prob_image->GetZaxis()->SetTitle("fit1 p value");
    fit1prob_image->GetZaxis()->SetTitleSize(0.02);
    fit1prob_image->GetZaxis()->SetLabelSize(0.02);

	// Fit probability map for second angle dist
	TH2F * fit2prob_image = new TH2F("fit2prob_image","fit2prob_image",numcol,umin,umax,numrow,vmin,vmax);
	fit2prob_image->SetStats(kFALSE);
    fit2prob_image->GetXaxis()->SetTitle("u [mm]");
    fit2prob_image->GetYaxis()->SetTitle("v [mm]");
    fit2prob_image->GetZaxis()->SetTitle("fit2 p value");
    fit2prob_image->GetZaxis()->SetTitleSize(0.02);
    fit2prob_image->GetZaxis()->SetLabelSize(0.02);

	// Fit probability map for merged angle dist
	TH2F * fitsumprob_image = new TH2F("fitsumprob_image","fitsumprob_image",numcol,umin,umax,numrow,vmin,vmax);
	fitsumprob_image->SetStats(kFALSE);
    fitsumprob_image->GetXaxis()->SetTitle("u [mm]");
    fitsumprob_image->GetYaxis()->SetTitle("v [mm]");
    fitsumprob_image->GetZaxis()->SetTitle("fitsum p value");
    fitsumprob_image->GetZaxis()->SetTitleSize(0.02);
    fitsumprob_image->GetZaxis()->SetLabelSize(0.02);


	// Fit mean value of first distribution
	TH2F * theta1mean_image = new TH2F("theta1mean_image","theta1mean_image",numcol,umin,umax,numrow,vmin,vmax);
	theta1mean_image->SetStats(kFALSE);
    theta1mean_image->GetXaxis()->SetTitle("u [mm]");
    theta1mean_image->GetYaxis()->SetTitle("v [mm]");
    theta1mean_image->GetZaxis()->SetTitle("theta1 mean value[rad]");
    theta1mean_image->GetZaxis()->SetTitleSize(0.02);
    theta1mean_image->GetZaxis()->SetLabelSize(0.02);


    // Fit mean value of second distribution
	TH2F * theta2mean_image = new TH2F("theta2mean_image","theta2mean_image",numcol,umin,umax,numrow,vmin,vmax);
	theta2mean_image->SetStats(kFALSE);
    theta2mean_image->GetXaxis()->SetTitle("u [mm]");
    theta2mean_image->GetYaxis()->SetTitle("v [mm]");
    theta2mean_image->GetZaxis()->SetTitle("theta2 mean value[rad]");
    theta2mean_image->GetZaxis()->SetTitleSize(0.02);
    theta2mean_image->GetZaxis()->SetLabelSize(0.02);

	// Fit mean value of first distribution
	TH2F * correctedtheta1mean_image = new TH2F("correctedtheta1mean_image","correctedtheta1mean_image",numcol,umin,umax,numrow,vmin,vmax);
	correctedtheta1mean_image->SetStats(kFALSE);
    correctedtheta1mean_image->GetXaxis()->SetTitle("u [mm]");
    correctedtheta1mean_image->GetYaxis()->SetTitle("v [mm]");
    correctedtheta1mean_image->GetZaxis()->SetTitle("theta1 mean value[rad]");
    correctedtheta1mean_image->GetZaxis()->SetTitleSize(0.02);
    correctedtheta1mean_image->GetZaxis()->SetLabelSize(0.02);

	// Fit mean value of second distribution
	TH2F * correctedtheta2mean_image = new TH2F("correctedtheta2mean_image","correctedtheta2mean_image",numcol,umin,umax,numrow,vmin,vmax);
	correctedtheta2mean_image->SetStats(kFALSE);
    correctedtheta2mean_image->GetXaxis()->SetTitle("u [mm]");
    correctedtheta2mean_image->GetYaxis()->SetTitle("v [mm]");
    correctedtheta2mean_image->GetZaxis()->SetTitle("theta2 mean value[rad]");
    correctedtheta2mean_image->GetZaxis()->SetTitleSize(0.02);
    correctedtheta2mean_image->GetZaxis()->SetLabelSize(0.02);


	// Fit mean value of u residuals
	TH2F * uresidualmean_image = new TH2F("uresidualmean_image","uresidualmean_image",numcol,umin,umax,numrow,vmin,vmax);
	uresidualmean_image->SetStats(kFALSE);
    uresidualmean_image->GetXaxis()->SetTitle("u [mm]");
    uresidualmean_image->GetYaxis()->SetTitle("v [mm]");
    uresidualmean_image->GetZaxis()->SetTitle("u residual[µm]");
    uresidualmean_image->GetZaxis()->SetTitleSize(0.02);
    uresidualmean_image->GetZaxis()->SetLabelSize(0.02);

	// Fit mean value of v residuals
	TH2F * vresidualmean_image = new TH2F("vresidualmean_image","vresidualmean_image",numcol,umin,umax,numrow,vmin,vmax);
	vresidualmean_image->SetStats(kFALSE);
    vresidualmean_image->GetXaxis()->SetTitle("u [mm]");
    vresidualmean_image->GetYaxis()->SetTitle("v [mm]");
    vresidualmean_image->GetZaxis()->SetTitle("v residual[µm]");
    vresidualmean_image->GetZaxis()->SetTitleSize(0.02);
    vresidualmean_image->GetZaxis()->SetLabelSize(0.02);
	
	// #Tracks map
	TH2F * beamspot = new TH2F("beamspot","beamspot",numcol,umin,umax,numrow,vmin,vmax);
    beamspot->SetStats(kFALSE);
    beamspot->GetXaxis()->SetTitle("u [mm]");
    beamspot->GetYaxis()->SetTitle("v [mm]");
    beamspot->GetZaxis()->SetTitle("number of tracks");
    beamspot->GetZaxis()->SetTitleOffset(1.4);
    beamspot->GetZaxis()->SetTitleSize(0.02);
    beamspot->GetZaxis()->SetLabelSize(0.02);

	// Momentum map
	TH2F * BE_image = new TH2F("BE_image","BE_image",numcol,umin,umax,numrow,vmin,vmax);
    BE_image->SetStats(kFALSE);
    BE_image->GetXaxis()->SetTitle("u [mm]");
    BE_image->GetYaxis()->SetTitle("v [mm]");
    BE_image->GetZaxis()->SetTitle("momentum [GeV/c]");
    BE_image->GetZaxis()->SetTitleOffset(1.4);
    BE_image->GetZaxis()->SetTitleSize(0.02);
    BE_image->GetZaxis()->SetLabelSize(0.02);

	// vertex w mean map
	TH2F * vertex_w_mean_image = new TH2F("vertex_w_mean_image","vertex_w_mean_image",numcol,umin,umax,numrow,vmin,vmax);
    vertex_w_mean_image->SetStats(kFALSE);
    vertex_w_mean_image->GetXaxis()->SetTitle("u [mm]");
    vertex_w_mean_image->GetYaxis()->SetTitle("v [mm]");
    vertex_w_mean_image->GetZaxis()->SetTitle("mean vertex w position [mm]");
    vertex_w_mean_image->GetZaxis()->SetTitleOffset(1.4);
    vertex_w_mean_image->GetZaxis()->SetTitleSize(0.02);
    vertex_w_mean_image->GetZaxis()->SetLabelSize(0.02);

	// vertex w rms map
	TH2F * vertex_w_rms_image = new TH2F("vertex_w_rms_image","vertex_w_rms_image",numcol,umin,umax,numrow,vmin,vmax);
    vertex_w_rms_image->SetStats(kFALSE);
    vertex_w_rms_image->GetXaxis()->SetTitle("u [mm]");
    vertex_w_rms_image->GetYaxis()->SetTitle("v [mm]");
    vertex_w_rms_image->GetZaxis()->SetTitle("vertex w position rms [mm]");
    vertex_w_rms_image->GetZaxis()->SetTitleOffset(1.4);
    vertex_w_rms_image->GetZaxis()->SetTitleSize(0.02);
    vertex_w_rms_image->GetZaxis()->SetLabelSize(0.02);

	// vertex chi2 map
	TH2F * vertex_chi2_image = new TH2F("vertex_chi2_image","vertex_chi2_image",numcol,umin,umax,numrow,vmin,vmax);
    vertex_chi2_image->SetStats(kFALSE);
    vertex_chi2_image->GetXaxis()->SetTitle("u [mm]");
    vertex_chi2_image->GetYaxis()->SetTitle("v [mm]");
    vertex_chi2_image->GetZaxis()->SetTitle("vertex chi2");
    vertex_chi2_image->GetZaxis()->SetTitleOffset(1.4);
    vertex_chi2_image->GetZaxis()->SetTitleSize(0.02);
    vertex_chi2_image->GetZaxis()->SetLabelSize(0.02);

	// vertex multiplicity map
	TH2F * vertex_multiplicity_image = new TH2F("vertex_multiplicity_image","vertex_multiplicity_image",numcol,umin,umax,numrow,vmin,vmax);
    vertex_multiplicity_image->SetStats(kFALSE);
    vertex_multiplicity_image->GetXaxis()->SetTitle("u [mm]");
    vertex_multiplicity_image->GetYaxis()->SetTitle("v [mm]");
    vertex_multiplicity_image->GetZaxis()->SetTitle("vertex multiplicity");
    vertex_multiplicity_image->GetZaxis()->SetTitleOffset(1.4);
    vertex_multiplicity_image->GetZaxis()->SetTitleSize(0.02);
    vertex_multiplicity_image->GetZaxis()->SetLabelSize(0.02);

	// vertex track mean u residual map
	TH2F * res_u_vtx_trk_image = new TH2F("res_u_vtx_trk_image","res_u_vtx_trk_image",numcol,umin,umax,numrow,vmin,vmax);
    res_u_vtx_trk_image->SetStats(kFALSE);
    res_u_vtx_trk_image->GetXaxis()->SetTitle("u [mm]");
    res_u_vtx_trk_image->GetYaxis()->SetTitle("v [mm]");
    res_u_vtx_trk_image->GetZaxis()->SetTitle("vtx trk mean u residual [mm]");
    res_u_vtx_trk_image->GetZaxis()->SetTitleOffset(1.4);
    res_u_vtx_trk_image->GetZaxis()->SetTitleSize(0.02);
    res_u_vtx_trk_image->GetZaxis()->SetLabelSize(0.02);

	// vertex track mean v residual map
	TH2F * res_v_vtx_trk_image = new TH2F("res_v_vtx_trk_image","res_v_vtx_trk_image",numcol,umin,umax,numrow,vmin,vmax);
    res_v_vtx_trk_image->SetStats(kFALSE);
    res_v_vtx_trk_image->GetXaxis()->SetTitle("u [mm]");
    res_v_vtx_trk_image->GetYaxis()->SetTitle("v [mm]");
    res_v_vtx_trk_image->GetZaxis()->SetTitle("vtx trk mean v residual [mm]");
    res_v_vtx_trk_image->GetZaxis()->SetTitleOffset(1.4);
    res_v_vtx_trk_image->GetZaxis()->SetTitleSize(0.02);
    res_v_vtx_trk_image->GetZaxis()->SetLabelSize(0.02);

	// vertex track u residual rms map
	TH2F * res_u_rms_vtx_trk_image = new TH2F("res_u_rms_vtx_trk_image","res_u_rms_vtx_trk_image",numcol,umin,umax,numrow,vmin,vmax);
    res_u_rms_vtx_trk_image->SetStats(kFALSE);
    res_u_rms_vtx_trk_image->GetXaxis()->SetTitle("u [mm]");
    res_u_rms_vtx_trk_image->GetYaxis()->SetTitle("v [mm]");
    res_u_rms_vtx_trk_image->GetZaxis()->SetTitle("vtx trk u residual rms [mm]");
    res_u_rms_vtx_trk_image->GetZaxis()->SetTitleOffset(1.4);
    res_u_rms_vtx_trk_image->GetZaxis()->SetTitleSize(0.02);
    res_u_rms_vtx_trk_image->GetZaxis()->SetLabelSize(0.02);

	// vertex track v residual rms map
	TH2F * res_v_rms_vtx_trk_image = new TH2F("res_v_rms_vtx_trk_image","res_v_rms_vtx_trk_image",numcol,umin,umax,numrow,vmin,vmax);
    res_v_rms_vtx_trk_image->SetStats(kFALSE);
    res_v_rms_vtx_trk_image->GetXaxis()->SetTitle("u [mm]");
    res_v_rms_vtx_trk_image->GetYaxis()->SetTitle("v [mm]");
    res_v_rms_vtx_trk_image->GetZaxis()->SetTitle("vtx trk v residual rms [mm]");
    res_v_rms_vtx_trk_image->GetZaxis()->SetTitleOffset(1.4);
    res_v_rms_vtx_trk_image->GetZaxis()->SetTitleSize(0.02);
    res_v_rms_vtx_trk_image->GetZaxis()->SetLabelSize(0.02);

	for(int col=0; col<numcol; col++)
	{
		for(int row=0; row<numrow; row++)
		{
			// Get the root file of this prtial image
			filenameB.Form("-part-%i-%i.root",int(col/max_u_pixels)+1,num_vsplits-int(row/max_v_pixels)); 
			//cout<<"Part "<<(int(col/max_u_pixels)+1)<<" "<<(num_vsplits-int(row/max_v_pixels))<<endl;

			filename=filenameA+filenameB;
			X0file = new TFile(filename, "READ");

			x0_image_aid=(TH2F*)X0file->Get("mapping/result/x0_image");
			x0err_image_aid=(TH2F*)X0file->Get("mapping/result/x0err_image");
			x0relerr_image_aid=(TH2F*)X0file->Get("mapping/result/x0relerr_image");
			fit1chi2ndof_image_aid=(TH2F*)X0file->Get("mapping/result/fit1chi2ndof_image");
			fit2chi2ndof_image_aid=(TH2F*)X0file->Get("mapping/result/fit2chi2ndof_image");
			fitsumchi2ndof_image_aid=(TH2F*)X0file->Get("mapping/result/fitsumchi2ndof_image");
			fit1prob_image_aid=(TH2F*)X0file->Get("mapping/result/fit1prob_image");
			fit2prob_image_aid=(TH2F*)X0file->Get("mapping/result/fit2prob_image");
			fitsumprob_image_aid=(TH2F*)X0file->Get("mapping/result/fitsumprob_image");
			theta1mean_image_aid=(TH2F*)X0file->Get("mapping/result/theta2mean_image");
			theta2mean_image_aid=(TH2F*)X0file->Get("mapping/result/theta2mean_image");
			correctedtheta1mean_image_aid=(TH2F*)X0file->Get("mapping/result/correctedtheta1mean_image");
			correctedtheta2mean_image_aid=(TH2F*)X0file->Get("mapping/result/correctedtheta2mean_image");
			uresidualmean_image_aid=(TH2F*)X0file->Get("mapping/result/uresidualmean_image");
			vresidualmean_image_aid=(TH2F*)X0file->Get("mapping/result/vresidualmean_image");
			beamspot_aid=(TH2F*)X0file->Get("mapping/result/beamspot");
			BE_image_aid=(TH2F*)X0file->Get("mapping/result/BE_image");

			vertex_w_mean_image_aid=(TH2F*)X0file->Get("mapping/result/vertex_w_image");				    
			vertex_w_rms_image_aid=(TH2F*)X0file->Get("mapping/result/vertex_w_rms_image");
			vertex_chi2_image_aid=(TH2F*)X0file->Get("mapping/result/vertex_chi2_image");
			vertex_multiplicity_image_aid=(TH2F*)X0file->Get("mapping/result/vertex_multiplicity_image");

			res_u_vtx_trk_image_aid=(TH2F*)X0file->Get("mapping/result/u_res_mean_vtx_trk_image");
			res_v_vtx_trk_image_aid=(TH2F*)X0file->Get("mapping/result/v_res_mean_vtx_trk_image");
			res_u_rms_vtx_trk_image_aid=(TH2F*)X0file->Get("mapping/result/u_res_rms_vtx_trk_image");
			res_v_rms_vtx_trk_image_aid=(TH2F*)X0file->Get("mapping/result/v_res_rms_vtx_trk_image");

			fit1prob_histo_aid=(TH1F*)X0file->Get("mapping/result/fit1prob_histo");
			fit2prob_histo_aid=(TH1F*)X0file->Get("mapping/result/fit2prob_histo");
			fitsumprob_histo_aid=(TH1F*)X0file->Get("mapping/result/fitsumprob_histo");

			fit1prob_histo->Add(fit1prob_histo_aid);
			fit2prob_histo->Add(fit2prob_histo_aid);
			fitsumprob_histo->Add(fitsumprob_histo_aid);

			// Copy entries of first images
			x0_image->SetBinContent(col+1,row+1,x0_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			x0_image->SetBinError(col+1,row+1,x0err_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			x0err_image->SetBinContent(col+1,row+1,x0err_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			x0relerr_image->SetBinContent(col+1,row+1,x0relerr_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			fit1chi2ndof_image->SetBinContent(col+1,row+1,fit1chi2ndof_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			fit2chi2ndof_image->SetBinContent(col+1,row+1,fit2chi2ndof_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			fitsumchi2ndof_image->SetBinContent(col+1,row+1,fitsumchi2ndof_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			fit1prob_image->SetBinContent(col+1,row+1,fit1prob_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			fit2prob_image->SetBinContent(col+1,row+1,fit2prob_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			fitsumprob_image->SetBinContent(col+1,row+1,fitsumprob_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			theta1mean_image->SetBinContent(col+1,row+1,theta1mean_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			theta2mean_image->SetBinContent(col+1,row+1,theta2mean_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			correctedtheta1mean_image->SetBinContent(col+1,row+1,correctedtheta1mean_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			correctedtheta2mean_image->SetBinContent(col+1,row+1,correctedtheta2mean_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			uresidualmean_image->SetBinContent(col+1,row+1,uresidualmean_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			vresidualmean_image->SetBinContent(col+1,row+1,vresidualmean_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			beamspot->SetBinContent(col+1,row+1,beamspot_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			BE_image->SetBinContent(col+1,row+1,BE_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));

			vertex_w_mean_image->SetBinContent(col+1,row+1,vertex_w_mean_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			vertex_w_rms_image->SetBinContent(col+1,row+1,vertex_w_rms_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			vertex_chi2_image->SetBinContent(col+1,row+1,vertex_chi2_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			vertex_multiplicity_image->SetBinContent(col+1,row+1,vertex_multiplicity_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
    		res_u_vtx_trk_image->SetBinContent(col+1,row+1,res_u_vtx_trk_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			res_v_vtx_trk_image->SetBinContent(col+1,row+1,res_v_vtx_trk_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
    		res_u_rms_vtx_trk_image->SetBinContent(col+1,row+1,res_u_rms_vtx_trk_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));
			res_v_rms_vtx_trk_image->SetBinContent(col+1,row+1,res_v_rms_vtx_trk_image_aid->GetBinContent(col%max_u_pixels+1,row%max_v_pixels+1));

			delete X0file;
			X0file = ((TFile *)0);


		}
	}

	Resultsfile->cd("");
	
	x0_image->Write();
	x0err_image->Write();
	x0relerr_image->Write();
	fit1chi2ndof_image->Write();
	fit2chi2ndof_image->Write();
	fitsumchi2ndof_image->Write();
	fit1prob_image->Write();
	fit2prob_image->Write();
	fitsumprob_image->Write();
	theta1mean_image->Write();
	theta2mean_image->Write();
	correctedtheta1mean_image->Write();
	correctedtheta2mean_image->Write();
	uresidualmean_image->Write();
	vresidualmean_image->Write();

	beamspot->Write();
	BE_image->Write();

	fit1prob_histo->Write();
	fit2prob_histo->Write();
	fitsumprob_histo->Write();

	hscatt_theta1_vs_resu->Write();
	hscatt_theta1_vs_resv->Write();
	hscatt_theta2_vs_resu->Write();
	hscatt_theta2_vs_resv->Write();

	vertex_w_mean_image->Write();
	vertex_w_rms_image->Write();
	vertex_chi2_image->Write();
	vertex_multiplicity_image->Write();

	res_u_vtx_trk_image->Write();
	res_v_vtx_trk_image->Write();
	res_u_rms_vtx_trk_image->Write();
	res_v_rms_vtx_trk_image->Write();

	Resultsfile->Close();

    return 0;

	
}
