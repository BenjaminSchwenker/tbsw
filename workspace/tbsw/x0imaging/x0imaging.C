//#include <iostream.h>
#include <fstream>
using namespace std ;



  // Highland model of a MSC angle distribution, the parameters are:

  /*
	* par[0]:  Expected beam energy at u,v=0;
	* par[1]:  Beam particle charge
	* par[2]:  Beam particle mass
	* par[3]:  Radiation length X/X0
	* par[4]:  Expected angle reconstruction error
	* par[5]:  Normalization
	* par[6]:  mean value

*/
  
  // Highland model of multiple scattering: Simple gaussian with a well defined standard deviation depending on X/X0 and the beam energy.
  // The overall function describing the kink angle distributions is the Highland function convoluted with a gaussian function due to the finite angle resolution on the target plane. 
  Double_t highlandfunction(Double_t *x, Double_t *par)
  { 

	// Other parameters
	double recoerror=par[4];  //expected reconstruction error

    // particle parameters

	// mass of beam particle
	double mass;   
	mass=par[2];  

	// charge of beam particle
	double charge;   
	charge=par[1];

	// beam energy
	double p=par[0];

	// calibrated momentum
	double E=TMath::Sqrt(p*p+mass*mass);  // energy in GeV

	double beta;  //relative velocity
	beta=p/E;

	// Radiation length computed from the other parameters
	double XX0=par[3];

	// Combination of Highland width and reconstruction error
	double sigma=TMath::Sqrt(pow(0.0136*charge/(p*beta)*TMath::Sqrt(XX0)*(1.0+0.038*TMath::Log(XX0)),2)+pow(recoerror,2));

	// function value at a certain theta value
	double value=par[5]*TMath::Gaus(x[0],par[6],sigma);

	return value;
  }// End definition of highland model


  // Returns the mean value of the angle reco error squared
  double getanglerecovar(TFile* file)
  {
	//TTree in input root file, that contains the MSC projected angle distributions and reconstruction error distribution
	file->cd("");
	TTree *msc_tree = (TTree*)file->Get("MSCTree");

	// Draw reconstruction error 1 histogram
	msc_tree->Draw("theta1_var", "", "P*");

	// Get mean value
	double recovar = msc_tree->GetHistogram()->GetMean();

	// Get maximum value
	double recovar_error = msc_tree->GetHistogram()->GetMeanError();	//max methode

	// The recovar error shouldn't be too large! 
	cout<<"Angle Reconstruction Variance is "<<recovar<<" +/- "<<recovar_error<<endl; 

	return recovar;
  }



  // This function fills histograms corresponding to certain u v values with msc angle distributions 
  void getcorrection(TFile* file1, TFile* file2, std::vector<double> means, std::vector<double> plotranges, int numberofbins, const int numcol, const int numrow, double umin, double vmin, double umax, double vmax, int vertex_multiplicity_min, int vertex_multiplicity_max)
  {
	// parameters which are read out from the root file
	Double_t theta1;
	Double_t theta2;
	Double_t u;
	Double_t v;
	Int_t vertex_multiplicity=1;
	
	//TTree in input root file, that contains the MSC projected angle distributions and reconstruction error distribution
	file1->cd("");

	TTree * msc_tree = (TTree*) file1->Get("MSCTree");

	msc_tree->SetBranchAddress("theta1",&theta1);
	msc_tree->SetBranchAddress("theta2",&theta2);
	msc_tree->SetBranchAddress("u",&u);
	msc_tree->SetBranchAddress("v",&v);

	int test=msc_tree->SetBranchAddress("vertex_multiplicity",&vertex_multiplicity);

	file2->cd("");
	file2->cd("mapping/raw");

	// arrays of msc angle histograms
	TH1F *histo_theta1[numcol][numrow];
	TH1F *histo_theta2[numcol][numrow];
	TH1F *histo_uresidual[numcol][numrow];
	TH1F *histo_vresidual[numcol][numrow];
	TH1F *histo_thetasum[numcol][numrow];
	TH2F *histo_2d[numcol][numrow];

	for (int i=0; i<numcol; i++)
	{
		for (int j=0; j<numrow; j++)
		{
		 	histo_theta1[i][j] = new TH1F("","",numberofbins,means.at(0)-1.0*plotranges.at(0),means.at(0)+plotranges.at(0));
		 	histo_theta2[i][j] = new TH1F("","",numberofbins,means.at(1)-1.0*plotranges.at(1),means.at(1)+plotranges.at(1));
		}
	}

	// Loop over all events
	for(int i=0; i< msc_tree->GetEntries(); i++)
	{

		if(i%1000000==0) cout<<"Mean Correction loop, Track No. "<<i<<endl;
		msc_tree->GetEntry(i);
		
		// position within the map area
		double u_pos=u-umin;
		double v_pos=v-vmin;

		// side lengths of the map area
		double u_length=umax-umin;
		double v_length=vmax-vmin;

		// skip this entry, when u or v is outside of the mapping area
	    if (u_pos<0) continue;
	    if (u_pos>=u_length) continue;
	    if (v_pos<0) continue;
	    if (v_pos>=v_length) continue;
        
		// Apply cut on vertex multiplicity
        if (vertex_multiplicity>vertex_multiplicity_max||vertex_multiplicity<vertex_multiplicity_min) continue;

		// Determine column and row number from the position within the map area and the number of rows and columns
		int col=floor(u_pos*numcol/u_length);
		int row=floor(v_pos*numrow/v_length);

		// Fill histograms
		histo_theta1[col][row]->Fill(theta1);
		histo_theta2[col][row]->Fill(theta2);
	}

	cout<<"Write angle histograms "<<endl;

	for (int i=0; i<numcol; i++)
	{
		for (int j=0; j<numrow; j++)
		{
			// Name of the histograms
			TString histoname;
			histoname.Form("area(%i,%i)",i,j);

			 histo_theta1[i][j]->Write("theta1_uncorrected_"+histoname);
			 histo_theta2[i][j]->Write("theta2_uncorrected_"+histoname);
		}
	}

  }


  // This function fills histograms corresponding to certain u v values with msc angle distributions 
  void savehistos(TFile* file1, TFile* file2, int numberofbins, double histo_range, const int numcol, const int numrow, double umin, double vmin, double umax, double vmax, int vertex_multiplicity_min, int vertex_multiplicity_max)
  {
	// parameters which are read out from the root file
	Double_t theta1;
	Double_t theta2;
	Double_t u;
	Double_t v;
	Double_t u_in;
	Double_t v_in;
	Double_t u_out;
	Double_t v_out;

	Int_t vertex_multiplicity=1;
	Double_t vertex_w;
	Double_t vertex_chi2;
	Double_t vertex_u;
	Double_t vertex_v;

	// Array of mean theta1 and theta2 values in each map pixel
	double mean1[numcol][numrow];
	double mean2[numcol][numrow];
	
	//TTree in input root file, that contains the MSC projected angle distributions and reconstruction error distribution
	file1->cd("");

	TTree * msc_tree = (TTree*) file1->Get("MSCTree");

	// Set branch adresses for parameters connected to the scattering angles
	msc_tree->SetBranchAddress("theta1",&theta1);
	msc_tree->SetBranchAddress("theta2",&theta2);
	msc_tree->SetBranchAddress("u",&u);
	msc_tree->SetBranchAddress("v",&v);
	msc_tree->SetBranchAddress("u_in",&u_in);
	msc_tree->SetBranchAddress("v_in",&v_in);
	msc_tree->SetBranchAddress("u_out",&u_out);
	msc_tree->SetBranchAddress("v_out",&v_out);

	// Set branch adresses for parameters connected to the vertex fit
	msc_tree->SetBranchAddress("vertex_w",&vertex_w);
	msc_tree->SetBranchAddress("vertex_chi2ndf",&vertex_chi2);
	int test=msc_tree->SetBranchAddress("vertex_multiplicity",&vertex_multiplicity);
	msc_tree->SetBranchAddress("vertex_u",&vertex_u);
	msc_tree->SetBranchAddress("vertex_v",&vertex_v);

	file2->cd("");
	file2->cd("mapping/raw");

	// arrays of msc angle histograms
	TH1F *histo_theta1[numcol][numrow];
	TH1F *histo_theta2[numcol][numrow];
	TH1F *histo_uresidual[numcol][numrow];
	TH1F *histo_vresidual[numcol][numrow];
	TH1F *histo_thetasum[numcol][numrow];
	TH2F *histo_2d[numcol][numrow];

	// arrays of vertex parameter histograms
	TH1F *histo_vertex_w[numcol][numrow];
	TH1F *histo_vertex_multiplicity[numcol][numrow];
	TH1F *histo_vertex_chi2[numcol][numrow];
	TH1F *histo_vtx_trk_u_res[numcol][numrow];
	TH1F *histo_vtx_trk_v_res[numcol][numrow];

	for (int i=0; i<numcol; i++)
	{
		for (int j=0; j<numrow; j++)
		{

			// Get the histograms generated in the getcorrection function
			// Name of the histograms
			TString aidhistoname;
			aidhistoname.Form("area(%i,%i)",i,j);

			// Get histogram
			TH1* histogram1=(TH1*)file2->Get("mapping/raw/theta1_uncorrected_"+aidhistoname);
			TH1* histogram2=(TH1*)file2->Get("mapping/raw/theta2_uncorrected_"+aidhistoname);
			
			// Determine plot range from uncorrected histograms
			double limits=histo_range/2.0*(histogram1->GetRMS()+histogram2->GetRMS());

		 	histo_theta1[i][j] = new TH1F("","",numberofbins,-limits,limits);
		 	histo_theta2[i][j] = new TH1F("","",numberofbins,-limits,limits);
		 	histo_uresidual[i][j] = new TH1F("","",1000,-1.0,1.0);
		 	histo_vresidual[i][j] = new TH1F("","",1000,-1.0,1.0);
		 	histo_thetasum[i][j] = new TH1F("","",numberofbins,-limits,limits);
			histo_2d[i][j] = new TH2F("","",numberofbins,-limits,limits,numberofbins,-limits,limits);

			// Set range and title of vertex histograms
		 	histo_vertex_w[i][j] = new TH1F("","",900,-15.0,15.0);
		 	histo_vertex_chi2[i][j] = new TH1F("","",200,0.0,11.0);
		 	histo_vertex_multiplicity[i][j] = new TH1F("","",10,0.0,10);
		 	histo_vtx_trk_u_res[i][j] = new TH1F("","",1000,-5.0,5.0);
		 	histo_vtx_trk_v_res[i][j] = new TH1F("","",1000,-5.0,5.0);

			// Save mean of both distributions in arrays
			mean1[i][j]=histogram1->GetMean();
			mean2[i][j]=histogram2->GetMean();

            histogram1->Delete();
            histogram2->Delete();

		}
	}

	file2->cd("");
	file2->cd("mapping/raw");

	// Loop over all events, find the corresponding image pixel from the u,v values and fill the histograms
	for(int i=0; i< msc_tree->GetEntries(); i++)
	{

		if(i%1000000==0) cout<<"Track No. "<<i<<endl;
		msc_tree->GetEntry(i);
		
		// position within the map area
		double u_pos=u-umin;
		double v_pos=v-vmin;

		// side lengths of the map area
		double u_length=umax-umin;
		double v_length=vmax-vmin;

		// skip this entry, when u or v is outside of the image area
	    if (u_pos<0) continue;
	    if (u_pos>=u_length) continue;
	    if (v_pos<0) continue;
	    if (v_pos>=v_length) continue;

		// Apply cut on vertex multiplicity
        if (vertex_multiplicity>vertex_multiplicity_max||vertex_multiplicity<vertex_multiplicity_min) continue;

		// Determine column and row number from the position within the map area and the number of rows and columns
		int col=floor(u_pos*numcol/u_length);
		int row=floor(v_pos*numrow/v_length);

		// mean correction of theta
		theta1=theta1-mean1[col][row];
		theta2=theta2-mean2[col][row];

		// Fill histograms
		histo_theta1[col][row]->Fill(theta1);
		histo_theta2[col][row]->Fill(theta2);
		histo_thetasum[col][row]->Fill(theta1);
		histo_thetasum[col][row]->Fill(theta2);
		histo_2d[col][row]->Fill(theta1,theta2);

		histo_uresidual[col][row]->Fill(u_in-u_out);
		histo_vresidual[col][row]->Fill(v_in-v_out);

		histo_vertex_w[col][row]->Fill(vertex_w);
		histo_vertex_chi2[col][row]->Fill(vertex_chi2);
		histo_vertex_multiplicity[col][row]->Fill(vertex_multiplicity);

		histo_vtx_trk_u_res[col][row]->Fill(vertex_u-u);
		histo_vtx_trk_v_res[col][row]->Fill(vertex_v-v);
	}

	cout<<"Write histograms "<<endl;

	for (int i=0; i<numcol; i++)
	{
		for (int j=0; j<numrow; j++)
		{
			// Name of the histograms
			TString histoname;
			histoname.Form("area(%i,%i)",i,j);

			histo_theta1[i][j]->Write("theta1_"+histoname);
			histo_theta1[i][j]->Delete();
			histo_theta2[i][j]->Write("theta2_"+histoname);
			histo_theta2[i][j]->Delete();
			histo_uresidual[i][j]->Write("uresidual_"+histoname);
			histo_uresidual[i][j]->Delete();
			histo_vresidual[i][j]->Write("vresidual_"+histoname);
			histo_vresidual[i][j]->Delete();
			histo_thetasum[i][j]->Write("sumhisto_"+histoname);
			histo_thetasum[i][j]->Delete();
			histo_2d[i][j]->Write("2Dhisto_"+histoname);
			histo_2d[i][j]->Delete();

			histo_vertex_w[i][j]->Write("vertex_w_"+histoname);
			histo_vertex_w[i][j]->Delete();
			histo_vertex_chi2[i][j]->Write("vertex_chi2_"+histoname);
			histo_vertex_chi2[i][j]->Delete();
			histo_vertex_multiplicity[i][j]->Write("vertex_multiplicity_"+histoname);
			histo_vertex_multiplicity[i][j]->Delete();

			histo_vtx_trk_u_res[i][j]->Write("res_u_vtx_trk_"+histoname);
			histo_vtx_trk_u_res[i][j]->Delete();
			histo_vtx_trk_v_res[i][j]->Write("res_v_vtx_trk_"+histoname);
			histo_vtx_trk_v_res[i][j]->Delete();

		}
		
	}

  }


// Determine fit range for a kink angle histogram
// This is done by finding the first and last bin above a certain threshold.
// The fit range is half of the distance between these two bins in rad
double DetermineFitrange(TH1* histo,double rangevalue)
{

    // Clone histo
	TH1F *h2 = (TH1F*) histo->Clone();

	cout<<"RMS value of distribution: "<<histo->GetRMS()<<endl;
	cout<<"Selected range parameter: "<<rangevalue<<" -> Fit range up to y=1/("<<rangevalue<<"*e)"<<endl;
	double fitrange = sqrt(2.0*rangevalue)*histo->GetRMS();

	// Use RMS value as a rough measure of the fit range for a gaussian fit
	TF1 *f1 = new TF1("f1","gaus(x)",-fitrange,fitrange);
	f1->SetLineStyle(2);
	TFitResultPtr fitr=h2->Fit("f1","RS");

	// Repeat fit in case it failed
	if(fitr!=0)
	{
		cout<<"Fit of angle distribution failed with status: "<<fitr<<endl;
		cout<<"Repeat fit "<<endl;
		h2->Fit("f1","RM");
	}

	// Use the determined sigma value to calculate the fit range
	double sigma = f1->GetParameter(2);
	fitrange=sqrt(2.0*rangevalue)*sigma;
	cout<<"Determined fit range: " << fitrange<<endl<<endl;

	return fitrange;
}

  // Function to fit the MSC angle histograms and fill the map histograms
  void fithisto( TFile* file, int fittype, double maxchi2ndof_fit, double rangevalue, int col, int numcol, int row,int numrow, double* parameters, TString fitoptions)
  { 
	// Calculate number of parameters
	int num_parameters = 7;
	//cout<<num_parameters<<endl;

	// Open the histograms
 
	// histogram name
	TString histoname;
	histoname.Form("area(%i,%i)",col,row);
	file->cd("");
	TH1* histogram1=(TH1*)file->Get("mapping/raw/theta1_"+histoname);
	TH1* histogram2=(TH1*)file->Get("mapping/raw/theta2_"+histoname);
	TH1* histogramu=(TH1*)file->Get("mapping/raw/uresidual_"+histoname);
	TH1* histogramv=(TH1*)file->Get("mapping/raw/vresidual_"+histoname);
	TH1* histogramsum=(TH1*)file->Get("mapping/raw/sumhisto_"+histoname);

	// Get vertex histos
	TH1* histogram_vertex_w=(TH1*)file->Get("mapping/raw/vertex_w_"+histoname);
	TH1* histogram_vertex_chi2=(TH1*)file->Get("mapping/raw/vertex_chi2_"+histoname);
	TH1* histogram_vertex_multiplicity=(TH1*)file->Get("mapping/raw/vertex_multiplicity_"+histoname);

	TH1* histogram_resu_vtx_trk=(TH1*)file->Get("mapping/raw/res_u_vtx_trk_"+histoname);
	TH1* histogram_resv_vtx_trk=(TH1*)file->Get("mapping/raw/res_v_vtx_trk_"+histoname);


	double uncorrected_mean1=histogram1->GetMean();
	double uncorrected_mean2=histogram1->GetMean();

    if((file->Get("mapping/raw/theta1_uncorrected_"+histoname)!=NULL)&&(file->Get("mapping/raw/theta2_uncorrected_"+histoname)!=NULL))
	{

		TH1* histogram_uncorrected1=(TH1*)file->Get("mapping/raw/theta1_uncorrected_"+histoname);
		TH1* histogram_uncorrected2=(TH1*)file->Get("mapping/raw/theta2_uncorrected_"+histoname);

		if((histogram_uncorrected1->GetEntries()>0)&&(histogram_uncorrected2->GetEntries()>0))
		{
			uncorrected_mean1=histogram_uncorrected1->GetMean();
			uncorrected_mean2=histogram_uncorrected2->GetMean();
		}

		else
		{

			uncorrected_mean1=0;
			uncorrected_mean2=0;

		}
		
	}

	file->cd("");
	file->cd("mapping/fit");

	// Make a copy of the histogram
	TH1* fithistogram1=(TH1*)histogram1->Clone("fithisto");
	TH1* fithistogram2=(TH1*)histogram2->Clone("fithisto2");
	TH1* fithistogramsum=(TH1*)histogramsum->Clone("fithistosum");

	// Fit result parameters of both angle distribution

	// mean of the gaussian and its error
	double mean1, mean2, meansum;
	double mean_error1,mean_error2, mean_errorsum;

	// quality parameters of the fit
	double chi2ndof1,chi2ndof2,chi2ndofsum;
	double prob1,prob2,probsum;

	// X/X0 values from the fit
	double XX01,XX02,XX0sum;
	double XX0err1,XX0err2,XX0errsum;

	// vertex parameters
	double vertex_chi2,vertex_w_mean,vertex_w_rms;
	double vertex_multiplicity;
	double vtx_trk_res_u_mean,vtx_trk_res_v_mean,vtx_trk_res_u_rms,vtx_trk_res_v_rms;


	// Variables used to calculate the fit range of the histograms
	int bin1;
	int bin2;
	double fitrange;

	double uresidual_mean;
	double vresidual_mean;

	double minvalue=1.0/(rangevalue*2.7);

	int NumberOfTracks=fithistogram1->GetEntries();

	// Get residual values for this image bin from histogram
	uresidual_mean=histogramu->GetMean();
	vresidual_mean=histogramv->GetMean();

	// Get vertex values for this image bin from histogram
	vertex_w_mean=histogram_vertex_w->GetMean();
	vertex_w_rms=histogram_vertex_w->GetRMS();
	vertex_chi2=histogram_vertex_chi2->GetMean();
	vertex_multiplicity=histogram_vertex_multiplicity->GetMean();
	vtx_trk_res_u_mean=histogram_resu_vtx_trk->GetMean();
	vtx_trk_res_v_mean=histogram_resv_vtx_trk->GetMean();
	vtx_trk_res_u_rms=histogram_resu_vtx_trk->GetRMS();
	vtx_trk_res_v_rms=histogram_resv_vtx_trk->GetRMS();

	if(NumberOfTracks>400)
	{
		//Fit histogram with Highland model

		// Fit functions definition
		// The fit range is determined from the RMS value of the histogram
		// The Highland model is validated up to the width sigma_1/e, where the 
		// distribution reaches the values max/e. The fit range should be therefore limited to 
		// The intervall [-sqrt(2)*RMS,+sqrt(2)*RMS]
		fitrange =DetermineFitrange(fithistogram1,rangevalue);
		TF1 *fit1 = new TF1("theta1_fit",highlandfunction,-fitrange,fitrange,num_parameters);

		fitrange =DetermineFitrange(fithistogram2,rangevalue);
		TF1 *fit2 = new TF1("theta2_fit",highlandfunction,-fitrange,fitrange,num_parameters);

		fitrange =DetermineFitrange(fithistogramsum,rangevalue);
		TF1 *fitsum = new TF1("thetasum_fit",highlandfunction,-fitrange,fitrange,num_parameters);



		// Set starting values
   		fit1->SetParameters(parameters);
   		fit2->SetParameters(parameters);
   		fitsum->SetParameters(parameters);

		for(int i=0; i<num_parameters;i++)
		{
			if(i!=3&&i!=5)
			{
   				fit1->FixParameter(i,parameters[i]);
   				fit2->FixParameter(i,parameters[i]);
   				fitsum->FixParameter(i,parameters[i]);
			}
		}	

		// The mean value of the distribution should not be shifted more than 100 µrad -> limit parameter 6
		// --> Just don't fit the mean, just fix it at 0.0
   		//fit1->SetParLimits(6,-0.0001,0.0001);
   		//fit2->SetParLimits(6,-0.0001,0.0001);
   		//fitsum->SetParLimits(6,-0.0001,0.0001);

   		fit1->SetParLimits(3,0.00001,200.0);
   		fit2->SetParLimits(3,0.00001,200.0);
   		fitsum->SetParLimits(3,0.00001,200.0);

		TFitResultPtr fitr=fithistogram1->Fit("theta1_fit",fitoptions);
		if(fitr!=0)
		{
			cout<<"Fit of first angle distribution failed with status: "<<fitr<<endl;
			cout<<"Repeat fit "<<endl;
			fithistogram1->Fit("theta1_fit",fitoptions);
		}
		fitr=fithistogram2->Fit("theta2_fit",fitoptions);
		if(fitr!=0)
		{
			cout<<"Fit of second angle distribution failed with status: "<<fitr<<endl;
			cout<<"Repeat fit "<<endl;
			fithistogram2->Fit("theta2_fit",fitoptions);
		}
		fitr=fithistogramsum->Fit("thetasum_fit",fitoptions);
		if(fitr!=0)
		{
			cout<<"Fit of combined angle distribution failed with status: "<<fitr<<endl;
			cout<<"Repeat fit "<<endl;
			fithistogramsum->Fit("thetasum_fit",fitoptions);
		}

		// Names of the 14 local parameters
		const int num_localparameters=7;
		TString name[num_localparameters];
		name[0]="E";
		name[1]="z[e]";
		name[2]="m";
		name[3]="X/X0";
		name[4]="reco err";
		name[5]="norm";
		name[6]="mean";

		gStyle->SetOptFit(1111);

		for(int iname=0;iname<num_localparameters;iname++) 
		{
			fit1->SetParName(iname,name[iname]);			
			fit2->SetParName(iname,name[iname]);
			fitsum->SetParName(iname,name[iname]);
		}

		// Chi2ndof of the fit
		chi2ndof1=fit1->GetChisquare()/(fit1->GetNDF()*1.0);
		chi2ndof2=fit2->GetChisquare()/(fit2->GetNDF()*1.0);
		chi2ndofsum=fitsum->GetChisquare()/(fitsum->GetNDF()*1.0);	

		// define chi2 cut
		double chi2_cut=maxchi2ndof_fit;

		// Use the chi2 values for quality cuts
		if((fittype==0)&&((chi2ndof1+chi2ndof2)>chi2_cut*2.0))
		{
			// The fit didn't work: Dont trust this X/X0 and theta mean value
			XX01=0.0;
			XX0err1=99.0;
			XX02=0.0;
			XX0err2=99.0;
			mean1=0;
			mean2=0;
			// save the fitted histograms
			file->cd("");
			file->cd("mapping/badfit");
			fithistogram1->GetListOfFunctions()->Add(fit1);
			fithistogram1->Write("theta1_"+histoname+"_fit");
			fithistogram2->GetListOfFunctions()->Add(fit2);
			fithistogram2->Write("theta2_"+histoname+"_fit");
			fithistogramsum->GetListOfFunctions()->Add(fitsum);
			fithistogramsum->Write("sumhisto_"+histoname+"_fit");

		}

		else if ((fittype==1)&&(chi2ndofsum>chi2_cut))
		{
			// The fit didn't work: Dont trust this X/X0 and theta mean value
			XX0sum=0.0;
			XX0errsum=99.0;
			mean1=0;
			mean2=0;

			// save the fitted histograms
			file->cd("");
			file->cd("mapping/badfit");
			fithistogram1->GetListOfFunctions()->Add(fit1);
			fithistogram1->Write("theta1_"+histoname+"_fit");
			fithistogram2->GetListOfFunctions()->Add(fit2);
			fithistogram2->Write("theta2_"+histoname+"_fit");
			fithistogramsum->GetListOfFunctions()->Add(fitsum);
			fithistogramsum->Write("sumhisto_"+histoname+"_fit");

		}

		else
		{

			// Extract Fit parameters of the first distribution
			XX01=fit1->GetParameter(3);
			XX0err1=fit1->GetParError(3);
			mean1=fit1->GetParameter(6);
			prob1=fit1->GetProb();

			// Extract Fit parameters of the second distribution
			XX02=fit2->GetParameter(3);
			XX0err2=fit2->GetParError(3);
			mean2=fit2->GetParameter(6);
			prob2=fit2->GetProb();

			// Extract Fit parameters of the merged distribution
			XX0sum=fitsum->GetParameter(3);
			XX0errsum=fitsum->GetParError(3);
			meansum=fitsum->GetParameter(6);
			probsum=fitsum->GetProb();

			// save the fitted histograms
			file->cd("");
			file->cd("mapping/fit");
			fithistogram1->GetListOfFunctions()->Add(fit1);
			fithistogram1->Write("theta1_"+histoname+"_fit");
			fithistogram2->GetListOfFunctions()->Add(fit2);
			fithistogram2->Write("theta2_"+histoname+"_fit");
			fithistogramsum->GetListOfFunctions()->Add(fitsum);
			fithistogramsum->Write("sumhisto_"+histoname+"_fit");
		}



	}

	// If the histogram has too few tracks for the fit, just set all relevant fit results to unrealistic values
	else
	{

		XX01=0.0;
		XX0err1=99.0;
		mean1=0.0;
		chi2ndof1=100.0;
		prob1=-1.0;

		XX02=0.0;
		XX0err2=99.0;
		mean2=0.0;
		chi2ndof2=100.0;
		prob2=-1.0;

		XX0sum=0.0;
		XX0errsum=99.0;
		meansum=0.0;
		chi2ndofsum=100.0;
		probsum=-1.0;


		cout<<"Too few tracks in this pixel... no fit done!"<<endl;
	}

	// Go to results directory
	file->cd("");
	file->cd("mapping/result");

	// Get the chi2 of the fit of the merged distribution
	TH2* chi2ndofsummap=(TH2*)file->Get("mapping/result/fitsumchi2ndof_image");

	// Get the chi2 of the fit of the theta1 distribution
	TH2* chi2ndof1map=(TH2*)file->Get("mapping/result/fit1chi2ndof_image");

	// Get the chi2 of the fit of the theta2 distribution
	TH2* chi2ndof2map=(TH2*)file->Get("mapping/result/fit2chi2ndof_image");

	// Get the mean1 histogram
	TH2* meanmap1=(TH2*)file->Get("mapping/result/theta1mean_image");

	// Get the mean2 histogram
	TH2* meanmap2=(TH2*)file->Get("mapping/result/theta2mean_image");

	// Get the mean1 histogram
	TH2* correctedmeanmap1=(TH2*)file->Get("mapping/result/correctedtheta1mean_image");

	// Get the mean2 histogram
	TH2* correctedmeanmap2=(TH2*)file->Get("mapping/result/correctedtheta2mean_image");

	// Get the mean u residual histogram
	TH2* uresidual_meanmap=(TH2*)file->Get("mapping/result/uresidualmean_image");

	// Get the mean v residual histogram
	TH2* vresidual_meanmap=(TH2*)file->Get("mapping/result/vresidualmean_image");

    // Get the momentum 2D histogram
	TH2* mommap=(TH2*)file->Get("mapping/result/BE_image");

    // Get the vertex w mean 2D histogram
	TH2* vertex_w_mean_map=(TH2*)file->Get("mapping/result/vertex_w_image");

    // Get the vertex w rms 2D histogram
	TH2* vertex_w_rms_map=(TH2*)file->Get("mapping/result/vertex_w_rms_image");

    // Get the vertex chi2 2D histogram
	TH2* vertex_chi2_map=(TH2*)file->Get("mapping/result/vertex_chi2_image");

    // Get the vertex multiplicity 2D histogram
	TH2* vertex_multiplicity_map=(TH2*)file->Get("mapping/result/vertex_multiplicity_image");

	// Get the mean u residual (vertex  vs track) histogram
	TH2* u_res_mean_vtx_trk_map=(TH2*)file->Get("mapping/result/u_res_mean_vtx_trk_image");

	// Get the mean u residual (vertex  vs track) histogram
	TH2* v_res_mean_vtx_trk_map=(TH2*)file->Get("mapping/result/v_res_mean_vtx_trk_image");

	// Get the u residual rms (vertex  vs track) histogram
	TH2* u_res_rms_vtx_trk_map=(TH2*)file->Get("mapping/result/u_res_rms_vtx_trk_image");

	// Get the u residual rms (vertex  vs track) histogram
	TH2* v_res_rms_vtx_trk_map=(TH2*)file->Get("mapping/result/v_res_rms_vtx_trk_image");

	// Fill both maps containing the theta means
	meanmap1->SetBinContent(col+1,row+1,uncorrected_mean1);
	meanmap2->SetBinContent(col+1,row+1,uncorrected_mean2);

	// Fill both maps containing the theta means
	correctedmeanmap1->SetBinContent(col+1,row+1,mean1);
	correctedmeanmap2->SetBinContent(col+1,row+1,mean2);

	// Fill both maps containing the residual means
	uresidual_meanmap->SetBinContent(col+1,row+1,uresidual_mean*1E3);
	vresidual_meanmap->SetBinContent(col+1,row+1,vresidual_mean*1E3);

    chi2ndof1map->SetBinContent(col+1,row+1,chi2ndof1);
	chi2ndof2map->SetBinContent(col+1,row+1,chi2ndof2);
	chi2ndofsummap->SetBinContent(col+1,row+1,chi2ndofsum);

    // Fill the momentum image
	mommap->SetBinContent(col+1,row+1,parameters[0]);

    // Fill the vertex images
	vertex_w_mean_map->SetBinContent(col+1,row+1,vertex_w_mean);
	vertex_w_rms_map->SetBinContent(col+1,row+1,vertex_w_rms);
	vertex_chi2_map->SetBinContent(col+1,row+1,vertex_chi2);
	vertex_multiplicity_map->SetBinContent(col+1,row+1,vertex_multiplicity);
	// Fill both maps containing the residual means and rms
	u_res_mean_vtx_trk_map->SetBinContent(col+1,row+1,vtx_trk_res_u_mean);
	v_res_mean_vtx_trk_map->SetBinContent(col+1,row+1,vtx_trk_res_v_mean);
	u_res_rms_vtx_trk_map->SetBinContent(col+1,row+1,vtx_trk_res_u_rms);
	v_res_rms_vtx_trk_map->SetBinContent(col+1,row+1,vtx_trk_res_v_rms);

	double X0;
	double X0err;

	// Fill histograms with mean X0 (in %) and mean absolute X0 error (in %) of both fits


	// Fitting type:
	// 0: gaussian fit function with cuts on the tails, use the two distributions seperately to fit X/X0 and calculate the sum
	// 1: gaussian fit function with cuts on the tails, use the merged distribution to fit X/X0
	//		
	if(fittype==0)
	{
		X0=100.0*(XX01+XX02)/2.0;
	    X0err=100.0*sqrt(pow(XX0err1,2)+pow(XX0err2,2))/2.0;

	}


    // Merged theta distribution fit
	else
	{
		X0=100.0*XX0sum;
	    X0err=100.0*XX0errsum;

	}

	// Get the X0 histogram
	TH2* X0map=(TH2*)file->Get("mapping/result/x0_image");
	cout<<"X0 in area "<<col<<","<<row<<" is: "<<X0<<"%"<<endl;
	X0map->SetBinContent(col+1,row+1,X0);

	// Get the X0 fit error histogram
	TH2* X0errmap=(TH2*)file->Get("mapping/result/x0err_image");
	cout<<"X0 error in area "<<col<<","<<row<<" is: "<<X0err<<"%"<<endl;
	X0errmap->SetBinContent(col+1,row+1,X0err);

	// Get the relative X0 fit error histogram
	TH2* X0relerrmap=(TH2*)file->Get("mapping/result/x0relerr_image");

	double X0relerr;

	// Fill it with mean relative X0 error (in %) of both fits
	if(X0!=0.0)
	{
		X0relerr=X0err/X0*100.0;
	}

	else
	{
		X0relerr=99.0;
	}
	X0relerrmap->SetBinContent(col+1,row+1,X0relerr);

	// Get the # tracks map
	TH2* nummap=(TH2*)file->Get("mapping/result/beamspot");

	// Fill it with mean prob of both fits
	nummap->SetBinContent(col+1,row+1,NumberOfTracks);

	// Get the prob map (2D) histogram (theta1)
	TH2* probmap1=(TH2*)file->Get("mapping/result/fit1prob_image");

	// Get the prob histogram (theta1)
	TH1* probhisto1=(TH1*)file->Get("mapping/result/fit1prob_histo");

	probmap1->SetBinContent(col+1,row+1,prob1);
	probhisto1->Fill(prob1);

	// Get the prob map (2D) histogram (theta2)
	TH2* probmap2=(TH2*)file->Get("mapping/result/fit2prob_image");

	// Get the prob histogram (theta2)
	TH1* probhisto2=(TH1*)file->Get("mapping/result/fit2prob_histo");

	probmap2->SetBinContent(col+1,row+1,prob2);
	probhisto2->Fill(prob2);

	// Get the prob map (2D) histogram (sum)
	TH2* probmapsum=(TH2*)file->Get("mapping/result/fitsumprob_image");

	// Get the prob histogram (sum)
	TH1* probhistosum=(TH1*)file->Get("mapping/result/fitsumprob_histo");

	probmapsum->SetBinContent(col+1,row+1,probsum);
	probhistosum->Fill(probsum);

	// Write the histogram to the rootfile
	if(col==numcol-1&&row==numrow-1)
	{

		chi2ndof1map->Write();
		chi2ndof2map->Write();
		chi2ndofsummap->Write();

		probmap1->Write();
		probmap2->Write();
		probmapsum->Write();

		probhisto1->Write();
		probhisto2->Write();
		probhistosum->Write();

		X0map->Write();
		X0errmap->Write();
        X0relerrmap->Write();

		meanmap1->Write();
		meanmap2->Write();
		correctedmeanmap1->Write();
		correctedmeanmap2->Write();

		uresidual_meanmap->Write();
		vresidual_meanmap->Write();

		nummap->Write();
		mommap->Write();

		vertex_w_mean_map->Write();
		vertex_w_rms_map->Write();
		vertex_chi2_map->Write();
		vertex_multiplicity_map->Write();
		u_res_mean_vtx_trk_map->Write();
		v_res_mean_vtx_trk_map->Write();
		u_res_rms_vtx_trk_map->Write();
		v_res_rms_vtx_trk_map->Write();

	}

	//delete all histograms from memory
	fithistogram1->SetDirectory(gROOT);
	delete fithistogram1;

	fithistogram2->SetDirectory(gROOT);
	delete fithistogram2;

	fithistogramsum->SetDirectory(gROOT);
	delete fithistogramsum;

	//delete all histograms from memory
	histogram1->SetDirectory(gROOT);
	delete histogram1;

	histogram2->SetDirectory(gROOT);
	delete histogram2;

	//delete all histograms from memory
	histogramu->SetDirectory(gROOT);
	delete histogramu;

	histogramv->SetDirectory(gROOT);
	delete histogramv;

	histogramsum->SetDirectory(gROOT);
	delete histogramsum;

    histogram_vertex_w->SetDirectory(gROOT);
    delete histogram_vertex_w;

    histogram_vertex_chi2->SetDirectory(gROOT);
    delete histogram_vertex_chi2;

    histogram_resu_vtx_trk->SetDirectory(gROOT);
    delete histogram_resu_vtx_trk;

    histogram_resv_vtx_trk->SetDirectory(gROOT);
    delete histogram_resv_vtx_trk;

  }

// Function, which returns beam momentum value for every point on the target plane
// This is necessary because the beam profile at DESY often have beam energy gradients in the order of a few MeV/mm
Double_t GetMomentum(double meanvalue,double ugrad,double vgrad, double u, double v)
{
	double p;
	p=meanvalue+u*ugrad+v*vgrad;
	return p;
}
  


// This script is used to create a map of a plane in a test beam telescope. The input is a TTree including 
// MSC projected scattering angle distributions and reconstruction errors.
int x0imaging()
{
	gSystem->Load("libProof.so");
	gSystem->Load("libTreePlayer.so");

	gROOT->Reset(); 

	// display mode
	gStyle->SetPalette(1);
	gStyle->SetOptStat(11111111);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);

    gROOT->ForceStyle();

    // Read config file
    //------------------
    TEnv mEnv("x0image-partial.cfg");

    // Read calibration results config file
    //------------------
    TEnv mEnv_res("x0cal_result.cfg");

	// TString for the input root file name
	TString histoname,range;
	TString filename=mEnv.GetValue("x0filename", "X0-merge");

	// Copy the X0 Analysis Root file 
	TFile *X0file = new TFile(filename, "READ");

	//Open the copied file
	filename=filename+mEnv.GetValue("x0fileidentifier", "-part-1-1");
	TFile *rootfile = new TFile(filename+".root", "RECREATE");

	// Create directories containing map histograms, fits and results
	rootfile->mkdir("mapping");
	rootfile->mkdir("mapping/raw");
	rootfile->mkdir("mapping/fit");
	rootfile->mkdir("mapping/badfit");
	rootfile->mkdir("mapping/result");

	// Number of Rows and Columns of the sensor map
	int numcol = mEnv.GetValue("maxupixels", 100);
	int numrow = mEnv.GetValue("maxvpixels", 50);

	// u minimum and v maximum value (in mm)
	double umin=mEnv.GetValue("umin", -10.0);
	double vmax=mEnv.GetValue("vmax", 5.0);

	// u and v length of the map (in mm)
	double ulength=mEnv.GetValue("ulength", 20.0);
	double vlength=mEnv.GetValue("vlength", -10.0);

	// u and v pitch (length of one pixel in the map, in mm)
	double upitch=ulength/(1.0*numcol);
	double vpitch=vlength/(1.0*numrow);

	// umax and vmin can be computed from the parameters already implemented 
	double umax=umin+ulength;
	double vmin=vmax-vlength;

    // calculate the u value of the center of the image
    double u_center=umin+0.5*ulength;

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

    // Vertex multiplicity cut (should be 1 for default X0 analysis)
	int vertex_multiplicity_min=mEnv.GetValue("vertexmultiplicitymin", 1);
	int vertex_multiplicity_max=mEnv.GetValue("vertexmultiplicitymax", 1);

	cout<<"Minimal vertex multiplicity:"<<vertex_multiplicity_min<<endl;
	cout<<"Maximal vertex multiplicity:"<<vertex_multiplicity_max<<endl;

    // Vertex multiplicity cut (should be 1 for default X0 analysis)
	double maxchi2ndof_fit=mEnv.GetValue("maxchi2ndof", 10.0);

    // Fit range parameter
	double rangevalue=mEnv.GetValue("fitrange_parameter", 2.0);

    // Fit options
	TString fitoptions=mEnv.GetValue("fit_options", "RMELS");

	// Choose the type of fit
	// 0: gaussian fit function with cuts on the tails, both kink distributions are used seperately
	// 1: gaussian fit function with cuts on the tails, use only 1 fit on the merged histogram consisting of both distributions
	int fittype=1;

	TTree * tree = (TTree*) X0file->Get("MSCTree");

 	std::vector<double> plotranges,means;

	// Draw theta1 histogram
	tree->Draw("theta1", "", "P*");

	// Get plot range and mean value for first iteration of theta1 histograms
	plotranges.push_back(4*tree->GetHistogram()->GetRMS());
	means.push_back(tree->GetHistogram()->GetMean());

	// Draw theta2 histogram
	tree->Draw("theta2", "", "P*");

	// Get plot range and mean value for first iteration of theta1 histograms
	plotranges.push_back(4*tree->GetHistogram()->GetRMS());
	means.push_back(tree->GetHistogram()->GetMean());
	
	// The number of bins of angle histograms
	int numberofbins=mEnv.GetValue("num_bins", 50);

	// Range parameter of angle histograms
	double histo_range=mEnv.GetValue("histo_range", 5.0);

	cout<<"The first Scattering angle distributions will be plotted with range "<<plotranges.at(0)<<" rad and "<<numberofbins<<" bins!"<<endl;
	cout<<"The mean value of the histogram is "<<means.at(0)<<" rad !"<<endl;

	cout<<"The second Scattering angle distributions will be plotted with range "<<plotranges.at(1)<<" rad and "<<numberofbins<<" bins!"<<endl;
	cout<<"The mean value of the histogram is "<<means.at(1)<<" rad !"<<endl;

	// Calibration factor lambda, used to change the reconstruction error to include systematical errors
	// The calibration factor is either taken from the x0 calibration results cfg file or
	// from the cfg, where it must be inserted manually
	double lambda_default=mEnv.GetValue("lambda", 1.0);
	double lambda=mEnv_res.GetValue("lambda_start", lambda_default);
	double recoerror=sqrt(getanglerecovar(X0file))*lambda;
	cout<<"The reconstruction error is "<<recoerror*1E6<<" µrad!"<<endl;
	cout<<"This includes the calibration factor of "<<lambda<<endl;

	// Beam energy in GeV
    // The particle momenta are distributed due to the beam generation at the desy facility
    // Particles of different momenta are sperated by a dipole magnet and a collimator. Due to
    // the finite size of the collimator opening particles with momenta in a certain range p0+/-delta_p
    // can traverse into the test beam area
    // We expect a linear distribution with slope corresponding to ~500MeV/20mm
	// Take either results from x0 calibration or manually inserted values from cfg file
	double mom0_default=mEnv.GetValue("momentumoffset", 4.0);           // in GeV
	double mom0=mEnv_res.GetValue("momentumoffset", mom0_default);
    double mom_uslope_default=mEnv.GetValue("momentumugradient", 0.0);    // in GeV/mm
	double mom_uslope=mEnv_res.GetValue("momentumugradient", mom_uslope_default);
    double mom_vslope_default=mEnv.GetValue("momentumvgradient", 0.0);    // in GeV/mm
	double mom_vslope=mEnv_res.GetValue("momentumvgradient", mom_vslope_default);
	cout<<"The beam energy is "<<mom0<<" GeV!"<<endl;
	cout<<"The beam energy gradient (u direction) is "<<mom_uslope<<" GeV/mm!"<<endl;
	cout<<"The beam energy gradient (v direction) is "<<mom_vslope<<" GeV/mm!"<<endl;


	//particle charge
	double charge=1;

	//particle mass
	double mass=0.000511;

	rootfile->cd("");
	rootfile->cd("mapping/result");

	// Scatter theta1 vs residual u
	TH2F * hscatt_theta1_vs_resu = new TH2F("hscatt_theta1_vs_resu","hscatt_theta1_vs_resu",100,-0.1,0.1,150,means.at(0)-plotranges.at(0),means.at(0)+plotranges.at(0));
    hscatt_theta1_vs_resu->SetStats(kFALSE);
    hscatt_theta1_vs_resu->GetXaxis()->SetTitle("u residual[mm]");
    hscatt_theta1_vs_resu->GetYaxis()->SetTitle("theta1[rad]");
    hscatt_theta1_vs_resu->GetZaxis()->SetTitle("Number of tracks");
    hscatt_theta1_vs_resu->GetZaxis()->SetTitleOffset(1.4);
    hscatt_theta1_vs_resu->GetZaxis()->SetTitleSize(0.02);
    hscatt_theta1_vs_resu->GetZaxis()->SetLabelSize(0.02);

	// Scatter theta2 vs residual u
	TH2F * hscatt_theta2_vs_resu = new TH2F("hscatt_theta2_vs_resu","hscatt_theta2_vs_resu",100,-0.1,0.1,150,means.at(1)-plotranges.at(1),means.at(1)+plotranges.at(1));
    hscatt_theta2_vs_resu->SetStats(kFALSE);
    hscatt_theta2_vs_resu->GetXaxis()->SetTitle("u residual[mm]");
    hscatt_theta2_vs_resu->GetYaxis()->SetTitle("theta2[rad]");
    hscatt_theta2_vs_resu->GetZaxis()->SetTitle("Number of tracks");
    hscatt_theta2_vs_resu->GetZaxis()->SetTitleOffset(1.4);
    hscatt_theta2_vs_resu->GetZaxis()->SetTitleSize(0.02);
    hscatt_theta2_vs_resu->GetZaxis()->SetLabelSize(0.02);

	// Scatter theta1 vs residual v
	TH2F * hscatt_theta1_vs_resv = new TH2F("hscatt_theta1_vs_resv","hscatt_theta1_vs_resv",100,-0.1,0.1,150,means.at(0)-plotranges.at(0),means.at(0)+plotranges.at(0));
    hscatt_theta1_vs_resv->SetStats(kFALSE);
    hscatt_theta1_vs_resv->GetXaxis()->SetTitle("v residual[mm]");
    hscatt_theta1_vs_resv->GetYaxis()->SetTitle("theta1[rad]");
    hscatt_theta1_vs_resv->GetZaxis()->SetTitle("Number of tracks");
    hscatt_theta1_vs_resv->GetZaxis()->SetTitleOffset(1.4);
    hscatt_theta1_vs_resv->GetZaxis()->SetTitleSize(0.02);
    hscatt_theta1_vs_resv->GetZaxis()->SetLabelSize(0.02);

	// Scatter theta2 vs residual v
	TH2F * hscatt_theta2_vs_resv = new TH2F("hscatt_theta2_vs_resv","hscatt_theta2_vs_resv",100,-0.1,0.1,150,means.at(1)-plotranges.at(1),means.at(1)+plotranges.at(1));
    hscatt_theta2_vs_resv->SetStats(kFALSE);
    hscatt_theta2_vs_resv->GetXaxis()->SetTitle("v residual[mm]");
    hscatt_theta2_vs_resv->GetYaxis()->SetTitle("theta2[rad]");
    hscatt_theta2_vs_resv->GetZaxis()->SetTitle("Number of tracks");
    hscatt_theta2_vs_resv->GetZaxis()->SetTitleOffset(1.4);
    hscatt_theta2_vs_resv->GetZaxis()->SetTitleSize(0.02);
    hscatt_theta2_vs_resv->GetZaxis()->SetLabelSize(0.02);



	tree->Draw("theta1:(u_out-u_in)>>hscatt_theta1_vs_resu","","colz");
	hscatt_theta1_vs_resu->Write();

	tree->Draw("theta2:(u_out-u_in)>>hscatt_theta2_vs_resu","","colz");
	hscatt_theta2_vs_resu->Write();

	tree->Draw("theta1:(v_out-v_in)>>hscatt_theta1_vs_resv","","colz");
	hscatt_theta1_vs_resv->Write();

	tree->Draw("theta2:(v_out-v_in)>>hscatt_theta2_vs_resv","","colz");
	hscatt_theta2_vs_resv->Write();

	rootfile->cd("");
	
	getcorrection(X0file, rootfile, means, plotranges, numberofbins, numcol, numrow, umin, vmin, umax, vmax, vertex_multiplicity_min, vertex_multiplicity_max);
	savehistos(X0file, rootfile, numberofbins, histo_range, numcol, numrow, umin, vmin, umax, vmax, vertex_multiplicity_min, vertex_multiplicity_max);

	X0file->Close();

	rootfile->cd("");
	rootfile->cd("mapping/result");

	// X0 map
	TH2F * x0_image = new TH2F("x0_image","x0_image",numcol,umin,umax,numrow,vmin,vmax);
    x0_image->SetStats(kFALSE);
    x0_image->SetMinimum(0);
    x0_image->GetXaxis()->SetTitle("u [mm]");
    x0_image->GetYaxis()->SetTitle("v [mm]");
    x0_image->GetZaxis()->SetTitle("X/X0 [%]");
    x0_image->GetZaxis()->SetTitleSize(0.02);
    x0_image->GetZaxis()->SetLabelSize(0.02);

	// X0 statistical error map (absolute value)
	TH2F * x0err_image = new TH2F("x0err_image","x0err_image",numcol,umin,umax,numrow,vmin,vmax);
    x0err_image->SetStats(kFALSE);
    x0err_image->SetMinimum(0);
    x0err_image->GetXaxis()->SetTitle("u [mm]");
    x0err_image->GetYaxis()->SetTitle("v [mm]");
    x0err_image->GetZaxis()->SetTitle("X/X0 [%]");
    x0err_image->GetZaxis()->SetTitleSize(0.02);
    x0err_image->GetZaxis()->SetLabelSize(0.02);

	// X0 statistical error map (relative value)
	TH2F * x0relerr_image = new TH2F("x0relerr_image","x0relerr_image",numcol,umin,umax,numrow,vmin,vmax);
    x0relerr_image->SetStats(kFALSE);
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

    // Fit Probability distribution for first angle dist
	TH1F * fit1prob_histo = new TH1F("fit1prob_histo","fit1prob_histo",50,0.0,1.0);
	fit1prob_histo->SetStats(kFALSE);
    fit1prob_histo->GetXaxis()->SetTitle("fit1 p value");
    fit1prob_histo->GetYaxis()->SetTitle("number of fits");

	// Fit probability map for first angle dist
	TH2F * fit1prob_image = new TH2F("fit1prob_image","fit1prob_image",numcol,umin,umax,numrow,vmin,vmax);
	fit1prob_image->SetStats(kFALSE);
    fit1prob_image->GetXaxis()->SetTitle("u [mm]");
    fit1prob_image->GetYaxis()->SetTitle("v [mm]");
    fit1prob_image->SetMinimum(-0.1);
    fit1prob_image->GetZaxis()->SetTitle("fit1 p value");
    fit1prob_image->GetZaxis()->SetTitleSize(0.02);
    fit1prob_image->GetZaxis()->SetLabelSize(0.02);


	// Fit Probability distribution for second angle dist
	TH1F * fit2prob_histo = new TH1F("fit2prob_histo","fit2prob_histo",50,0.0,1.0);
	fit2prob_histo->SetStats(kFALSE);
    fit2prob_histo->GetXaxis()->SetTitle("fit2 p value");
    fit2prob_histo->GetYaxis()->SetTitle("number of fits");

	// Fit probability map for second angle dist
	TH2F * fit2prob_image = new TH2F("fit2prob_image","fit2prob_image",numcol,umin,umax,numrow,vmin,vmax);
	fit2prob_image->SetStats(kFALSE);
    fit2prob_image->GetXaxis()->SetTitle("u [mm]");
    fit2prob_image->GetYaxis()->SetTitle("v [mm]");
    fit2prob_image->SetMinimum(-0.1);
    fit2prob_image->GetZaxis()->SetTitle("fit2 p value");
    fit2prob_image->GetZaxis()->SetTitleSize(0.02);
    fit2prob_image->GetZaxis()->SetLabelSize(0.02);

	// Fit Probability distribution for merged angle dist
	TH1F * fitsumprob_histo = new TH1F("fitsumprob_histo","fitsumprob_histo",50,0.0,1.0);
	fitsumprob_histo->SetStats(kFALSE);
    fitsumprob_histo->GetXaxis()->SetTitle("fitsum p value");
    fitsumprob_histo->GetYaxis()->SetTitle("number of fits");

	// Fit probability map for merged angle dist
	TH2F * fitsumprob_image = new TH2F("fitsumprob_image","fitsumprob_image",numcol,umin,umax,numrow,vmin,vmax);
	fitsumprob_image->SetStats(kFALSE);
    fitsumprob_image->GetXaxis()->SetTitle("u [mm]");
    fitsumprob_image->GetYaxis()->SetTitle("v [mm]");
    fitsumprob_image->SetMinimum(-0.1);
    fitsumprob_image->GetZaxis()->SetTitle("fitsum p value");
    fitsumprob_image->GetZaxis()->SetTitleSize(0.02);
    fitsumprob_image->GetZaxis()->SetLabelSize(0.02);


	// Fit mean value of first scattering angle distribution
	TH2F * theta1mean_image = new TH2F("theta1mean_image","theta1mean_image",numcol,umin,umax,numrow,vmin,vmax);
	theta1mean_image->SetStats(kFALSE);
    theta1mean_image->GetXaxis()->SetTitle("u [mm]");
    theta1mean_image->GetYaxis()->SetTitle("v [mm]");
    theta1mean_image->GetZaxis()->SetTitle("theta1 mean value[rad]");
    theta1mean_image->GetZaxis()->SetTitleSize(0.02);
    theta1mean_image->GetZaxis()->SetLabelSize(0.02);


    // Fit mean value of second scattering angle distribution
	TH2F * theta2mean_image = new TH2F("theta2mean_image","theta2mean_image",numcol,umin,umax,numrow,vmin,vmax);
	theta2mean_image->SetStats(kFALSE);
    theta2mean_image->GetXaxis()->SetTitle("u [mm]");
    theta2mean_image->GetYaxis()->SetTitle("v [mm]");
    theta2mean_image->GetZaxis()->SetTitle("theta2 mean value[rad]");
    theta2mean_image->GetZaxis()->SetTitleSize(0.02);
    theta2mean_image->GetZaxis()->SetLabelSize(0.02);

	// Fit mean value of centered first scattering angle distribution
	TH2F * correctedtheta1mean_image = new TH2F("correctedtheta1mean_image","correctedtheta1mean_image",numcol,umin,umax,numrow,vmin,vmax);
	correctedtheta1mean_image->SetStats(kFALSE);
    correctedtheta1mean_image->GetXaxis()->SetTitle("u [mm]");
    correctedtheta1mean_image->GetYaxis()->SetTitle("v [mm]");
    correctedtheta1mean_image->GetZaxis()->SetTitle("theta1 mean value[rad]");
    correctedtheta1mean_image->GetZaxis()->SetTitleSize(0.02);
    correctedtheta1mean_image->GetZaxis()->SetLabelSize(0.02);

	// Fit mean value of centered second scattering angle distribution
	TH2F * correctedtheta2mean_image = new TH2F("correctedtheta2mean_image","correctedtheta2mean_image",numcol,umin,umax,numrow,vmin,vmax);
	correctedtheta2mean_image->SetStats(kFALSE);
    correctedtheta2mean_image->GetXaxis()->SetTitle("u [mm]");
    correctedtheta2mean_image->GetYaxis()->SetTitle("v [mm]");
    correctedtheta2mean_image->GetZaxis()->SetTitle("theta2 mean value[rad]");
    correctedtheta2mean_image->GetZaxis()->SetTitleSize(0.02);
    correctedtheta2mean_image->GetZaxis()->SetLabelSize(0.02);


	// u residuals of downstream and upstream estimation
	TH2F * uresidualmean_image = new TH2F("uresidualmean_image","uresidualmean_image",numcol,umin,umax,numrow,vmin,vmax);
	uresidualmean_image->SetStats(kFALSE);
    uresidualmean_image->GetXaxis()->SetTitle("u [mm]");
    uresidualmean_image->GetYaxis()->SetTitle("v [mm]");
    uresidualmean_image->GetZaxis()->SetTitle("u residual[µm]");
    uresidualmean_image->GetZaxis()->SetTitleSize(0.02);
    uresidualmean_image->GetZaxis()->SetLabelSize(0.02);

	// v residuals of downstream and upstream estimation
	TH2F * vresidualmean_image = new TH2F("vresidualmean_image","vresidualmean_image",numcol,umin,umax,numrow,vmin,vmax);
	vresidualmean_image->SetStats(kFALSE);
    vresidualmean_image->GetXaxis()->SetTitle("u [mm]");
    vresidualmean_image->GetYaxis()->SetTitle("v [mm]");
    vresidualmean_image->GetZaxis()->SetTitle("v residual[µm]");
    vresidualmean_image->GetZaxis()->SetTitleSize(0.02);
    vresidualmean_image->GetZaxis()->SetLabelSize(0.02);
	
	// Beam spot image from track intersections
	TH2F * beamspot = new TH2F("beamspot","beamspot",numcol,umin,umax,numrow,vmin,vmax);
    beamspot->SetStats(kFALSE);
    beamspot->GetXaxis()->SetTitle("u [mm]");
    beamspot->GetYaxis()->SetTitle("v [mm]");
    beamspot->GetZaxis()->SetTitle("number of tracks");
    beamspot->GetZaxis()->SetTitleOffset(1.4);
    beamspot->GetZaxis()->SetTitleSize(0.02);
    beamspot->GetZaxis()->SetLabelSize(0.02);

	// Momentum image
	TH2F * BE_image = new TH2F("BE_image","BE_image",numcol,umin,umax,numrow,vmin,vmax);
    BE_image->SetStats(kFALSE);
    BE_image->GetXaxis()->SetTitle("u [mm]");
    BE_image->GetYaxis()->SetTitle("v [mm]");
    BE_image->GetZaxis()->SetTitle("momentum [GeV/c]");
    BE_image->GetZaxis()->SetTitleOffset(1.4);
    BE_image->GetZaxis()->SetTitleSize(0.02);
    BE_image->GetZaxis()->SetLabelSize(0.02);

	// Mean Vertex w position image
	TH2F * vertex_w_image = new TH2F("vertex_w_image","vertex_w_image",numcol,umin,umax,numrow,vmin,vmax);
    vertex_w_image->SetStats(kFALSE);
    vertex_w_image->GetXaxis()->SetTitle("u [mm]");
    vertex_w_image->GetYaxis()->SetTitle("v [mm]");
    vertex_w_image->GetZaxis()->SetTitle("vertex w [mm]");
    vertex_w_image->GetZaxis()->SetTitleOffset(1.4);
    vertex_w_image->GetZaxis()->SetTitleSize(0.02);
    vertex_w_image->GetZaxis()->SetLabelSize(0.02);

	// Vertex w position RMS image
	TH2F * vertex_w_rms_image = new TH2F("vertex_w_rms_image","vertex_w_rms_image",numcol,umin,umax,numrow,vmin,vmax);
    vertex_w_rms_image->SetStats(kFALSE);
    vertex_w_rms_image->GetXaxis()->SetTitle("u [mm]");
    vertex_w_rms_image->GetYaxis()->SetTitle("v [mm]");
    vertex_w_rms_image->GetZaxis()->SetTitle("vertex w RMS [mm]");
    vertex_w_rms_image->GetZaxis()->SetTitleOffset(1.4);
    vertex_w_rms_image->GetZaxis()->SetTitleSize(0.02);
    vertex_w_rms_image->GetZaxis()->SetLabelSize(0.02);

	// Mean Vertex multiplicity image
	TH2F * vertex_multiplicity_image = new TH2F("vertex_multiplicity_image","vertex_multiplicity_image",numcol,umin,umax,numrow,vmin,vmax);
    vertex_multiplicity_image->SetStats(kFALSE);
    vertex_multiplicity_image->GetXaxis()->SetTitle("u [mm]");
    vertex_multiplicity_image->GetYaxis()->SetTitle("v [mm]");
    vertex_multiplicity_image->GetZaxis()->SetTitle("mean vertex multiplicity");
    vertex_multiplicity_image->GetZaxis()->SetTitleOffset(1.4);
    vertex_multiplicity_image->GetZaxis()->SetTitleSize(0.02);
    vertex_multiplicity_image->GetZaxis()->SetLabelSize(0.02);

	// Mean Vertex chi2 image
	TH2F * vertex_chi2_image = new TH2F("vertex_chi2_image","vertex_chi2_image",numcol,umin,umax,numrow,vmin,vmax);
    vertex_chi2_image->SetStats(kFALSE);
    vertex_chi2_image->GetXaxis()->SetTitle("u [mm]");
    vertex_chi2_image->GetYaxis()->SetTitle("v [mm]");
    vertex_chi2_image->GetZaxis()->SetTitle("mean vertex chi2");
    vertex_chi2_image->GetZaxis()->SetTitleOffset(1.4);
    vertex_chi2_image->GetZaxis()->SetTitleSize(0.02);
    vertex_chi2_image->GetZaxis()->SetLabelSize(0.02);

	// Residual mean of Vertex position u and weighted means intersection from down and upstream track
	TH2F * u_res_mean_vtx_trk_image = new TH2F("u_res_mean_vtx_trk_image","u_res_mean_vtx_trk_image",numcol,umin,umax,numrow,vmin,vmax);
    u_res_mean_vtx_trk_image->SetStats(kFALSE);
    u_res_mean_vtx_trk_image->GetXaxis()->SetTitle("u [mm]");
    u_res_mean_vtx_trk_image->GetYaxis()->SetTitle("v [mm]");
    u_res_mean_vtx_trk_image->GetZaxis()->SetTitle("mean u residual vtx-trk [mm]");
    u_res_mean_vtx_trk_image->GetZaxis()->SetTitleOffset(1.4);
    u_res_mean_vtx_trk_image->GetZaxis()->SetTitleSize(0.02);
    u_res_mean_vtx_trk_image->GetZaxis()->SetLabelSize(0.02);

	// Residual mean of Vertex position v and weighted means intersection from down and upstream track
	TH2F * v_res_mean_vtx_trk_image = new TH2F("v_res_mean_vtx_trk_image","v_res_mean_vtx_trk_image",numcol,umin,umax,numrow,vmin,vmax);
    v_res_mean_vtx_trk_image->SetStats(kFALSE);
    v_res_mean_vtx_trk_image->GetXaxis()->SetTitle("u [mm]");
    v_res_mean_vtx_trk_image->GetYaxis()->SetTitle("v [mm]");
    v_res_mean_vtx_trk_image->GetZaxis()->SetTitle("mean v residual vtx-trk [mm]");
    v_res_mean_vtx_trk_image->GetZaxis()->SetTitleOffset(1.4);
    v_res_mean_vtx_trk_image->GetZaxis()->SetTitleSize(0.02);
    v_res_mean_vtx_trk_image->GetZaxis()->SetLabelSize(0.02);

	// Residual rms of Vertex position u and weighted means intersection from down and upstream track
	TH2F * u_res_rms_vtx_trk_image = new TH2F("u_res_rms_vtx_trk_image","u_res_rms_vtx_trk_image",numcol,umin,umax,numrow,vmin,vmax);
    u_res_rms_vtx_trk_image->SetStats(kFALSE);
    u_res_rms_vtx_trk_image->GetXaxis()->SetTitle("u [mm]");
    u_res_rms_vtx_trk_image->GetYaxis()->SetTitle("v [mm]");
    u_res_rms_vtx_trk_image->GetZaxis()->SetTitle("u residual rms vtx-trk [mm]");
    u_res_rms_vtx_trk_image->GetZaxis()->SetTitleOffset(1.4);
    u_res_rms_vtx_trk_image->GetZaxis()->SetTitleSize(0.02);
    u_res_rms_vtx_trk_image->GetZaxis()->SetLabelSize(0.02);

	// Residual rms of Vertex position v and weighted means intersection from down and upstream track
	TH2F * v_res_rms_vtx_trk_image = new TH2F("v_res_rms_vtx_trk_image","v_res_rms_vtx_trk_image",numcol,umin,umax,numrow,vmin,vmax);
    v_res_rms_vtx_trk_image->SetStats(kFALSE);
    v_res_rms_vtx_trk_image->GetXaxis()->SetTitle("u [mm]");
    v_res_rms_vtx_trk_image->GetYaxis()->SetTitle("v [mm]");
    v_res_rms_vtx_trk_image->GetZaxis()->SetTitle("v residual rms vtx-trk [mm]");
    v_res_rms_vtx_trk_image->GetZaxis()->SetTitleOffset(1.4);
    v_res_rms_vtx_trk_image->GetZaxis()->SetTitleSize(0.02);
    v_res_rms_vtx_trk_image->GetZaxis()->SetLabelSize(0.02);


	for(int col=0; col<numcol; col++)
	{
		for(int row=0; row<numrow; row++)
		{

			cout<<"fit histogram in (col,row): ("<<col<<","<<row<<")"<<endl;
			cout<<"Fit range value: "<<rangevalue<<endl;

            // Calculate the u position of this bin
            double u=umin+col*upitch;
            double v=vmin+row*vpitch;

           	// Determine the momentum value from the u position and the p distribution parameters
           	double mom=GetMomentum(mom0, mom_uslope, mom_vslope, u, v);

			/*	* par[0]:  Expected beam energy
				* par[1]:  Beam particle charge
				* par[2]:  Beam particle mass
				* par[3]:  Radiation length X/X0
				* par[4]:  Calibrated angle reconstruction error
			*/

			double parameters[7]={mom,charge,mass,0.01,recoerror,300,0.0};

			// fit the histograms
			fithisto(rootfile, fittype, maxchi2ndof_fit, rangevalue, col,numcol,row,numrow,parameters,fitoptions);
		}
	}

	rootfile->Close();

    return 0;

	
}

