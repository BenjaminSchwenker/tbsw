//#include <iostream.h>
#include <fstream>
using namespace std ;


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
  void getcorrection(TFile* file1, TFile* file2, double plotrange, int numberofbins, const int numcol, const int numrow, double umin, double vmin, double umax, double vmax)
  {
	// parameters which are read out from the root file
	Double_t theta1;
	Double_t theta2;
	Double_t u;
	Double_t v;
	
	//TTree in input root file, that contains the MSC projected angle distributions and reconstruction error distribution
	file1->cd("");

	TTree * msc_tree = (TTree*) file1->Get("MSCTree");

	msc_tree->SetBranchAddress("theta1_val",&theta1);
	msc_tree->SetBranchAddress("theta2_val",&theta2);
	msc_tree->SetBranchAddress("u",&u);
	msc_tree->SetBranchAddress("v",&v);

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
		 	histo_theta1[i][j] = new TH1F("","",numberofbins,-plotrange,plotrange);
		 	histo_theta2[i][j] = new TH1F("","",numberofbins,-plotrange,plotrange);
		}
	}

	// Loop over all events
	for(int i=0; i< msc_tree->GetEntries(); i++)
	{

		if(i%500000==0) cout<<"Mean Correction loop, Track No. "<<i<<endl;
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

		// Determine column and row number from the position within the map area and the number of rows and columns
		int col=floor(u_pos*numcol/u_length);
		int row=floor(v_pos*numrow/v_length);

		// Fill histograms
		histo_theta1[col][row]->Fill(theta1);
		histo_theta2[col][row]->Fill(theta2);
	}

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
  void savehistos(TFile* file1, TFile* file2, double plotrange, int numberofbins, const int numcol, const int numrow, double umin, double vmin, double umax, double vmax, int correctmean)
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

	// Array of mean theta1 and theta2 values in each map pixel
	double mean1[numcol][numrow];
	double mean2[numcol][numrow];
	
	//TTree in input root file, that contains the MSC projected angle distributions and reconstruction error distribution
	file1->cd("");

	TTree * msc_tree = (TTree*) file1->Get("MSCTree");

	msc_tree->SetBranchAddress("theta1_val",&theta1);
	msc_tree->SetBranchAddress("theta2_val",&theta2);
	msc_tree->SetBranchAddress("u",&u);
	msc_tree->SetBranchAddress("v",&v);
	msc_tree->SetBranchAddress("u_in",&u_in);
	msc_tree->SetBranchAddress("v_in",&v_in);
	msc_tree->SetBranchAddress("u_out",&u_out);
	msc_tree->SetBranchAddress("v_out",&v_out);

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
		 	histo_theta1[i][j] = new TH1F("","",numberofbins,-plotrange,plotrange);
		 	histo_theta2[i][j] = new TH1F("","",numberofbins,-plotrange,plotrange);
		 	histo_uresidual[i][j] = new TH1F("","",1000,-1.0,1.0);
		 	histo_vresidual[i][j] = new TH1F("","",1000,-1.0,1.0);
		 	histo_thetasum[i][j] = new TH1F("","",numberofbins,-plotrange,plotrange);
			histo_2d[i][j] = new TH2F("","",numberofbins,-plotrange,plotrange,numberofbins,-plotrange,plotrange);

			// Get the histograms generated in the getcorrection function
			// Name of the histograms
			TString aidhistoname;
			aidhistoname.Form("area(%i,%i)",i,j);


			if(correctmean==1)
			{
				// Get histogram
				TH1* histogram1=(TH1*)file2->Get("mapping/raw/theta1_uncorrected_"+aidhistoname);
				TH1* histogram2=(TH1*)file2->Get("mapping/raw/theta2_uncorrected_"+aidhistoname);

				// Save mean of both distributions in arrays
				mean1[i][j]=histogram1->GetMean();
				mean2[i][j]=histogram2->GetMean();

                                histogram1->Delete();
                                histogram2->Delete();
			}

		}
	}

	// Loop over all events
	for(int i=0; i< msc_tree->GetEntries(); i++)
	{

		if(i%500000==0) cout<<"Track No. "<<i<<endl;
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

		// Determine column and row number from the position within the map area and the number of rows and columns
		int col=floor(u_pos*numcol/u_length);
		int row=floor(v_pos*numrow/v_length);

		if(correctmean==1)
		{
			// if the mean correction is enabled change theta accordingly
			theta1=theta1-mean1[col][row];
			theta2=theta2-mean2[col][row];
		}


		// Fill histograms
		histo_theta1[col][row]->Fill(theta1);
		histo_theta2[col][row]->Fill(theta2);

		histo_uresidual[col][row]->Fill(u_in-u_out);
		histo_vresidual[col][row]->Fill(v_in-v_out);

		histo_thetasum[col][row]->Fill(theta1);
		histo_thetasum[col][row]->Fill(theta2);

		histo_2d[col][row]->Fill(theta1,theta2);
	}

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
		}
		
	}

  }


  double calculateX0(double sigma, double recoerror, double mom, int model)
  {
	// Initialize X0min and X0max values
	double X0min=0.0001;
	double X0max=0.5;

	// Iniatialize XX0 value as mean between starting values of lower and upper XX0 bound
	double X0=(X0min+X0max)/2.0;

	double theta_stddev=sqrt(sigma*sigma-recoerror*recoerror);

	// Calculate difference between measured stddev and stddev(X0) with this starting X/X0 value
	double difference;
	if(model==0)
	{
		difference=theta_stddev-0.0136/mom*sqrt(X0)*(1+0.038*log(X0));
	}

	if(model==1)
	{
		difference=theta_stddev-0.0150/mom*sqrt(X0)*(0.9184+0.0361*log(X0));
	}

	if(model==2)
	{
		double Z=13;
		double hs=(Z+1)/Z*log(287*pow(Z,-0.5))/log(159*pow(Z,-0.333));
		double norm=0.000225*X0/(hs*pow(mom,2));
		difference=theta_stddev-pow(norm,0.5)*sqrt(0.851+0.00314*log(X0/hs)-0.001825*pow(log(X0/hs),2));
	}

	if(model==3)

	{
		double Z=13;
		double hs=(Z+1)/Z*log(287*pow(Z,-0.5))/log(159*pow(Z,-0.333));
		double norm=0.000225*X0/(hs*pow(mom,2));
		difference=theta_stddev-pow(norm,0.5);
	}

	// For loop used to minimize difference
	for(Int_t i=1; sqrt(pow(difference,2))>1E-8; i++)
	{
		// If the difference is larger than 0 the current X/X0 value is too small, therefore the lower bound of X/X0 (X0min)

		// is now set to be the current X/X0 value
		if(difference>0){X0min=X0;}

		// If the difference is smaller than 0 the current X/X0 value is too large, therefore the upper bound of X/X0 (X0max)
		// is now set to be the current X/X0 value
		else{X0max=X0;}

		// Calculate the new X/X0 value from the updated a and b values. Once again its the mean of a and b.

		X0=(X0min+X0max)/2.0;

		// Recalculate F with the updated X/X0 value
		if(model==0)
		{
			difference=theta_stddev-0.0136/mom*sqrt(X0)*(1+0.038*log(X0));

		}

		if(model==1)
		{
			difference=theta_stddev-0.0150/mom*sqrt(X0)*(0.9184+0.0361*log(X0));
		}

		if(model==2)
		{
			double Z=13;
			double hs=(Z+1)/Z*log(287*pow(Z,-0.5))/log(159*pow(Z,-0.333));
			double norm=0.000225*X0/(hs*pow(mom,2));
			difference=theta_stddev-pow(norm,0.5)*sqrt(0.851+0.00314*log(X0/hs)-0.001825*pow(log(X0/hs),2));
		}

		if(model==3)

		{
			double Z=13;
			double hs=(Z+1)/Z*log(287*pow(Z,-0.5))/log(159*pow(Z,-0.333));
			double norm=0.000225*X0/(hs*pow(mom,2));
			difference=theta_stddev-pow(norm,0.5);
		}

		if(i==100)
		{
			X0=0.0;
			break;
		}
	}

	// Return X0 in percent
	X0*=100;

	return X0;
  }

  double getX0err(double mom, double sigma, double sigma_error, double recoerror, int model)
  {
	// There are two cases to consider: 
	// 1) if the fitted standard deviation minus the fitting error is smaller than the reconstruction error only the 
	//    max value can be used to calculate the X0 statistical error
	// 2) if the reconstruction error is smaller than the fitted standard deviation minus the fitting error the statistical
	//    X0 error can be calculated from the min and max X0 values


	// Output of this function
	double X0_error;

	// Definition of the 4 X0 values, which can be calculated at the limits of the statistical error intervalls
	double X0_sigma_max,X0_sigma_min;


	// Calculate the max and min std dev values within the statistical errors
	double sigma_max=sigma+sigma_error;

	double sigma_min;


	// Compute max and min X0 values (within the statistical fit uncertainties) for the angle distribution
	X0_sigma_max=calculateX0(sigma_max, recoerror, mom, model);
	cout<<"X0max: "<<X0_sigma_max<<endl;


	if(sigma-sigma_error>recoerror)
	{
		sigma_min=sigma-sigma_error;

		X0_sigma_min=calculateX0(sigma_min, recoerror, mom, model);
		cout<<"X0min1: "<<X0_sigma_min<<endl;


		// The statistical fit error is given by half the difference between the max and min value
		X0_error=(X0_sigma_max-X0_sigma_min)/2.0;

	}


	else
	{

		// The statistical fit error is given by the difference between the max and mean value
		X0_error=(X0_sigma_max-calculateX0(sigma, recoerror, mom, model));

	}

	return X0_error;
  }


  // Function to fit the MSC angle histograms and fill the map histograms
  void fithisto( TFile* file, int fittype, int model,int col, int numcol, int row,int numrow, double mom, double recoerror)
  {
	// Open the histograms
 
	// histogram name
	TString histoname;
	histoname.Form("area(%i,%i)",col,row);
	file->cd("");
	histogram1=(TH1*)file->Get("mapping/raw/theta1_"+histoname);
	histogram2=(TH1*)file->Get("mapping/raw/theta2_"+histoname);
	histogramu=(TH1*)file->Get("mapping/raw/uresidual_"+histoname);
	histogramv=(TH1*)file->Get("mapping/raw/vresidual_"+histoname);
	histogramsum=(TH1*)file->Get("mapping/raw/sumhisto_"+histoname);

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

	// std deviation of the gaussian and its error
	double sigma1,sigma2,sigmasum;
	double sigma_error1,sigma_error2,sigma_errorsum;

	// quality parameters of the fit
	double chi2ndof_1,chi2ndof_2,chi2ndofsum;
	double prob1,prob2,probsum;

	int bin1;
	int bin2;
	double fitrange;

	double uresidual_mean;
	double vresidual_mean;

	double minvalue=1.0/(2.0*2.7);

	int NumberOfTracks=fithistogram1->GetEntries();

	uresidual_mean=histogramu->GetMean();
	vresidual_mean=histogramv->GetMean();

	if(NumberOfTracks>100)
	{
		//Fit histogram with standard gaussian
		//fithistogram1->Draw("");

		// Fitting type:
		// 0: gaussian fit function with cuts on the tails, use the two distributions seperately
		// 1: gaussian fit function with cuts on the tails, use the merged distribution
		// 2: No fit, instead read out the RMS value of the histogram
                // 3: gaussian fit function with cuts on the tails, use only theta1
                // 4: gaussian fit function with cuts on the tails, use only theta2
		//
		if(fittype==0||fittype==3||fittype==4)
		{
			// Fit histograms
			bin1 = fithistogram1->FindFirstBinAbove(fithistogram1->GetMaximum()*minvalue);
			bin2 = fithistogram1->FindLastBinAbove(fithistogram1->GetMaximum()*minvalue);
			fitrange = (fithistogram1->GetBinCenter(bin2) - fithistogram1->GetBinCenter(bin1))/2.0;
			fithistogram1->Fit("gaus","","",-fitrange,fitrange);

			bin1 = fithistogram2->FindFirstBinAbove(fithistogram2->GetMaximum()*minvalue);
			bin2 = fithistogram2->FindLastBinAbove(fithistogram2->GetMaximum()*minvalue);
			fitrange = (fithistogram2->GetBinCenter(bin2) - fithistogram2->GetBinCenter(bin1))/2.0;
			fithistogram2->Fit("gaus","","",-fitrange,fitrange);
	
			bin1 = fithistogramsum->FindFirstBinAbove(fithistogramsum->GetMaximum()*minvalue);
			bin2 = fithistogramsum->FindLastBinAbove(fithistogramsum->GetMaximum()*minvalue);
			fitrange = (fithistogramsum->GetBinCenter(bin2) - fithistogramsum->GetBinCenter(bin1))/2.0;
			fithistogramsum->Fit("gaus","","",-fitrange,fitrange);

			// Get the fit function	
			TF1 *fit1 = fithistogram1->GetFunction("gaus");
			TF1 *fit2 = fithistogram2->GetFunction("gaus");
			TF1 *fitsum = fithistogramsum->GetFunction("gaus");

			// Extract Fit parameters of the first distribution
			mean1=fit1->GetParameter(1);
			sigma1=fit1->GetParameter(2);
			sigma_error1=fit1->GetParError(2);

			// Extract Quality parameters of the first distribution
			chi2ndof_1=fit1->GetChisquare()/(fit1->GetNDF()*1.0);
			prob1=fit1->GetProb();

			// Extract Fit parameters of the second distribution
			mean2=fit2->GetParameter(1);
			sigma2=fit2->GetParameter(2);
			sigma_error2=fit2->GetParError(2);

			// Extract Quality parameters of the second distribution
			chi2ndof_2=fit2->GetChisquare()/(fit2->GetNDF()*1.0);
			prob2=fit2->GetProb();

			probsum=fitsum->GetProb();
		}

		if(fittype==1)
		{

			// Fit histograms
			bin1 = fithistogram1->FindFirstBinAbove(fithistogram1->GetMaximum()*minvalue);
			bin2 = fithistogram1->FindLastBinAbove(fithistogram1->GetMaximum()*minvalue);
			fitrange = (fithistogram1->GetBinCenter(bin2) - fithistogram1->GetBinCenter(bin1))/2.0;
			fithistogram1->Fit("gaus","","",-fitrange,fitrange);

			bin1 = fithistogram2->FindFirstBinAbove(fithistogram2->GetMaximum()*minvalue);
			bin2 = fithistogram2->FindLastBinAbove(fithistogram2->GetMaximum()*minvalue);
			fitrange = (fithistogram2->GetBinCenter(bin2) - fithistogram2->GetBinCenter(bin1))/2.0;
			fithistogram2->Fit("gaus","","",-fitrange,fitrange);
	
			bin1 = fithistogramsum->FindFirstBinAbove(fithistogramsum->GetMaximum()*minvalue);
			bin2 = fithistogramsum->FindLastBinAbove(fithistogramsum->GetMaximum()*minvalue);
			fitrange = (fithistogramsum->GetBinCenter(bin2) - fithistogramsum->GetBinCenter(bin1))/2.0;
			fithistogramsum->Fit("gaus","","",-fitrange,fitrange);

			// Get the fit function	
			TF1 *fit1 = fithistogram1->GetFunction("gaus");
			TF1 *fit2 = fithistogram2->GetFunction("gaus");
			TF1 *fitsum = fithistogramsum->GetFunction("gaus");

			// Extract Fit parameters of the merged distribution
			meansum=fitsum->GetParameter(1);
			sigmasum=fitsum->GetParameter(2);
			sigma_errorsum=fitsum->GetParError(2);

			// Extract the mean of the fit from the two projected angle distributions
			mean1=fit1->GetParameter(1);
			mean2=fit2->GetParameter(1);

			// Extract Quality parameters of the merged distribution
			chi2ndofsum=fitsum->GetChisquare()/(fitsum->GetNDF()*1.0);
			probsum=fitsum->GetProb();

			prob1=fit1->GetProb();
			prob2=fit2->GetProb();
		}

		if(fittype==2)
		{
			// Read out RMS values
			sigma1=fithistogram1->GetRMS();
			sigma2=fithistogram2->GetRMS();	

			// Read out errors of RMS values
			sigma_error1=fithistogram1->GetRMSError();
			sigma_error2=fithistogram2->GetRMSError();	

			// Read out RMS values
			mean1=fithistogram1->GetMean();
			mean2=fithistogram2->GetMean();

			// No fit -> fit quality parameters are 0
			chi2ndof_1=0;
			chi2ndof_2=0;
			prob1=0;
			prob2=0;
			probsum=0;

		}

		// save the fitted histograms
		file->cd("");
		file->cd("mapping/fit");
		fithistogram1->Write("theta1_"+histoname+"_fit");
		fithistogram2->Write("theta2_"+histoname+"_fit");
		fithistogramsum->Write("sumhisto_"+histoname+"_fit");
	}

	else
	{
		mean1=0.0;
		sigma1=0.0;
		// Set the error of the std dev to 0, elsewise there will be problems later (division by 0)
		sigma_error1=0.0;
		chi2ndof_1=0.0;
		prob1=0.0;

		mean2=0.0;
		sigma2=0.0;
		// Set the error of the std dev to 0, elsewise there will be problems later (division by 0)
		sigma_error2=0.0;
		chi2ndof_2=0.0;
		prob2=0.0;

		sigmasum=0.0;

		cout<<"Histograms (nearly) empty... no fit possible!"<<endl;
	}

	// Go to results directory
	file->cd("");
	file->cd("mapping/result");

	// Get the chi2 histogram
	chi2map=(TH2*)file->Get("mapping/result/hchi2map");

	// Get the mean1 histogram
	meanmap1=(TH2*)file->Get("mapping/result/hmeanmap1");

	// Get the mean2 histogram
	meanmap2=(TH2*)file->Get("mapping/result/hmeanmap2");

	// Get the mean u residual histogram
	uresidual_meanmap=(TH2*)file->Get("mapping/result/huresidualmeanmap");

	// Get the mean v residual histogram
	vresidual_meanmap=(TH2*)file->Get("mapping/result/hvresidualmeanmap");
        // Get the momentum 2D histogram
	mommap=(TH2*)file->Get("mapping/result/hmommap");

	// Fill both maps containing the theta means
	meanmap1->SetBinContent(col+1,row+1,mean1);
	meanmap2->SetBinContent(col+1,row+1,mean2);

	// Fill both maps containing the residual means
	uresidual_meanmap->SetBinContent(col+1,row+1,uresidual_mean*1E3);
	vresidual_meanmap->SetBinContent(col+1,row+1,vresidual_mean*1E3);

        // Fill the momentum image
	mommap->SetBinContent(col+1,row+1,mom);

	double X0;
	double X0err;

	// Fill histograms with mean X0 (in %) and mean absolute X0 error (in %) of both fits
	if(fittype==0||fittype==2)
	{
		if(sigma1>recoerror&&sigma2>recoerror)
		{
			X0=(calculateX0(sigma1, recoerror, mom, model)+calculateX0(sigma2, recoerror, mom, model))/2.0;
	                X0err=sqrt(pow(getX0err(mom,sigma1,sigma_error1,recoerror,model),2)+pow(getX0err(mom,sigma2,sigma_error2,recoerror,model),2))/2.0;
		}

		else
		{
			X0=0.0;
	                X0err=99.0;
		}
                chi2map->SetBinContent(col+1,row+1,(chi2ndof_1+chi2ndof_2)/2.0);
	}

        // Fill it with X0 value and X0 error from theta1 fits
	else if(fittype==3)
	{
		if(sigma1>recoerror)
		{
			X0=calculateX0(sigma1, recoerror, mom, model);
                        X0err=getX0err(mom,sigma1,sigma_error1,recoerror,model);
		}

		else
		{
			X0=0.0;
	                X0err=99.0;
		}
                chi2map->SetBinContent(col+1,row+1,chi2ndof_1);
	}

        // Fill it with X0 value and X0 error from theta2 fits
	else if(fittype==4)
	{
		if(sigma2>recoerror)
		{
			X0=calculateX0(sigma2, recoerror, mom, model);
                        X0err=getX0err(mom,sigma2,sigma_error2,recoerror,model);
		}

		else
		{
			X0=0.0;
		}
		chi2map->SetBinContent(col+1,row+1,chi2ndof_2);
	}

        // Merged theta distribution fit
	else
	{
		if(sigmasum>recoerror)
		{
			X0=calculateX0(sigmasum, recoerror, mom, model);
	                X0err=getX0err(mom,sigmasum,sigma_errorsum,recoerror,model);

		}

		else
		{
			X0=0.0;
	                X0err=99.0;
		}

		chi2map->SetBinContent(col+1,row+1,chi2ndofsum);
	}

	// Get the X0 histogram
	X0map=(TH2*)file->Get("mapping/result/hX0map");
	cout<<"X0 in area "<<col<<","<<row<<" is: "<<X0<<"%"<<endl;
	X0map->SetBinContent(col+1,row+1,X0);

	// Get the X0 fit error histogram
	X0errmap=(TH2*)file->Get("mapping/result/hX0errmap");
	cout<<"X0 error in area "<<col<<","<<row<<" is: "<<X0err<<"%"<<endl;
	X0errmap->SetBinContent(col+1,row+1,X0err);

	// Get the relative X0 fit error histogram
	X0relerrmap=(TH2*)file->Get("mapping/result/hX0relerrmap");

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
	TH2* nummap=(TH2*)file->Get("mapping/result/hnummap");

	// Fill it with mean prob of both fits
	nummap->SetBinContent(col+1,row+1,NumberOfTracks);

	// Get the prob map (2D) histogram (theta1)
	TH2* probmap1=(TH2*)file->Get("mapping/result/hprobmap1");

	// Get the prob histogram (theta1)
	TH1* probhisto1=(TH1*)file->Get("mapping/result/hprobhisto1");

	probmap1->SetBinContent(col+1,row+1,prob1);
	probhisto1->Fill(prob1);

	// Get the prob map (2D) histogram (theta2)
	probmap2=(TH2*)file->Get("mapping/result/hprobmap2");

	// Get the prob histogram (theta2)
	TH1* probhisto2=(TH1*)file->Get("mapping/result/hprobhisto2");

	probmap2->SetBinContent(col+1,row+1,prob2);
	probhisto2->Fill(prob2);

	// Get the prob map (2D) histogram (sum)
	TH2* probmapsum=(TH2*)file->Get("mapping/result/hprobmapsum");

	// Get the prob histogram (sum)
	TH1* probhistosum=(TH1*)file->Get("mapping/result/hprobhistosum");

	probmapsum->SetBinContent(col+1,row+1,probsum);
	probhistosum->Fill(probsum);

	// Write the histogram to the rootfile
	if(col==numcol-1&&row==numrow-1)
	{

		X0map->Write();
		X0errmap->Write();

		chi2map->Write();
		meanmap1->Write();
		meanmap2->Write();

		uresidual_meanmap->Write();
		vresidual_meanmap->Write();

		nummap->Write();
                X0relerrmap->Write();

		mommap->Write();

		probmap1->Write();
		probmap2->Write();
		probmapsum->Write();

		probhisto1->Write();
		probhisto2->Write();
		probhistosum->Write();
	}

	//delete all histograms from memory
	fithistogram1->SetDirectory(gROOT);
	delete fithistogram1;

	fithistogram2->SetDirectory(gROOT);
	delete fithistogram2;

	fithistogramsum->SetDirectory(gROOT);
	delete fithistogramsum;
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

	// TString for the input root file name
	TString histoname,range;
	TString filename=mEnv.GetValue("x0filename", "X0-merge");

	// Set preprocessing parameter
	// 0: Use standard procedure
	// 1: Additional preprocessing loop to correct possible shifts in the angle distributions
	int correctmean=1;

	// Copy the X0 Analysis Root file 
	TFile *X0file = new TFile(filename+".root", "READ");

	//Open the copied file
	filename=filename+mEnv.GetValue("x0fileidentifier", "-part-1-1");
	TFile *rootfile = new TFile(filename+".root", "RECREATE");

	// Create directories containing map histograms, fits and results
	rootfile->mkdir("mapping");
	rootfile->mkdir("mapping/raw");
	rootfile->mkdir("mapping/fit");
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

	// Binning and range of the histograms
	int nbins=100;
	double range1=-0.0015;
	double range2=0.0015;

	// Choose the fitting function
	// 0: gaussian fit function with cuts on the tails, both kink distributions are used seperately
	// 1: gaussian fit function with cuts on the tails, use only 1 fit on the merged histogram consisting of both distributions
	// 2: No fit, instead read out the RMS value of the histogram
        // 3: gaussian fit function with cuts on the tails, use only theta1
        // 4: gaussian fit function with cuts on the tails, use only theta2
	int fittype=1;

	// Choose MSC model to calculate X/X0 from the standard deviation
	// 0: Highland model
	// 1: Frühwirth parametrization for the Highland model
	// 2: Frühwirth core width model
	// 3: Frühwirth RMS model
	int model=0;

	//Fit range should be quite small to reduce the influence of the MSC tails
	double plotrange=0.003;
	int numberofbins=200;

	// Calibration factor lambda, used to change the reconstruction error to include systematical errors
	double lambda=mEnv.GetValue("calibrationfactor", 1.0);
	double recoerror=sqrt(getanglerecovar(X0file))*lambda;
	cout<<"The reconstruction error is "<<recoerror*1E6<<" µrad!"<<endl;

	// Beam energy in GeV
        // The particle momenta are distributed due to the beam generation at the desy facility
        // Particles of different momenta are sperated by a dipole magnet and a collimator. Due to
        // the finite size of the collimator opening particles with momenta in a certain range p0+/-delta_p
        // can traverse into the test beam area
        // We expect a linear distribution with slope corresponding to ~500MeV/20mm
	double mom0=mEnv.GetValue("momentumoffset", 4.0);           // in GeV
    double mom_uslope=mEnv.GetValue("momentumugradient", 0.0);    // in GeV/mm
    double mom_vslope=mEnv.GetValue("momentumvgradient", 0.0);    // in GeV/mm
	cout<<"The beam energy is "<<mom0<<" GeV!"<<endl;
	cout<<"The beam energy gradient (u direction) is "<<mom_uslope<<" GeV/mm!"<<endl;
	cout<<"The beam energy gradient (v direction) is "<<mom_vslope<<" GeV/mm!"<<endl;
	
	

	// Save the theta1 histogram of this pixel to the root file
	if(correctmean==1)
	{
		getcorrection(X0file, rootfile, plotrange, numberofbins, numcol, numrow, umin, vmin, umax, vmax);
	}
	savehistos(X0file, rootfile, plotrange, numberofbins, numcol, numrow, umin, vmin, umax, vmax, correctmean);

	X0file->Close();

	rootfile->cd("");
	rootfile->cd("mapping/result");

	// X0 map
	TH2F * hX0map = new TH2F("hX0map","hX0map",numcol,umin,umax,numrow,vmin,vmax);
        hX0map->SetStats(kFALSE);
        hX0map->SetMaximum(10);
        hX0map->SetMinimum(0);
        hX0map->GetXaxis()->SetTitle("u [mm]");
        hX0map->GetYaxis()->SetTitle("v [mm]");
        hX0map->GetZaxis()->SetTitle("X/X0 [%]");
        hX0map->GetZaxis()->SetTitleSize(0.02);
        hX0map->GetZaxis()->SetLabelSize(0.02);

	// X0 statistical error map (absolute value)
	TH2F * hX0errmap = new TH2F("hX0errmap","hX0errmap",numcol,umin,umax,numrow,vmin,vmax);
        hX0errmap->SetStats(kFALSE);
        hX0errmap->SetMaximum(10);
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

			cout<<"fit histogram in (col,row): ("<<col<<","<<row<<")"<<endl;

                        // Calculate the u position of this bin
                        double u=umin+col*upitch;
                        double v=vmin+row*vpitch;

                        // Determine the momentum value from the u position and the p distribution parameters
                        double mom=GetMomentum(mom0, mom_uslope, mom_vslope, u, v);


			// fit the histograms
			fithisto(rootfile, fittype,model,col,numcol,row,numrow,mom, recoerror);
		}
	}

	rootfile->Close();

        return 0;

	
}
