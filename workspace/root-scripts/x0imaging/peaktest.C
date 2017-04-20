//#include <iostream.h>
#include <fstream>
using namespace std ;

  // This function fills histograms corresponding to certain u v values with msc angle distributions 
  void savehistos(TFile* file1, TFile* file2, double plotrange, int numberofbins, const int numcol, const int numrow, double umin, double vmin, double umax, double vmax)
  {
	// parameters which are read out from the root file
	Double_t theta1;
	Double_t theta2;
	Double_t u;
	Double_t v;
	
	//TTree in input root file, that contains the MSC projected angle distributions and reconstruction error distribution
	file1->cd("");

	TTree * msc_tree = (TTree*) file1->Get("MSCTree");

	msc_tree->SetBranchAddress("theta1",&theta1);
	msc_tree->SetBranchAddress("theta2",&theta2);
	msc_tree->SetBranchAddress("u",&u);
	msc_tree->SetBranchAddress("v",&v);

	file2->cd("");
	file2->cd("raw");

	// arrays of msc angle histograms
	TH1F *histo_theta1[numcol][numrow];
	TH1F *histo_theta2[numcol][numrow];
	TH1F *histo_thetasum[numcol][numrow];
	TH2F *histo_2d[numcol][numrow];

	TH1F *histo_trackchi2[numcol][numrow];

	for (int i=0; i<numcol; i++)
	{
		for (int j=0; j<numrow; j++)
		{
		 	histo_theta1[i][j] = new TH1F("","",numberofbins,-plotrange,plotrange);
		 	histo_theta2[i][j] = new TH1F("","",numberofbins,-plotrange,plotrange);
		 	histo_thetasum[i][j] = new TH1F("","",numberofbins,-plotrange,plotrange);
			histo_2d[i][j] = new TH2F("","",numberofbins,-plotrange,plotrange,numberofbins,-plotrange,plotrange);

			// Get the histograms generated in the getcorrection function
			// Name of the histograms
			TString aidhistoname;
			aidhistoname.Form("area(%i,%i)",i,j);
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

		// skip this entry, when u or v is outside of the testing area
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
			histo_thetasum[i][j]->Write("sumhisto_"+histoname);
			histo_thetasum[i][j]->Delete();
			histo_2d[i][j]->Write("2Dhisto_"+histoname);
			histo_2d[i][j]->Delete();

		}
		
	}

  }

  // This function fills histograms corresponding to certain u v values with msc angle distributions 
  void testdistribution(TFile* rootfile, int col, int numcol, int row, int numrow, double rms_criterium, int binrange, double peak_height)
  {

	// histogram name
	TString histoname;
	histoname.Form("area(%i,%i)",col,row);
	rootfile->cd("");
	histo=(TH1*)rootfile->Get("raw/sumhisto_"+histoname);

	// Clone the current histogram
	TH1* testedhistogram=(TH1*)histo->Clone("testedhistogram");

	// Get maximum of histogram, the test of the distributions starts there
	double maxvalue=histo->GetMaximum();
	int maxbin=histo->GetMaximumBin();

	cout<<"Histogram in (col,row): ("<<col<<","<<row<<")"<<endl;
	cout<<"Maximal value: "<<maxvalue<<endl;
	cout<<"at bin: "<<maxbin<<endl;

	int numbins=histo->GetNbinsX();

	int currentbin=maxbin+binrange;

	cout<<"Total number of bins: "<<numbins<<endl;

	// Get RMS value of histogram
	double RMS=histo->GetRMS();
	cout<<"RMS value of this histograms is "<<RMS<<" rad!"<<endl;
	double RMS_in_bins=RMS/(histo->GetBinCenter(1)-histo->GetBinCenter(0));
	cout<<"RMS value of this histograms is "<<RMS_in_bins<<" bins!"<<endl;

	int teststatistic=0;

        // Loop over bins on the right side of maximal value
	while((teststatistic==0)&&(currentbin<=numbins))
	{

		double currentbin_value=histo->GetBinContent(currentbin);	
		double left_mean=0;
		double right_mean=0;

		for(int i=currentbin+1;i<currentbin+binrange+1;i++) right_mean+=histo->GetBinContent(i)/binrange;
		for(int i=currentbin-binrange;i<currentbin;i++) left_mean+=histo->GetBinContent(i)/binrange;

		// Print out current bin
		cout<<"Current bin: "<<currentbin<<endl;

		//  value of current pixel
		cout<<"Value of current pixel: "<<currentbin_value<<endl;

		// mean value left of current pixel
		cout<<"Mean value left of current pixel: "<<left_mean<<endl;

		// mean value right of current pixel
		cout<<"Mean value right of current pixel: "<<right_mean<<endl;

		//calculate width of peak if currentbin lies in a local minimum:
		// Calculate the mean value of some bins left and right (depending on binrange), if the currentbin value * peakheight is smaller than the mean values
		// a peak was found
		if((peak_height*currentbin_value<left_mean)&&(peak_height*currentbin_value<right_mean)) 
		{
			// Calculate RMS of peak
			double peak_rms=currentbin-maxbin;

			// RMS value of histogram
			cout<<"RMS value of histogram: "<<RMS_in_bins<<" bins"<<endl;

			// RMS value of peak
			cout<<"RMS value of peak: "<<peak_rms<<" bins"<<endl;

			// Compare peak RMS to overall RMS of histogram
			if(peak_rms<rms_criterium*RMS_in_bins) 
			{
				
				teststatistic=1;

				// Print out comparison result
				cout<<"The RMS of the peak is smaller than  "<<rms_criterium*RMS_in_bins<<" bins"<<endl;
				cout<<"This angle distribution shows the characteristical peak structure of digital effects!"<<endl;
			}


		}

		// teststatitsic
		cout<<"Teststatistic: "<<teststatistic<<endl;

		// Go to next bin
		currentbin++;
	}

	// Save the histogram according to the teststatistic value
	if(teststatistic==1) 
	{
		rootfile->cd("");
		rootfile->cd("peaky");
		testedhistogram->Write("sumhisto_"+histoname+"_tested");
	}

	else 
	{
		rootfile->cd("");
		rootfile->cd("smooth");
		testedhistogram->Write("sumhisto_"+histoname+"_tested");
	}

	// Get teststatistic histogram
	image=(TH2*)rootfile->Get("result/hteststatistic");
	counter=(TH1F*)rootfile->Get("result/hteststatistic_counter");
	rootfile->cd("");
	rootfile->cd("result");
	image->SetBinContent(col+1,row+1,teststatistic);
	counter->Fill(teststatistic);
	if(col==numcol-1&&row==numrow-1) 
	{
		image->Write();
		counter->Write();
	}
  }


// This script is used to create a map of a plane in a test beam telescope. The input is a TTree including 
// MSC projected scattering angle distributions and reconstruction errors.
int peaktest()
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

	// TString for the input root file name
	TString filename,histoname,range;
	filename="X0-merge-calice-coghits";

	// Copy the X0 Analysis Root file 
	TFile *X0file = new TFile(filename+".root", "READ");

	//Open the copied file
	filename=filename+"-test_statistic-1mmbins";
	TFile *rootfile = new TFile(filename+".root", "RECREATE");

	// Create directories containing map histograms, fits and results

	rootfile->mkdir("raw");
	rootfile->mkdir("peaky");
	rootfile->mkdir("smooth");
	rootfile->mkdir("result");

	rootfile->cd("");
	rootfile->cd("result");

	// Number of Rows and Columns of the sensor map
	int numcol = 20;
	int numrow = 10;

	// u minimum and v maximum value (in mm)
	double umin=-10.0;
	double vmax=5.0;

	// u and v length of the map (in mm)
	double ulength=20.0;
	double vlength=10.0;

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

	//Fit range should be quite small to reduce the influence of the MSC tails
	double plotrange=0.003;
	int numberofbins=200;

	int binrange=3;
	double rms_criterium=0.5;
	double peak_height=1.1;

	// hteststatistic
	TH2F * hteststatistic = new TH2F("hteststatistic","hteststatistic",numcol,umin,umax,numrow,vmin,vmax);
        hteststatistic->SetStats(kFALSE);
        hteststatistic->SetMaximum(2.0);
        hteststatistic->SetMinimum(-1.0);
        hteststatistic->GetXaxis()->SetTitle("u [mm]");
        hteststatistic->GetYaxis()->SetTitle("v [mm]");
        hteststatistic->GetZaxis()->SetTitle("teststatistic");

	// hteststatistic
	TH1F * hteststatistic_counter = new TH1F("hteststatistic_counter","hteststatistic_counter",5,0,5);
        hteststatistic->SetStats(kFALSE);
        hteststatistic->GetXaxis()->SetTitle("teststatistic");

	savehistos(X0file, rootfile, plotrange, numberofbins, numcol, numrow, umin, vmin, umax, vmax);

	for(int col=0; col<numcol; col++)
	{
		for(int row=0; row<numrow; row++)
		{

			cout<<"fit histogram in (col,row): ("<<col<<","<<row<<")"<<endl;

			// fit the histograms
			testdistribution(rootfile, col, numcol, row, numrow, rms_criterium, binrange, peak_height);
		}
	}

	X0file->Close();

	rootfile->Close();

        return 0;

	
}
