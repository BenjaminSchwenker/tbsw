#include <iostream>
#include <cmath>
#include <fstream>
#include "TSystem.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TString.h"
#include "TMath.h"
#include "TROOT.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <cstdlib>
#include <map>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include "TChain.h"
#include "TCanvas.h"
#include "TObjString.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TKey.h"
#include "TFractionFitter.h"
#include "TLegend.h"

// The purpose of this script is to do a simultaneous fit of a set of angle distributions, which are taken from the X0_merge.root input file. Each of these distributions 
// should come from an area with a homogenious thickness and material composition. The hole class defined in this script can be used to select the area and give them some attributes
// like for example the thickness, position and the X0 value in this area. The angle distributions are then fitted by a function which corresponds to the convolution 
// between a function given by a theoretical model for the multiple scattering (highland or moliere model) and a gaussian error function, which depends on the expected telescope angular resolution
// and a calibration factor. The fit estimates an optimal value for the calibration factor.

// This script can be used by starting root and afterwards typing the command .x calibration.C+.

using namespace std;

// MeasurementArea Class describing a measurement area on the calibration target plane. Parameters of the measurement areas are position and side length
// as well as the thickness, density, atomic number and atomic mass of the measurement region
class MeasurementArea
{
	private:

	// Parameter declarations
	double center_u,center_v;  	// Position of the center of the MA in mm
	double length_u,length_v;  	// Side length of the MA in mm
	double thickness;          	// Thickness of the material in the region of the MA in mm
	double density;				// Density of the target material in the area (in g/cm³)
	double Z;					// Atomic number of the target material in the area
	double A; 					// Atomic mass of the target material in the area

	public:

	// Constructors

	MeasurementArea(double, double, double, double, double, double, double, double);	// Constructor

	// Reading parameters of the measurement area

	double Get_u_min() { return center_u-0.5*length_u; }	// Return min u value of MA (mm)
	double Get_u_max() { return center_u+0.5*length_u; }	// Return max u value of MA (mm)
	double Get_v_min() { return center_v-0.5*length_v; }	// Return min v value of MA (mm)
	double Get_v_max() { return center_v+0.5*length_v; }	// Return max v value of MA (mm)
	double Get_u_center() { return center_u; }				// Return center u value of MA (mm)
	double Get_v_center() { return center_v; }				// Return center v value of MA (mm)
	double Get_u_length() { return length_u; }				// Return the side length (in u direction) of MA (mm)
	double Get_v_length() { return length_v; }				// Return the side length (in v direction) of MA (mm)
	double Get_thickness()	{ return thickness; }			// Return thickness (mm)
	double Get_density() { return density; }				// Return density (g/cm³)
	double Get_Z() { return Z; }							// Return atomic number Z
	double Get_A() { return A; }							// Return atomic mass A

	

	void PrintParameters() 		// Print all parameters
	{
		std::cout<<"u center: "<<center_u<<" mm"<<std::endl;
		std::cout<<"v center: "<<center_v<<" mm"<<std::endl;

		std::cout<<"u length: "<<length_u<<" mm"<<std::endl;
		std::cout<<"v length: "<<length_v<<" mm"<<std::endl;

		std::cout<<"thickness: "<<thickness<<" mm"<<std::endl;

		std::cout<<"density: "<<density<<" g/cm³"<<std::endl;
		std::cout<<"atomic number Z: "<<Z<<std::endl;
		std::cout<<"atomic mass A: "<<A<<std::endl;
	}
};	
// Constructor definition
MeasurementArea::MeasurementArea(double ucenter, double vcenter, double ulength, double vlength, double thick, double dens, double atom_num, double atom_mass )
{
	center_u=ucenter;	// mm
	center_v=vcenter;	// mm
	length_u=ulength;	// mm
	length_v=vlength;	// mm
	thickness=thick;	// mm
	density=dens;		// g/cm³
	Z=atom_num;
	A=atom_mass;
}

// Class describing the complete calibration grid. The grid consists of a set of holes with a specific radiation length.
// This class can be used to described the aluminium target which was used during the March and November 2015 PXD test beams. 
// The second constructor can be used to model this exact aluminium target
class Grid
{
	private:

	// Declaration of parameters

	std::vector<MeasurementArea> m_MeasurementAreas;		// Declaration the holes in the grid				

	public:

	// Constructor
	Grid(TEnv*);			// Create set of measurement areas from cfg file		

	const std::vector<MeasurementArea>& GetMeasurementAreas() const {return  m_MeasurementAreas;}
    std::vector<MeasurementArea>& GetMeasurementAreas() {return  m_MeasurementAreas;}

	void PrintGridParameters()
	{
		for(int i=0;i<m_MeasurementAreas.size();i++) 
		{
			cout<<"-----------------------"<<endl;
			cout<<"Measurement Area "<<i<<endl;
			cout<<"-----------------------"<<endl;			
			m_MeasurementAreas.at(i).PrintParameters();
		}
	}

	// Add another measurement area
	//void AddMeasurementArea(MeasurementArea MA)	{ MeasurementAreas.push_back(MA); }
};	

Grid::Grid(TEnv* mEnv)
{


	// Determine the number of fit functions from the cfg file
	int num_MA=0; // Number of additional measurement areas
	int num_line=0; // Number of Lines, which are used for BE gradient calibration

	

	TString MAname;
	MAname.Form("MA%i",num_MA+1);

	TString linename;
	linename.Form("line%i",num_line+1);

	while(mEnv->GetValue(MAname+".exist",0)!=0)
	{

		// Read out measurement area parameters
		double ucenter=mEnv->GetValue(MAname+".ucenter", 1.0);
		double vcenter=mEnv->GetValue(MAname+".vcenter", 1.0);
		double ulength=mEnv->GetValue(MAname+".ulength", 1.0);
		double vlength=mEnv->GetValue(MAname+".vlength", 1.0);
		double thickness=mEnv->GetValue(MAname+".thickness", 1.8);
		double Z=mEnv->GetValue(MAname+".atomicnumber", 13.0);
		double A=mEnv->GetValue(MAname+".atomicmassnumber", 27.0);
		double density=mEnv->GetValue(MAname+".density", 2.7);

		// Define measurement area based on these parameters
		MeasurementArea MA(ucenter,vcenter,ulength,vlength,thickness,density,Z,A);
  
		// Add measurement area to predefined grid
		m_MeasurementAreas.push_back(MA);

		num_MA++;
		MAname.Form("MA%i",num_MA+1);
	}

	while(mEnv->GetValue(linename+".exist",0)!=0)
	{
		if(mEnv->GetValue(linename+".usteplength",0.0)!=0.0) 
		{
			for(double d=mEnv->GetValue(linename+".startu",0.0);d<(mEnv->GetValue(linename+".startu",0.0)+mEnv->GetValue(linename+".ulength",-1.0));d+=mEnv->GetValue(linename+".usteplength",10.0))
			{
				// Compute/Read out measurement area parameters
				double ucenter=d;
				double vcenter=mEnv->GetValue(linename+".startv", 1000.0);
				double ulength=mEnv->GetValue(linename+".usteplength", 0.1);
				double vlength=mEnv->GetValue(linename+".vlength", 0.1);
				double thickness=mEnv->GetValue(linename+".thickness", 1.8);
				double Z=mEnv->GetValue(linename+".atomicnumber", 13.0);
				double A=mEnv->GetValue(linename+".atomicmassnumber", 27.0);
				double density=mEnv->GetValue(linename+".density", 2.7);

				// Define measurement area based on these parameters
				MeasurementArea MA(ucenter,vcenter,ulength,vlength,thickness,density,Z,A);
			  
				// Add measurement area to predefined grid
				m_MeasurementAreas.push_back(MA);
			}
		}


		else if(mEnv->GetValue(linename+".vsteplength",0.0)!=0.0) 
		{

			for(double d=mEnv->GetValue(linename+".startv",0.0);d<(mEnv->GetValue(linename+".startv",0.0)+mEnv->GetValue(linename+".vlength",-1.0));d+=mEnv->GetValue(linename+".vsteplength",10.0))
			{
				// Compute/Read out measurement area parameters
				double ucenter=mEnv->GetValue(linename+".startu", 1000.0);
				double vcenter=d;
				double ulength=mEnv->GetValue(linename+".ulength", 0.1);
				double vlength=mEnv->GetValue(linename+".vsteplength", 0.1);
				double thickness=mEnv->GetValue(linename+".thickness", 1.8);
				double Z=mEnv->GetValue(linename+".atomicnumber", 13.0);
				double A=mEnv->GetValue(linename+".atomicmassnumber", 27.0);
				double density=mEnv->GetValue(linename+".density", 2.7);

				// Define measurement area based on these parameters
				MeasurementArea MA(ucenter,vcenter,ulength,vlength,thickness,density,Z,A);
		  
				// Add measurement area to predefined grid
				m_MeasurementAreas.push_back(MA);
			}

		}
		// in this case there is nothing to do
		else continue;

		num_line++;
		linename.Form("line%i",num_line+1);
	}	

}

// Functions used in this script
double DetermineFitrange(TH1*,double);
double calculateB(double);
double getrecoerror(TFile*);
void savehisto(TFile*,TFile*,TString,TString, double, double, double, double, int );
void correcthisto(TFile*,TFile*, TString , double , double , double , double );
double* fit(TFile*,TString, Grid, std::vector<double>, double, TString, std::vector<bool>);
void shiftbins(TH1*, double);
void calibrationfit();
int** GetParameterMapping(int);


  // Highland model of a MSC angle distribution, the parameters are:

  /*
	* par[0]:  Expected beam energy at u,v=0;
	* par[1]:  Beam particle charge
	* par[2]:  Beam particle mass
	* par[3]:  Target material density
	* par[4]:  Target material atomic number
	* par[5]:  Target material atomic weight
	* par[6]:  Thickness
	* par[7]:  Expected angle reconstruction error
	* par[8]:  reco error calibration factor
	* par[9]: Normalization
	* par[10]:  u coordinate
	* par[11]: v coordinate
	* par[12]: u BE gradient
	* par[13]: v BE gradient

*/

// Function, which returns beam momentum value for every point on the target plane
// This is necessary because the beam profile at DESY often have beam energy gradients in the order of a few MeV/mm
Double_t GetMomentum(double meanvalue,double ugrad,double vgrad, double u, double v)
{
	double p;
	p=meanvalue+u*ugrad+v*vgrad;
	return p;
}
  
// Highland model of multiple scattering: Simple gaussian with a well defined standard deviation depending on X/X0 and the beam energy.
// The overall function describing the kink angle distributions is the Highland function convoluted with a gaussian function due to the finite angle resolution on the target plane. 
Double_t highlandfunction(Double_t *x, Double_t *par)
{  

	// atomic number of target material
	double Z; 
	Z=par[4];  

	// atomic weight of target material
	double A;
	A=par[5];

	//density of the target material
	double density;  
	density=par[3]; 

	// thickness of the target material
	double dm1=par[6]; // in mm
	double d1;	// in cm
	d1=dm1/10.0;

	// Other parameters
	double exp_recoerror=par[7];  //expected reconstruction error
	double lambda=par[8];	//calibration factor of the reconstrcution error

	double recoerror=exp_recoerror*lambda;   // calibrated reco error

        // particle parameters

	// mass of beam particle
	double mass;   
	mass=par[2];  

	// charge of beam particle
	double charge;   
	charge=par[1];

	// beam energy
	double p=GetMomentum(par[0],par[10],par[11],par[12],par[13]);

	// calibrated momentum
	double E=TMath::Sqrt(p*p+mass*mass);  // energy in GeV

	double beta;  //relative velocity
	beta=p/E;

	// Radiation length computed from the other parameters
	double X0=716.4*A/((Z+1)*Z*density*TMath::Log(287.0/TMath::Sqrt(Z)));

	// Combination of Highland width and reconstruction error
	double sigma=TMath::Sqrt(pow(0.0136*charge/(p*beta)*TMath::Sqrt(d1/X0)*(1.0+0.038*TMath::Log(d1/X0)),2)+pow(recoerror,2));

	// function value at a certain theta value
	double value=par[9]*TMath::Gaus(x[0],0.0,sigma);

	return value;
}// End definition of highland model


  // Moliere model of a MSC angle distribution, the parameters are:

  /* 
	* par[0]:   Expected beam energy at u,v=0;
	* par[1]:   Beam particle charge
	* par[2]:   Beam particle mass
	* par[3]:   Target material density
	* par[4]:   Target material atomic number
	* par[5]:   Target material atomic weight
	* par[6]:   Thickness
	* par[7]:   Expected angle reconstruction error
	* par[8]:   reco error calibration factor
	* par[9]:   Normalization
	* par[10]:  u coordinate
	* par[11]:  v coordinate
	* par[12]:  u BE gradient
	* par[13]:  v BE gradient

*/
  
// Moliere model of multiple scattering: Function also describing the tails of multiple scattering distributions. The function depends on material properties, which
// can be reduced to X0, the material thickness and the beam energy.
// The overall function describing the kink angle distributions is the Moliere function convoluted with a gaussian function due to the finite angle resolution on the target plane. 
Double_t molierefunction(Double_t *x, Double_t *par)
{  
	// atomic number of target material
	double Z; 
	Z=par[4];  

	// atomic weight of target material
	double A;
	A=par[5];

	//density of the target material
	double density;  
	density=par[3]; 

	// thickness of the target material
	double dm1=par[6]; // in mm
	double d1;	// in cm
	d1=dm1/10.0;

	// Other parameters
	double exp_recoerror=par[7];  //expected reconstruction error
	double lambda=par[8];	//calibration factor of the reconstrcution error

	double recoerror=exp_recoerror*lambda;   // calibrated reco error

    // particle parameters

	// mass of beam particle
	double mass;   
	mass=par[2];  

	// charge of beam particle
	double charge;   
	charge=par[1];

	// beam energy
	double p=GetMomentum(par[0],par[10],par[11],par[12],par[13]);   // Energy in GeV

	double E=TMath::Sqrt(p*p+mass*mass);

	double beta;  //relative velocity
	beta=p/E;

	//Moliere model parameters

	// areal density (proportional to the thickness) for thickness 1
	double arealdensity1=density*d1; 

	// em coupling constant
	double alpha=1.0/(137.0*beta)*Z; 

	// aid variable for thickness 1
	double log_omega_b1=8.215+log(pow(Z,(-0.6667))*(arealdensity1/A)*pow(alpha,2)/(1.13+3.76*pow(alpha,2)))/log(10.0); 
	//cout<<"log Omega 1 is "<<log_omega_b1<<endl;

	// parameter that will be used in the masterformula of the overall angle distribution (thickness1)
	double B1=calculateB(log_omega_b1);

	// chi_C parameter (dependend on thickness, momentum etc), given in µrad
	double chi_C=22.9*charge*Z/(1000.0*p)*TMath::Sqrt(arealdensity1/(pow(beta,2)*A))*TMath::Pi()/180.0;

	//Histogram Definitions
	double low,high;

	// number of points
	const int numpoints_aid=19;
	const int numpoints=2*numpoints_aid-1;

	// number of bins in the final histogram
	const int numbins=81;

	// upper limit of histogram
	high=4.05*chi_C*TMath::Sqrt(B1);
	
	// lower limit of histogram
	low=-high;

	// values taken from table 2
	double phi_values[numbins];
	double phi_values_rad[numbins];

	// Calculate phi(rad) from the normalized values
	for(int i=0; i<numbins; i++)
	{
		phi_values[i]=-4.0+i*0.1;
		phi_values_rad[i]=phi_values[i]*chi_C*TMath::Sqrt(B1);
	}
	
	// Get f1 and f2 values for the the distribution (only positive angle values)
	double f1_values_aid[numpoints_aid] ={0.0206,-0.0246,-0.1336,-0.2440,-0.2953,-0.2630,-0.1622,-0.0423,
				    0.0609,0.1274,0.147,0.142,0.1225,0.100,0.078,0.059,0.045,0.0316,0.0194};
	double f2_values_aid[numpoints_aid] ={0.416,0.299,0.019,-0.229,-0.292,-0.174,0.010,0.138,0.146,0.094,
				    0.045,-0.049,-0.071,-0.064,-0.043,-0.024,-0.010,0.001,0.006};

	//  f1 and f2 values for the total distribution (also negative angles, but some missing values)
	double f1_values_aid2[numpoints];
	double f2_values_aid2[numpoints];

	//  f1 and f2 values for the total distribution (also negative angles) between phi=-4.0 and 4.0 with an 0.1 increment
	double f1_values[numbins];
	double f2_values[numbins];

	// Get f1 and f2 values for the whole (also negative angles) distribution 
	for(int i=0; i<numpoints;i++)
	{
		if(i<numpoints_aid-1)
		{
			f1_values_aid2[i]=f1_values_aid[numpoints_aid-(i+1)];
			f2_values_aid2[i]=f2_values_aid[numpoints_aid-(i+1)];
		}

		else
		{
			f1_values_aid2[i]=f1_values_aid[i-(numpoints_aid-1)];
			f2_values_aid2[i]=f2_values_aid[i-(numpoints_aid-1)];
		}

	}

	// add temporary values to the bins which will be filled with the extrapolated values later
	for(int i=0; i<numbins;i++)
	{
			if(i<3)
			{			
				f1_values[i]=f1_values_aid2[0];
				f2_values[i]=f2_values_aid2[0];
			}

			else if(i<7)	
			{
				f1_values[i]=f1_values_aid2[1];
				f2_values[i]=f2_values_aid2[1];
			}

			else if(i>numbins-4)	
			{
				f1_values[i]=f1_values_aid2[numpoints-1];
				f2_values[i]=f2_values_aid2[numpoints-1];
			}

			else if(i>numbins-8)	
			{
				f1_values[i]=f1_values_aid2[numpoints-2];
				f2_values[i]=f2_values_aid2[numpoints-2];
			}



			else if(i<41)
			{
				if(i%2==1)
				{
					f1_values[i]=f1_values_aid2[i/2-1];
					f2_values[i]=f2_values_aid2[i/2-1];
				}
				else
				{
					f1_values[i]=f1_values_aid2[i/2-2];
					f2_values[i]=f2_values_aid2[i/2-2];
				}
			}

			else
			{
				if(i%2==1)
				{
					f1_values[i]=f1_values_aid2[i/2-2];
					f2_values[i]=f2_values_aid2[i/2-2];
				}
				else
				{
					f1_values[i]=f1_values_aid2[i/2-2];
					f2_values[i]=f2_values_aid2[i/2-2];

				}
			}

	}

	// Declaration of the histogram containing the values of the first function given in table 2
	TH1F * h_f1_table = new TH1F("f1_table","f1_table",numbins,low,high);
        h_f1_table->SetStats(kFALSE);

	// Declaration of the histogram containing the values of the first function given in table 2
	TH1F * h_f2_table = new TH1F("f2_table","f2_table",numbins,low,high);
        h_f2_table->SetStats(kFALSE);

	// Declaration of the histogram containing the overall moliere angle distibution (including reconstruction errors)
	TH1F * h_moliere_conv = new TH1F("moliere_conv","moliere_conv",numbins,low,high);
        h_moliere_conv->SetStats(kFALSE);

	// Declaration of the histogram containing the gaussian distribution
	TH1F * h_gaus = new TH1F("gausfunc","gausfunc",numbins,low,high);
        h_gaus->SetStats(kFALSE);

	// Declaration of the histogram containing the overall angle distibution
	TH1F * h_total1 = new TH1F("total1","total1",numbins,low,high);
        h_total1->SetStats(kFALSE);

	// fill the f1 and f2 histograms with the corresponding values
	for(int i=0; i<numbins;i++)
	{
			h_f1_table->SetBinContent(i+1,f1_values[i]);
			h_f2_table->SetBinContent(i+1,f2_values[i]);		
	}

	// fill the gaussian function with the corresponding values
	for(int i=0; i<numbins;i++)
	{
		// value of the gaussian function at this bin
		double value=2.0/TMath::Sqrt(TMath::Pi())*TMath::Exp(-1.0*pow((phi_values_rad[i]/(chi_C*TMath::Sqrt(B1))),2));

		// Bin number corresponding to the current phi value
		int bin_No=i+1;

		// Set histogram values
		h_gaus->SetBinContent(bin_No,value);	
	}

	// Interpolation of additional points of the f1 and f2 histograms
	for(int i=1; i<numbins-1;i+=2)
	{

	double interpolation1=0;
	double interpolation2=0;

		if(i<40)
		{
			// f1 interpolation
			double newpoint=0.5*(phi_values_rad[i]+phi_values_rad[i-1]);

			interpolation1=h_f1_table->Interpolate(newpoint);
			interpolation2=h_f2_table->Interpolate(newpoint);
		}

		else 
		{
			// f1 interpolation
			double newpoint=0.5*(phi_values_rad[i]+phi_values_rad[i+1]);

			interpolation1=h_f1_table->Interpolate(newpoint);
			interpolation2=h_f2_table->Interpolate(newpoint);
		}

			// Set histogram values
			h_f1_table->SetBinContent(i+1,interpolation1);
			h_f2_table->SetBinContent(i+1,interpolation2);	
	

	} 


	// fill the overall distribution with the values calculated from gaussian f1 and f2
	for(int i=0; i<numbins;i++)
	{
		// Bin number corresponding to the current phi value
		int bin_No=i+1;

		// value of the overall function computed from the other functions (thickness 1)
		double total_value1=h_gaus->GetBinContent(bin_No)+h_f1_table->GetBinContent(bin_No)/B1+h_f2_table->GetBinContent(bin_No)/(B1*B1);

		h_total1->SetBinContent(bin_No,total_value1);
	}

	// Do a normalization of the histograms
	// h_total1->Scale(1.0/h_total1->Integral("width"));

	double convolutionintegral=0.0;

	// Numerical convolution of moliere density and reco error gaussian function (sum over all bins of the histogram) 
	// at a certain point defined by j loop
	for(int i=0; i<numbins;i++)
	{
		// current value of the sum
		convolutionintegral+=h_total1->GetBinContent(i+1)*TMath::Gaus(phi_values_rad[i]-x[0],0.0,recoerror);
	}

	//delete all histograms from memory
	h_f1_table->SetDirectory(gROOT);
	delete h_f1_table;

	h_f2_table->SetDirectory(gROOT);
	delete h_f2_table;

	h_moliere_conv->SetDirectory(gROOT);
	delete h_moliere_conv;

	h_gaus->SetDirectory(gROOT);
	delete h_gaus;

	h_total1->SetDirectory(gROOT);
	delete h_total1;

	// function value at a certain theta value
	return par[9]*convolutionintegral;
} // End definition moliere model


// definition of shared parameters like the parameter numbers of each fit function and 
// a globalChi2 structure

// Number of parameters per fit function
const int num_localparameters=14;

// Global chi2 object, will be used in fit
struct GlobalChi2 { 
    std::vector< ROOT::Math::IMultiGenFunction * >& getChi2Vec() {return fChi2Vec;}
	void SetArray() {ipar=GetParameterMapping(fChi2Vec.size());}

    double operator() (const double *par) const {

      double ret =0;

	  
	  for(int j=0;j<fChi2Vec.size();j++)
	  {
        // read function args
        double p1[num_localparameters];
        for (int i = 0; i < num_localparameters; ++i) 
		{
				p1[i] = par[ipar[j][i] ];
		}
        // evaluate function and sum up
        ROOT::Math::IMultiGenFunction * func  =  fChi2Vec.at(j);

		// Add up all local chi2 values of the single fit functions
		ret += (*func)(p1);
       
      }
	  // return global chi2 value
	  return ret;
}


 
    std::vector< ROOT::Math::IMultiGenFunction * > fChi2Vec;
	int ** ipar;
  };


// Get Relation between global fit parameters and local fit parameters
// Parameter mapping
int **GetParameterMapping(int numfuncs)
{

	  int **ipar=0;

	  int newparsperfunction=7;

      ipar = new int*[numfuncs];

      for (int i = 0; i < numfuncs; i++)
      {
            ipar[i] = new int[num_localparameters];

            for (int j = 0; j < num_localparameters; j++)
            {
                  
				// There are basically 2 cases: The first function has the Parameters 0-14,
				//								the other function have a new Parameter number at the 3rd parameter (density), the 4th parameter (Z), the 4th parameter (A)
				//								the 6th Parameter (thickness of target material), the 9th Parameter (~#tracks), the 10th parameter (u coordinate) and 11th parameter (v coordinate)
						
				if(i==0) 
				{
					ipar[i][j] = j;
				}
				else
				{
					if((j!=3)&&(j!=4)&&(j!=5)&&(j!=6)&&(j!=9)&&(j!=10)&&(j!=11)) ipar[i][j]=j;
					else if(j==3) ipar[i][j]=14+(i-1)*newparsperfunction;
					else if(j==4) ipar[i][j]=15+(i-1)*newparsperfunction;
					else if(j==5) ipar[i][j]=16+(i-1)*newparsperfunction;
					else if(j==6) ipar[i][j]=17+(i-1)*newparsperfunction;
					else if(j==9) ipar[i][j]=18+(i-1)*newparsperfunction;
					else if(j==10) ipar[i][j]=19+(i-1)*newparsperfunction;
					else ipar[i][j]=20+(i-1)*newparsperfunction;
				}
				//cout<<"Parameter mapping: Global parameter number of local parameter "<<j<<" in fit function "<<i<<" is "<<ipar[i][j]<<endl;
            }
      }

      return ipar;
}


// Function used in Moliere model to calculate a parameter
double calculateB(double log_omega_b)
{
	// Initialize Bmin and Bmax values
	double Bmin=1.0;
	double Bmax=30.0;

	// Iniatialize B value as mean between starting values of lower and upper B bound
	double B=(Bmin+Bmax)/2.0;

	// Calculate difference between measured omega b and omega b(B) with this starting B value
	double difference;
	difference=pow(10,(log_omega_b))-1.167*exp(B)/B;

	// For loop used to minimize difference between measured omega b and omega b(B)
	for(Int_t i=1; sqrt(pow(difference,2))>1E-8; i++)
	{
		// If the difference is larger than 0 the current B value is too small, therefore the lower bound of B (Bmin)
		// is now set to be the current X/X0 value
		if(difference>0){Bmin=B;}

		// If the difference is smaller than 0 the current B value is too large, therefore the upper bound of B (Bmax)
		// is now set to be the current X/X0 value
		else{Bmax=B;}

		// Calculate the new B value from the updated limits
		B=(Bmin+Bmax)/2.0;
		difference=pow(10,(log_omega_b))-1.167*exp(B)/B;
		if(i==100)
		{
			B=1.0;
			break;
		}
	}

	return B;
}

// Determine fit range for a kink angle histogram
// This is done by finding the first and last bin above a certain threshold.
// The fit range is half of the distance between these two bins in rad
double DetermineFitrange(TH1* histo,double minvalue)
{

	int bin1 = histo->FindFirstBinAbove(histo->GetMaximum()*minvalue);
	int bin2 = histo->FindLastBinAbove(histo->GetMaximum()*minvalue);
	double fitrange = (histo->GetBinCenter(bin2) - histo->GetBinCenter(bin1))/2.0;

	return fitrange;
}


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
	cout<<"Angle Reconstruction Variance: 	"<<recovar<<" +/- "<<recovar_error<<"rad^2"<<endl; 

	return recovar;
}


// Shift a decentralized angle distribution in order to get a mean value of 0
void shiftbins(TH1* histogram, double mean1)
{
	// The mean offset has to be corrected before merging the two histograms
	// First the bin value of the offset has to be found out
	TAxis *xaxis =  histogram->GetXaxis();
	Int_t binx = xaxis->FindBin(mean1)-xaxis->FindBin(0.0);

	Double_t stats[5]={0,0,0,0,0};
	histogram->PutStats(stats); // reset mean value, etc
	Int_t nbins=histogram->GetNbinsX();

	//Now the histogram has to be shifted accordingly
	if(binx>=0)
	{
		for (Int_t j=1;j<=nbins-binx;j++) 
		{
		   	 histogram->SetBinContent(j, histogram->GetBinContent(j+binx));
		}
		for (Int_t j=nbins-binx;j<=nbins;j++)  
		{
			histogram->SetBinContent(j,0);
		}
	}
	else
	{
		for (Int_t j=nbins;j>=2-binx;j--)  
		{
			histogram->SetBinContent(j, histogram->GetBinContent(j+binx));
		}
		for (Int_t j=1;j<=1-binx;j++) 
		{
		   	 histogram->SetBinContent(j, histogram->GetBinContent(0));
		}
	}
}

// This function fills histograms corresponding to certain u v values with msc angle distributions 
void correcthisto(TFile* file,TFile* file2, TString histoname, TString range, double umin, double umax, double vmin, double vmax )
{
	//TTree in input root file, that contains the MSC projected angle distributions
	file->cd("");
	TTree *msc_tree = (TTree*)file->Get("MSCTree");

	// Definition of the four cuts used to define the measurement region
	TString cutcondition;

	cutcondition.Form("u>%f",umin);
	TString cutall=cutcondition;

	cutcondition.Form("u<%f",umax);
	cutall=cutall+"&&"+cutcondition;

	cutcondition.Form("v>%f",vmin);
	cutall=cutall+"&&"+cutcondition;

	cutcondition.Form("v<%f",vmax);
	cutall=cutall+"&&"+cutcondition;


	// Draw histogram of the first scattering angle in the given u and v range and save it
	msc_tree->Draw("theta1>>histo1("+range+")",cutall,"");

	// Get the histogram and save it in the raw folder
	TH1 * aidhistogram=msc_tree->GetHistogram();
	aidhistogram->SetTitle("#theta_{1} distribution");
	aidhistogram->GetXaxis()->SetTitle("#theta_{1} [rad]");

	// Go to raw directory
	file2->cd("");
	file2->cd("grid/raw/");

	aidhistogram->Write("theta1_uncorrected_"+histoname);

	// Delete histogram1 from memory
	aidhistogram->SetDirectory(gROOT);
	delete aidhistogram;

	// Draw histogram of second scattering angle in the given u and v range and save it
	msc_tree->Draw("theta2>>histo2("+range+")",cutall,"");

	// Get the histogram and save it in the raw folder
	TH1 * aidhistogram2=msc_tree->GetHistogram();
	aidhistogram2->SetTitle("#theta_{2} distribution");
	aidhistogram2->GetXaxis()->SetTitle("#theta_{2} [rad]");

	// Go to raw directory
	file2->cd("");
	file2->cd("grid/raw/");

	aidhistogram2->Write("theta2_uncorrected_"+histoname);

	// Delete histogram2 from memory
	aidhistogram2->SetDirectory(gROOT);
	delete aidhistogram2;

	// Draw sum histogram of the scattering angles in the given u and v range and save it
	TH1 * hsum1;
	msc_tree->Draw("theta1>>hsum1("+range+")",cutall,"");
	hsum1=msc_tree->GetHistogram();

	TH1 * hsum2;
	msc_tree->Draw("theta2>>hsum2("+range+")",cutall,"");
	hsum2=msc_tree->GetHistogram();

	// Give the two histograms to a listOperator
	TList *list = new TList;
        list->Add(hsum1);
        list->Add(hsum2);

	// Merge all histograms in the list
        TH1 *h = (TH1*)hsum1->Clone("h");
        h->Reset();
        h->Merge(list);
	h->SetTitle("merged projected angle distribution");
	h->GetXaxis()->SetTitle("#theta_{proj.} [rad]");

	// Go to raw directory
	file2->cd("");
	file2->cd("grid/raw/");

	h->Write("sumhisto_uncorrected_"+histoname);

	// Delete sum histogram from memory
	h->SetDirectory(gROOT);
	delete h;

	// Delete histogram 1 from memory
	hsum1->SetDirectory(gROOT);
	delete hsum1;

	// Delete histogram 2 from memory
	hsum2->SetDirectory(gROOT);
	delete hsum2;
}


// Function to save the projected angle histograms of different regions in the u-v plane, the region lies within the given 
// u and v min and max values.
void savehisto(TFile* file, TFile* file2, TString histoname, TString range, double umin, double umax, double vmin, double vmax, int correctmean)
{
	//TTree in input root file, that contains the MSC projected angle distributions
	file->cd("");
	TTree *msc_tree = (TTree*)file->Get("MSCTree");

	// Array of mean theta1 and theta2 values in each map pixel
	double mean1=0.0;
	double mean2=0.0;

	// Definition of the four cuts used to define the measurement region
	TString cutcondition;

	cutcondition.Form("u>%f",umin);
	TString cutall=cutcondition;

	cutcondition.Form("u<%f",umax);
	cutall=cutall+"&&"+cutcondition;

	cutcondition.Form("v>%f",vmin);
	cutall=cutall+"&&"+cutcondition;

	cutcondition.Form("v<%f",vmax);
	cutall=cutall+"&&"+cutcondition;

	if(correctmean==1)
	{
		// Get histogram
		TH1* histogram1=(TH1*)file2->Get("grid/raw/theta1_uncorrected_"+histoname);
		TH1* histogram2=(TH1*)file2->Get("grid/raw/theta2_uncorrected_"+histoname);

		// Save mean of both distributions in arrays
		mean1=histogram1->GetMean();
		mean2=histogram2->GetMean();
	}

	// Draw histogram of the first scattering angle in the given u and v range and save it
	msc_tree->Draw("theta1>>h("+range+")",cutall,"");

	// Get the histogram and save it in the raw folder
	TH1 * aidhistogram=msc_tree->GetHistogram();
	
	// Shift the bins of the histogram to reduce the offset
	if(correctmean==1) shiftbins(aidhistogram,mean1);

	aidhistogram->SetTitle("#theta_{1} distribution");
	aidhistogram->GetXaxis()->SetTitle("#theta_{1} [rad]");

	// Go to raw directory
	file2->cd("");
	file2->cd("grid/raw/");

	aidhistogram->Write("theta1_"+histoname);

	// Delete histogram1 from memory
	aidhistogram->SetDirectory(gROOT);
	delete aidhistogram;

	// Draw histogram of second scattering angle in the given u and v range and save it
	msc_tree->Draw("theta2>>h("+range+")",cutall,"");

	// Get the histogram and save it in the raw folder
	TH1 * aidhistogram2=msc_tree->GetHistogram();

	// Shift the bins of the histogram to reduce the offset
	if(correctmean==1) shiftbins(aidhistogram2,mean2);

	aidhistogram2->SetTitle("#theta_{2} distribution");
	aidhistogram2->GetXaxis()->SetTitle("#theta_{2} [rad]");

	// Go to raw directory
	file2->cd("");
	file2->cd("grid/raw/");

	aidhistogram2->Write("theta2_"+histoname);

	// Delete histogram2 from memory
	aidhistogram2->SetDirectory(gROOT);
	delete aidhistogram2;

	// Draw sum histogram of the scattering angles in the given u and v range and save it
	TH1 * hsum1;
	msc_tree->Draw("theta1>>hsum1("+range+")",cutall,"");
	hsum1=msc_tree->GetHistogram();

	// Shift the bins of the histogram to reduce the offset
	if(correctmean==1) shiftbins(hsum1,mean1);

	TH1 * hsum2;
	msc_tree->Draw("theta2>>hsum2("+range+")",cutall,"");
	hsum2=msc_tree->GetHistogram();

	// Shift the bins of the histogram to reduce the offset
	if(correctmean==1) shiftbins(hsum2,mean2);

	// Give the two histograms to a list Operator
	TList *list = new TList;
        list->Add(hsum1);
        list->Add(hsum2);

	// Merge all histograms in the list
        TH1 *h = (TH1*)hsum1->Clone("h");
        h->Reset();
        h->Merge(list);
	h->SetTitle("merged projected angle distribution");
	h->GetXaxis()->SetTitle("#theta_{proj.} [rad]");

	// Go to raw directory
	file2->cd("");
	file2->cd("grid/raw/");
	h->Write("sumhisto_"+histoname);

	// Delete sum histogram from memory
	h->SetDirectory(gROOT);
	delete h;

	// Delete histogram 1 from memory
	hsum1->SetDirectory(gROOT);
	delete hsum1;

	// Delete histogram 2 from memory
	hsum2->SetDirectory(gROOT);
	delete hsum2;

}

// Function to fit the MSC angle histograms simultaneously, it returns a pointer to the fit results
double* fit( TFile* file, Grid grid, std::vector<double> beamoptions, double recoerr, TString model, std::vector<bool> fitoptions)
{ 

	// Read out the parameters from the beamoptions vector
	double z=beamoptions.at(0);
	double mass=beamoptions.at(1);

	double BE_mean=beamoptions.at(2);
	double BE_ugrad=beamoptions.at(3);
	double BE_vgrad=beamoptions.at(4);


	// Read out the parameters from the fitoptions vector
	bool fixlambda=fitoptions.at(0);
	bool fixmomentumoffset=fitoptions.at(1);
	bool fixmomentumugradient=fitoptions.at(2);
	bool fixmomentumvgradient=fitoptions.at(3);


	// Some parameter definitions

	double lambda_startvalue=beamoptions.at(5);

	// Number of new parameters for every new measurement area
	int newparsperfunction=7;

	// histogram name
	TString histoname;

	// Fit range
	double fitrange;

	const int num_fitfunctions=grid.GetMeasurementAreas().size();

	// fitresults Array
	double *fitresults = new double[8];

	// Fit and raw angle histograms
	TH1 * fithistogramsum[num_fitfunctions];
	TH1 * histogramsum;

	// Copy raw histograms
	for(int i=0; i< num_fitfunctions; i++)
	{
		// histogram name
		histoname.Form("measurementarea%i",i+1);
		file->cd("");

		histogramsum=(TH1*)file->Get("grid/raw/sumhisto_"+histoname);

		// Make a copy of the histogram
		fithistogramsum[i]=(TH1*)histogramsum->Clone("sumhisto_"+histoname+"_fit");

		file->cd("grid/fit/");
	}
	
	// Minvalue and fit range depend on the model: In Case of the moliere model the fitrange can be a little larger
	double minvalue;
	// Declaration of fit functions
	TF1 *fitFcn[num_fitfunctions];
	TString fctname;

	ROOT::Fit::DataOptions opt; 
	ROOT::Fit::DataRange range[num_fitfunctions];

	// loop for definition of the fit functions, the fitrange is determined for every one of them
	for(int i=0;i<num_fitfunctions;i++)
	{
		// Fitrange depends on model
		if(model=="moliere") minvalue=1.0/(10*2.71);
		else minvalue=1.0/(1.0*2.71);

		// And fit range is also smaller, when there is only air
		if(grid.GetMeasurementAreas().at(i).Get_thickness()<0.0001) minvalue=1.0/(1.5*2.71);

		fitrange=DetermineFitrange(fithistogramsum[i],minvalue);
		fctname.Form("fitFcn%i",i);

		if(model=="moliere")
		{
			// Use Gaussian function with width corresponding to the calibrated angle resolution, if material is only air or extremely thin
			// Also in this case the fit range is set to be a little smaller
			if(grid.GetMeasurementAreas().at(i).Get_thickness()<0.0001) fitFcn[i] = new TF1(fctname,"[9]*TMath::Gaus(x,0.0*[0]*[1]*[2]*[3]*[4]*[5]*[6]*[10]*[11]*[12]*[13],[7]*[8])",-fitrange,fitrange);


			// Use Moliere model in case the material is not just air
			else fitFcn[i] = new TF1(fctname,molierefunction,-fitrange,fitrange,num_localparameters);
		}

		else
		{
			// Use Gaussian function with width corresponding to the calibrated angle resolution, if material is only air or extremely thin
			if(grid.GetMeasurementAreas().at(i).Get_thickness()<0.0001) fitFcn[i] = new TF1(fctname,"[9]*TMath::Gaus(x,0.0*[0]*[1]*[2]*[3]*[4]*[5]*[6]*[10]*[11]*[12]*[13],[7]*[8])",-fitrange,fitrange);
			// Use Highland model in case the material is not just air
			else fitFcn[i] = new TF1(fctname,highlandfunction,-fitrange,fitrange,num_localparameters);
		}

		// set the data range	
		cout<<"Fitrange at d="<<grid.GetMeasurementAreas().at(i).Get_thickness()<<" mm: "<<fitrange<<" rad"<<endl;
		cout<<"Min value at d="<<grid.GetMeasurementAreas().at(i).Get_thickness()<<" mm: "<<minvalue<<endl;
		range[i].SetRange(-fitrange,fitrange);
	}

    // Vector of wrapped multi tf1
    std::vector<ROOT::Math::WrappedMultiTF1> wf;

	// Vector of bin data
	std::vector<ROOT::Fit::BinData> data;

	// Vector of chi2 functions
    std::vector<ROOT::Fit::Chi2Function> chi2;
	// The stored chi2 functions mustn't be reallocated, else there will be a crash.
	// Therefore reserve enough space for the vector
	chi2.reserve(3*num_fitfunctions);

	// data size definition 
	// needed during the fit
	int datasize;

	cout<<"1"<<endl;
	
    // Create globalChi2 object, which will be used during the fitting to return fit chi2 values and update the global parameters
	GlobalChi2 globalChi2;

    for(int i=0;i<num_fitfunctions;i++) 
	{
		// Create wrapped multi function entry
		ROOT::Math::WrappedMultiTF1 wf_entry(*fitFcn[i],1);
		wf.push_back(wf_entry);
		
		// Create data entry and fill data store and increase data size variable
		ROOT::Fit::BinData data_entry(opt,range[i]);
		data.push_back(data_entry);
		ROOT::Fit::FillData(data.at(i), fithistogramsum[i]);
		datasize+=data.at(i).Size();
	}


	for(int i=0;i<num_fitfunctions;i++)
	{
		// Create chi2 function entry and fill global chi2
		chi2.push_back(ROOT::Fit::Chi2Function(data[i], wf[i]));

		globalChi2.getChi2Vec().push_back(&(chi2.at(i))); 
	}

	// Set parameter mapping array in the global chi2 object
	globalChi2.SetArray();


	cout<<"2"<<endl;
    
	ROOT::Fit::Fitter fitter;

	// Number of global parameters: There are num_fitfunctions fit functions with three new parameters each: Scale, u coordinate, vcoordinate, thickness
	// Then there are 8 parameters, which are used in all fit functions (beam energy calibration factor, angle calibration factor, material properties, etc)
	static const int num_globalparameters = num_fitfunctions*newparsperfunction+(num_localparameters-newparsperfunction);

	double par0[num_globalparameters];

	// Set the entries of the globalparameters array

	// First fit function has num_localparameters new parameters:
	double aid_array[num_localparameters]={ BE_mean,z,mass,grid.GetMeasurementAreas().at(0).Get_density(),grid.GetMeasurementAreas().at(0).Get_Z(),grid.GetMeasurementAreas().at(0).Get_A(),
											grid.GetMeasurementAreas().at(0).Get_thickness(),recoerr,lambda_startvalue,700.0,grid.GetMeasurementAreas().at(0).Get_u_center(),
											grid.GetMeasurementAreas().at(0).Get_v_center(),BE_ugrad,BE_vgrad};
	for(int i=0;i<num_localparameters;i++) par0[i]=aid_array[i];

	// Afterwards for each fit functions we get 3 new parameters
	for(int i=1;i<num_fitfunctions;i++)
	{
		par0[14+(i-1)*newparsperfunction]=grid.GetMeasurementAreas().at(i).Get_density();
		par0[15+(i-1)*newparsperfunction]=grid.GetMeasurementAreas().at(i).Get_Z();
		par0[16+(i-1)*newparsperfunction]=grid.GetMeasurementAreas().at(i).Get_A();
		par0[17+(i-1)*newparsperfunction]=grid.GetMeasurementAreas().at(i).Get_thickness();
		par0[18+(i-1)*newparsperfunction]=700;
		par0[19+(i-1)*newparsperfunction]=grid.GetMeasurementAreas().at(i).Get_u_center();
		par0[20+(i-1)*newparsperfunction]=grid.GetMeasurementAreas().at(i).Get_v_center();
	}

	// create before the parameter settings in order to fix or set range on them
	fitter.Config().SetParamsSettings(num_globalparameters,par0);

	// fix constant parameters 1 to 7
	for(int i=1;i<8;i++) fitter.Config().ParSettings(i).Fix();

	// fix u and v position of first measurement area
	fitter.Config().ParSettings(10).Fix();
	fitter.Config().ParSettings(11).Fix();

	// Fix momentum offset?
	if(fixmomentumoffset) fitter.Config().ParSettings(0).Fix();

	// Fix the beam momentum gradient parameters?
	if((fixmomentumugradient)) fitter.Config().ParSettings(12).Fix();
	if((fixmomentumvgradient)) fitter.Config().ParSettings(13).Fix();

	// Limit the beam momentum gradient parameters, they should be small anyway (< 10MeV/mm)
	// and we dont want to see negative beam energies at any position on the target plane
	fitter.Config().ParSettings(12).SetLimits(-0.035,0.035);
	fitter.Config().ParSettings(13).SetLimits(-0.035,0.035);

	// The beam energy mean value should also always be positive
	fitter.Config().ParSettings(0).SetLimits(0.5*BE_mean,1.5*BE_mean);

	// fix density, A, Z, thickness and coordinate parameters for all fitfunctions
	for(int i=1;i<num_fitfunctions;i++)
	{
		fitter.Config().ParSettings(14+(i-1)*newparsperfunction).Fix();
		fitter.Config().ParSettings(15+(i-1)*newparsperfunction).Fix();
		fitter.Config().ParSettings(16+(i-1)*newparsperfunction).Fix();
		fitter.Config().ParSettings(17+(i-1)*newparsperfunction).Fix();
		fitter.Config().ParSettings(19+(i-1)*newparsperfunction).Fix();
		fitter.Config().ParSettings(20+(i-1)*newparsperfunction).Fix();
	}

	if(fixlambda) fitter.Config().ParSettings(8).Fix();				 //fix lambda?

	// set limits on the calibration factor parameter
	fitter.Config().ParSettings(8).SetLimits(0.6,1.4);

	//Set limits on fit function Normalizations
	fitter.Config().ParSettings(9).SetLimits(10.0,200000.0);
	for(int i=1;i<num_fitfunctions;i++) fitter.Config().ParSettings(18+(i-1)*newparsperfunction).SetLimits(10.0,200000.0);

	fitter.Config().MinimizerOptions().SetPrintLevel(10);
	fitter.Config().SetMinimizer("Minuit2","Migrad"); 

	// fit FCN function directly 
	// (specify optionally data size and flag to indicate that is a chi2 fit)
	fitter.FitFCN(num_globalparameters,globalChi2,0,datasize,true);
	ROOT::Fit::FitResult result = fitter.Result();
	result.Print(std::cout);

	int ** parameter_mapping=GetParameterMapping(num_fitfunctions);

	for(int i=0;i<num_fitfunctions;i++)
	{	

		int parameters[num_localparameters];
		for(int j=0;j<num_localparameters;j++)
		{
			parameters[j]=parameter_mapping[i][j];
		}

		fitFcn[i]->SetFitResult( result, parameters);

		fitrange=DetermineFitrange(fithistogramsum[i],minvalue);
		fitFcn[i]->SetRange(-fitrange,fitrange);  
		fitFcn[i]->SetLineColor(kRed);
		fithistogramsum[i]->GetListOfFunctions()->Add(fitFcn[i]);
		histoname.Form("gridpoint%i",i+1);
		// display mode
		gStyle->SetOptFit(1111111);
		fithistogramsum[i]->Write("thetasum_"+histoname+"_fit");
	}

	// Plot the fit results, print them in terminal and save them to histogram
		
	// Save lambda value and error as result
	fitresults[0]=result.Parameter(8);
	fitresults[1]=result.Error(8);

	// Save BE mean value (at u,v=0) as result
	fitresults[2]=result.Parameter(0);
	fitresults[3]=result.Error(0);

	// Save BE u gradient as result
	fitresults[4]=result.Parameter(12);
	fitresults[5]=result.Error(12);

	// Save BE v gradient as result
	fitresults[6]=result.Parameter(13);
	fitresults[7]=result.Error(13);

	
	// output of the single fit function chi2 values and the quadratic sum of them
	double chi2_summation=0.0;

	for(int i=0; i<num_fitfunctions;i++)
	{
		fctname.Form("fitFcn%i",i);

		cout<<"--------------------------"<<endl;
		cout<<"Fit function "<<i<<endl;
		//cout<<"Chi2 value for this individual fit is: "<<fitFcn[i]->GetChisquare()<<endl;
		cout<<"Chi2 value for this individual fit is: "<<fithistogramsum[i]->Chisquare(fitFcn[i],"R")<<endl;
		cout<<"Number of degrees of freedom is: "<<fitFcn[i]->GetNDF()<<endl;
		cout<<"Alternatively: Chi2 is: "<<chi2.at(i)(fitFcn[i]->GetParameters())<<endl;
		cout<<"--------------------------"<<endl;

		chi2_summation+=pow(fithistogramsum[i]->Chisquare(fitFcn[i],"R"),2);
	}
	
	cout<<"The quadratic sum of these chi2 values is: "<<sqrt(chi2_summation)<<endl;
	TString pdfname,Title;

	if(num_fitfunctions>4)
	{
		TCanvas *c1 = new TCanvas("c1","multipads",900,1000);
		TPad *pad1 = new TPad("pad1","pad1",0.01,0.51,0.49,0.99);
		pad1->Draw();
		TPad *pad2 = new TPad("pad2","pad2",0.51,0.51,0.99,0.99);
		pad2->Draw();
		TPad *pad3 = new TPad("pad3","pad3",0.01,0.01,0.49,0.49);
		pad3->Draw();
		TPad *pad4 = new TPad("pad4","pad4",0.51,0.01,0.99,0.49);
		pad4->Draw();

		int aid=num_fitfunctions/4;

	   	pad1->cd();
		Title.Form("Area %i: d=%fmm",0,grid.GetMeasurementAreas().at(0).Get_thickness());
		fithistogramsum[0]->SetTitle(Title);
		fithistogramsum[0]->Draw();

		pad2->cd();
		Title.Form("Area %i: d=%fmm",aid,grid.GetMeasurementAreas().at(aid).Get_thickness());
		fithistogramsum[aid]->SetTitle(Title);
		fithistogramsum[aid]->Draw();

		pad3->cd();
		Title.Form("Area %i: d=%fmm",num_fitfunctions-(aid+1),grid.GetMeasurementAreas().at(num_fitfunctions-(aid+1)).Get_thickness());
		fithistogramsum[num_fitfunctions-(aid+1)]->SetTitle(Title);
		fithistogramsum[num_fitfunctions-(aid+1)]->Draw();

		pad4->cd();
		Title.Form("Area %i: d=%fmm",num_fitfunctions-1,grid.GetMeasurementAreas().at(num_fitfunctions-1).Get_thickness());
		fithistogramsum[num_fitfunctions-1]->SetTitle(Title);
		fithistogramsum[num_fitfunctions-1]->Draw();

		pdfname=model+"_results1.pdf";

		c1->SaveAs(pdfname); 
	}


	if(num_fitfunctions>11)
	{
		TCanvas *c2 = new TCanvas("c2","multipads2",900,1000);
		TPad *pads1 = new TPad("pads1","pad1",0.01,0.76,0.32,0.99);
		pads1->Draw();
		TPad *pads2 = new TPad("pads2","pad2",0.34,0.76,0.65,0.99);
		pads2->Draw();
		TPad *pads3 = new TPad("pads3","pad3",0.67,0.76,0.99,0.99);
		pads3->Draw();
		TPad *pads4 = new TPad("pads4","pad4",0.01,0.51,0.32,0.74);
		pads4->Draw();
		TPad *pads5 = new TPad("pads5","pad5",0.34,0.51,0.65,0.74);
		pads5->Draw();
		TPad *pads6 = new TPad("pads6","pad6",0.67,0.51,0.99,0.74);
		pads6->Draw();
		TPad *pads7 = new TPad("pads7","pad7",0.01,0.26,0.32,0.49);
		pads7->Draw();
		TPad *pads8 = new TPad("pads8","pad8",0.34,0.26,0.65,0.49);
		pads8->Draw();
		TPad *pads9 = new TPad("pads9","pad9",0.67,0.26,0.99,0.49);
		pads9->Draw();
		TPad *pads10 = new TPad("pads10","pad10",0.01,0.01,0.32,0.24);
		pads10->Draw();
		TPad *pads11 = new TPad("pads11","pad11",0.34,0.01,0.65,0.24);
		pads11->Draw();
		TPad *pads12 = new TPad("pads12","pad12",0.67,0.01,0.99,0.24);
		pads12->Draw();

		int aid=num_fitfunctions/12;

	   	pads1->cd();
		Title.Form("Measurement area %i: d=%fmm",0,grid.GetMeasurementAreas().at(0).Get_thickness());
		fithistogramsum[0]->SetTitle(Title);
		fithistogramsum[0]->Draw();

		pads2->cd();
		Title.Form("Measurement area %i: d=%fmm",aid,grid.GetMeasurementAreas().at(aid).Get_thickness());
		fithistogramsum[aid]->SetTitle(Title);
		fithistogramsum[aid]->Draw();

		pads3->cd();
		Title.Form("Measurement area %i: d=%fmm",2*aid,grid.GetMeasurementAreas().at(2*aid).Get_thickness());
		fithistogramsum[2*aid]->SetTitle(Title);
		fithistogramsum[2*aid]->Draw();

		pads4->cd();
		Title.Form("Measurement area %i: d=%fmm",3*aid,grid.GetMeasurementAreas().at(3*aid).Get_thickness());
		fithistogramsum[3*aid]->SetTitle(Title);
		fithistogramsum[3*aid]->Draw();

	   	pads5->cd();
		Title.Form("Measurement area %i: d=%fmm",4*aid,grid.GetMeasurementAreas().at(4*aid).Get_thickness());
		fithistogramsum[4*aid]->SetTitle(Title);
		fithistogramsum[4*aid]->Draw();

		pads6->cd();
		Title.Form("Measurement area %i: d=%fmm",5*aid,grid.GetMeasurementAreas().at(5*aid).Get_thickness());
		fithistogramsum[5*aid]->SetTitle(Title);
		fithistogramsum[5*aid]->Draw();

		pads7->cd();
		Title.Form("Measurement area %i: d=%fmm",num_fitfunctions-(5*aid+1),grid.GetMeasurementAreas().at(num_fitfunctions-(5*aid+1)).Get_thickness());
		fithistogramsum[num_fitfunctions-(5*aid+1)]->SetTitle(Title);
		fithistogramsum[num_fitfunctions-(5*aid+1)]->Draw();

		pads8->cd();
		Title.Form("Measurement area %i: d=%fmm",num_fitfunctions-(4*aid+1),grid.GetMeasurementAreas().at(num_fitfunctions-(4*aid+1)).Get_thickness());
		fithistogramsum[num_fitfunctions-(4*aid+1)]->SetTitle(Title);
		fithistogramsum[num_fitfunctions-(4*aid+1)]->Draw();

	   	pads9->cd();
		Title.Form("Measurement area %i: d=%fmm",num_fitfunctions-(3*aid+1),grid.GetMeasurementAreas().at(num_fitfunctions-(3*aid+1)).Get_thickness());
		fithistogramsum[num_fitfunctions-(3*aid+1)]->SetTitle(Title);
		fithistogramsum[num_fitfunctions-(3*aid+1)]->Draw();

		pads10->cd();
		Title.Form("Measurement area %i: d=%fmm",num_fitfunctions-(2*aid+1),grid.GetMeasurementAreas().at(num_fitfunctions-(2*aid+1)).Get_thickness());
		fithistogramsum[num_fitfunctions-(2*aid+1)]->SetTitle(Title);
		fithistogramsum[num_fitfunctions-(2*aid+1)]->Draw();

		pads11->cd();
		Title.Form("Measurement area %i: d=%fmm",num_fitfunctions-(aid+1),grid.GetMeasurementAreas().at(num_fitfunctions-(aid+1)).Get_thickness());
		fithistogramsum[num_fitfunctions-(aid+1)]->SetTitle(Title);
		fithistogramsum[num_fitfunctions-(aid+1)]->Draw();

		pads12->cd();
		Title.Form("Measurement area %i: d=%fmm",num_fitfunctions-1,grid.GetMeasurementAreas().at(num_fitfunctions-1).Get_thickness());
		fithistogramsum[num_fitfunctions-1]->SetTitle(Title);
		fithistogramsum[num_fitfunctions-1]->Draw();

		pdfname=model+"_results2.pdf";
		c2->SaveAs(pdfname); 	
	}


	// Copy the X0 Analysis Root file and rename it
	TFile *resultsfile = new TFile("X0calibration_results.root", "RECREATE");

	TH1F * resultshist=new TH1F("resultshist","results of calibration",4,1,4);

    resultshist->GetXaxis()->SetBinLabel( 1, "lambda" );
    resultshist->GetXaxis()->SetBinLabel( 2, "BE_mean[GeV]" );
    resultshist->GetXaxis()->SetBinLabel( 3, "BE_u_grad[GeV/mm]" );
    resultshist->GetXaxis()->SetBinLabel( 3, "BE_u_grad[GeV/mm]" );

	// Save the results to histogram
	resultshist->SetBinContent(1,fitresults[0]);
	resultshist->SetBinError(1,fitresults[1]);

	resultshist->SetBinContent(2,fitresults[2]);
	resultshist->SetBinError(2,fitresults[3]);

	resultshist->SetBinContent(3,fitresults[4]);
	resultshist->SetBinError(3,fitresults[5]);

	resultshist->SetBinContent(4,fitresults[6]);
	resultshist->SetBinError(4,fitresults[7]);

	resultshist->Write("results");

	resultsfile->Close();
	

	// Return pointer to fitresults
	return fitresults;
} // End of fit function


  // Function used to get the angle distributions of tracks crossing certain regions in the u-v measurement plane,
  // perform a moliere fit on the distributions to estimate the calibrationfactors mu and lambda. Usually this script is used on 
  // measurement data of a plane with a precisely known material distribution ( for example a aluminum grid with a set 
  // of holes with different thicknesses.
  void calibrationfit()
  { 

	#if defined(__CINT__) && !defined(__MAKECINT__) 
   	cout << "ERROR: This script can run only using ACliC, you must run it by doing: " << endl;
   	cout << "\t .x $ROOTSYS/tutorials/fit/combinedFit.C+" << endl;
   	return;
	#endif

	gSystem->Load("libProof.so");
	gSystem->Load("libTreePlayer.so");

	gROOT->Reset(); 
	gROOT->SetStyle("Plain"); 	

    // Read config file
    //------------------
    TEnv *mEnv=new TEnv("x0calibration.cfg");

	// Choose the multiple scattering model
	// "moliere": Moliere model
	// "highland": Highland model
	TString model=mEnv->GetValue("model", "highland");

	// Fix lambda or mu or slope parameter during fit?
	bool fixlambda=mEnv->GetValue("fixlambda", 0);
    bool fix_offset=mEnv->GetValue("fix_momentumoffset", 0);
	bool fix_u_gradient=mEnv->GetValue("fix_momentumugradient", 0);
	bool fix_v_gradient=mEnv->GetValue("fix_momentumvgradient", 0);

	// Integer determining, whether a offset correction is applied to the angular distributions
	// 1: correction
	// everything else, no correction
	int correctmean=mEnv->GetValue("correctmean", 1 );

	// Get the lambda starting value
	double 	lambda_start_default=mEnv->GetValue("lambda_start", 1.0);

	// Define and set the parameters used in the Moliere fit
	// BE: beam energy (GeV), z: charge of beam particle (e), mass: mass of beam particle (GeV),
	// d_layer: thickness per layer(mm), num_layers: total number of layers, recoerr: angle reconstruction error (rad)
	double z,p,mass,recoerr;
	
	cout<<"---------------------------------------------------"<<endl;
	cout<<"-----------------Parameter settings----------------"<<endl;
	cout<<"---------------------------------------------------"<<endl;
	// Set mean beam energy (GeV) and slope
	double BE_mean_default=mEnv->GetValue("momentumumoffset", 4.0);
	double BE_ugrad_default=mEnv->GetValue("momentumugradient", 0.0);			// energy slope in GeV/mm in u direction
	double BE_vgrad_default=mEnv->GetValue("momentumvgradient", 0.0);			// energy slope in GeV/mm in u direction

	// Set beam particle charge (e)	
	z= mEnv->GetValue("particle.charge", 1);

	// Set beam particle mass (GeV)	
	mass= mEnv->GetValue("particle.mass", 4.0);
	cout<<"Mass: 					"<<mass<<" GeV"<<endl;
	cout<<"Charge: 				"<<z<<" e"<<endl;

	// Create empty grid object that will be filled and used in the fit
	Grid grid(mEnv);

	// Total number of measurement areas
	const int num_fitfunctions=grid.GetMeasurementAreas().size();

	cout<<"Total number of measurement areas: 	"<<num_fitfunctions<<endl;

	cout<<"---------------------------------------------------"<<endl;
	cout<<"---------------------------------------------------"<<endl;

	// Definition of the measurement areas can be done in two ways. Either a fixed 3 x 3 alu grid with a regular thickness increase between grid points and fixed positions or 
    // 12 completely independent measurement areas with completely unreleated thicknesses and positions.

    double u_MA_center[num_fitfunctions];
	double v_MA_center[num_fitfunctions];

	// u minimum and v maximum values (in mm)
	double umin[num_fitfunctions];
	double vmin[num_fitfunctions];
	double umax[num_fitfunctions];
	double vmax[num_fitfunctions];
	
	// Print out the measurement areas, which will be used for the fit
	grid.PrintGridParameters();

	// Readout max and min values of all measurement areas of the calibration grid
	for(int i=0;i<num_fitfunctions;i++)
	{
		umin[i]=grid.GetMeasurementAreas().at(i).Get_u_min();
		umax[i]=grid.GetMeasurementAreas().at(i).Get_u_max();
		vmin[i]=grid.GetMeasurementAreas().at(i).Get_v_min();
		vmax[i]=grid.GetMeasurementAreas().at(i).Get_v_max();
	}

	// TString for the input root file name
	TString filename,histoname,range;
	filename.Form("X0");

	// Copy the X0 Analysis Root file and rename it
	TFile *X0file = new TFile(filename+".root", "READ");

	//Open the copied file
	filename=filename+"CalibrationDQM_"+model+"fit";
	TFile *rootfile = new TFile(filename+".root", "RECREATE");

	// Set number of layers of target	
	recoerr=sqrt(getanglerecovar(X0file));
	cout<<"Angle reconstruction error: 	"<<recoerr<<" rad"<<endl;

	// Create directories containing angle histograms and fits
	rootfile->mkdir("grid");
	rootfile->mkdir("grid/raw");
	rootfile->mkdir("grid/fit");

	rootfile->cd("");

	// Binning and range of the histograms

	int nbins=500;
	double range1=-0.0025;
	double range2=0.0025;

	// Set the range and number of bins of the histogram
	range.Form("%i,%f,%f",nbins,range1,range2);

	rootfile->mkdir("grid/raw/");
	rootfile->mkdir("grid/fit/");

	// Loop over number of measurement areas
	for(Int_t i=0; i<num_fitfunctions; i++)
	{
			// Set the histogram name as a string
			histoname.Form("measurementarea%i",i+1);

			cout<<"save histogram of measurement area "<<i+1<<endl;

			// Save the angle histograms of the current measurement area to the root file
			correcthisto(X0file,rootfile, histoname, range, umin[i], umax[i], vmin[i], vmax[i]);
			savehisto(X0file,rootfile, histoname, range, umin[i], umax[i], vmin[i], vmax[i], correctmean);

	}// end of first loop over measurement areas

	X0file->Close();

	// Open results config file
    //------------------
    TEnv *mEnv_res=new TEnv("x0cal_result.cfg");
	// Check whether there are entries, of this is the case use these entries as starting values

	double lambda_start=mEnv_res->GetValue("lambda", lambda_start_default);
	double BE_mean=mEnv_res->GetValue("momentumumoffset", BE_mean_default);
	double BE_ugrad=mEnv_res->GetValue("momentumugradient", BE_ugrad_default);	
	double BE_vgrad=mEnv_res->GetValue("momentumvgradient", BE_ugrad_default);


	std::vector<double> beamoptions;
	// Use abs, because z,mass and BE_mean should be positive
	beamoptions.push_back(abs(z));
	beamoptions.push_back(abs(mass));
	beamoptions.push_back(abs(BE_mean));
	beamoptions.push_back(BE_ugrad);
	beamoptions.push_back(BE_vgrad);
	// lambda is not really a beam parameter, but this is the only place this fits in
	beamoptions.push_back(lambda_start);
	

	std::vector<bool> fitoptions;
	fitoptions.push_back(fixlambda);
	fitoptions.push_back(fix_offset);
	fitoptions.push_back(fix_u_gradient);
	fitoptions.push_back(fix_v_gradient);

	cout<<"Beam particle mean energy start value is				"<<BE_mean<<" GeV"<<endl;
	cout<<"Beam particle energy u gradient start value is				"<<BE_ugrad<<" GeV/mm"<<endl;
	cout<<"Beam particle energy v gradient start value is				"<<BE_vgrad<<" GeV/mm"<<endl;
	cout<<"Beam particle lambda		 start value is				"<<lambda_start<<endl;

	double* iresults=fit(rootfile, grid, beamoptions, recoerr, model, fitoptions);

	cout<<" The lambda calibration factor is: "<<iresults[0]<<" +/- "<<iresults[1]<<endl;
	cout<<" The mean beam energy at (0,0) is: "<<iresults[2]<<" +/- "<<iresults[3]<<"GeV"<<endl;
	cout<<" The BE gradient in u direction is: "<<iresults[4]<<" +/- "<<iresults[5]<<"GeV/mm"<<endl;
	cout<<" The BE gradient in v direction is: "<<iresults[6]<<" +/- "<<iresults[7]<<"GeV/mm"<<endl;

	TString value;
	// Fill cfg file with result values
	value.Form("%f",iresults[0]);
	mEnv_res->SetValue("lambda_start",value);
	value.Form("%f",iresults[2]);
	mEnv_res->SetValue("momentumoffset",value);
	value.Form("%f",iresults[4]);
	mEnv_res->SetValue("momentumugradient",value);
	value.Form("%f",iresults[6]);
	mEnv_res->SetValue("momentumvgradient",value);
	mEnv_res->Print();
	mEnv_res->SaveLevel(kEnvLocal);

	rootfile->Close();

  gApplication->Terminate();
  }

