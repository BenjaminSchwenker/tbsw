/// The purpose of this script is to do a simultaneous fit of a set of angle distributions, which are taken from the X0_merge.root input file. Each of these distributions 
/// should come from an area with a homogenious thickness and material composition. The hole class defined in this script can be used to select the area and give them some attributes
/// like for example the thickness, position and the X0 value in this area. The angle distributions are then fitted by a function which corresponds to the convolution 
/// between a function given by a theoretical model for the multiple scattering (highland or moliere model) and a gaussian error function, which depends on the expected telescope angular resolution
/// and a calibration factor. The fit estimates an optimal value for the calibration factor.
///
/// This script can be used by starting root and afterwards typing the command .x calibrationfit.C+.
///
/// \author Ulf Stolzenberg 
/// \author Benjamin Schwenker

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include "TROOT.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/PoissonLikelihoodFCN.h"
#include "HFitInterface.h"
#include "TH1.h"
#include "TH1F.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TEnv.h"
#include "TFitResult.h"




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
	double X0; 					// Radiation length constant of the target material in the area
	int run_min; 				// Minimal run number of data, which contains the specified material parameters
	int run_max;				// Maximal run number of data, which contains the specified material parameters
	int max_angles;				// Maximal number of scattering angles in the scattering angle distribution

	public:

	// Constructors

	MeasurementArea(double, double, double, double, double, double, double, double, double, int, int, int);	// Constructor

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
	double Get_X0() { return X0; }							// Return X0 (mm)
	double Get_Z() { return Z; }							// Return atomic number Z
	double Get_A() { return A; }							// Return atomic mass A

	int Get_run_min() { return run_min; }					// Return run_min
	int Get_run_max() { return run_max; }					// Return run_max

	int Get_max_angles() { return max_angles; }			    // Return max_angles

	

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
		std::cout<<"X0: "<<X0<<" mm"<<std::endl;

		std::cout<<"Min. run number: "<<run_min<<std::endl;
		std::cout<<"Max run number: "<<run_max<<std::endl;

		std::cout<<"Max number of angles: "<<max_angles<<std::endl;
	}
};	
// Constructor definition
MeasurementArea::MeasurementArea(double ucenter, double vcenter, double ulength, double vlength, double thick, double X0constant, double dens, double atom_num, double atom_mass, int run_minimum, int run_maximum, int angle_maximum)
{
	center_u=ucenter;	// mm
	center_v=vcenter;	// mm
	length_u=ulength;	// mm
	length_v=vlength;	// mm
	thickness=thick;	// mm
	density=dens;		// g/cm³
	X0=X0constant;		// mm
	Z=atom_num;
	A=atom_mass;
	run_min=run_minimum;
	run_max=run_maximum;
	max_angles=angle_maximum;
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
		for(size_t i=0;i<m_MeasurementAreas.size();i++) 
		{
			cout<<endl<<"-----------------------"<<endl;
			cout<<"Measurement Area "<<i<<endl;
			cout<<"-----------------------"<<endl;			
			m_MeasurementAreas.at(i).PrintParameters();
			if(i==m_MeasurementAreas.size()-1) cout<<endl;

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
	int num_rectangle=0; // Number of rectangular areasle, which include a whole grid of measurement areas

	

	TString MAname;
	MAname.Form("MA%i",num_MA+1);

	TString linename;
	linename.Form("line%i",num_line+1);

	TString rectanglename;
	rectanglename.Form("rectangle%i",num_rectangle+1);

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
		double X0=mEnv->GetValue(MAname+".X0", 88.97);
		int run_min=mEnv->GetValue(MAname+".minrunnumber", -1);
		int run_max=mEnv->GetValue(MAname+".maxrunnumber", -1);
		int max_angle=mEnv->GetValue(MAname+".maxanglenumber", -1);

		// Define measurement area based on these parameters
		MeasurementArea MA(ucenter,vcenter,ulength,vlength,thickness,X0,density,Z,A,run_min,run_max,max_angle);
  
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
				double X0=mEnv->GetValue(linename+".X0", 88.97);
				int run_min=mEnv->GetValue(linename+".minrunnumber", -1);
				int run_max=mEnv->GetValue(linename+".maxrunnumber", -1);
				int max_angle=mEnv->GetValue(linename+".maxanglenumber", -1);

				// Define measurement area based on these parameters
				MeasurementArea MA(ucenter,vcenter,ulength,vlength,thickness,X0,density,Z,A,run_min,run_max,max_angle);
			  
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
				double X0=mEnv->GetValue(linename+".X0", 88.97);
				int run_min=mEnv->GetValue(linename+".minrunnumber", -1);
				int run_max=mEnv->GetValue(linename+".maxrunnumber", -1);
				int max_angle=mEnv->GetValue(linename+".maxanglenumber", -1);

				// Define measurement area based on these parameters
				MeasurementArea MA(ucenter,vcenter,ulength,vlength,thickness,X0,density,Z,A,run_min,run_max,max_angle);
		  
				// Add measurement area to predefined grid
				m_MeasurementAreas.push_back(MA);
			}

		}
		// in this case there is nothing to do
		else continue;

		num_line++;
		linename.Form("line%i",num_line+1);
	}	

	while(mEnv->GetValue(rectanglename+".exist",0)!=0)
	{

		cout<<"rectanglename: "<<rectanglename<<endl;
		if((mEnv->GetValue(rectanglename+".usteplength",0.0)!=0.0)&&(mEnv->GetValue(rectanglename+".vsteplength",0.0)!=0.0)) 
		{

			cout<<"inside if statement"<<endl;
			for(double d_u=mEnv->GetValue(rectanglename+".startu",0.0);d_u<(mEnv->GetValue(rectanglename+".startu",0.0)+mEnv->GetValue(rectanglename+".ulength",-1.0));d_u+=mEnv->GetValue(rectanglename+".usteplength",10.0))
			{
				for(double d_v=mEnv->GetValue(rectanglename+".startv",0.0);d_v<(mEnv->GetValue(rectanglename+".startv",0.0)+mEnv->GetValue(rectanglename+".vlength",-1.0));d_v+=mEnv->GetValue(rectanglename+".vsteplength",10.0))
				{
					// Compute/Read out measurement area parameters
					double ucenter=d_u;
					double vcenter=d_v;
					double ulength=mEnv->GetValue(rectanglename+".usteplength", 0.1);
					double vlength=mEnv->GetValue(rectanglename+".vsteplength", 0.1);
					double thickness=mEnv->GetValue(rectanglename+".thickness", 1.8);
					double Z=mEnv->GetValue(rectanglename+".atomicnumber", 13.0);
					double A=mEnv->GetValue(rectanglename+".atomicmassnumber", 27.0);
					double density=mEnv->GetValue(rectanglename+".density", 2.7);
					double X0=mEnv->GetValue(rectanglename+".X0", 88.97);
					int run_min=mEnv->GetValue(rectanglename+".minrunnumber", -1);
					int run_max=mEnv->GetValue(rectanglename+".maxrunnumber", -1);
					int max_angle=mEnv->GetValue(rectanglename+".maxanglenumber", -1);

					// Define measurement area based on these parameters
					MeasurementArea MA(ucenter,vcenter,ulength,vlength,thickness,X0,density,Z,A,run_min,run_max,max_angle);
				  
					// Add measurement area to predefined grid
					m_MeasurementAreas.push_back(MA);
				}
			}
		}

		// in this case there is nothing to do
		else continue;

		num_rectangle++;
		rectanglename.Form("rectangle%i",num_rectangle+1);
	}

}

// Functions used in this script
double DetermineFitrange(TH1F*,double);
double calculateB(double);
double getrecoerror(TFile*);
void savehisto(TFile*,TFile*,TString,TString, double, double, double, double, int );
void correcthisto(TFile*,TFile*, TString , double , double , double , double );
double* fit(TFile*,TString, Grid, std::vector<double>, double, TString, std::vector<bool>);
void shiftbins(TH1F*, double);
//void calibrationfit();
int** GetParameterMapping(int);

  // Highland model of a MSC angle distribution, the parameters are:

  /*
	* par[0]:  Expected beam energy at u,v=0;
	* par[1]:  Beam particle charge
	* par[2]:  Beam particle mass
	* par[3]:  Target material density (not used here)
	* par[4]:  Target material atomic number (not used here)
	* par[5]:  Target material atomic weight (not used here)
	* par[6]:  Thickness
	* par[7]:  Expected angle reconstruction error
	* par[8]:  reco error calibration factor
	* par[9]:  Normalization
	* par[10]: u coordinate
	* par[11]: v coordinate
	* par[12]: u BE gradient
	* par[13]: v BE gradient
	* par[14]: mean of angle distribution
	* par[15]: Target material radiation length 

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
    //double Z;
    //Z=par[4];

	// atomic weight of target material
    //double A;
    //A=par[5];

	//density of the target material
    //double density;
    //density=par[3];

	// thickness of the target material
	double d1=par[6]; // in mm

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
	//double X0=716.4*A/((Z+1)*Z*density*TMath::Log(287.0/TMath::Sqrt(Z)));  	// This formula should only be used in case X0 is unknown (1-2% deviation from PDG value)
	double X0=par[15];

	// Combination of Highland width and reconstruction error
	double sigma=TMath::Sqrt(pow(recoerror,2)+pow(0.0136*charge/(p*beta)*TMath::Sqrt(d1/X0)*(1.0+0.038*TMath::Log(d1/X0)),2));

	// function value at a certain theta value
	double value=par[9]*TMath::Gaus(x[0],par[14],sigma);

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
	* par[14]:  mean of angle distribution
	* par[15]:  Target material radiation length (not used here)

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

	// Radiation length (not used here)
    //double X0=par[15];

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
		convolutionintegral+=h_total1->GetBinContent(i+1)*TMath::Gaus(phi_values_rad[i]-x[0],2*par[14],recoerror);
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
// a globalEstimator structure

// Number of parameters per fit function
const int num_localparameters=16;

// Number of new parameters per fit function
const int newparsperfunction=9;

// Create the GlobalCHi2 structure
struct GlobalEstimator { 
  std::vector< ROOT::Math::IMultiGenFunction * >& getEstimatorVec() {return fEstimatorVec;}
  void SetArray() {ipar=GetParameterMapping(fEstimatorVec.size());}
  
  double operator() (const double *par) const {
    double ret =0;
      
	for(size_t j=0;j<fEstimatorVec.size();j++)
	{
      // read function args
      double p1[num_localparameters];
      for (int i = 0; i < num_localparameters; ++i) 
      {
        p1[i] = par[ipar[j][i] ];
	  }
      // evaluate function and sum up
      ROOT::Math::IMultiGenFunction * func  =  fEstimatorVec.at(j);
      
      // Add up all local estimator values of the single fit functions
	  ret += (*func)(p1); 
    }
	// return global estimator value
    return ret;
  }
  
  std::vector< ROOT::Math::IMultiGenFunction * > fEstimatorVec;
  int ** ipar;
};

// Get Relation between global fit parameters and local fit parameters
// Parameter mapping
int **GetParameterMapping(int numfuncs)
{

	  int **ipar=0;

      ipar = new int*[numfuncs];

      for (int i = 0; i < numfuncs; i++)
      {
            ipar[i] = new int[num_localparameters];

            for (int j = 0; j < num_localparameters; j++)
            {
                  
				// There are basically 2 cases: The first function has the Parameters 0-15,
				//								the other function have a new Parameter number at the 3rd parameter (density), the 4th parameter (Z), the 5th parameter (A)
				//								the 6th Parameter (thickness of target material), the 9th Parameter (~#tracks), the 10th parameter (u coordinate), the 11th parameter (v coordinate),
				//								the 14th parameter (mean value) and the 15th parameter (X0 value)
						
				if(i==0) 
				{
					ipar[i][j] = j;
				}
				else
				{
					if((j!=3)&&(j!=4)&&(j!=5)&&(j!=6)&&(j!=9)&&(j!=10)&&(j!=11)&&(j!=14)&&(j!=15)) ipar[i][j]=j;
					else if(j==3) ipar[i][j]=num_localparameters+(i-1)*newparsperfunction;
					else if(j==4) ipar[i][j]=num_localparameters+1+(i-1)*newparsperfunction;
					else if(j==5) ipar[i][j]=num_localparameters+2+(i-1)*newparsperfunction;
					else if(j==6) ipar[i][j]=num_localparameters+3+(i-1)*newparsperfunction;
					else if(j==9) ipar[i][j]=num_localparameters+4+(i-1)*newparsperfunction;
					else if(j==10) ipar[i][j]=num_localparameters+5+(i-1)*newparsperfunction;
					else if(j==11) ipar[i][j]=num_localparameters+6+(i-1)*newparsperfunction;
					else if(j==14) ipar[i][j]=num_localparameters+7+(i-1)*newparsperfunction;
					else ipar[i][j]=num_localparameters+8+(i-1)*newparsperfunction;
				}
				cout<<"Parameter mapping: Global parameter number of local parameter "<<j<<" in fit function "<<i<<" is "<<ipar[i][j]<<endl;
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
double DetermineFitrange(TH1F* histo,double rangevalue)
{
  // Clone histo
  TH1F *h2 = (TH1F*) histo->Clone();

  cout<<"RMS value of distribution: "<<histo->GetRMS()<<endl;
  cout<<"Selected range parameter: "<<rangevalue<<" -> Fit range up to y=1/("<<rangevalue<<"*e)"<<endl;
  double fitrange = sqrt(2.0*rangevalue)*histo->GetRMS();
  
  // Use RMS value as a rough measure of the fit range for a gaussian fit
  TF1 *f1 = new TF1("f1","gaus(x)",-fitrange,fitrange);
  f1->SetParameter(2,histo->GetRMS());
  f1->FixParameter(1,0);
  f1->SetParLimits(2,0.0,2*histo->GetRMS());
  f1->SetLineStyle(2);
  TFitResultPtr fitr=h2->Fit("f1","RS");

  cout<<endl<<"Fit results:  "<<endl;
  cout<<"Fit is valid: "<<fitr->IsValid()<<endl;
  cout<<"Fit status: "<<fitr->Status()<<endl;
  
  // Repeat fit in case it failed
  if(!fitr->IsValid())
  {
    cout<<"Fit of angle distribution failed with status: "<<fitr<<endl;
    cout<<"Repeat fit "<<endl;
    h2->Fit("f1","RM");

	cout<<endl<<"Fall back Fit results: "<<endl;
	cout<<"Fit is valid: "<<fitr->IsValid()<<endl;
	cout<<"Fit status: "<<fitr->Status()<<endl;
  }

  // Use the determined sigma value to calculate the fit range
  double sigma = f1->GetParameter(2);

  // In case the second fit failed as well, simply use the RMS of the histo in the 3 sigma range
  if(!fitr->IsValid()) 
  {
	double corerange=3*histo->GetRMS();

	// Limit histo range to core of distribution
	histo->GetXaxis()->SetRangeUser(-corerange,corerange);

	// Get RMS of core distribution without outliers
	sigma = histo->GetRMS();
    cout<<"Fall back fit of angle distribution failed with status: "<<fitr<<"! Just use RMS of histogram."<<endl;
  }

  fitrange=sqrt(2.0*rangevalue)*sigma;
  cout<<"Determined fit range: " << fitrange<<endl<<endl;

  delete h2;
  delete f1;
  return fitrange;
}


// Returns the mean value of the angle reco error squared
double getanglerecovar(TFile* file)
{
	//TTree in input root file, that contains the MSC projected angle distributions and reconstruction error distribution
	file->cd("");
	TTree *msc_tree = (TTree*)file->Get("MSCTree");

	// Draw reconstruction error 1 histogram
	msc_tree->Draw("theta1_var>>h_help_theta1_var", "", "P*");
    TH1F *h_help_theta1_var = (TH1F*)gDirectory->Get("h_help_theta1_var");
    
	// Get mean value 1
	double recovar1 = h_help_theta1_var->GetMean();

	// Get mean error 1
	double recovar1_error = h_help_theta1_var->GetMeanError();

	// Draw reconstruction error 1 histogram
	msc_tree->Draw("theta2_var>>h_help_theta2_var", "", "P*");
    TH1F *h_help_theta2_var = (TH1F*)gDirectory->Get("h_help_theta2_var");

	// Get mean value 2
	double recovar2 = h_help_theta2_var->GetMean();

	// Get mean error 2
	double recovar2_error = h_help_theta2_var->GetMeanError();

	// Calculate mean value
	double recovar=(recovar1+recovar2)/2.0;
	double recovar_error=sqrt(recovar1_error*recovar1_error+recovar2_error*recovar2_error)/2.0;

	// The recovar error shouldn't be too large! 
	cout<<endl<<"Angle Reconstruction Variance: 	"<<recovar<<" +/- "<<recovar_error<<"rad^2"<<endl; 

	return recovar;
}


// Shift a decentralized angle distribution in order to get a mean value of 0
void shiftbins(TH1F* histogram, double mean1)
{
	// The mean offset has to be corrected before merging the two histograms
	// First the bin value of the offset has to be found out
	TAxis *xaxis =  histogram->GetXaxis();
	Int_t binx = xaxis->FindBin(mean1)-xaxis->FindBin(0.0);

	Double_t stats[5]={0,0,0,0,0};
	histogram->PutStats(stats); // reset mean value, etc

	// Get correct number of histogram entries
	int n_entries=histogram->GetEntries();

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
		   	 histogram->SetBinContent(j,0);
		}
	}

	histogram->SetEntries(n_entries);
}

// Function to save the projected angle histograms of different regions in the u-v plane, the region lies within the given 
// u and v min and max values.
void savehisto(std::vector<TString> inputfiles, TFile* file2, TString histoname, int numbins, double histo_range, std::vector <double> cutparameters_position, std::vector <int> cutparameters_run, std::vector <int> cutparameters_vertex_multiplicity, int max_angles, int correctmean, TString nametag)
{

    // Find relevant input files for the given run number range
    std::vector<TString> relevantfiles;
    for(size_t ifile=0;ifile<inputfiles.size();ifile++)
	{
        TString tmpstring;
        tmpstring=inputfiles.at(ifile);
		tmpstring.Remove(0,3);
		int runnumber=tmpstring.Atoi();

		if((cutparameters_run.at(0)<0 && cutparameters_run.at(1)<0 ) || (runnumber>=cutparameters_run.at(0)&&runnumber<=cutparameters_run.at(1)))
		{
			relevantfiles.push_back(inputfiles.at(ifile));
		}
	}

	TFile* tmpfile=new TFile(relevantfiles.at(0),"READ");

	//TTree in input root file, that contains the MSC projected angle distributions
	tmpfile->cd("");
	TTree *tmp_tree = (TTree*)tmpfile->Get("MSCTree");

	TString areacuts;
    areacuts.Form("(u>%.2f)&&(u<%.2f)&&(v>%.2f)&&(v<%.2f)",cutparameters_position.at(0),cutparameters_position.at(1),cutparameters_position.at(2),cutparameters_position.at(3));
    
	// Draw histogram of the first scattering angle in the given u and v range and save it
	tmp_tree->Draw("theta1>>h_help",areacuts,"");
    TH1F *h_help = (TH1F*)gDirectory->Get("h_help");
    
	double limits=histo_range*(h_help->GetRMS());
    
    cout << "Helper histogram RMS: " << h_help->GetRMS()  << endl; 
    
	tmpfile->Close();
    
	// Temporary histos with angle distributions fulfilling the cut conditions
	TH1F * tmp_anglehisto[2];
	tmp_anglehisto[0]=new TH1F("theta1_histo","theta1_histo",numbins,-limits,limits);
	tmp_anglehisto[1]=new TH1F("theta2_histo","theta2_histo",numbins,-limits,limits);

    for(size_t ifile=0;ifile<relevantfiles.size();ifile++)
	{

		TFile* inputfile=new TFile(relevantfiles.at(ifile),"READ");

		//TTree in input root file, that contains the MSC projected angle distributions
		inputfile->cd("");
		TTree *msc_tree = (TTree*)inputfile->Get("MSCTree");

		// Array of mean theta1 and theta2 values in each map pixel
		double mean1=0.0;
		double mean2=0.0;

		// parameters which are read out from the root file
		Double_t theta1;
		Double_t theta2;
		Double_t u;
		Double_t v;

		Int_t RunNo;
		Int_t vertex_multiplicity=1;

		// Set branch adresses for parameters connected to the scattering angles
		msc_tree->SetBranchAddress("theta1",&theta1);
		msc_tree->SetBranchAddress("theta2",&theta2);
		msc_tree->SetBranchAddress("u",&u);
		msc_tree->SetBranchAddress("v",&v);

		msc_tree->SetBranchAddress("iRun",&RunNo);
		msc_tree->SetBranchAddress("vertex_multiplicity",&vertex_multiplicity);

		if(correctmean==1)
		{
			// Get histogram
			TH1F* histogram1=(TH1F*)file2->Get("grid/raw/theta1_uncorrected_"+histoname);
			TH1F* histogram2=(TH1F*)file2->Get("grid/raw/theta2_uncorrected_"+histoname);

			// Save mean of both distributions in arrays
			mean1=histogram1->GetMean();
			mean2=histogram2->GetMean();

			limits=histo_range/2.0*(histogram1->GetRMS()+histogram2->GetRMS());
		}

		for(int ientry=0;ientry<msc_tree->GetEntries();ientry++)
		{

			if(ientry%100000==0) cout<<"Tree entry "<<ientry<<endl;
			msc_tree->GetEntry(ientry);

			bool u_condition,v_condition;	
			bool run_condition=true;
			bool vertex_multiplicity_condition=true;

			u_condition=(u>cutparameters_position.at(0))&&(u<cutparameters_position.at(1));
			v_condition=(v>cutparameters_position.at(2))&&(v<cutparameters_position.at(3));

			if(cutparameters_run.at(0)>-1 && cutparameters_run.at(1)>-1 ) run_condition=(RunNo>=cutparameters_run.at(0))&&(RunNo<=cutparameters_run.at(1));
			else if (cutparameters_run.at(0)>-1) run_condition=(RunNo>=cutparameters_run.at(0));
			else if (cutparameters_run.at(1)>-1) run_condition=(RunNo<=cutparameters_run.at(1));

			if(msc_tree->GetBranchStatus("vertex_multiplicity")) vertex_multiplicity_condition=(vertex_multiplicity>=cutparameters_vertex_multiplicity.at(0))&&(vertex_multiplicity<=cutparameters_vertex_multiplicity.at(1));

			if(u_condition&&v_condition&&run_condition&&vertex_multiplicity_condition)
			{
					tmp_anglehisto[0]->Fill(theta1-mean1);
					tmp_anglehisto[1]->Fill(theta2-mean2);
			}

			if((tmp_anglehisto[0]->GetEntries()>=max_angles)&&(max_angles>-1)) break;
		}
	}

	// Give the two histograms to a list Operator
	TList *list = new TList;
	TH1F* hsum1=(TH1F*)tmp_anglehisto[0]->Clone("hsum1");
	TH1F* hsum2=(TH1F*)tmp_anglehisto[1]->Clone("hsum2");
    list->Add(hsum1);
    list->Add(hsum2);

	// Merge all histograms in the list
    TH1F *h = (TH1F*)hsum1->Clone("h");
    h->Reset();
    h->Merge(list);
	h->SetTitle("merged projected angle distribution");
	h->GetXaxis()->SetTitle("#theta_{proj.} [rad]");

	// Set histo axis titels etc
	tmp_anglehisto[0]->SetTitle("#theta_{1} distribution");
	tmp_anglehisto[0]->GetXaxis()->SetTitle("#theta_{1} [rad]");

	tmp_anglehisto[1]->SetTitle("#theta_{2} distribution");
	tmp_anglehisto[1]->GetXaxis()->SetTitle("#theta_{2} [rad]");
/*
	// Shift the bins of the histogram to reduce the offset, but only in case the offset is larger than the histogram bin size
	if((correctmean==1)&&(abs(mean1)>(2*limits/numbins))) 
	{
		cout<<endl<<"Correct histogram 1 offset"<<endl;
		shiftbins(tmp_anglehisto[0],mean1);
	}

	// Shift the bins of the histogram to reduce the offset, but only in case the offset is larger than the histogram bin size
	if((correctmean==1)&&(abs(mean2)>(2*limits/numbins))) 
	{
		cout<<endl<<"Correct histogram 2 offset"<<endl;
		shiftbins(tmp_anglehisto[1],mean2);
	}

	// Shift the bins of the histogram to reduce the offset, but only in case the offset is larger than the histogram bin size
	if((correctmean==1)&&(abs(mean1+mean2)>(2*limits/numbins))) 
	{
		cout<<endl<<"Correct histogram sum offset"<<endl;
		shiftbins(h,0.5*(mean1+mean2));
	}
*/
	// Go to raw directory
	file2->cd("");
	file2->cd("grid/raw/");

	// Write histogram to disk
	if(correctmean==1) tmp_anglehisto[0]->Write("theta1_"+nametag+histoname);
	else tmp_anglehisto[0]->Write("theta1_"+nametag+histoname);

	// Write histogram to disk
	if(correctmean==1) tmp_anglehisto[1]->Write("theta2_"+nametag+histoname);
	else tmp_anglehisto[1]->Write("theta2_"+nametag+histoname);

	// Write histogram to disk
	if(correctmean==1) h->Write("sumhisto_"+nametag+histoname);
	else h->Write("sumhisto_"+nametag+histoname);

	// Delete sum histogram from memory
	h->SetDirectory(gROOT);
	delete h;

	// Delete sum histogram from memory
	hsum1->SetDirectory(gROOT);
	delete hsum1;

	// Delete sum histogram from memory
	hsum2->SetDirectory(gROOT);
	delete hsum2;

	// Delete histogram1 from memory
	tmp_anglehisto[0]->SetDirectory(gROOT);
	delete tmp_anglehisto[0];

	// Delete histogram1 from memory
	tmp_anglehisto[1]->SetDirectory(gROOT);
	delete tmp_anglehisto[1];		

}

// Function to fit the MSC angle histograms simultaneously, it returns a pointer to the fit results
double* fit( TFile* file, Grid grid, std::vector<double> beamoptions, double recoerr, TString model, std::vector<bool> fitoptions, double rangevalue)
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
	bool Use_LogLikelihoodfit=fitoptions.at(4);


	// Some parameter definitions

	double lambda_startvalue=beamoptions.at(5);

	// histogram name
	TString histoname;

	// Fit range
	double fitrange;

	const int num_fitfunctions=grid.GetMeasurementAreas().size();

	// fitresults Array
	double *fitresults = new double[8];

	// Fit and raw angle histograms
	std::vector<TH1F *> histo_vec;
	TH1F * histogramsum;

	// Copy raw histograms
	for(int i=0; i< num_fitfunctions; i++)
	{
		// histogram name
		histoname.Form("measurementarea%i",i+1);
		file->cd("");

		histogramsum=(TH1F*)file->Get("grid/raw/sumhisto_"+histoname);

		// Make a copy of the histogram
		TH1F* fithistogramsum;
		fithistogramsum=(TH1F*)histogramsum->Clone("sumhisto_"+histoname+"_fit");
		histo_vec.push_back(fithistogramsum);

		file->cd("grid/fit/");
	}
	
	// The fit range depend on the model: In Case of the moliere model the fitrange can be a little larger
	// Declaration of fit functions
	std::vector<TF1 *> fitFcn_vec;
	TF1* fitFcn;
	TString fctname;
    
	std::vector<double> range_vec;

	// loop for definition of the fit functions, the fitrange is determined for every one of them
	for(int i=0;i<num_fitfunctions;i++)
	{

		fitrange=DetermineFitrange(histo_vec.at(i),rangevalue);
		fctname.Form("fitFcn%i",i);

		if(model=="moliere")
		{
			// Use Gaussian function with width corresponding to the calibrated angle resolution, if material is only air or extremely thin
			// Also in this case the fit range is set to be a little smaller
			if(grid.GetMeasurementAreas().at(i).Get_thickness()<0.0002) fitFcn = new TF1(fctname,"[9]*TMath::Gaus(x,0.0*[0]*[1]*[2]*[3]*[4]*[5]*[6]*[10]*[11]*[12]*[13]*[15]+[14],[7]*[8])",-fitrange,fitrange);


			// Use Moliere model in case the material is not just air
			else fitFcn = new TF1(fctname,molierefunction,-fitrange,fitrange,num_localparameters);
		}

		else
		{
			// Use Gaussian function with width corresponding to the calibrated angle resolution, if material is only air or extremely thin
			if(grid.GetMeasurementAreas().at(i).Get_thickness()<0.0002) fitFcn = new TF1(fctname,"[9]*TMath::Gaus(x,0.0*[0]*[1]*[2]*[3]*[4]*[5]*[6]*[10]*[11]*[12]*[13]*[15]+[14],[7]*[8])",-fitrange,fitrange);
			// Use Highland model in case the material is not just air
			else fitFcn = new TF1(fctname,highlandfunction,-fitrange,fitrange,num_localparameters);
		}

		// Fill vector with pointers to fit functions
		fitFcn_vec.push_back(fitFcn);
        
		// set the data range
        cout<<"Fit function number: "<<i<<endl;
		cout<<"The range value from the cfg file is "<<rangevalue<<endl;	
		cout<<"Fitrange at d="<<grid.GetMeasurementAreas().at(i).Get_thickness()<<" mm: "<<fitrange<<" rad"<<endl<<endl;
		
        range_vec.push_back(fitrange);
	}
    
    cout<<"Done reading definition of fitfunctions"<< endl;
    
    // Create globalChi2 object
	GlobalEstimator globalestimator;
    // Total datasize
	int datasize = 0;
     
    for(int i=0;i<num_fitfunctions;i++) 
	{ 
      // Create wrapped multi function entry
	  ROOT::Math::WrappedMultiTF1 * wf_entry = new ROOT::Math::WrappedMultiTF1(*fitFcn_vec.at(i),1);
      
      // Create data entry and fill fit data 
      ROOT::Fit::DataOptions opt;
      ROOT::Fit::DataRange range;
      range.SetRange(-range_vec.at(i), range_vec.at(i));
        
      TH1F * histoptr =  histo_vec.at(i);
        
      ROOT::Fit::BinData * data_entry = new ROOT::Fit::BinData(opt,range); 
      ROOT::Fit::FillData(*data_entry, histoptr);
      datasize+=data_entry->Size();
       
      // Create estimator function entry and fill global estimator 
      if(Use_LogLikelihoodfit) 
      {
        ROOT::Fit::PoissonLLFunction* log_entry = new ROOT::Fit::PoissonLLFunction(*data_entry, *wf_entry);
        globalestimator.getEstimatorVec().push_back(log_entry); 
      }
      else 
      {
        ROOT::Fit::Chi2Function* chi2_entry = new ROOT::Fit::Chi2Function(*data_entry, *wf_entry);
        globalestimator.getEstimatorVec().push_back(chi2_entry); 
      }       
    } 
     
	// Set parameter mapping array in the global chi2 object
	globalestimator.SetArray();
    
    cout<<"Prepare the fitter for the calibration fit"<< endl;
    
	ROOT::Fit::Fitter fitter;

	// Number of global parameters: There are num_fitfunctions fit functions with three new parameters each: Scale, u coordinate, vcoordinate, thickness
	// Then there are 8 parameters, which are used in all fit functions (beam energy calibration factor, angle calibration factor, material properties, etc)
	static const int num_globalparameters = num_fitfunctions*newparsperfunction+(num_localparameters-newparsperfunction);

    std::vector<double> par0(num_globalparameters,0);

	// Set the entries of the globalparameters array

	// First fit function has num_localparameters new parameters:
	double aid_array[num_localparameters]={ BE_mean,z,mass,grid.GetMeasurementAreas().at(0).Get_density(),grid.GetMeasurementAreas().at(0).Get_Z(),grid.GetMeasurementAreas().at(0).Get_A(),
											grid.GetMeasurementAreas().at(0).Get_thickness(),recoerr,lambda_startvalue,700.0,grid.GetMeasurementAreas().at(0).Get_u_center(),
											grid.GetMeasurementAreas().at(0).Get_v_center(),BE_ugrad,BE_vgrad,0.0, grid.GetMeasurementAreas().at(0).Get_X0()};
	for(int i=0;i<num_localparameters;i++) par0[i]=aid_array[i];

	// Afterwards for each fit functions we get several new parameters
	for(int i=1;i<num_fitfunctions;i++)
	{
		par0[num_localparameters+(i-1)*newparsperfunction]=grid.GetMeasurementAreas().at(i).Get_density();
		par0[num_localparameters+1+(i-1)*newparsperfunction]=grid.GetMeasurementAreas().at(i).Get_Z();
		par0[num_localparameters+2+(i-1)*newparsperfunction]=grid.GetMeasurementAreas().at(i).Get_A();
		par0[num_localparameters+3+(i-1)*newparsperfunction]=grid.GetMeasurementAreas().at(i).Get_thickness();
		par0[num_localparameters+4+(i-1)*newparsperfunction]=700;
		par0[num_localparameters+5+(i-1)*newparsperfunction]=grid.GetMeasurementAreas().at(i).Get_u_center();
		par0[num_localparameters+6+(i-1)*newparsperfunction]=grid.GetMeasurementAreas().at(i).Get_v_center();
		par0[num_localparameters+7+(i-1)*newparsperfunction]=0.0;
		par0[num_localparameters+8+(i-1)*newparsperfunction]=grid.GetMeasurementAreas().at(i).Get_X0();
	}

	// create before the parameter settings in order to fix or set range on them
    fitter.Config().SetParamsSettings(num_globalparameters,&par0[0]);

	// fix constant parameters 1 to 7
	for(int i=1;i<8;i++) fitter.Config().ParSettings(i).Fix();

	// fix X0 constant
	fitter.Config().ParSettings(15).Fix();

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

	// fix X0, density, A, Z, thickness and coordinate parameters for all fitfunctions
	for(int i=1;i<num_fitfunctions;i++)
	{
		fitter.Config().ParSettings(num_localparameters+(i-1)*newparsperfunction).Fix();
		fitter.Config().ParSettings(num_localparameters+1+(i-1)*newparsperfunction).Fix();
		fitter.Config().ParSettings(num_localparameters+2+(i-1)*newparsperfunction).Fix();
		fitter.Config().ParSettings(num_localparameters+3+(i-1)*newparsperfunction).Fix();
		fitter.Config().ParSettings(num_localparameters+5+(i-1)*newparsperfunction).Fix();
		fitter.Config().ParSettings(num_localparameters+6+(i-1)*newparsperfunction).Fix();
		fitter.Config().ParSettings(num_localparameters+8+(i-1)*newparsperfunction).Fix();
	}

	if(fixlambda) fitter.Config().ParSettings(8).Fix();				 //fix lambda?

	// set limits on the calibration factor parameter
	fitter.Config().ParSettings(8).SetLimits(0.6,1.4);

	//Set limits on fit function Normalizations and mean angle values
	fitter.Config().ParSettings(9).SetLimits(10.0,200000.0);
	fitter.Config().ParSettings(14).SetLimits(-0.0001,+0.0001);
	for(int i=1;i<num_fitfunctions;i++) 
	{
		fitter.Config().ParSettings(num_localparameters+4+(i-1)*newparsperfunction).SetLimits(10.0,200000.0);
		fitter.Config().ParSettings(num_localparameters+7+(i-1)*newparsperfunction).SetLimits(-0.0001,+0.0001);
	}

	fitter.Config().MinimizerOptions().SetPrintLevel(1);
	fitter.Config().SetMinimizer("Minuit2","Migrad"); 
    
    cout<<"Perform the actual fit of calibration data"<< endl;
    
	// fit FCN function directly 
	// (specify optionally data size and flag to indicate that is a chi2 fit)
	fitter.FitFCN(num_globalparameters,globalestimator,0,datasize,!Use_LogLikelihoodfit);
    
    cout<<"Done with fitting. Read back the fit results"<< endl;
    
	ROOT::Fit::FitResult result = fitter.Result();
	result.Print(std::cout);

	int ** parameter_mapping=GetParameterMapping(num_fitfunctions);


	// Names of the 14 local parameters
	TString name[num_localparameters];
	name[0]="E[GeV]";
	name[1]="z[e]";
	name[2]="m[GeV/c^2]";
	name[3]="#rho[g/cm^3]";
	name[4]="Z";
	name[5]="A";
	name[6]="X[mm]";
	name[7]="#sigma_{err}[rad]";
	name[8]="#lambda";
	name[9]="norm";
	name[10]="u[mm]";
	name[11]="v[mm]";
	name[12]="#nablaE_{u}[GeV/mm]";
	name[13]="#nablaE_{v}[GeV/mm]";
	name[14]="#theta_{mean}[rad]";
	name[15]="X_{0}[mm]";

	for(int i=0;i<num_fitfunctions;i++)
	{	

		int parameters[num_localparameters];
		for(int j=0;j<num_localparameters;j++)
		{
			parameters[j]=parameter_mapping[i][j];
		}

		fitFcn_vec.at(i)->SetFitResult( result, parameters);
	
		fitrange=range_vec.at(i);  
		TF1 * fitfunc=fitFcn_vec.at(i);
		fitfunc->SetRange(-fitrange,fitrange);  
		fitfunc->SetLineColor(kRed);

		histo_vec.at(i)->GetListOfFunctions()->Add(fitFcn_vec.at(i));
		histoname.Form("gridpoint%i",i+1);
		// display mode
		gStyle->SetOptFit(1111);
		histo_vec.at(i)->Write("thetasum_"+histoname+"_fit");

		for(int iname=0;iname<num_localparameters;iname++) fitfunc->SetParName(iname,name[iname]);
	}
    
    cout<<"Plot the fit results, print them in terminal and save them to histogram"<< endl;
    
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
    //double chi2_summation=0.0;

    for(size_t i=0; i<num_fitfunctions;i++)
	{
        fctname.Form("fitFcn%lu",i);

		cout<<"--------------------------"<<endl;
		cout<<"Fit function "<<i<<endl;
		cout<<"Chi2 value for this individual fit is: "<<histo_vec.at(i)->Chisquare(fitFcn_vec.at(i),"R")<<endl;
		cout<<"Number of degrees of freedom is: "<<fitFcn_vec.at(i)->GetNDF()<<endl;
		cout<<"--------------------------"<<endl;
	}
	
	TString pdfname,Title;

	for(int i=0;i<TMath::Ceil(double(num_fitfunctions)/4.0);i++)
	{
		TString canvasname;
		canvasname.Form("cfits_%i",i);
		TCanvas *c = new TCanvas(canvasname,canvasname,900,1000);
		std::vector<TPad*> pads;
		TPad *pad1 = new TPad("pad1","pad1",0.01,0.51,0.49,0.99);
		pad1->Draw();
		pads.push_back(pad1);
		TPad *pad2 = new TPad("pad2","pad2",0.51,0.51,0.99,0.99);
		pad2->Draw();
		pads.push_back(pad2);
		TPad *pad3 = new TPad("pad3","pad3",0.01,0.01,0.49,0.49);
		pad3->Draw();
		pads.push_back(pad3);
		TPad *pad4 = new TPad("pad4","pad4",0.51,0.01,0.99,0.49);
		pad4->Draw();
		pads.push_back(pad4);

		for(int j=0;j<4;j++)
		{
		   	pads.at(j)->cd();
			if(((4*i)+j)<num_fitfunctions)
			{
				Title.Form("Area %i: d=%fmm",(4*i)+j,grid.GetMeasurementAreas().at((4*i)+j).Get_thickness());
				histo_vec.at((4*i)+j)->SetTitle(Title);
				histo_vec.at((4*i)+j)->Draw();
                cout<<"fitfunction "<<(4*i)+j<<" of "<<num_fitfunctions<<endl;
            }
		}

		pdfname=model+"_results_"+canvasname+".pdf";

		c->SaveAs(pdfname); 
	}
      
	// Do a self consistency check:
	// Check whether a fit with the determined calibration parameters yields the correct radiation length value for all the
    // angle distributions used in the calibration fit
    cout <<"Start self consistency fit"<< endl;

	TString fitoption;
	if(model=="moliere") fitoption="RMELS";
	else fitoption="RMES";

	// Declaration of another set of fit functions, these will be used for the self consistency check
	std::vector<TF1 *> fitFcn_vec2;

	TH1F* h_d = new TH1F("h_d","h_d",num_fitfunctions,0.5,num_fitfunctions+0.5);
	TH1F* h_d_true = new TH1F("h_d_true","h_d_true",num_fitfunctions,0.5,num_fitfunctions+0.5);
	h_d_true->SetMinimum(-0.2);
	h_d_true->GetXaxis()->SetTitle("Measurement area");
	if(num_fitfunctions<15) 
	{
		h_d_true->GetXaxis()->SetNdivisions(num_fitfunctions);
	}
	else
	{
		h_d_true->GetXaxis()->SetNdivisions(num_fitfunctions/4);
	}
	h_d_true->GetYaxis()->SetTitle("Thickness [mm]");
	h_d_true->SetTitle("Self-consistency check");
    
	// loop for definition of the fit functions, the fitrange is determined for every one of them
    for(size_t i=0;i<num_fitfunctions;i++)
	{

		fitrange=range_vec.at(i);
        fctname.Form("fitFcn%lu",i);
        
		if(model=="moliere") {
          fitFcn = new TF1(fctname,molierefunction,-fitrange,fitrange,num_localparameters);
        } else {
          fitFcn = new TF1(fctname,highlandfunction,-fitrange,fitrange,num_localparameters);              
        }
         
		// Set starting values
		double lambda=fitresults[0];
		double BE_mean=fitresults[2];
		double BE_ugrad=fitresults[4];
		double BE_vgrad=fitresults[6];
        std::vector<double> parameters_temp{ BE_mean,z,mass,grid.GetMeasurementAreas().at(i).Get_density(),grid.GetMeasurementAreas().at(i).Get_Z(),grid.GetMeasurementAreas().at(i).Get_A(),
                    1.0,recoerr,lambda,700.0,grid.GetMeasurementAreas().at(i).Get_u_center(),
                    grid.GetMeasurementAreas().at(i).Get_v_center(),BE_ugrad,BE_vgrad,0.0,grid.GetMeasurementAreas().at(i).Get_X0()};
        parameters_temp.resize(num_localparameters);
        fitFcn->SetParameters(&parameters_temp[0]);

        for(int ii=0; ii<num_localparameters;ii++)
		{
            if(ii!=6&&ii!=9&&ii!=14)
			{
                fitFcn->FixParameter(ii,parameters_temp[ii]);
			}
		}	
 
   		fitFcn->SetParLimits(6,1E-8,200.0);
   		fitFcn->SetParLimits(14,-1E-4,+1E-4);

		// Fill vector with pointers to fit functions
		fitFcn_vec2.push_back(fitFcn);
	}

    

    for(size_t i=0;i<TMath::Ceil(double(num_fitfunctions)/4.0);i++)
	{
		TString canvasname;
        canvasname.Form("c%lu",i);
		TCanvas *c = new TCanvas(canvasname,canvasname,900,1000);
		std::vector<TPad*> pads;
		TPad *pad1 = new TPad("pad1","pad1",0.01,0.51,0.49,0.99);
		pad1->Draw();
		pads.push_back(pad1);
		TPad *pad2 = new TPad("pad2","pad2",0.51,0.51,0.99,0.99);
		pad2->Draw();
		pads.push_back(pad2);
		TPad *pad3 = new TPad("pad3","pad3",0.01,0.01,0.49,0.49);
		pad3->Draw();
		pads.push_back(pad3);
		TPad *pad4 = new TPad("pad4","pad4",0.51,0.01,0.99,0.49);
		pad4->Draw();
		pads.push_back(pad4);

        for(size_t j=0;j<4;j++)
		{
		   	pads.at(j)->cd();
			if(((4*i)+j)<num_fitfunctions)
			{
                Title.Form("Area %lu: d=%fmm",(4*i)+j,grid.GetMeasurementAreas().at((4*i)+j).Get_thickness());
				histo_vec.at((4*i)+j)->SetTitle(Title);


				TFitResultPtr fitr=histo_vec.at((4*i)+j)->Fit(fitFcn_vec2.at((4*i)+j),fitoption);
				if(fitr!=0)
				{
					cout<<"Fit of angle distribution failed with status: "<<fitr<<endl;
					cout<<"Repeat fit "<<endl;
					histo_vec.at((4*i)+j)->Fit(fitFcn_vec2.at((4*i)+j),fitoption);
				}
				
				
				histo_vec.at((4*i)+j)->Draw();
                cout<<"fitfunction "<<(4*i)+j+1<<" of "<<num_fitfunctions<<endl;

				h_d->SetBinContent((4*i)+j+1,fitFcn_vec2.at((4*i)+j)->GetParameter(6));
				h_d->SetBinError((4*i)+j+1,fitFcn_vec2.at((4*i)+j)->GetParError(6));
				h_d_true->SetBinContent((4*i)+j+1,grid.GetMeasurementAreas().at((4*i)+j).Get_thickness());
            }
		}

		pdfname="x0fit_"+model+"_"+canvasname+".pdf";

		c->SaveAs(pdfname); 
	}

	TCanvas *c = new TCanvas("c1","c1",900,1000);
	gStyle->SetOptStat(0);
	h_d->SetLineColor(2);
	h_d_true->Draw("");
	h_d->Draw("same");
   	TLegend* legend = new TLegend(0.55,0.15,0.85,0.25);
   	legend->AddEntry(h_d,"Measured values","l");
   	legend->AddEntry(h_d_true,"Truth values","l");
   	legend->Draw();
	c->SaveAs("selfconsistency.pdf");

	file->cd("selfconsistency/");
	h_d->Write("d_measured");
	h_d_true->Write("d_truth");

	// Create a results root file and save the fit results in a histogram
	TFile *resultsfile = new TFile("X0calibration_results.root", "RECREATE");
	TH1F * resultshist=new TH1F("resultshist","results of calibration",4,1,4);

    resultshist->GetXaxis()->SetBinLabel( 1, "#lambda" );
    resultshist->GetXaxis()->SetBinLabel( 2, "E_{mean}[GeV]" );
    resultshist->GetXaxis()->SetBinLabel( 3, "#nablaE_{u}[GeV/mm]" );
    resultshist->GetXaxis()->SetBinLabel( 4, "#nablaE_{v}[GeV/mm]" );
    
	// Save the results to histogram
	for(int i=0;i<4;i++)
	{
		resultshist->SetBinContent(i+1,fitresults[2*i]);
		resultshist->SetBinError(i+1,fitresults[2*i+1]);
	}

	resultshist->Write("results");
	resultsfile->Close();
	
	// Return pointer to fitresults
	return fitresults;
} // End of fit function

void GetInputFiles(std::vector<TString>& inputfiles, const char *dirname=".", const char *ext="run") 
{ 
	TSystemDirectory dir(dirname, dirname); 
	TList *files = dir.GetListOfFiles(); 
	if (files) 
	{ 
		TSystemFile *file; 
		TString fname; 
		TIter next(files); 
		while ((file=(TSystemFile*)next())) 
		{ 
			fname = file->GetName(); 
			if (!file->IsDirectory() && fname.Contains(ext)) 
			{ 	
                inputfiles.push_back(fname); 
			} 
		} 
	} 
}


  // Function used to get the angle distributions of tracks crossing certain regions in the u-v measurement plane,
  // perform a moliere fit on the distributions to estimate the calibrationfactors mu and lambda. Usually this script is used on 
  // measurement data of a plane with a precisely known material distribution ( for example a aluminum grid with a set 
  // of holes with different thicknesses.
  int main(int , char **)
  //void calibrationfit()
  {     
    // Read config file
    //------------------
    TEnv *mEnv=new TEnv("x0.cfg");

	// Choose the multiple scattering model
	// "moliere": Moliere model
	// "highland": Highland model
	TString model=mEnv->GetValue("model", "highland");

	// Fix lambda or mu or slope parameter during fit?
	bool fixlambda=mEnv->GetValue("fixlambda", 0);
    bool fix_offset=mEnv->GetValue("fix_momentumoffset", 0);
	bool fix_u_gradient=mEnv->GetValue("fix_momentumugradient", 0);
	bool fix_v_gradient=mEnv->GetValue("fix_momentumvgradient", 0);

	// Use log likelihood estimator?
	// True: Use log likelihood estimator
	// False: Use Chi2 estimator
	bool Use_loglikelihood_estimator=mEnv->GetValue("use_loglikelihood", 1);

	// The number of bins of angle histograms
	int numbins=mEnv->GetValue("cali_num_bins", 50);

	// Range parameter of angle histograms
	double histo_range=mEnv->GetValue("cali_histo_range", 5.0);

	// Read out the parameter that determines the range of the fit
	double rangevalue=mEnv->GetValue("fitrange_parameter", 1.0);

	// Integer determining, whether a offset correction is applied to the angular distributions
	// 1: correction
	// everything else, no correction
	int correctmean=mEnv->GetValue("correctmean", 1 );

	// Get the lambda starting value
	double 	lambda_start_default=mEnv->GetValue("lambda", 99.0);

	// Define and set the parameters used in the Moliere fit
	// BE: beam energy (GeV), z: charge of beam particle (e), mass: mass of beam particle (GeV),
	// d_layer: thickness per layer(mm), num_layers: total number of layers, recoerr: angle reconstruction error (rad)
    double z,mass,recoerr;
	
	cout<<"---------------------------------------------------"<<endl;
	cout<<"-----------------Parameter settings----------------"<<endl;
	cout<<"---------------------------------------------------"<<endl;
	// Set mean beam energy (GeV) and slope
	double BE_mean_default=mEnv->GetValue("momentumoffset", 99.0);
	double BE_ugrad_default=mEnv->GetValue("momentumugradient", 99.0);			// energy slope in GeV/mm in u direction
    double BE_vgrad_default=mEnv->GetValue("momentumvgradient", 99.0);			// energy slope in GeV/mm in u direction

	int vertexmultiplicitymin=mEnv->GetValue("vertexmultiplicitymin", 1);			// energy slope in GeV/mm in u direction
	int vertexmultiplicitymax=mEnv->GetValue("vertexmultiplicitymax", 1);			// energy slope in GeV/mm in u direction

	// Set beam particle charge (e)	
	z= mEnv->GetValue("particle.charge", 1);

	// Set beam particle mass (GeV)	
	mass= mEnv->GetValue("particle.mass", 99.0);
	cout<<"Mass: 					"<<mass<<" GeV"<<endl;
	cout<<"Charge: 				"<<z<<" e"<<endl;

	// Create empty grid object that will be filled and used in the fit
	Grid grid(mEnv);

	// Total number of measurement areas
	const size_t num_fitfunctions=grid.GetMeasurementAreas().size();

	cout<<"Total number of measurement areas: 	"<<num_fitfunctions<<endl;

	cout<<"---------------------------------------------------"<<endl;
	cout<<"---------------------------------------------------"<<endl;

	// Definition of the measurement areas can be done in two ways. Either a fixed 3 x 3 alu grid with a regular thickness increase between grid points and fixed positions or 
    // 12 completely independent measurement areas with completely unreleated thicknesses and positions.

	// u minimum and v maximum values (in mm)
	//double umin;
	//double vmin;
	//double umax;
	//double vmax;
	
	// Print out the measurement areas, which will be used for the fit
	grid.PrintGridParameters();

	// TString with input root file name
	TString histoname;
    std::vector<TString> inputfiles;

    GetInputFiles(inputfiles);
	for(size_t iinput=0;iinput<inputfiles.size();iinput++)
	{
		cout<<"input file:"<<inputfiles.at(iinput)<<endl;
	}
	// Open the first input file
	// This is needed to determine the angle resolution
	TFile *Inputfile = new TFile(inputfiles.at(0), "READ");
	recoerr=sqrt(getanglerecovar(Inputfile));
	cout<<"Angle reconstruction error: 	"<<recoerr<<" rad"<<endl;
	Inputfile->Close();
    
	//Open the x0 calibration file
	TString filename="X0Calibration-"+model+"fit";
	TFile *rootfile = new TFile(filename+".root", "RECREATE");
     
	// Create directories containing angle histograms and fits
	rootfile->mkdir("grid");
	rootfile->mkdir("grid/raw");
	rootfile->mkdir("grid/fit");
	rootfile->mkdir("selfconsistency");

	rootfile->cd("");

	// Binning and range of the histograms
	cout<<"histo_range: "<<histo_range<<endl;
	cout<<"numbins: "<<numbins<<endl<<endl;


	rootfile->mkdir("grid/raw/");
	rootfile->mkdir("grid/fit/");

	// Loop over number of measurement areas
	for(size_t i=0; i<num_fitfunctions; i++)
	{
			// Set the histogram name as a string
            histoname.Form("measurementarea%lu",i+1);

			cout<<endl<<endl<<"Measurement area "<<i+1<<endl;

			std::vector <double> cutparameters_position;
			std::vector <int> cutparameters_run;
			std::vector <int> cutparameters_vertex_multiplicity;
			int maxangles;

		    cutparameters_position.push_back(grid.GetMeasurementAreas().at(i).Get_u_min());
		    cutparameters_position.push_back(grid.GetMeasurementAreas().at(i).Get_u_max());
		    cutparameters_position.push_back(grid.GetMeasurementAreas().at(i).Get_v_min());
		    cutparameters_position.push_back(grid.GetMeasurementAreas().at(i).Get_v_max());

			cutparameters_run.push_back(grid.GetMeasurementAreas().at(i).Get_run_min());
			cutparameters_run.push_back(grid.GetMeasurementAreas().at(i).Get_run_max());

			maxangles=grid.GetMeasurementAreas().at(i).Get_max_angles();

			cutparameters_vertex_multiplicity.push_back(vertexmultiplicitymin);
			cutparameters_vertex_multiplicity.push_back(vertexmultiplicitymax);

			// Save the angle histograms of the current measurement area to the root file
			cout<<endl<<"Correct angle distribution offsets..."<<endl;
			TString nametag="uncorrected_";
			savehisto(inputfiles,rootfile, histoname, numbins, histo_range, cutparameters_position, cutparameters_run, cutparameters_vertex_multiplicity, maxangles, 0, nametag);
			cout<<endl<<"Write histos..."<<endl;
			nametag="";
			savehisto(inputfiles,rootfile, histoname, numbins, histo_range, cutparameters_position, cutparameters_run, cutparameters_vertex_multiplicity, maxangles, correctmean, nametag);


	}// end of first loop over measurement areas


	// Open results config file
    //------------------
    TEnv *mEnv_res=new TEnv("x0cal_result.cfg");
	// Check whether there are entries, of this is the case use these entries as starting values

	double lambda_start=mEnv_res->GetValue("lambda_start", lambda_start_default);
	double BE_mean=mEnv_res->GetValue("momentumoffset", BE_mean_default);
	double BE_ugrad=mEnv_res->GetValue("momentumugradient", BE_ugrad_default);	
    double BE_vgrad=mEnv_res->GetValue("momentumvgradient", BE_vgrad_default);


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
	fitoptions.push_back(Use_loglikelihood_estimator);

	cout<<endl<<"Beam particle mean energy start value is				"<<BE_mean<<" GeV"<<endl;
	cout<<"Beam particle energy u gradient start value is		"<<BE_ugrad<<" GeV/mm"<<endl;
	cout<<"Beam particle energy v gradient start value is		"<<BE_vgrad<<" GeV/mm"<<endl;
	cout<<"Beam particle lambda		 start value is				"<<lambda_start<<endl<<endl;

	if(fixlambda) 	cout<<"Calibration factor Lambda will not be calibrated"<<endl;
	else cout<<"Calibration factor Lambda will be calibrated"<<endl;

	if(fix_offset) cout<<"Beam particle mean energy will not be calibrated"<<endl<<endl;
	else cout<<"Beam particle mean energy will be calibrated"<<endl<<endl;

	cout<<"Using "<<model<<" fit model during calibration and self-consistency check!"<<endl;
	if(Use_loglikelihood_estimator) cout<<"Using log-likelihood estimator fit during calibration and self-consistency check!"<<endl<<endl;
	else cout<<"Using chi2 estimator fit during calibration and self-consistency check!"<<endl<<endl;

	double* iresults=fit(rootfile, grid, beamoptions, recoerr, model, fitoptions, rangevalue);

	cout<<endl<<" The lambda calibration factor is: "<<iresults[0]<<" +/- "<<iresults[1]<<endl;
	cout<<" The mean beam energy at (0,0) is: "<<iresults[2]<<" +/- "<<iresults[3]<<"GeV"<<endl;
	cout<<" The kappa calibration factor is: "<<BE_mean/iresults[2]<<" +/- "<<iresults[3]*BE_mean/(iresults[2]*iresults[2])<<endl;
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
    return 0;
  
  }

