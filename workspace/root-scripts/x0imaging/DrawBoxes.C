#include <fstream>
using namespace std ;

// Hole Class describing a hole in the calibration grid. Parameters of the hole are position and side length
// as well as the radiation length in the hole region
class Hole
{
	private:

	// Parameter declarations

	double center_u,center_v;  	// Position of the center of the hole
	double length_u,length_v;  	// Side length of the hole
	double thickness;          	// Radiation length of the material in the region of the hole
	double BE;			// Beam energy of particle beam in the region of the hole

	public:

	// Constructors

	Hole(double, double, double, double, double, double);	// Constructor
	Hole();						// Default Constructor

	// Reading hole parameters

	double Get_u_min() { return center_u-0.5*length_u; }	// Return min u value of hole
	double Get_u_max() { return center_u+0.5*length_u; }	// Return max u value of hole
	double Get_v_min() { return center_v-0.5*length_v; }	// Return min v value of hole
	double Get_v_max() { return center_v+0.5*length_v; }	// Return max v value of hole
	double Get_u_center() { return center_u; }		// Return center u value of hole
	double Get_v_center() { return center_v; }		// Return center v value of hole
	double Get_u_length() { return length_u; }		// Return the side length (in u direction) of hole
	double Get_v_length() { return length_v; }		// Return the side length (in v direction) of hole
	double Get_thickness()	{ return thickness; }		// Return thickness
	double Get_BE()		{ return BE; }			// Return Beam energy


	// Setting the hole parameters

	void SetHoleCenter(double a, double b)			// Reset the center Position of the hole
	{
		center_u=a;
		center_v=b;
	}

	void SetHoleLength(double a, double b)			// Reset the side length of the hole
	{
		length_u=a;
		length_v=b;
	}

	void SetHolethickness(double a) { thickness=a; } 	// Reset the thickness of the hole

	void SetBeamEnergy(double a) { BE=a; } 		// Reset the Beamenergy of the particles in the hole area

	void SetHoleParameters(double ucenter, double vcenter, double ulength, double vlength, double thick, double beamenergy)	// Reset all parameters of the hole
	{
		center_u=ucenter;
		center_v=vcenter;
		length_u=ulength;
		length_v=vlength;
		thickness=thick;
		BE=beamenergy;
	}

	void PrintHoleParameters()
	{
		cout<<"u center: "<<center_u<<endl;
		cout<<"v center: "<<center_v<<endl;

		cout<<"u length: "<<length_u<<endl;
		cout<<"v length: "<<length_v<<endl;

		cout<<"thickness: "<<thickness<<endl;
		cout<<"beam energy: "<<BE<<endl;
	}
};	

// Constructor definition
Hole::Hole(double ucenter, double vcenter, double ulength, double vlength, double thick, double beamenergy)
{
	center_u=ucenter;
	center_v=vcenter;
	length_u=ulength;	// mm
	length_v=vlength;	// mm
	thickness=thick;	// mm
	BE=beamenergy;		// GeV
}

// Default Constructor
Hole::Hole()
{
	center_u=0.0;
	center_v=0.0;
	length_u=1.0;	// mm
	length_v=1.0;	// mm
	thickness=0.0;	// mm
	BE=4.0;			// GeV
}

// Class describing the complete calibration grid. The grid consists of a set of holes with a specific radiation length.
class Grid
{
	private:

	// Declaration of parameters

	double center_u,center_v;  					// Position of the center of the central hole
	double lengthfactor;
	double thickness_inc;        				// Radiation length increment from hole to hole
	double density;								// Density of the target material
	double Z;									// Atomic number of the target material
	double A; 									// Atomic mass of the target material

	double mirrored_u;
	double mirrored_v;

	std::vector<Hole> holes;					// Declaration the holes in the grid

	public:

	// Constructors

	Grid(double, double, double);
	Grid(double, double, double, double, double, double, double, double, bool, bool, bool);			

	// Read out Grid Parameters

	// Return min u value of i-th hole
	double Get_u_min(int i) { return holes[i].Get_u_center()-0.5*holes[i].Get_u_length(); }
	
	// Return max u value of i-th hole
	double Get_u_max(int i)	{ return holes[i].Get_u_center()+0.5*holes[i].Get_u_length(); }	

	// Return min v value of i-th hole
	double Get_v_min(int i)	{ return holes[i].Get_v_center()-0.5*holes[i].Get_v_length(); }

	// Return max v value of i-th hole
	double Get_v_max(int i) { return holes[i].Get_v_center()+0.5*holes[i].Get_v_length(); }

	// Return center u value of i-th hole
	double Get_u_center(int i) { return holes[i].Get_u_center(); }

	// Return center v value of i-th hole
	double Get_v_center(int i) { return holes[i].Get_v_center(); }

	// Return thickness of i-th hole
	double Get_thickness(int i) { return holes[i].Get_thickness(); }

	// Return density of i-th hole
	double Get_density() { return density; }

	// Return atomic number of i-th hole
	double Get_Z() { return Z; }

	// Return atomic mass of i-th hole
	double Get_A() { return A; }

	// Return Beam Energy in area of the i-th hole
	double Get_BE(int i) { return holes[i].Get_BE(); }

	// Return number of measurement areas/holes
	int GetNumHoles() { return holes.size(); }

	// Return holes vector
	std::vector<Hole> Get_holes() { return holes; }

	void PrintGridParameters()
	{
		for(int i=0;i<holes.size();i++) holes.at(i).PrintHoleParameters();
		cout<<"Density: "<<density<<endl;
		cout<<"Z: "<<Z<<endl;
		cout<<"A: "<<A<<endl;
	}

	// Set grid parameters

	// Change hole parameters
	void SetHoleParameters(int ihole, double u_center, double v_center, double u_length, double v_length, double thick, double dens, double atom_num, double atom_mass, double energy)
	{
		holes.at(ihole).SetHoleParameters(u_center, v_center, u_length, v_length, thick, energy);
	}

	// Change hole beam energy
	void SetHoleBE(int ihole, double energy)
	{
		holes.at(ihole).SetBeamEnergy(energy);
	}

	// Add a hole
	void AddHole(Hole hole)
	{
		holes.push_back(hole);
	}

	// Set grid to be equal to another grid
	void SetGrid(Grid grid)
	{
		holes=grid.Get_holes();
		density=grid.Get_density();
		Z=grid.Get_Z();
		A=grid.Get_A();
	}
};	


Grid::Grid(double ucenter, double vcenter, double length, double thick_increase, double dens, double atomicnumber, double atomicmass, double BE,  bool i, bool j, bool uv)
{
	center_u=ucenter;
	center_v=vcenter;
	double lengthfactor=length;

	thickness_inc=thick_increase;
	density=dens;
	Z=atomicnumber;
	A=atomicmass;	

	int help=0;

	// Must the u axis be mirrored?
	if(i==true) mirrored_u=-1.0;
	else mirrored_u=1.0;

	// Must the v axis be mirrored?
	if(j==true) mirrored_v=-1.0;
	else mirrored_v=1.0;

	// Set the parameters of all holes, which are part of the 3x3 grid
	for(int m=0; m<9; m++)
	{
		if(m==0)
		{
            // Case: u=first coordinate and v=second coordinate
		    if(uv==false) 
		    {
				Hole hole(center_u+mirrored_u*2.8, center_v+mirrored_v*2.1, 1.0*lengthfactor, 1.0*lengthfactor, 0.0, BE);
				holes.push_back(hole);
			}
            // Case: v=first coordinate and u=second coordinate
            else 
			{
				Hole hole(center_v+mirrored_v*2.8, center_u+mirrored_u*2.1, 1.0*lengthfactor, 1.0*lengthfactor, 0.0, BE);
				holes.push_back(hole);
			}
		}


		else if(m==8)
		{
           	// Case: u=first coordinate and v=second coordinate
		   	if(uv==false) 
			{
				Hole hole(center_u-mirrored_u*2.0, center_v-mirrored_v*2.5, 1.0*lengthfactor, 2.0*lengthfactor, m*1.0*thickness_inc, BE);
				holes.push_back(hole);
			}
           	// Case: v=first coordinate and u=second coordinate
		   	else 
			{
				Hole hole(center_v-mirrored_v*2.5, center_u-mirrored_u*2.0, 2.0*lengthfactor, 1.0*lengthfactor, m*1.0*thickness_inc, BE);
				holes.push_back(hole);
			}
		}

		else
		{
		   if((m==3)||(m==6)) help++;

		   	// Case: u=first coordinate and v=second coordinate
		   	if(uv==false) 
			{

				Hole hole(center_u+mirrored_u*(1-help)*2.0, center_v+mirrored_v*(1+3*help-m)*2.0, 1.0*lengthfactor, 1.0*lengthfactor, m*1.0*thickness_inc, BE);
				holes.push_back(hole);
			}
          	// Case: v=first coordinate and u=second coordinate
           	else 
			{
				Hole hole(center_v+mirrored_v*(1+3*help-m)*2.0, center_u+mirrored_u*(1-help)*2.0, 1.0*lengthfactor, 1.0*lengthfactor, m*1.0*thickness_inc, BE);
				holes.push_back(hole);
			}
		}
	}

	// Leave the default values for the last three holes
	
}

Grid::Grid(double dens, double atomicnumber, double atomicmass)
{
	density=dens;
	Z=atomicnumber;
	A=atomicmass;
}





// This script is used to create a map of a plane in a test beam telescope. The input is a TTree including 
// MSC projected scattering angle distributions and reconstruction errors.
void DrawBoxes()
{
	gSystem->Load("libProof.so");
	gSystem->Load("libTreePlayer.so");
	gSystem->Load("libRooFit") ;

	gROOT->Reset(); 
	gROOT->SetStyle("Plain"); 

	// display mode
	gStyle->SetPalette(1);
	gStyle->SetOptStat(11111111);

    // Read config file
    //------------------
    TEnv mEnv("x0calibration.cfg");


	//number of fit functions required for the number of measurement areas given by the cfg file

	// Determine the number of fit functions from the cfg file
	int num_grid=0; // Number of measurement grids (corresponds to 9 holes)
	int num_hole=0; // Number of additional holes

	TString gridname;
	gridname.Form("grid");
	TString holename;
	holename.Form("hole%i",num_hole+1);
	
	if(mEnv.GetValue(gridname+".exist",0)!=0) num_grid=1;

	cout<<"Number of grids: "<<num_grid<<endl;

	while(mEnv.GetValue(holename+".exist",0)!=0)
	{
		num_hole++;
		holename.Form("hole%i",num_hole+1);
	}

	cout<<"Number of additional holes: "<<num_hole<<endl;
	
	const int num_fitfunctions=num_grid*9+num_hole;

	cout<<"Total number of measurement areas: "<<num_fitfunctions<<endl;

	// Choose the multiple scattering model
	// "moliere": Moliere model
	// "highland": Highland model
	TString model=mEnv.GetValue("model", "highland");

	// Integer determining, whether a offset correction is applied to the angular distributions
	// 1: correction
	// everything else, no correction
	int correctmean=mEnv.GetValue("correctmean", 1 );

	// Define and set the parameters used in the Moliere fit
	// BE: beam energy (GeV), z: charge of beam particle (e), mass: mass of beam particle (GeV),
	// density: density of the target material (g/cm³), Z: atomic number of target material, A: atomic mass (amu)
	// d_layer: thickness per layer(mm), num_layers: total number of layers, recoerr: angle reconstruction error (rad)
	double z,Z,A,p,mass,density,recoerr;
	
	cout<<"---------------------------------------------------"<<endl;
	cout<<"-----------------Parameter settings----------------"<<endl;
	cout<<"---------------------------------------------------"<<endl;
	// Set mean beam energy (GeV) and slope
	p=mEnv.GetValue("particle.momentummean", 4.0);
	cout<<"Beam particle momentum: 	"<<p<<" GeV/c"<<endl;

	double BE=TMath::Sqrt(p*p+mass*mass);   // energy in GeV
	double BE_Slope=mEnv.GetValue("particle.momentumslope", 0.0);			// energy slope in GeV/mm
	cout<<"Mean Beam energy: 		"<<BE<<" GeV"<<endl;

	// Set beam particle charge (e)	
	z= mEnv.GetValue("particle.charge", 1);

	// Set beam particle mass (GeV)	
	mass= mEnv.GetValue("particle.mass", 4.0);
	cout<<"Mass: 				"<<mass<<" GeV"<<endl;
	cout<<"Charge: 		"<<z<<" e"<<endl;

	// Set target material density (g/cm³)	
	density= mEnv.GetValue("material.density", 2.7);

	// Set atomic number of target material	
	Z=mEnv.GetValue("material.atomicnumber", 13);

	// Set atomic mass of target material (amu)
	A=mEnv.GetValue("material.atomicmass", 27);

	cout<<"Density: 			"<<density<<" g/cm^3"<<endl;
	cout<<"A: 				"<<A<<" "<<endl;
	cout<<"Z: 				"<<Z<<" "<<endl;


	// Print number of measurement areas on target plane
	cout<<"Number of measurement areas: 		"<<num_fitfunctions<<" "<<endl;

	cout<<"---------------------------------------------------"<<endl;
	cout<<"---------------------------------------------------"<<endl;

	// Definition of the measurement areas can be done in two ways. Either a fixed 3 x 3 alu grid with a regular thickness increase between grid points and fixed positions or 
    // 12 completely independent measurement areas with completely unreleated thicknesses and positions.

    double u_hole_center[num_fitfunctions];
	double v_hole_center[num_fitfunctions];

	// u minimum and v maximum values (in mm)
	double umin[num_fitfunctions];
	double vmin[num_fitfunctions];
	double umax[num_fitfunctions];
	double vmax[num_fitfunctions];

	double sidelength;

	// Definition of grid class object
	// Center of the 3x3 grid
	double center_u=mEnv.GetValue("grid.ucenter", 0.0);
	double center_v=mEnv.GetValue("grid.vcenter", 0.0);

    sidelength=mEnv.GetValue("grid.hole_side_length", 0.8);	

	// Thickness increment for the 9 holes of the 3x3 grid
	double thickness_inc=mEnv.GetValue("grid.thicknessincrease", 0.3);

	// These parameters determine the orientation and mirroring of the 3x3 grid
	bool mirroring_u=mEnv.GetValue("grid.mirror_u", 0);
	bool mirroring_v=mEnv.GetValue("grid.mirror_v", 0);
	bool switch_uv=mEnv.GetValue("grid.switch_uv", 0);

	// Create empty grid object that will be filled and used in the fit
	Grid grid(density,Z,A);

	// The Grid class creates the 3x3 grid, the three additional holes outside of the grid have to be added manually
	if(num_grid==1)
	{
		Grid grid1(center_u, center_v, sidelength, thickness_inc, density, Z, A, BE, mirroring_u, mirroring_v, switch_uv);
		grid.SetGrid(grid1);
	}

	for(int i=0; i<num_hole;i++)
	{
		TString holename;
		holename.Form("hole%i",i+1);

		double ucenter=mEnv.GetValue(holename+".ucenter", 1.0);
		double vcenter=mEnv.GetValue(holename+".vcenter", 1.0);
		double ulength=mEnv.GetValue(holename+".hole_u_length", 1.0);
		double vlength=mEnv.GetValue(holename+".hole_v_length", 1.0);
		double thickness=mEnv.GetValue(holename+".thickness", 1.8);

		Hole hole(ucenter,vcenter,ulength,vlength,thickness,BE);

		//cout<<"Hole info: "<<endl;
		//hole.PrintHoleParameters();
  
		grid.AddHole(hole);
		
	}

	// Readout max and min values of all measurement areas of the calibration grid
	for(int i=0;i<num_fitfunctions;i++)
	{
		umin[i]=grid.Get_u_min(i);
		umax[i]=grid.Get_u_max(i);
		vmin[i]=grid.Get_v_min(i);
		vmax[i]=grid.Get_v_max(i);
	}


    // TString for the input root file name
    TString filename,histoname,range;
	filename="X0image";

	// Copy the X0 Analysis Root file 
	TFile *rootFile = new TFile(filename+".root", "READ");

	// Create mapping histopgrams
	rootFile->cd("");
    gStyle->SetPalette(1,0);

	// X0 map
	hX0map = (TH2*)rootFile->Get("mapping/result/hX0map");
	rootFile->cd("mapping/result");

	// Create a new canvas
	TCanvas * c = new TCanvas("c", "c", 1400, 1000);
    hX0map->SetMaximum(2.3);
	hX0map->SetStats(kFALSE);
    hX0map->Draw("colz");

    hX0map->GetZaxis()->SetTitle("X/X_{0}[%]");
    hX0map->GetZaxis()->SetLabelSize(0.07);
    hX0map->GetZaxis()->SetTitleSize(0.07);
    hX0map->GetZaxis()->SetTitleOffset(0.10);

    //X0map->GetXaxis()->SetTitleSize(3);
    hX0map->GetXaxis()->SetTitle("u[mm]");
    hX0map->GetXaxis()->SetLabelSize(0.07);
    hX0map->GetXaxis()->SetTitleSize(0.07);
    hX0map->GetXaxis()->SetTitleOffset(0.6);

	hX0map->GetYaxis()->SetTitle("v[mm]");
    hX0map->GetYaxis()->SetLabelSize(0.07);
    hX0map->GetYaxis()->SetTitleSize(0.07);
    hX0map->GetYaxis()->SetTitleOffset(0.5);

    hX0map->SetTitle("");

	for(int i=0;i<num_fitfunctions;i++)
	{
	 	TBox* box=new TBox(umin[i], vmin[i], umax[i], vmax[i]); 
		TString number;
		number.Form("%i",i+1);
	 	TText* text=new TText(0.5*(umax[i]+umin[i])-0.25*(umax[i]-umin[i]),0.5*(vmax[i]+vmin[i])-0.25*(vmax[i]-vmin[i]),number);
		if(i>2) text->SetTextColor(1);
		else text->SetTextColor(0);


		if(i<9)
		{
			text->DrawText(0.5*(umax[i]+umin[i])-0.25*(umax[i]-umin[i]),0.5*(vmax[i]+vmin[i])-0.25*(vmax[i]-vmin[i]),number);
		}

		else
		{
			text->DrawText(0.5*(umax[i]+umin[i])-0.42*(umax[i]-umin[i]),0.5*(vmax[i]+vmin[i])-0.25*(vmax[i]-vmin[i]),number);
		}

		box->SetFillStyle(0);
		box->SetLineColor(1);
		box->SetLineWidth(2);
		box->Draw();
	}

    c->Modified();
    c->Update();

    TPaletteAxis *palette = (TPaletteAxis*)hX0map->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.9025);
    palette->SetX2NDC(0.925);

	c->SaveAs(filename+"_Boxes.pdf");
}


