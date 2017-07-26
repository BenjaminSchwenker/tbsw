// TBVertexFitter implementation file 
// 
//  @Author E. Heinrich, University of Göttingen
//  <mailto:erik.heinrich@stud.uni-goettingen.de>

// Marlin includes 
#include <marlin/Global.h>
#include <streamlog/streamlog.h>
  
// DEPFETTrackTool includes
#include "TBVertexFitter.h"

// C++ STL includes 
#include <cmath>
#include <limits>
#include <vector>

// Namespaces
using namespace CLHEP;
using namespace marlin;
using namespace std; 

namespace depfet {

/** Constructor 
 */
TBVertexFitter::TBVertexFitter(int ipl)
{                   
  
  ndim = 4; // dimension of measurement vector p = (du,dv,u,v)

  plnr = ipl; //plane number of DUT

}

TBVertex TBVertexFitter::FitVertex(TBTrackState& InState, TBTrackState& OutState)
{
  //Take relevant matrices from trackstates
  HepMatrix In_pars = InState.GetPars();
  HepMatrix Out_pars = OutState.GetPars();
  HepMatrix In_covs = InState.GetCov();
  HepMatrix Out_covs = OutState.GetCov();

  //Arrays containing the relevant 4D submatrices as 5th Dimension q/p is irrelevant for vertex reconstruction 
  std::vector<HepMatrix> p = {In_pars.sub(1,ndim,1,1), Out_pars.sub(1,ndim,1,1)};
  std::vector<HepMatrix> V = {In_covs.sub(1,ndim,1,ndim), Out_covs.sub(1,ndim,1,ndim)};

  //Initial vertex guess at (0,0,0) and large covariance
  HepMatrix r(3,1,0);
  HepMatrix C(3,3,0);
	C[0][0] = 100*100;
	C[1][1] = 100*100;
	C[2][2] = 100*100;
  //Initial chi2
  double chi2 = 0;

  //Loop over trackstates
  for (int i=0; i<2; i++) {

    //Jacobian B = dh/dq for slope states q = (a,b)
    HepMatrix B(ndim,2,0);
	B[0][0] = 1;	//dh1/da
	B[1][1] = 1;	//dh2/db

    //Calculate filter matrices
    HepMatrix G = V[i].inverse();
    HepMatrix W = (B.T()*G*B).inverse();
    HepMatrix GB = G - G * B * W * B.T() * G;

    //Jacobian A = dh/dr for vertex state r = (x,y,z)
    HepMatrix A(ndim,3,0);
	A[2][0] = 1;	//dh3/dx
	A[3][0] = 1;	//dh4/dy
	A[2][1] = -p[i][0][0];	//dh3/dz
	A[3][2] = -p[i][1][0]; 	//dh4/dz

    //Update vertex, covariance
    HepMatrix Ctemp = C;
    C = (C.inverse() + A.T()*GB*A).inverse();

    HepMatrix rtemp = r;
    r = C * (Ctemp.inverse()*r + A.T()*GB*p[i]);

    //calculate slope state
    HepMatrix q = W * B.T() * G * (p[i] - A * r);

    //compute residual and chi2 increment
    HepMatrix zeta = p[i] - A * r - B * q;
    chi2 += (zeta.T() * G * zeta + (r - rtemp).T() * Ctemp.inverse() * (r - rtemp))[0][0];

  }

  TBVertex vertex(r,C,chi2);
  
  return vertex;
}

} // Namespace;

