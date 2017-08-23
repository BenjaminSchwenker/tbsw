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
TBVertexFitter::TBVertexFitter(int ipl) : A(5,3,0), B(5,2,0)
{                   

  plnr = ipl; //plane number of DUT

  //Define Matrices used for vertex fit
  //Jacobian B = dh/dq for slope states q = (a,b)
	B[0][0] = 1;	//dh1/da
	B[1][1] = 1;	//dh2/db

  //Jacobian A = dh/dr for vertex state r = (x,y,z)
	A[2][0] = 1;	//dh3/dx
	A[3][1] = 1;	//dh4/dy

}

bool TBVertexFitter::FitVertex(TBVertex& Vertex)
{
  int ierr = 0; 

  //Initial vertex guess at (0,0,0) and large covariance
  HepMatrix r(3,1,0);
  HepMatrix C(3,3,0);
	C[0][0] = 100*100;
	C[1][1] = 100*100;
	C[2][2] = 100*100;
  //Initial chi2 and degrees of freedom
  double chi2 = 0;
  int ndf = -3;
  //Initialisation of residual vector
  HepMatrix zeta(4,1,0);

  //Loop over trackstates
  for (int i=0; i < Vertex.GetStates().size(); i++) {

    //Copy current trackstate
    HepMatrix p = Vertex.GetStates()[i].GetPars();
    HepMatrix V = Vertex.GetStates()[i].GetCov();

    //Calculate filter matrices
    HepMatrix G = V.inverse(ierr);
    if (ierr != 0) {
      streamlog_out(ERROR) << "ERR: Matrix inversion failed. Quit fitting!"
                           << std::endl;
      return true;
    }	
    HepMatrix W = (B.T()*G*B).inverse();
    HepMatrix GB = G - G * B * W * B.T() * G;

    //Modify Jacobian A to represent current trackstate
    A[2][2] = -p[0][0];
    A[3][2] = -p[1][0];

    //Update vertex, covariance
    HepMatrix Ctemp = C;
    C = (C.inverse() + A.T()*GB*A).inverse();

    HepMatrix rtemp = r;
    r = C * (Ctemp.inverse()*r + A.T()*GB*p);

    //calculate slope state
    HepMatrix q = W * B.T() * G * (p - A * r);

    //compute residual and chi2 increment
    zeta = p - A * r - B * q;
    chi2 += (zeta.T() * G * zeta + (r - rtemp).T() * Ctemp.inverse() * (r - rtemp))[0][0];

    //increase Number degrees of freedom
    ndf += 2;

  }

  Vertex.SetRes(zeta);
  Vertex.SetPos(r);
  Vertex.SetCov(C);
  Vertex.SetChi2(chi2);
  Vertex.SetNdf(ndf);
  
  return false;
}

} // Namespace;

