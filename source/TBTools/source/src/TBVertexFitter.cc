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
TBVertexFitter::TBVertexFitter(int ipl, TBDetector det) : A(5,3,0), B(5,3,0)
{                   

  plnr = ipl; //plane number of DUT
  detector = det; // TBDetector object

  //Define Matrices used for vertex fit
  //Jacobian B = dh/dq for slope states q = (a,b,(q/p))
	B[0][0] = 1;	//dh1/da
	B[1][1] = 1;	//dh2/db
	B[4][2] = 1;	//dh5/d(q/p)

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
  C[0][0] = 10000*10000;
  C[1][1] = 10000*10000;
  C[2][2] = 10000*10000;

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

  // Transform vertex position and vertex position covariance into global coordinates

  // First get the reference frame of the DUT
  ReferenceFrame localframe = detector.GetDet(plnr).GetNominal();

  // Get rotation matrix and position offset of local coordinate system of DUT plane
  HepMatrix Rotation = localframe.GetRotRef().T();
  HepVector Offset = localframe.GetPosRef();

  /*cout<<"Rotation Matrix"<<endl;
  cout<<Rotation<<endl;

  cout<<"Offset"<<endl;
  cout<<Offset<<endl;

  cout<<"position"<<endl;
  cout<<r<<endl;*/

  // Now transform the local vertex position and its covariance
  HepVector r_global=Rotation*r+Offset;
  HepMatrix C_global=Rotation.T()*C*Rotation;

  /*cout<<"Vertex position in local coords:"<<endl;
  cout<<r<<endl;

  cout<<"Vertex position in global coords"<<endl;
  cout<<r_global<<endl;

  cout<<"Vertex position cov in local coords:"<<endl;
  cout<<C<<endl;

  cout<<"Vertex position cov in global coords"<<endl;
  cout<<C_global<<endl;*/

  Vertex.SetRes(zeta);
  Vertex.SetPos(r);
  Vertex.SetGlobalPos(r_global);
  Vertex.SetCov(C);
  Vertex.SetGlobalCov(C_global);
  Vertex.SetChi2(chi2);
  Vertex.SetNdf(ndf);
  
  return false;
}

} // Namespace;

