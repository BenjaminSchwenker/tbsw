// C/C++ includes
#include <iostream>
 	
// ROOT includes
#include <TRandom.h>
#include <TMath.h>
 	
// local includes
#include "ThreeDModel.h"
 	
using namespace std;
using namespace CLHEP;

namespace depfet {
 	
/** 
 *
 *  Definition of functions and procedures for 3D manipulation with matrices.
 *
 *  Rotations are always defined to transform global coordinates to local
 *  coordinates. (GEANT convention).
 *
 */
	
/** Fill a rotation matrix with the Karimaki angles alpha, beta, gamma. */
void FillRotMatrixKarimaki(HepMatrix & m, double alpha, double beta, double gamma)
{
  // Rotation matrix definition from Veikko Karimaki
  HepMatrix rot_alpha(3, 3, 1);
  rot_alpha[1][1] = cos(alpha);
  rot_alpha[1][2] = -sin(alpha);
  rot_alpha[2][2] = rot_alpha[1][1];
  rot_alpha[2][1] = -rot_alpha[1][2];
  HepMatrix rot_beta(3, 3, 1);
  rot_beta[0][0] = cos(beta);
  rot_beta[0][2] = sin(beta);
  rot_beta[2][2] = rot_beta[0][0];
  rot_beta[2][0] = -rot_beta[0][2];
  HepMatrix rot_gamma(3, 3, 1);
  rot_gamma[0][0] = cos(gamma);
  rot_gamma[0][1] = -sin(gamma);
  rot_gamma[1][1] = rot_gamma[0][0];
  rot_gamma[1][0] = -rot_gamma[0][1];
  //   cout << "Rotation Matrix A " << rot_alpha << endl;
  //   cout << "Rotation Matrix B " << rot_beta << endl;
  //   cout << "Rotation Matric C " << rot_gamma << endl;
  m = rot_gamma * rot_beta * rot_alpha;
  //   cout << "Combined Rotation " << m << endl;
}
 	
/** Extract from a rotation matrix the Karimaki angles alpha, beta, gamma. */
void GetAnglesKarimaki(const HepMatrix & m, double & alpha, double & beta, double & gamma)
{
  // Compute from the given matrix m the Karimaki angles
  // this only works if for all angles abs(angle) < Pi/2
  //   cout << m << endl << m[2][0] << endl << m[2][1] << endl << m[1][0] << endl;
  beta = -asin(m[2][0]);
  double cbeta = cos(beta);
  alpha = asin(m[2][1]/cbeta);
  gamma = asin(m[1][0]/cbeta);
}
 	
/** Test the Karimaki rotation matrices. */
void testAnglesKarimaki()
{
  double Pi = TMath::Pi();
  double alpha = gRandom->Uniform(Pi)-Pi/2.;
  double beta  = gRandom->Uniform(Pi)-Pi/2.;
  double gamma = gRandom->Uniform(Pi)-Pi/2.; 
  //   cout << "alpha = " << alpha << endl
  //        << "beta  = " << beta  << endl
  //        << "gamma = " << gamma << endl;
  HepMatrix m;
  FillRotMatrixKarimaki(m, alpha, beta, gamma);
  double aa, bb, cc;
  GetAnglesKarimaki(m, aa, bb, cc);
  //   cout << "after FillRotMatrixKarimaki and GetAnglesKarimaki" << endl;
  //   cout << "alpha = " << aa << endl
  //        << "beta  = " << b << endl
  //        << "gamma = " << c << endl;
  if (TMath::Abs(aa-alpha) > 1E-7) {
    cout << "testAnglesKarimaki failed with angle alpha" << endl;
  }
  if (TMath::Abs(bb-beta) > 1E-7) {
    cout << "testAnglesKarimaki failed with angle beta" << endl;
  }
  if (TMath::Abs(cc-gamma) > 1E-7) {
    cout << "testAnglesKarimaki failed with angle gamma" << endl;
  }
}

} // Namespace;

