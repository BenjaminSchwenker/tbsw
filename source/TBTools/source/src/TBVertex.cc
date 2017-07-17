

// DEPFETTrackTools includes
#include "TBVertex.h"

using namespace std;
using namespace CLHEP;

namespace depfet {

/**
 * File TBVertex.cc
 * 
 * Definition of class TBVertex.
 * 
 *  @Author E. Heinrich, University of Göttingen
 *  <mailto:erik.heinrich@stud.uni-goettingen.de>
 */

// Constructors

TBVertex::TBVertex() : Pos(3, 1, 0), Cov(3, 3, 0)  , chi2(0)
{;}

TBVertex::TBVertex(HepMatrix aPos, HepMatrix aCov, double achi2) : Pos(3, 1, 0), Cov(3,3,0)  , chi2(0)
{ 
  // Set vertex position 
  Pos = aPos;
  // Set vertex covariance
  Cov = aCov;
  // Set fit chi²
  chi2 = achi2;
}

} // Namespace
