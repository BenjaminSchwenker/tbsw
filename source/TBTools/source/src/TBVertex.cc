

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

TBVertex::TBVertex() : Pos(3, 1, 0), GlobalPos(3, 1, 0), Cov(3,3,0), GlobalCov(3,3,0), chi2(0)
{;}

TBVertex::TBVertex(HepMatrix aPos, HepMatrix aGlobalPos, HepMatrix aCov, HepMatrix aGlobalCov, double achi2) : Pos(3, 1, 0), GlobalPos(3, 1, 0), Cov(3,3,0), GlobalCov(3,3,0), chi2(0)
{ 
  // Set vertex position 
  Pos = aPos;
  // Set global vertex position 
  GlobalPos = aGlobalPos;
  // Set vertex covariance
  Cov = aCov;
  // Set global vertex covariance
  GlobalCov = aGlobalCov;
  // Set fit chi²
  chi2 = achi2;
}

} // Namespace
