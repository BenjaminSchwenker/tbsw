

// DEPFETTrackTools includes
#include "TBTrackState.h"

using namespace std;
using namespace CLHEP;

namespace depfet {

/**
 * File TBTrackState.cc
 * 
 * Definition of class TBTrackState.
 * 
 *  @Author B. Schwenker, University of Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

// Constructors

TBTrackState::TBTrackState() : Pars(4, 1, 0), Cov(4, 0)  , Plane(0)
{;}

TBTrackState::TBTrackState(HepMatrix aPars, HepSymMatrix aCov) : Pars(4, 1, 0), Cov(4,0)  , Plane(0)
{ 
  // Set track parameters 
  Pars = aPars;
  // Set track covariance
  Cov = aCov;
}

// Get intersection coordinates 
HepMatrix TBTrackState::GetXCoord()  
{ 
  return Pars.sub(3,4,1,1); 
}



} // Namespace
