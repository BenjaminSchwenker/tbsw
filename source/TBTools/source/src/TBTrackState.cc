

// TBTools includes
#include "TBTrackState.h"

using namespace std;

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

TBTrackState::TBTrackState() :  Plane(0)
{
    Pars=TrackState::Zero();
    Cov=TrackStateCov::Zero();

}

TBTrackState::TBTrackState(TrackState aPars, TrackStateCov aCov, int aPlane) : Plane(aPlane)
{ 
  // Set the plane number
  Plane = aPlane
  // Set track state vector  
  Pars = aPars;
  // Set track state covariance matrix
  Cov = aCov;
}

// Get intersection coordinates 
Eigen::Vector2d TBTrackState::GetXCoord()
{ 
  return Pars.block<2,1>(2,0);

}



} // Namespace
