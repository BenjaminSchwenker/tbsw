// TBTools includes
#include "TBTrackState.h"


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

TBTrackState::TBTrackState() 
{
    Plane = 0;
    Pars=TrackState::Zero();
    Cov=TrackStateCovariance::Zero();
}

TBTrackState::TBTrackState(const TrackState& aPars, const TrackStateCovariance& aCov, int aPlane) 
{ 
  // Set the plane number
  Plane = aPlane;
  // Set track state vector  
  Pars = aPars;
  // Set track state covariance matrix
  Cov = aCov;
}

// Get intersection coordinates 
Eigen::Vector2d TBTrackState::GetXCoord() const
{ 
  return Pars.block<2,1>(2,0);
}

Eigen::Matrix2d TBTrackState::GetXCoordCovariance() const
{
  return Cov.block<2,2>(2,2);
}

} // Namespace
