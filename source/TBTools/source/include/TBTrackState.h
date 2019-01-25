#ifndef TBTRACKSTATE_H
#define TBTRACKSTATE_H 1

#include<Eigen/Core>
typedef Eigen::Vector5d TrackState;
typedef Eigen::Matrix5d TrackStateCov;

namespace depfet { 


/** Class TBTrackState
 *  
 *  The class TBTrackState represents a local track state on surface. 
 *  Without magnetic field, a particle trajectory can locally be 
 *  represented as a straight line. We use the following 5 variables:
 *  u', v', u, v, q/p
 *  where u'=du/dw and v'=dv/dw are direction tangents wrt. the w=0 
 *  plane in a local UVW coordinate system. The number q/p is the ration 
 *  of charge in units of e and momentum in units of GeV. 
 *  The reference frame is uniquely determined by a reference plane 
 *  number.
 *  
 *  @Author B. Schwenker, University of Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */


class TBTrackState {
   
 public: // members 
  
  // Parameter vector
  TrackState Pars;
  // Parameter covariance (5x5 matrix)
  TrackStateCov Cov;
  // Reference plane number
  int Plane; 
    
  
 public: // functions 
   
  // Constructors
  TBTrackState(); 
  TBTrackState(TrackState aPars, TrackStateCov aCov=TrackStateCov::Zero(), int aPlane=0);
  
  // Dimension of track state 
  int GetDim() { return 5;}
  
  // Get/Set track parameters  
  void SetPars(const TrackState& aPars) { Pars= aPars; }
  TrackState&  GetPars() { return Pars; }
   
  // Get/Set parameter covariance 
  void SetCov(const TrackStateCov& aCov ) { Cov = aCov; }
  TrackStateCov&  GetCov() { return Cov; }

  // Get/Set plane number 
  void SetPlane(int aPlane ) { Plane = aPlane; }
  int  GetPlane() { return Plane; }
  
  // Get intersection coordinates 
  Eigen::Vector2d GetXCoord();
  
};

} // Namespace

#endif 

