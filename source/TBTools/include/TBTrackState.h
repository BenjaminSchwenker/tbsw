#ifndef TBTRACKSTATE_H
#define TBTRACKSTATE_H 1

#include<Eigen/Dense>
typedef Eigen::Matrix<double,5,1> TrackState;
typedef Eigen::Matrix<double,5,5> TrackStateCovariance;



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
  TrackStateCovariance Cov;
  // Reference plane number
  int Plane; 
    
  
 public: // functions 
   
  // Constructors
  TBTrackState(); 
  TBTrackState(const TrackState& aPars, const TrackStateCovariance& aCov=TrackStateCovariance::Zero(), int aPlane=0);
  
  // Dimension of track state 
  int GetDim() const { return 5;}
  
  // Get/Set track parameters  
  void SetPars(const TrackState& aPars) { Pars= aPars; }
  TrackState&  GetPars() { return Pars; }
  const TrackState&  GetPars() const { return Pars; }
   
  // Get/Set parameter covariance 
  void SetCov(const TrackStateCovariance& aCov ) { Cov = aCov; }
  const TrackStateCovariance&  GetCov() const { return Cov; }
  TrackStateCovariance&  GetCov() { return Cov; }

  // Get/Set plane number 
  void SetPlane(int aPlane ) { Plane = aPlane; }
  int  GetPlane() const { return Plane; }
  
  // Get intersection coordinates 
  Eigen::Vector2d GetXCoord() const;

  // Get covariance for intersection coordinates 
  Eigen::Matrix2d GetXCoordCovariance() const;
  
};

} // Namespace

#endif 

