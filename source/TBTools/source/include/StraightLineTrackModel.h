#ifndef StraightLineTrackModel_H
#define StraightLineTrackModel_H 1

// DEPFETTrackTools includes
#include "GenericTrackModel.h"

namespace depfet {

 
//! Class StraightLineTrackModel 
/* 
 *  The class StraightLineTrackModel implements track dynamics inside a volume 
 *  with zero magnetic field. 
 * 
 *  Relative to the W=0 plane of a local UVW coordinate frame, the local 
 *  track state is defined by 4 variables:
 *  x = (u', v', u, v)
 *  where u'=du/dw and v'=dv/dw in the local detector UVW coordinate system.
 *  The local surface is spanned by U, V unit vectors.
 *
 *  The routines below cover tracking without a magnetic field. 
 *
 *  A) 3D trajectory model = straight line
 *  
 *  q(s) = q0 + s*t0
 *     
 *  where s the fligth length and q(s) is a point on the trajectory in local UVW 
 *  coordinates. he initial track parameters at s=0 are: 
 * 
 *  q0 = (u,v,0)   and  t0 = (u',v',1)  
 * 
 *  B) Track extrapolation 
 * 
 *  Track extrapolation is just a (nonlinear) mapping of track parameters
 *  
 *  x' = f(x) 
 *  
 *  from one local basis UVW to another one U'V'W'. The mapping f() will be
 *  nonlinear if the planes W=0 and W'=0 are not parallel to each other. One 
 *  way to linearize the problem is to linearize the mapping around a 
 *  reference trajectory 
 *
 *  x(i) = x::ref(i) + x::dev(i)
 *  
 *  Then, the mapping between the reference states between surfaces i and j 
 *  is given by 
 * 
 *  x::ref(j) = f[j,i](x::ref(i)) 
 *  
 * On the other hand, the mapping of the deviations is given by 
 * 
 *
 * x_k+1 = J(k+1,k) * x_k + G_k * w_k                    (*)
 * 
 *  J(j,k) : Transport matrix from k to j -> linearized around ref state  
 *  w_k      : Vector of all scatter angles for scatterers between k->k+1
 *  Q_k      : MSC covariance matrix for w_k -> w_k are zero mean variates 
 *  G_k      : Scatter gain matrix -> linearized around ref state   
 * 
 * The first part in equation (*) is just a straight line extrapolation, 
 * but the second part includes the influence of random scatterings w_k
 * along the way. In particular, w_k includes the scatterings at sensor 
 * k itself. From a physics perspective, x_k estimates 'in' states, i.e. 
 * the track parameters before the scattering in sensor k happens.
 * 
 * 
 * @Author B. Schwenker, University of Göttingen
 * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */


class StraightLineTrackModel : public GenericTrackModel {
  
 public:
  
  /** Default Constructor  
   */
  StraightLineTrackModel() { /* NOOP */ ; }
  
  //! Destructor
  ~StraightLineTrackModel() { /* NOOP */ ; }

  /** Track propagation in magnetic field
   */
  bool BfieldON() { return false;}
  
  /* Returns signed fligth length (mm) to surface fSurf from Surf. Track 
  * state given at Surf. 
  */
  double GetSignedStepLength(const CLHEP::HepMatrix& State, const depfet::ReferenceFrame& Surf, const depfet::ReferenceFrame& fSurf);

  /* Returns true if track hits the surface fSurf.
  */
  bool CheckHitsSurface(const CLHEP::HepMatrix& State, const depfet::ReferenceFrame& Surf, const depfet::ReferenceFrame& fSurf);
  
  /** Returns track state at surface fSurf for track model (State/Surf) 
  */
  CLHEP::HepMatrix Extrapolate(const CLHEP::HepMatrix& State, const depfet::ReferenceFrame& Surf, const depfet::ReferenceFrame& fSurf, bool& error);
  
  /** Extrapolate track model (State/Surf) by a flight length of L. Afterwards, Surf is
  *  normal to Z axis (beam axis). 
  */
  void Extrapolate(CLHEP::HepMatrix& State, depfet::ReferenceFrame& Surf, double length);
  
  /** Compute track derivatives for extrapolation from Surf to fSurf.
  *  Linearization point is State at Surf.  
  */
  int TrackJacobian( const CLHEP::HepMatrix& State, const depfet::ReferenceFrame& Surf, const depfet::ReferenceFrame& fSurf, CLHEP::HepMatrix& J);
  
  /** Get local scatter gain matrix
  *  It calculates the derivates of track parameters State on 
  *  scattering angles theta1 and theta2. 
  */
  void GetScatterGain(const CLHEP::HepMatrix& State, CLHEP::HepMatrix& G);
        
};
 

} // Namespace

#endif 
