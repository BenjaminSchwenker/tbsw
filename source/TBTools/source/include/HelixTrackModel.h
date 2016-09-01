#ifndef HelixTrackModel_H
#define HelixTrackModel_H 1

// DEPFETTrackTools includes
#include "GenericTrackModel.h"

namespace depfet {

 
//! Class HelixTrackModel 
/* 
 *  The class HelixTrackModel implements track dynamics inside a constant 
 *  magnetic field. 
 *  
 *  Relative to the W=0 surface of a local UVW coordinate frame, the local 
 *  track state is defined by 5 variables:
 *  x = (tu, tv, u, v, q/p)
 *  where tu=du/dw and tv=dv/dw in the local detector UVW coordinate system.
 *  
 *  In a (locally) homegeneous magnetic field, a charged particle trajectory 
 *  is a helix: 
 *    
 *  q(s) = q0 + rho*( 1-Cos(theta))*n + Sin(theta)*b )+ theta*Tan(lambda)*h;
 *  
 *  Here, s is the flight length, theta(s) is the turning angle and and q(s)
 *  is a point on the trajectory in local UVW coordinates. The tangent to the
 *  helix is given as
 * 
 *  t(s) = Cos(lambda)*( Sin(theta)*n + Cos(theta)*b + Tan(lambda)*h );
 *     
 *  For s=0, the initial position and flight direction of the particle in UVW 
 *  coordinates are
 * 
 *  q(s=0) = (u_0,v_0,0)   and  t(s=0) = (tu_0,tv=0,1)  
 *  
 *  Track extrapolation is just a mapping of track parameters
 *  
 *  x[i+1] = f_i( x[i] ) 
 *  
 *  from one local basis UVW[i] to track parameters x[i+1] in another basis 
 *  UVW[i+1]. The index i typically numbers the planes along the beam axis. 
 *  
 *  One way to linearize the problem is to linearize the mapping around a 
 *  reference trajectory 
 *
 *  x[i] = x::ref[i] + x::dev[i]
 *  
 *  Then, the mapping between the reference states is given by
 *  
 *  x::ref[i+1] = f_i(x::ref[i]) 
 *  
 *  and the mapping of the deviations is given by 
 * 
 *  x::dev[i+1] = J_i*x::dev[i] + G_i * w_i  
 *                 
 *  
 *  J(i)     : Transport matrix from i to i+1 -> linearized around refernce state  
 *  w_i      : Vector of all scatter angles for scatterers between k->k+1
 *  Q_i      : MSC covariance matrix for w_i -> w_i are zero mean Guassian random numbers 
 *  G_i      : Scatter gain matrix -> linearized around reference state   
 *  
 *  The first part in equation (*) is just a helix extrapolation. 
 *  The second part includes the influence of random scatterings w_i
 *  along the way. In particular, w_i includes the scatterings at the 
 *  w_i=0 plane itself.
 *  
 *  From a physics perspective, x[i] estimates 'in' states, i.e. 
 *  the track parameters before the scattering in sensor i happens.
 * 
 * 
 * @Author B. Schwenker, University of Göttingen
 * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */


class HelixTrackModel : public GenericTrackModel {
  
 public:
  
  /** Default Constructor  
   */
  HelixTrackModel(const CLHEP::HepVector& field);
  
  /** Track propagation in magnetic field
   */
  bool BfieldON() { return true;}
  
  //! Destructor
  ~HelixTrackModel() { /* NOOP */ ; }
  
  /* Returns signed fligth length (mm) to surface fSurf. Track starts at surface Surf and has state State.
  */
  double GetSignedStepLength(const CLHEP::HepMatrix& State, const depfet::ReferenceFrame& Surf, const depfet::ReferenceFrame& fSurf);
  
  /* Returns true if track hits the surface fSurf.
  */
  bool CheckHitsSurface(const CLHEP::HepMatrix& State, const depfet::ReferenceFrame& Surf, const depfet::ReferenceFrame& fSurf);
  
  /** Returns track state at surface fSurf. 
   */
  CLHEP::HepMatrix Extrapolate(const CLHEP::HepMatrix& State, const depfet::ReferenceFrame& Surf, const depfet::ReferenceFrame& fSurf, bool& error);
  
  /** Extrapolate track along helix for given flight length. Track parameters (State/Surf) are overwritten. 
  */
  void Extrapolate(CLHEP::HepMatrix& State, depfet::ReferenceFrame& Surf,  double length);
  
  /** Compute track derivatives for extrapolation from Surf to fSurf.
  *  Linearization point is State at Surf.  
  */
  int TrackJacobian( const CLHEP::HepMatrix& State, const depfet::ReferenceFrame& Surf, const depfet::ReferenceFrame& fSurf ,  CLHEP::HepMatrix& J);
  
  /** Get local scatter gain matrix
  *  It calculates the derivates of track parameters State on 
  *  scattering angles theta1 and theta2. 
  */
  void GetScatterGain(const CLHEP::HepMatrix& State, CLHEP::HepMatrix& G);
  
 // Private Functions --------------
 private: 
  
  /* Bisection method to find the root of the function GetDistanceToPlane. 
  */
  double BisectionMethod( 
   const CLHEP::HepMatrix& State, const depfet::ReferenceFrame& Surf, const depfet::ReferenceFrame& fSurf, 
   double maxlength);   
  
  /* Returns distance (mm) to surface fSurf after extrapolating track 
  * from Surf for some flight length 
  */
  double GetDistanceToPlane(const CLHEP::HepMatrix& State, const depfet::ReferenceFrame& Surf, const depfet::ReferenceFrame& fSurf, double length);
  
 // Private Members -----------------
 private:
  
  // Magentic field in Tesla
  CLHEP::HepVector Bfield;
  int ndim; 
  
};
 

} // Namespace

#endif 
