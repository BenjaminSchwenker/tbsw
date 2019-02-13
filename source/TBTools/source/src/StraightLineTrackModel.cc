// StraightLineTrackModel implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// TBTool includes
#include "StraightLineTrackModel.h"


#include <cmath>
#include <iostream>

// Namespaces
using namespace std; 

namespace depfet {

/* Returns signed fligth length (mm) to surface fSurf from Surf. Track 
 * state given at Surf. 
 */
double StraightLineTrackModel::GetSignedStepLength(const TrackState& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf)
{
  
  // Compute the intersection point of the trajectory with final 
  // surface fSurf
  
  // 3D representation of track: r(l) = xPoint + l*NDir
  // Note that NDir is normalized to unit length and 
  // l is the signed track length in mm. 
  
  // Intersection point on initial surface 
  Vector3d xPoint;
  xPoint<< State[2], State[3], 0;
  
  // Track direction on initial surface
  Vector3d Dir;
  Dir<< State[0], State[1], 1;
  Dir.normalize(); // Unit lenght vector

  
  // Get coord trafo from Surf to fSurf
  // perhaps use auto for the next two. Could help.
  Matrix3d Rot = fSurf.GetRotation()*Surf.GetRotation().transpose();
  Vector3d Origin = Surf.GetRotation()*(fSurf.GetPosition()-Surf.GetPosition());

  // Track direction on final surface   
  auto fDir = Rot*Dir;
  if (fDir[2] == 0) { 
    // No intersection with final plane  
    return 0;   
  }
   
  // Surface normal vector in beam direction 
  auto W = Vector3d::UnitZ();

  
  // Signed flight length between surfaces  
  double length = (Rot*(Origin-xPoint)).dot(W) / fDir[2];
  
  return length;  
}

/* Returns true if track hits the surface fSurf.
 */
bool StraightLineTrackModel::CheckHitsSurface(const TrackState& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf)
{
   
  // Track direction on initial surface
  Vector3d Dir;
  Dir<< State[0], State[1], 1;

  // Get coord trafo from Surf to fSurf
  Matrix3d Rot = fSurf.GetRotation()*Surf.GetRotation().transpose();
  
  // Track direction on final surface   
  auto fDir = Rot*Dir;
  if (fDir[2] == 0) { 
    // No intersection with final plane  
    return false;   
  }
  
  return true; 
}


/** Returns track state at surface fSurf from state at surface Surf 
 */
TrackState StraightLineTrackModel::Extrapolate(const TrackState& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf, bool& error)
{

  // No problem so far
  error = false; 
  
  // Compute 3D representation of track: r(s) = xPoint + s*Dir
  // All variables refer to initial surface coord.   

  // Intersection point on reference surface 
  Vector3d xPoint;
  xPoint<< State[2], State[3], 0;
  
  // Track direction on reference surface
  Vector3d Dir;
  Dir<< State[0], State[1], 1;
  
  // Get coord trafo from Surf to fSurf
  Matrix3d Rot = fSurf.GetRotation()*Surf.GetRotation().transpose();
  Vector3d Origin = Surf.GetRotation()*(fSurf.GetPosition()-Surf.GetPosition());
  
  // Track direction on final surface   
  auto fDir = Rot*Dir;
  if (fDir[2] == 0) { 
    // No intersection with fSurf
    error = true;
    return TrackState::Zero();
  }
   
  // Surface normal vector in beam direction 
  auto W=Vector3d::UnitZ();

  
  // This is like step lenght 
  double SX = (Rot*(Origin-xPoint)).dot(W) / fDir[2];
  
  // Intersection point with fSurf
  Vector3d fPoint = Rot*(xPoint+SX*Dir-Origin);
  
  // Track parameters on fSurf
  TrackState fState;
  fState << fDir[0]/fDir[2], fDir[1]/fDir[2], fPoint[0], fPoint[1], State[4];
  
  return fState; 
}




/** Extrapolate track model (State/Surf) by a flight length of L. Afterwards, Surf is normal 
 *  to Z axis (beam axis). 
 */
void StraightLineTrackModel::Extrapolate(TrackState& State, ReferenceFrame& Surf, double Length)
{
  
  // Get trafo from global coord to reference surface coord  
  const Matrix3d& Rot = Surf.GetRotRef();
  const Vector3d& Origin = Surf.GetPosRef();
  
  // Compute 3D representation of track: r(l) = xPoint + l*Dir
  

  // Intersection point on Surf
  Vector3d xPoint;
  xPoint<< State[2], State[3], 0;

   // Track direction on Surf
  Vector3d Dir;
  Dir<< State[0], State[1], 1;
  
  // Endpoint of extrapolation step in global coord.   
  Vector3d gPoint = Rot.transpose()*(xPoint+Length*Dir) + Origin ;
  
  // Track direction in global coord. 
  Vector3d gDir = Rot.transpose()*Dir;
  
  // Overwrite Surf  
  Surf.SetPosition(gPoint);              
  Surf.SetRotation(Matrix3d::Identity());  // Surface normal to Z axis
   
  // Overwrite State
  State[0] = gDir[0]/gDir[2];
  State[1] = gDir[1]/gDir[2];
  State[2] = 0;                       // Track crosses origin of Surf
  State[3] = 0;
 
  return;   
}


/** Compute track derivatives for extrapolation from initial surface Surf to 
 *  final surface fSurf. The expansion is around state p.  
 */
int StraightLineTrackModel::TrackJacobian( const TrackState& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf, TrackStateJacobian& J)
{
  
  // Get coord trafo from Surf to fSurf
  Matrix3d Rot = fSurf.GetRotation()*Surf.GetRotation().transpose();
  Vector3d Origin = Surf.GetRotation()*(fSurf.GetPosition()-Surf.GetPosition());
  
  // Intersection point on initial surface 
  Vector3d xPoint;
  xPoint<< State[2], State[3], 0;

  // Trajectory direction on initial surface
  Vector3d Dir;
  Dir<< State[0], State[1], 1;

  
  // Trajectory direction on final surface   
  auto fDir = Rot*Dir;
  if (fDir[2] == 0) { // No intersection with final plane 
    return 1;
  }
   
  // Surface normal vector in beam direction 
  auto W=Vector3d::UnitZ();
  
  double SX = (Rot*(Origin-xPoint)).dot( W) / fDir[2];
  
  // DEFINE TRANSPORTATION MATRIX OF TRAJ. PARAMETRTS FROM 
  // CURRENT SURFACE TO TARGET SURFACE 

  // For details of derivation, please refer to V. Karimäki "Straight Line 
  // Fit for Pixel and Strip Detectors with Arbitrary Plane Orientations", CMS Note.
  // Parameters at target surface:   U', V', U, V, Q/P
  // Parameters at current surface:  u', v', u, v, q/p  
  
  // Some abbreviations to save computation time
  
  double Up = fDir[0]/fDir[2];
  double Vp = fDir[1]/fDir[2];     
  J=TrackStateJacobian::Zero();
  // *** dU'/du'
  J(0,0) = (Rot(0,0)*fDir[2]-fDir[0]*Rot(2,0)) / (fDir[2]*fDir[2]) ;
  // *** dU'/dv'
  J(0,1) = (Rot(0,1)*fDir[2]-fDir[0]*Rot(2,1)) / (fDir[2]*fDir[2]) ;
  // *** dU'/du
  J(0,2) = 0;
  // *** dU'/dv
  J(0,3) = 0;
  
  // *** dV'/du' 
  J(1,0) = (Rot(1,0)*fDir[2]-fDir[1]*Rot(2,0)) / (fDir[2]*fDir[2]) ;
  // *** dV'/dv' 
  J(1,1) = (Rot(1,1)*fDir[2]-fDir[1]*Rot(2,1)) / (fDir[2]*fDir[2]) ;
  // *** dV'/du
  J(1,2) = 0;
  // *** dV'/dV
  J(1,3) = 0;
    
  // *** dU/du
  J(2,2) = Rot(0,0) - Rot(2,0) * Up;
  // *** dU/dv
  J(2,3) = Rot(0,1) - Rot(2,1) * Up;
  // *** dU/du'
  J(2,0) = SX * J(2,2);
  // *** dU/dv'
  J(2,1) = SX * J(2,3);
  
  // *** dV/du
  J(3,2) = Rot(1,0) - Rot(2,0) *  Vp;
  // *** dV/dv
  J(3,3) = Rot(1,1) - Rot(2,1) * Vp;
  // *** dV/du'
  J(3,0) = SX * J(3,2);
  // *** dV/dv'
  J(3,1) = SX * J(3,3);

  // *** d(Q/P)/d(q/p)
  J(4,4) = 1;

    
  // Everything ok
  return 0;
}


/** Numerical computation of track derivatives for extrapolation from Surf to fSurf. 
 *  Linearization point is State at Surf.  
 */
int StraightLineTrackModel::TrackJacobianNumerical( const TrackState& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf, TrackStateJacobian& J)
{
 

  TrackState  tmpState = State;
  bool error = false; 
   
  // Extrapolate track along helix to next detector plane 
  TrackState statePred = Extrapolate(tmpState, Surf, fSurf,  error);
  if (error) {
    std::cout << "ERR: Propagation to next detector failed. Quit fitting!" << std::endl;
    return 1;   
  }   

  TrackState difPred = statePred;

  J= TrackStateJacobian::Zero();

  // do column wise differntiation:
  for(int icol=0;icol<5;++icol){
    // choose step
    double h=std::fabs(tmpState[icol])*1.e-4;
    if(h<1e-8)h=1.e-8;
    // vary the state
    tmpState[icol]+=h;
    // extrap state 
    difPred = Extrapolate(tmpState, Surf, fSurf,  error); 
    if (error) {
      std::cout  << "ERR: Propagation to next detector failed. Quit fitting!" << std::endl; 
      return 1;  
    }       
    // difference
    difPred-=statePred;
    // remove variation from state
    tmpState[icol]-=h;
    // fill jacobian with difference quotient
    for(int irow=0;irow<5;++irow)
        J(irow,icol)=difPred[irow]/h;
  }
  
      
  // Everything ok
  return 0;
}
 

/** Get local scatter gain matrix
 *  It calculates the derivates of track parameters State on 
 *  scattering angles theta1 and theta2. 
 */
TrackStateGain StraightLineTrackModel::GetScatterGain(const TrackState& State)
{
   
  // Calculation of scatter gain matrix following Wolin and Ho (NIM A329 (1993) 493-500)
  // Transverse displacement of the track is ignored. A thin scatter is assumed. 
   
  // P3 and P4 are local direction tangents 
  double p3 = State[0];
  double p4 = State[1];
     
  // Needed is the rotation matrix comoving frame -> detector frame.  
  // We construct a basis for the comoving frame:
   
  // n_trk is parallel to unscattered track direction
  Vector3d n_trk;
  n_trk<< p3, p4, 1;
  n_trk.normalize();

  // v_trk is orthogonal to track dir and detector u axis  
  auto u_hat=Vector3d::UnitX();
  Vector3d v_trk = n_trk.cross(u_hat).normalized();

  // u_trk completes rigth handed system  
  Vector3d u_trk = v_trk.cross(n_trk);
  
  // Now, we can directly read off the direction cosines 
  double a1 = u_trk[0];  
  double a2 = v_trk[0];   
  double a3 = n_trk[0];
      
  double b1 = u_trk[1];  
  double b2 = v_trk[1];   
  double b3 = n_trk[1];   
  
  double g1 = u_trk[2];  
  double g2 = v_trk[2];   
  double g3 = n_trk[2];    
  
  // The two independent scattering angles are called 
  // theta1 and theta2. A change in these is affecting
  // the scattered track parameters ...
  TrackStateGain G;    
  
  // *** dU'/dtheta1
  G(0,0) = ( a1*g3 - a3*g1 ) / ( g3*g3 );
  // *** dU'/dtheta2
  G(0,1) = ( a2*g3 - a3*g2 ) / ( g3*g3 );
  // *** dV'/dtheta1
  G(1,0) = ( b1*g3 - b3*g1 ) / ( g3*g3 );
  // *** dV'/dtheta2
  G(1,1) = ( b2*g3 - b3*g2 ) / ( g3*g3 );

  // Scattering angles affect the track
  // do not affect impact point
  // *** dU/dtheta1
  G(2,0) = 0;
  // *** dU/dtheta2
  G(2,1) = 0;
  // *** dV/dtheta1
  G(3,0) = 0;
  // *** dV/dtheta2
  G(3,1) = 0;
    
  return G;
}




} // Namespace;

