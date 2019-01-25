// HelixTrackModel implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// TBTool includes
#include "HelixTrackModel.h"
#include "MaterialEffect.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "TMath.h"

// Namespaces
using namespace std; 

namespace depfet {

/** Constructor 
 */
HelixTrackModel::HelixTrackModel(const Vector3d& field) :  GenericTrackModel()
{ 
  // Magnetic field in Tesla
  Bfield = field;
  m_kappa = 1; 
  ndim = 5; // dimension of state vector
}


/* Returns signed fligth length (mm) of track to surface fSurf. Track starts at surface Surf and has state State.
 */
double HelixTrackModel::GetSignedStepLength(const TrackState& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf)
{
  // The signed length from Surf to fSurf must be smaller than that
  double maxlength = 2*(fSurf.GetZPosition() - Surf.GetZPosition());
  
  // Try to find signed flight length by a simpel numberical bisection
  double length = BisectionMethod(State, Surf, fSurf, maxlength);
  
  return length;
}



/** Returns track state at surface fSurf from state at surface Surf 
 */
TrackState HelixTrackModel::Extrapolate(const TrackState& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf, bool& error)
{

  error = false;

  // Check that surfaces are different
  if ( Surf == fSurf  ) { 
    TrackState newState = State;
    return newState; 
  }
  
  // Check that the helix really intersects surface fSurf
  if ( !CheckHitsSurface(State, Surf, fSurf) ) {
    error = true;
    return TrackState::Zero();
  }
    
  // Compute the signed step length in mm 
  double length = GetSignedStepLength(State, Surf, fSurf);
      
  // Extrapolate a copy of the current state and surface along 
  // the helix
  
  TrackState tmpState = State;
  ReferenceFrame tmpSurf = Surf;  
  Extrapolate(tmpState, tmpSurf,  length);
  
  // Track position in global coordinates 
  auto x = tmpSurf.GetPosition();
  
  // Track direction in global coord
  // Important: We assume that we see only the part of the helix where
  // the detector surface is traversed in +z direction.  
  Vector3d t;
  t<< tmpState[0], tmpState[1], 1;
    
  // Finally rotate the temporary state into the 
  // local coordinates of fSurf 
   
  auto xloc = fSurf.TransformPointToLocal(x);
  auto tloc = fSurf.TransformVecToLocal(t);
    
  // Overwrite State
  tmpState[0] = tloc[0]/tloc[2];
  tmpState[1] = tloc[1]/tloc[2];
  tmpState[2] = xloc[0];
  tmpState[3] = xloc[1];
  
  return tmpState; 
}


/** Extrapolate track along helix for flight length of L. Track parameters (State/Surf) are overwritten. 
 */
void HelixTrackModel::Extrapolate(TrackState& State, ReferenceFrame& Surf,  double length)
{
  
  // Track position in Surf coordinates
  Vector3d xloc;
  xloc<< State[2], State[3], 0;

  
  // Track direction in Surf coordinates
  Vector3d tloc;
  tloc<< State[0], State[1], 1;
  
  // Transform track state into global coordintes  
  auto x = Surf.TransformPointToGlobal(xloc);
  auto t = Surf.TransformVecToGlobal(tloc);
  
  // Construct basis vectors (n, b, h) of helix frame 
  
  // Important: t points in the direction of fligth
  // and has unit length 
  t.normalize();
  
  // Important: basic assumption tha the track intersects
  // surface in +z direction. 
  if (t[2] < 0  ) t*=-1;
                   
  double Bnorm = Bfield.norm(); 
  Vector3d h = Bfield/Bnorm;
    
  // n = t x h
  Vector3d n=t.cross(n).normalized();

  
  // b = h x n
  Vector3d b=h.cross(n).normalized();

  
  // Extrapolate along helix in global coordinates 
    
  double coslambda = t.dot(b);
  double sinlambda = t.dot(h);
  double tanlambda =  sinlambda / coslambda;
  double kappa = State[4];
  double theta = 0.29979*Bnorm*length*kappa/1000; 
  double rho = 1000*coslambda/kappa/0.29979/Bnorm;
  
  auto xnew = x + rho*(1-TMath::Cos(theta))*n + rho*TMath::Sin(theta)*b + rho*theta*tanlambda*h;
  auto tnew = coslambda*( TMath::Sin(theta)*n + TMath::Cos(theta)*b + tanlambda*h );
  
  // Overwrite Surf  
  Surf.SetPosition(xnew);                // Track crosses origin of Surf 
  Surf.SetRotation(Matrix3d::Identity());  // Surface normal to Z axis
   
  // Overwrite State
  State[0] = tnew[0]/tnew[2];
  State[1] = tnew[1]/tnew[2];
  State[2] = 0;
  State[3] = 0;
  
  return;   
}




/* Returns true if track hits the surface fSurf.
 */
bool HelixTrackModel::CheckHitsSurface(const TrackState& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf)
{
  // This is a simplified test to see if the helix hits surface fSurf. 
  // From a geometrical point of view, the flight distance should be 
  // smaller than two times the z distance between surfaces Surf and
  // fSurf.  
  double maxlength = 2*(fSurf.GetZPosition() - Surf.GetZPosition());

  // This is the distance from start point on helix (l=0) to surface fSurf
  double d1 = GetDistanceToPlane(State, Surf, fSurf,  0);

  // This is the distance from anotther point (l=maxlenght) on helix to surface fSurf
  double d2 = GetDistanceToPlane(State, Surf, fSurf,  maxlength);
  
  // An intersection needs a sign swap in the flight length intervall l=0 and l=maxlength
  if ( d1*d2 > 0.0) {
    return false; 
  } else {
    return true; 
  }
  
}



/* Returns distance (mm) to surface fSurf after extrapolating track 
 * some flight length 
 */
double HelixTrackModel::GetDistanceToPlane(const TrackState& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf,  double length)
{

  // Extrapolate a copy of the current state and surface along 
  // the helix. 
  auto tmpState = State;
  ReferenceFrame tmpSurf = Surf;  
  Extrapolate(tmpState, tmpSurf, length);

  // Track position in global coordinates 
  Vector3d xnew = tmpSurf.GetPosition();
    
  // Transform track state into loccal coordintes of fSurf  
  auto D = fSurf.TransformPointToLocal(xnew);

  // Return signed distance of xnew to surface fSurf
  return D[2];
}

/* Bisection method to find the root of the function GetDistanceToPlane
 */
double HelixTrackModel::BisectionMethod( 
  const TrackState& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf,
  double maxlength)
{
  
  // Numerical precision in mm
  //static const double about_zero_mag = 1E-3;
  static const double about_zero_mag = 1E-9;
  
  double pos_pt = maxlength;
  double neg_pt = 0;
   
  if ( GetDistanceToPlane(State,Surf, fSurf, pos_pt)*GetDistanceToPlane(State,Surf, fSurf, neg_pt) > 0  ) {
    cout << "NO ZERO IN BISECTION: THIS SHOULD NEVER HAPPEN " << endl;
    return 0;     
  }
  
  if ( GetDistanceToPlane(State,Surf, fSurf, pos_pt) < 0.0) {
    pos_pt = 0;  
    neg_pt = maxlength;
  }
  
  // main loop
  for (;;)
  {    
    double mid_pt = (pos_pt + neg_pt)/2.0;
    double f_mid_pt = GetDistanceToPlane(State,Surf, fSurf,  mid_pt);   

    if (fabs(f_mid_pt) < about_zero_mag)
      return mid_pt;
    
    if (f_mid_pt >= 0.0)
      pos_pt = mid_pt;
    else
    neg_pt = mid_pt;
  }
    
}



/** Numerical computation of track derivatives for extrapolation from Surf to fSurf. 
 *  Linearization point is State at Surf.  
 */
int HelixTrackModel::TrackJacobian( const TrackState& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf, TrackStateJacobian& J)
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
  for(int icol=0;icol<ndim;++icol){
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
    for(int irow=0;irow<ndim;++irow)
        J(irow,icol)=difPred[irow]/h;
  }
  
      
  // Everything ok
  return 0;
}
 

/** Get local scatter gain matrix
 *  It calculates the derivates of track parameters State on 
 *  scattering angles theta1 and theta2. 
 */

/*
void HelixTrackModel::GetScatterGain(const HepMatrix& State, CLHEP::HepMatrix& G)
{
   
  // The two independent scattering angles are called 
  // theta1 and theta2. A change in these is affecting
  // the scattered track parameters ...
   
  // To hold the scattered state 
  HepMatrix difPred = State;
   
  // The scattering angle in rad
  double h=1.e-5;
   
  {
    // scatter state 
    materialeffect::ScatterTrack(difPred, h, 0);      
    // difference
    difPred-=State;
    // fill gain matrix with difference quotient
    for(int irow=0;irow<ndim;++irow) G[irow][0]=difPred[irow][0]/h;
  }
  
  // reset variation
  difPred = State; 
   
  {
    // scatter state 
    materialeffect::ScatterTrack(difPred, 0, h);      
    // difference
    difPred-=State;
    // fill gain matrix with difference quotient
    for(int irow=0;irow<ndim;++irow) G[irow][1]=difPred[irow][0]/h;
  }
     
  // Everything ok
  return ;
}
*/


void HelixTrackModel::GetScatterGain(const TrackState& State, TrackStateGain& G)
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
  n_trk.normalize(); // Change vector to unit size
   
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
    
  return ;
}



} // Namespace;

