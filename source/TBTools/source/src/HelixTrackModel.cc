// HelixTrackModel implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// DEPFETTrackTool includes
#include "HelixTrackModel.h"
#include "MaterialEffect.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "TMath.h"

// Namespaces
using namespace CLHEP;
using namespace std; 

namespace depfet {

/** Constructor 
 */
HelixTrackModel::HelixTrackModel(const HepVector& field) :  GenericTrackModel()
{ 
  // Magnetic field in Tesla
  Bfield = field;
  m_kappa = 1; 
  ndim = 5; // dimension of state vector
}


/* Returns signed fligth length (mm) of track to surface fSurf. Track starts at surface Surf and has state State.
 */
double HelixTrackModel::GetSignedStepLength(const HepMatrix& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf)
{
  // The signed length from Surf to fSurf must be smaller than that
  double maxlength = 2*(fSurf.GetZPosition() - Surf.GetZPosition());
  
  // Try to find signed flight length by a simpel numberical bisection
  double length = BisectionMethod(State, Surf, fSurf, maxlength);
  
  return length;
}



/** Returns track state at surface fSurf from state at surface Surf 
 */
HepMatrix HelixTrackModel::Extrapolate(const HepMatrix& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf, bool& error)
{

  error = false;

  // Check that surfaces are different
  if ( Surf == fSurf  ) { 
    HepMatrix newState = State;
    return newState; 
  }
  
  // Check that the helix really intersects surface fSurf
  if ( !CheckHitsSurface(State, Surf, fSurf) ) {
    error = true;
    return HepMatrix(ndim,1,0);
  }
    
  // Compute the signed step length in mm 
  double length = GetSignedStepLength(State, Surf, fSurf);
      
  // Extrapolate a copy of the current state and surface along 
  // the helix
  
  HepMatrix tmpState = State; 
  ReferenceFrame tmpSurf = Surf;  
  Extrapolate(tmpState, tmpSurf,  length);
  
  // Track position in global coordinates 
  HepVector x = tmpSurf.GetPosition();
  
  // Track direction in global coord
  // Important: We assume that we see only the part of the helix where
  // the detector surface is traversed in +z direction.  
  HepVector t(3,0);    
  t[0] = tmpState[0][0]; 
  t[1] = tmpState[1][0];
  t[2] = 1;                    
    
  // Finally rotate the temporary state into the 
  // local coordinates of fSurf 
   
  HepVector xloc = fSurf.TransformPointToLocal(x); 
  HepVector tloc = fSurf.TransformVecToLocal(t); 
    
  // Overwrite State
  tmpState[0][0] = tloc[0]/tloc[2]; 
  tmpState[1][0] = tloc[1]/tloc[2];
  tmpState[2][0] = xloc[0]; 
  tmpState[3][0] = xloc[1]; 
  
  return tmpState; 
}


/** Extrapolate track along helix for flight length of L. Track parameters (State/Surf) are overwritten. 
 */
void HelixTrackModel::Extrapolate(HepMatrix& State, ReferenceFrame& Surf,  double length)
{
  
  // Track position in Surf coordinates
  HepVector xloc(3,0);
  xloc[0] = State[2][0];
  xloc[1] = State[3][0];
  xloc[2] = 0;
  
  // Track direction in Surf coordinates
  HepVector tloc(3,0);    
  tloc[0] = State[0][0]; 
  tloc[1] = State[1][0];
  tloc[2] = 1; 
  
  // Transform track state into global coordintes  
  HepVector x = Surf.TransformPointToGlobal(xloc); 
  HepVector t = Surf.TransformVecToGlobal(tloc);
  
  // Construct basis vectors (n, b, h) of helix frame 
  
  // Important: t points in the direction of fligth
  // and has unit length 
  t /= t.norm();
  
  // Important: basic assumption tha the track intersects
  // surface in +z direction. 
  if (t[2] < 0  ) t*=-1;
                   
  double Bnorm = Bfield.norm(); 
  HepVector h = Bfield/Bnorm;
    
  // n = t x h
  HepVector n(3,0);  
  n[0] = t[1]*h[2] - t[2]*h[1]; 
  n[1] = t[2]*h[0] - t[0]*h[2]; 
  n[2] = t[0]*h[1] - t[1]*h[0]; 
  n /= n.norm();
  
  // b = h x n
  HepVector b(3,0);  
  b[0] = h[1]*n[2] - h[2]*n[1]; 
  b[1] = h[2]*n[0] - h[0]*n[2]; 
  b[2] = h[0]*n[1] - h[1]*n[0];   
  b /= b.norm();
  
  // Extrapolate along helix in global coordinates 
    
  double coslambda = t[0]*b[0] + t[1]*b[1] + t[2]*b[2];  
  double sinlambda = t[0]*h[0] + t[1]*h[1] + t[2]*h[2];
  double tanlambda =  sinlambda / coslambda;
  double kappa = State[4][0];
  double theta = 0.29979*Bnorm*length*kappa/1000; 
  double rho = 1000*coslambda/kappa/0.29979/Bnorm;
  
  HepVector xnew = x + rho*(1-TMath::Cos(theta))*n + rho*TMath::Sin(theta)*b + rho*theta*tanlambda*h; 
  HepVector tnew = coslambda*( TMath::Sin(theta)*n + TMath::Cos(theta)*b + tanlambda*h );
  
  // Overwrite Surf  
  Surf.SetPosition(xnew);                // Track crosses origin of Surf 
  Surf.SetRotation(HepMatrix(3, 3, 1));  // Surface normal to Z axis 
   
  // Overwrite State
  State[0][0] = tnew[0]/tnew[2]; 
  State[1][0] = tnew[1]/tnew[2];
  State[2][0] = 0;                          
  State[3][0] = 0; 
  
  return;   
}




/* Returns true if track hits the surface fSurf.
 */
bool HelixTrackModel::CheckHitsSurface(const HepMatrix& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf)
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
double HelixTrackModel::GetDistanceToPlane(const HepMatrix& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf,  double length)
{

  // Extrapolate a copy of the current state and surface along 
  // the helix. 
  HepMatrix tmpState = State; 
  ReferenceFrame tmpSurf = Surf;  
  Extrapolate(tmpState, tmpSurf, length);

  // Track position in global coordinates 
  HepVector xnew = tmpSurf.GetPosition();
    
  // Transform track state into loccal coordintes of fSurf  
  HepVector D = fSurf.TransformPointToLocal(xnew); 

  // Return signed distance of xnew to surface fSurf
  return D[2];
}

/* Bisection method to find the root of the function GetDistanceToPlane
 */
double HelixTrackModel::BisectionMethod( 
  const HepMatrix& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf, 
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
int HelixTrackModel::TrackJacobian( const HepMatrix& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf, HepMatrix& J)
{
 
  HepMatrix j(ndim,ndim); 
  HepMatrix tmpState = State;
  bool error = false; 
   
  // Extrapolate track along helix to next detector plane 
  HepMatrix statePred = Extrapolate(tmpState, Surf, fSurf,  error);
  if (error) {
    std::cout << "ERR: Propagation to next detector failed. Quit fitting!" << std::endl;
    return 1;   
  }   

  HepMatrix difPred = statePred;
   
  // do column wise differntiation:
  for(int icol=0;icol<ndim;++icol){
    // choose step
    double h=std::fabs(tmpState[icol][0])*1.e-4;
    if(h<1e-8)h=1.e-8;
    // vary the state
    tmpState[icol][0]+=h;
    // extrap state 
    difPred = Extrapolate(tmpState, Surf, fSurf,  error); 
    if (error) {
      std::cout  << "ERR: Propagation to next detector failed. Quit fitting!" << std::endl; 
      return 1;  
    }       
    // difference
    difPred-=statePred;
    // remove variation from state
    tmpState[icol][0]-=h;
    // fill jacobian with difference quotient
    for(int irow=0;irow<ndim;++irow)j[irow][icol]=difPred[irow][0]/h;
  }
  
  // The final track jacobian 
  J = j; 

  
      
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


void HelixTrackModel::GetScatterGain(const HepMatrix& State, CLHEP::HepMatrix& G)
{
   
  // Calculation of scatter gain matrix following Wolin and Ho (NIM A329 (1993) 493-500)
  // Transverse displacement of the track is ignored. A thin scatter is assumed. 
   
  // P3 and P4 are local direction tangents 
  double p3 = State[0][0];    
  double p4 = State[1][0]; 
     
  // Needed is the rotation matrix comoving frame -> detector frame.  
  // We construct a basis for the comoving frame:
   
  // n_trk is parallel to unscattered track direction
  Hep3Vector n_trk;    
  n_trk[0] = p3;    
  n_trk[1] = p4;    
  n_trk[2] = 1;  
  n_trk /= n_trk.mag(); // Change vector to unit size
   
  // v_trk is orthogonal to track dir and detector u axis  
  Hep3Vector u_hat(1,0,0);
  Hep3Vector v_trk = n_trk.cross(u_hat);
  v_trk /= v_trk.mag();
  
  // u_trk completes rigth handed system  
  Hep3Vector u_trk = v_trk.cross(n_trk);   
  
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
  G[0][0] = ( a1*g3 - a3*g1 ) / ( g3*g3 );
  // *** dU'/dtheta2 
  G[0][1] = ( a2*g3 - a3*g2 ) / ( g3*g3 );
  // *** dV'/dtheta1     
  G[1][0] = ( b1*g3 - b3*g1 ) / ( g3*g3 );
  // *** dV'/dtheta2 
  G[1][1] = ( b2*g3 - b3*g2 ) / ( g3*g3 );
  
  // Scattering angles affect the track 
  // do not affect impact point 
  // *** dU/dtheta1     
  G[2][0] = 0;
  // *** dU/dtheta2 
  G[2][1] = 0;
  // *** dV/dtheta1     
  G[3][0] = 0;
  // *** dV/dtheta2   
  G[3][1] = 0;    
    
  return ;
}



} // Namespace;

