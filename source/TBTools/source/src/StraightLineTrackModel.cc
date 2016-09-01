// StraightLineTrackModel implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// DEPFETTrackTool includes
#include "StraightLineTrackModel.h"

#include "CLHEP/Vector/ThreeVector.h"

#include <cmath>

// Namespaces
using namespace CLHEP;
using namespace std; 

namespace depfet {

/* Returns signed fligth length (mm) to surface fSurf from Surf. Track 
 * state given at Surf. 
 */
double StraightLineTrackModel::GetSignedStepLength(const HepMatrix& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf)
{
  
  // Compute the intersection point of the trajectory with final 
  // surface fSurf
  
  // 3D representation of track: r(l) = xPoint + l*NDir
  // Note that NDir is normalized to unit length and 
  // l is the signed track length in mm. 
  
  // Intersection point on initial surface 
  HepVector xPoint(3,0);
  xPoint[0] = State[2][0];
  xPoint[1] = State[3][0];
  xPoint[2] = 0;
  
  // Track direction on initial surface
  HepVector Dir(3,0);
  Dir[0] = State[0][0];
  Dir[1] = State[1][0];
  Dir[2] = 1;
  Dir /= Dir.norm();  // Unit lenght vector  
  
  // Get coord trafo from Surf to fSurf
  HepMatrix Rot = fSurf.GetRotation()*Surf.GetRotation().T();  
  HepVector Origin = Surf.GetRotation()*(fSurf.GetPosition()-Surf.GetPosition());
  
  // Track direction on final surface   
  HepVector fDir = Rot*Dir;
  if (fDir[2] == 0) { 
    // No intersection with final plane  
    return 0;   
  }
   
  // Surface normal vector in beam direction 
  HepVector W(3,0); 
  W[2] = 1;  
  
  // Signed flight length between surfaces  
  double length = dot(Rot*(Origin-xPoint), W) / fDir[2];  
  
  return length;  
}

/* Returns true if track hits the surface fSurf.
 */
bool StraightLineTrackModel::CheckHitsSurface(const HepMatrix& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf)
{
   
  // Track direction on initial surface
  HepVector Dir(3,0);
  Dir[0] = State[0][0];
  Dir[1] = State[1][0];
  Dir[2] = 1;
   
  // Get coord trafo from Surf to fSurf
  HepMatrix Rot = fSurf.GetRotation()*Surf.GetRotation().T();  
  
  // Track direction on final surface   
  HepVector fDir = Rot*Dir;
  if (fDir[2] == 0) { 
    // No intersection with final plane  
    return false;   
  }
  
  return true; 
}


/** Returns track state at surface fSurf from state at surface Surf 
 */
HepMatrix StraightLineTrackModel::Extrapolate(const HepMatrix& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf, bool& error)
{

  // No problem so far
  error = false; 
  
  // Compute 3D representation of track: r(s) = xPoint + s*Dir
  // All variables refer to initial surface coord.   

  // Intersection point on reference surface 
  HepVector xPoint(3,0);
  xPoint[0] = State[2][0];
  xPoint[1] = State[3][0];
  xPoint[2] = 0;
  
  // Track direction on reference surface
  HepVector Dir(3,0);
  Dir[0] = State[0][0];
  Dir[1] = State[1][0];
  Dir[2] = 1;
  
  // Get coord trafo from Surf to fSurf
  HepMatrix Rot = fSurf.GetRotation()*Surf.GetRotation().T();  
  HepVector Origin = Surf.GetRotation()*(fSurf.GetPosition()-Surf.GetPosition());
  
  // Track direction on final surface   
  HepVector fDir = Rot*Dir;
  if (fDir[2] == 0) { 
    // No intersection with fSurf
    error = true;
    return HepMatrix();
  }
   
  // Surface normal vector in beam direction 
  HepVector W(3,0); 
  W[2] = 1;  
  
  // This is like step lenght 
  double SX = dot(Rot*(Origin-xPoint), W) / fDir[2];  
  
  // Intersection point with fSurf
  HepVector fPoint = Rot*(xPoint+SX*Dir-Origin);  
  
  // Track parameters on fSurf
  HepMatrix fState = State;  
  fState[0][0] = fDir[0]/fDir[2]; 
  fState[1][0] = fDir[1]/fDir[2];
  fState[2][0] = fPoint[0]; 
  fState[3][0] = fPoint[1]; 
  
  return fState; 
}




/** Extrapolate track model (State/Surf) by a flight length of L. Afterwards, Surf is normal 
 *  to Z axis (beam axis). 
 */
void StraightLineTrackModel::Extrapolate(HepMatrix& State, ReferenceFrame& Surf, double Length)
{
  
  // Get trafo from global coord to reference surface coord  
  HepMatrix& Rot = Surf.GetRotRef(); 
  HepVector& Origin = Surf.GetPosRef();
  
  // Compute 3D representation of track: r(l) = xPoint + l*Dir
  
  // Intersection point on Surf 
  HepVector xPoint(3,0);
  xPoint[0] = State[2][0];
  xPoint[1] = State[3][0];
  xPoint[2] = 0;
  
  // Track direction on Surf
  HepVector Dir(3,0);
  Dir[0] = State[0][0];
  Dir[1] = State[1][0];
  Dir[2] = 1;
  Dir /= Dir.norm(); 
  
  // Endpoint of extrapolation step in global coord.   
  HepVector gPoint = Rot.T()*(xPoint+Length*Dir) + Origin ;
  
  // Track direction in global coord. 
  HepVector gDir = Rot.T()*Dir;
  
  // Overwrite Surf  
  Surf.SetPosition(gPoint);              
  Surf.SetRotation(HepMatrix(3, 3, 1));  // Surface normal to Z axis 
   
  // Overwrite State
  State[0][0] = gDir[0]/gDir[2]; 
  State[1][0] = gDir[1]/gDir[2];
  State[2][0] = 0;                       // Track crosses origin of Surf    
  State[3][0] = 0; 
  
  return;   
}


/** Compute track derivatives for extrapolation from initial surface Surf to 
 *  final surface fSurf. The expansion is around state p.  
 */
int StraightLineTrackModel::TrackJacobian( const HepMatrix& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf, HepMatrix& J) 
{
  
  // Get coord trafo from Surf to fSurf
  HepMatrix Rot = fSurf.GetRotation()*Surf.GetRotation().T();  
  HepVector Origin = Surf.GetRotation()*(fSurf.GetPosition()-Surf.GetPosition());
  
  // Intersection point on initial surface 
  HepVector xPoint(3,0);
  xPoint[0] = State[2][0];
  xPoint[1] = State[3][0];
  xPoint[2] = 0;
  
  // Trajectory direction on initial surface
  HepVector Dir(3,0);
  Dir[0] = State[0][0];
  Dir[1] = State[1][0];
  Dir[2] = 1;
  
  // Trajectory direction on final surface   
  HepVector fDir = Rot*Dir;
  if (fDir[2] == 0) { // No intersection with final plane 
    return 1;
  }
   
  // Surface normal vector in beam direction 
  HepVector W(3,0); 
  W[2] = 1;  
  
  double SX = dot(Rot*(Origin-xPoint), W) / fDir[2];
  
  // DEFINE TRANSPORTATION MATRIX OF TRAJ. PARAMETRTS FROM 
  // CURRENT SURFACE TO TARGET SURFACE 
   
  HepMatrix j(5,5);
  
  // For details of derivation, please refer to V. Karimäki "Straight Line 
  // Fit for Pixel and Strip Detectors with Arbitrary Plane Orientations", CMS Note.
  // Parameters at target surface:   U', V', U, V, Q/P
  // Parameters at current surface:  u', v', u, v, q/p  
  
  // Some abbreviations to save computation time
  
  double Up = fDir[0]/fDir[2];
  double Vp = fDir[1]/fDir[2];     
  
  // *** dU'/du'
  j[0][0] = (Rot[0][0]*fDir[2]-fDir[0]*Rot[2][0]) / (fDir[2]*fDir[2]) ; 
  // *** dU'/dv'
  j[0][1] = (Rot[0][1]*fDir[2]-fDir[0]*Rot[2][1]) / (fDir[2]*fDir[2]) ; 
  // *** dU'/du
  j[0][2] = 0;
  // *** dU'/dv
  j[0][3] = 0;
  
  // *** dV'/du' 
  j[1][0] = (Rot[1][0]*fDir[2]-fDir[1]*Rot[2][0]) / (fDir[2]*fDir[2]) ; 
  // *** dV'/dv' 
  j[1][1] = (Rot[1][1]*fDir[2]-fDir[1]*Rot[2][1]) / (fDir[2]*fDir[2]) ; 
  // *** dV'/du
  j[1][2] = 0;
  // *** dV'/dV
  j[1][3] = 0;
    
  // *** dU/du
  j[2][2] = Rot[0][0] - Rot[2][0] * Up;
  // *** dU/dv
  j[2][3] = Rot[0][1] - Rot[2][1] * Up;
  // *** dU/du'
  j[2][0] = SX * j[2][2];
  // *** dU/dv'
  j[2][1] = SX * j[2][3];
  
  // *** dV/du
  j[3][2] = Rot[1][0] - Rot[2][0] * Vp;
  // *** dV/dv
  j[3][3] = Rot[1][1] - Rot[2][1] * Vp;
  // *** dV/du'
  j[3][0] = SX * j[3][2];
  // *** dV/dv'
  j[3][1] = SX * j[3][3];

  // *** d(Q/P)/d(q/p)
  j[4][4] = 1; 


  // The final track jacobian 
  J = j; 
    
  // Everything ok
  return 0;
}


/** Numerical computation of track derivatives for extrapolation from Surf to fSurf. 
 *  Linearization point is State at Surf.  
 */

/*

int StraightLineTrackModel::TrackJacobian( const HepMatrix& State, const ReferenceFrame& Surf, const ReferenceFrame& fSurf, HepMatrix& J)
{
 
  int ndim = 5; 
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

*/

/** Get local scatter gain matrix
 *  It calculates the derivates of track parameters State on 
 *  scattering angles theta1 and theta2. 
 */
void StraightLineTrackModel::GetScatterGain(const HepMatrix& State, CLHEP::HepMatrix& G)
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

