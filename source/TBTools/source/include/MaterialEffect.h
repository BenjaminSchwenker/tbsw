#ifndef MaterialEffect_H
#define MaterialEffect_H 1


// CLHEP includes 
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/SymMatrix.h>
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"


namespace materialeffect {

  /** Calculate variance of the projected angular deflection due to multiple
  *  scattering of a particle for x/x0 a la Highland
  *   
  *  Projected means that the angle is not the angle in space, but rather a
  *  one-dimensional projection on one axis of a plane that is perpendicular
  *  to the particle direction before scattering.  
  */ 
  double GetScatterTheta2(double x, double x0, double mass, double charge, double mom );  

  /** Scatter track at thin scatterer 
   *
   * The track is locally scattered by two scatter kink angles. The kink angles 
   * are defined as projections int the comoving frame, i.e. relative to the 
   * direction of the unscattered track. The offset through scattering is 
   * neglected (thin scatterer). 
   * 
   * Note: Initially, the track state is before scattering, and gets overwritten 
   * by state after scattering.
   */ 
  void ScatterTrack(CLHEP::HepMatrix& State, double kink_u, double kink_v);  
  
  // Highland radiation length for air [mm] (NTP: 20 deg, 1bar)   
  static const float X0_air = 305000;

} // Namespace


#endif
