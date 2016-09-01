#ifndef TBTRACKSTATE_H
#define TBTRACKSTATE_H 1

// CLHEP includes
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/SymMatrix.h>

namespace depfet { 


/** Class TBTrackState
 *  
 *  The class TBTrackState represents a local track state on surface. 
 *  Without magnetic field, a particle trajectory can locally be 
 *  represented as a straight line. We use the following 4 variables:
 *  u', v', u, v
 *  where u'=du/dw and v'=dv/dw are direction tangents wrt. the w=0 
 *  plane in a local UVW coordinate system.
 *  The reference frame is uniquely determined by a reference plane 
 *  number.
 *  
 *  @Author B. Schwenker, University of Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

class TBTrackState {
   
 public: // members 
  
  // Parameter vector (4x1 matrix)
  CLHEP::HepMatrix Pars;  
  // Parameter covariance (4x4 matrix)
  CLHEP::HepSymMatrix Cov;
  // Reference plane number
  int Plane; 
    
  
 public: // functions 
   
  // Constructors
  TBTrackState(); 
  TBTrackState(CLHEP::HepMatrix aPars, CLHEP::HepSymMatrix aCov=CLHEP::HepSymMatrix(4,0));  
  
  // Dimension of track state 
  int GetDim() { return 5;};
  
  // Get/Set track parameters  
  void SetPars(const CLHEP::HepMatrix& aPars) { Pars= aPars; }; 
  CLHEP::HepMatrix&  GetPars() { return Pars; };
   
  // Get/Set parameter covariance 
  void SetCov(const CLHEP::HepSymMatrix& aCov ) { Cov = aCov; }; 
  CLHEP::HepSymMatrix&  GetCov() { return Cov; };

  // Get/Set plane number 
  void SetPlane(int aPlane ) { Plane = aPlane; }; 
  int  GetPlane() { return Plane; };
  
  // Get intersection coordinates 
  CLHEP::HepMatrix  GetXCoord();
  
};

} // Namespace

#endif 

