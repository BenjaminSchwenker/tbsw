#ifndef GenericTrackModel_H
#define GenericTrackModel_H 1

// DEPFETTrackTools includes
#include "ReferenceFrame.h"

// CLHEP includes 
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/SymMatrix.h>

namespace depfet {

 
//! Class GenericTrackModel
/* 
 * The GenericTrackModel class defines an interface for the dynamics
 * of track propagation. The GenericTrackmodel is designed to play 
 * along with the GenericTrackFitter class. 
 * 
 * There are basically two different cases: A straight line model for 
 * zero magnetic field and a helix model for a constant magnetic field. 
 *  
 * @Author B. Schwenker, University of Göttingen
 * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */
  

class GenericTrackModel {
  
 public:
  
  double m_kappa; 

 /* Default Constructor  
  */
  GenericTrackModel() {;}
 
 //! Default destructor
  virtual ~GenericTrackModel() {;}

  /** Track propagation in magnetic field
   */
  virtual bool BfieldON() = 0; 
 
 /* Returns signed fligth length (mm) to surface fSurf. Track starts at surface Surf and has state State.
  */
  virtual double GetSignedStepLength(const CLHEP::HepMatrix& State, const depfet::ReferenceFrame& Surf, const depfet::ReferenceFrame& fSurf) = 0;
  
 /* Returns true if track hits the surface fSurf. Track starts at surface Surf and has state State.
  */
  virtual bool CheckHitsSurface(const CLHEP::HepMatrix& State, const depfet::ReferenceFrame& Surf, const depfet::ReferenceFrame& fSurf) = 0;

 /** Returns track state at surface fSurf. 
  */
  virtual CLHEP::HepMatrix Extrapolate(const CLHEP::HepMatrix& State, const depfet::ReferenceFrame& Surf, const depfet::ReferenceFrame& fSurf, bool& error) =0;
 
 /** Extrapolate track along helix for given flight length. Track parameters (State/Surf) are overwritten. 
 */
  virtual void Extrapolate(CLHEP::HepMatrix& State, depfet::ReferenceFrame& Surf,  double length) = 0;
  
 /** Compute track derivatives for extrapolation from Surf to fSurf.
  *  Linearization point is State at Surf.  
  */
  virtual int TrackJacobian( const CLHEP::HepMatrix& State, const depfet::ReferenceFrame& Surf, const depfet::ReferenceFrame& fSurf ,  CLHEP::HepMatrix& J) =0;
  
 /** Get local scatter gain matrix
  *  It calculates the derivates of track parameters State on 
  *  scattering angles theta1 and theta2. 
  */
 virtual void GetScatterGain(const CLHEP::HepMatrix& State, CLHEP::HepMatrix& G) =0;
 
  
};
 
} // Namespace

#endif 
