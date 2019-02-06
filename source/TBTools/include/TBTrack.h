#ifndef TBTrack_H
#define TBTrack_H 1

// DEPFETTrackTools includes
#include "TBTrackElement.h"
#include "TBDetector.h"
#include "TBTrackState.h" 

// C++ includes
#include <vector>

namespace depfet {

 
//! Class TBTrack
/** 
 * The TBTrack class represents a charged particle trajectory in a 
 * test beam detector. The trajectory is represented as a vector of 
 * intersections with position sensitive detectors ordered by their 
 * position along the beam line. These intersections are represented 
 * as instances of the class TBTrackElement (TE). 
 * 
 * For track fitting, the user must pass particle mass, charge and 
 * momentum to the TBTrack. An independent measurment of the momentum 
 * is not possible without a magnetic field and PID detectors.  
 *    
 * In addition, the TBTrack class holds a straight line reference 
 * trajectory. Track Fitters may use the reference trajectory for 
 * pattern recognition (->TrackFinder) and in the track fitter to 
 * initialize and linearize fitting. The reference trajectory is 
 * parametrized at the global Z=0 surface.
 * 
 * @Author B. Schwenker, University of Göttingen
 * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

class TBTrack {
  
 public: 
  
  /** Constructor    
   */
  TBTrack(TBDetector& detector); 
  
  /* Get particle charge in e
   */ 
  double GetCharge() 
    { return Charge; };
  
  /* Set particle charge in e
   */ 
  void SetCharge( double aCharge) 
    { Charge = aCharge; };
  
  /* Get particle mass in GeV 
   */ 
  double GetMass() 
    { return Mass; };
  
  /* Set particle mass in GeV 
   */ 
  void SetMass( double aMass) 
    { Mass = aMass; };
  
  /* Set particle abs momentum in GeV 
   */ 
  void SetMomentum(double aMom)
    { Mom = aMom; };
  
  /* Get particle 3-momentum in GeV 
   */ 
  double GetMomentum() 
    { return Mom; }; 
   
  /** Get track element by plane number 
   */  
  TBTrackElement& GetTE(int ipl);
  
  /** Get all track elements  
   */  
  std::vector<TBTrackElement>& GetTEs();
  
  /** Get number of active hits in track   
  */  
  int GetNumHits();
  
  /** Get number of track elements    
  */  
  int GetNumTEs();
   
  /** Get Chisqu of track fit 
  */
  double GetChiSqu()  
    { return ChiSqu; }
  
  /** Set fit chisqu
  */
  void SetChiSqu(double aChiSqu) 
    { ChiSqu = aChiSqu; }
   
  /** Get number degrees of freedom   
   */
  int GetNDF() 
    {return Ndof; }
  
  /** Set fit chisqu
  */
  void SetNDF(int aNdof) 
    { Ndof = aNdof; }
  
  /** Get reference state at Z=0
  */  
  TBTrackState& GetReferenceState() { return RefStateZ0; }; 
 
  /** Set reference state at Z=0
  */  
  void SetReferenceState(TBTrackState& State) {RefStateZ0 = State;}; 
  
 private:
  
  // Reference trajectory (at Z=0) 
  TBTrackState RefStateZ0;   

  // Particel hypothesis 
  double Mom; 
  double Charge;
  double Mass; 
   
  // ChiSqu of track fit
  double ChiSqu;
   
  // Number of degrees of freedom
  int Ndof; 

  // Vector of track elements
  std::vector<TBTrackElement> TEVec;     
};
 

} // Namespace

#endif 
