#ifndef TBTrack_H
#define TBTrack_H 1

// TBTools includes
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
 * trajectory. Track Fitters may use the reference track state 
 * (TBTrackState) for pattern recognition (->TrackFinder) and in 
 * the track fitter for linearization. 
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
  double GetCharge() const { return Charge; }
  
  /* Set particle charge in e
   */ 
  void SetCharge( double aCharge) { Charge = aCharge; }
  
  /* Get particle mass in GeV 
   */ 
  double GetMass() const { return Mass; }
  
  /* Set particle mass in GeV 
   */ 
  void SetMass( double aMass) { Mass = aMass; }
  
  /* Set particle abs momentum in GeV 
   */ 
  void SetMomentum(double aMom) { Mom = aMom; }
  
  /* Get particle 3-momentum in GeV 
   */ 
  double GetMomentum() const { return Mom; } 
   
  /** Get track element by plane number 
   */  
  TBTrackElement& GetTE(int ipl);
  
  /** Get track element by plane number 
   */  
  const TBTrackElement& GetTE(int ipl) const;
  
  /** Get all track elements  
   */  
  std::vector<TBTrackElement>& GetTEs();
  
  /** Get all track elements  
   */  
  const std::vector<TBTrackElement>& GetTEs() const;
  
  /** Get number of active hits in track   
  */  
  int GetNumHits() const;
  
  /** Get number of track elements    
  */  
  int GetNumTEs() const;
   
  /** Get Chisqu of track fit 
  */
  double GetChiSqu() const { return ChiSqu; }
  
  /** Set fit chisqu
  */
  void SetChiSqu(double aChiSqu) { ChiSqu = aChiSqu; }
   
  /** Get number degrees of freedom   
   */
  int GetNDF() const {return Ndof; }
  
  /** Set fit chisqu
  */
  void SetNDF(int aNdof) { Ndof = aNdof; }
  
  /** Get reference state 
  */  
  TBTrackState& GetReferenceState()  { return RefState; }
  
  /** Set reference state
  */  
  void SetReferenceState(const TBTrackState& State) {RefState = State;} 
  
 private:
  
  // Reference trajectory
  TBTrackState RefState;   
  
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
