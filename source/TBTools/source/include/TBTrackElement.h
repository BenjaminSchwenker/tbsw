#ifndef TBTrackElement_H
#define TBTrackElement_H 1

// DEPFETTrackTools includes
#include "Det.h"
#include "TBHit.h"
#include "TBTrackState.h"

namespace depfet { 

/** Class TBTrackElement 
 *  
 *  A TBTrackElement (TE) describes the intersection of a charged particle
 *  with a position sensitive detector. The current implementation focuses 
 *  on a (generic) pixel module providind 2D measurments of the intersection
 *  coordinates in the sensor plane. The TE holds information about
 *  
 *  A) Sensitive detector
 *  B) Measured hit  
 *  C) Local track state 
 *  
 *  A Det object represents a position sensitive detector. It provides all 
 *  information about material budget and sensor boundaries, position and 
 *  orientation. 
 *  
 *  A TBTrackElement can be assigned one TBHit object. The hit can be removed 
 *  if it proves incompatible to the track fit. The function HasHit() should 
 *  be queried to see if a hit is assigned.  
 *  
 *  The TBTrackElement owns a TBTrackState object. The track state provides a 
 *  local representation of the track state wrt. to the detector plane. Please
 *  call function IsFit() to see if state is valid. 
 *   
 *  @Author B. Schwenker, University of Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */
 	

class TBTrackElement 
{
  
 public:  
  
  /* Constructor
   */
  TBTrackElement(Det& aDetUnit);
    
  /* Destructor 
   */
  ~TBTrackElement(){;}; 
  
  /* Get detector 
   */
  Det& GetDet() { return DetUnit; } ;
  
  /** True if hit is set  
   */
  bool HasHit(); 
  
  /** Set measured hit  
   */
  void SetHit(TBHit& aHit); 
  
  /** Get hit - query HasHit() before 
   */
  TBHit& GetHit();
  
  /** Remove (bad) hit - Outlier rejection  
   */
  void RemoveHit();
  
  /** Get track state - query IsCrossed() before 
   */
  TBTrackState& GetState() {return State;}; 

  /** Set track state
   */
  void SetState(TBTrackState& aState); 
  
  /** Get fit flag 
   */
  bool IsCrossed() {return CrossedFlag; }; 
  
  /** Set fit flag 
   */
  void SetCrossed(bool Flag) { CrossedFlag = Flag; }; 
  
  /** Get ChiSqu 
   */
  double GetChiSqu() {return LocalChiSqu; }; 
  
  /** Set ChiSqu 
   */
  void SetChiSqu(double Chi2) { LocalChiSqu = Chi2; }; 
   
 private: 
  
  /* Detector unit   
   */ 
  Det DetUnit;
  
  /* Container for registered hits 
   */ 
  std::vector<TBHit> HitStore;   
     
  /* Track state on surface
   */
  TBTrackState State; 
  
  /* Flags detector crossed by track  
   */ 
  bool CrossedFlag; 
  
  /* Local ChiSqu - consistency of hit assignment 
   */  
  double LocalChiSqu;   
};

} // Namespace

#endif
