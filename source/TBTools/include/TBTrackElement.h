#ifndef TBTrackElement_H
#define TBTrackElement_H 1

// TBTools includes
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
 *  A TBTrackElement can be assigned one matched TBHit. The hit can be removed 
 *  if it proves incompatible to the track fit. The function HasHit() should 
 *  be queried to see if a hit is assigned.  
 *  
 *  A TBTrackState object. The track state provides a local track parameters 
 *  and tracl paremater covariance matrix. 
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
  ~TBTrackElement(){;}
  
  /* Get detector 
   */
  Det& GetDet() { return m_det; } 

  /* Get detector 
   */
  const Det& GetDet() const { return m_det; } 
  
  /** True if hit is set  
   */
  bool HasHit() const { return m_hasHit; }
  
  /** Set measured hit  
   */
  void SetHit(const TBHit& aHit) 
  { 
    m_hit = aHit;
    m_hasHit = true; 
  }
  
  /** Get hit - query HasHit() before 
   */
  TBHit& GetHit() {return m_hit;}

  /** Get hit - query HasHit() before 
   */
  const TBHit& GetHit() const {return m_hit;}
  
  /** Remove (bad) hit - Outlier rejection  
   */
  void RemoveHit() { m_hasHit=false;  }
  
  /** Get track state - query IsCrossed() before 
   */
  const TBTrackState& GetState() const {return m_state;}

  /** Get track state - query IsCrossed() before 
   */
  TBTrackState& GetState() {return m_state;}
  
  /** Set track state
   */
  void SetState(const TBTrackState& aState) {m_state = aState;}
  
  /** Get sensor crossed flag
   */
  bool IsCrossed() const {return m_isCrossed; }
  
  /** Set sensor crossed flag 
   */
  void SetCrossed(bool Flag) { m_isCrossed = Flag; }
  
  /** Get ChiSqu 
   */
  double GetChiSqu() const {return m_localChiSqu; }
  
  /** Set ChiSqu 
   */
  void SetChiSqu(double Chi2) { m_localChiSqu = Chi2; }
   
 private: 
  
  /* Reference to detector instance 
   */ 
  Det& m_det;
  
  /* Hit matched to track  
   */ 
  TBHit m_hit;   
     
  /* Track state on surface
   */
  TBTrackState m_state; 
  
  /* Flags detector crossed by track - defaults to false
   */ 
  bool m_isCrossed{false}; 
  
  /* Flags detector has a matched hit - defaults to false
   */ 
  bool m_hasHit{false}; 
  
  /* Local ChiSqu - defaults to 0
   */  
  double m_localChiSqu{0};   
};

} // Namespace

#endif
