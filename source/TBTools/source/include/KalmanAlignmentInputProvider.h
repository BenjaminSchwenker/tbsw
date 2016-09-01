#ifndef _KalmanAlignmentInputProvider_h
#define _KalmanAlignmentInputProvider_h

// DEPFETTrackTools includes
#include "AlignEvent.h"
#include "TBTrack.h"
#include "TBDetector.h"

namespace depfet { 


/** Class KalmanAlignmentInputProvider
 *  
 *  The class converts TBTrack<->AlignEvent. The AlignEvent class
 *  is used to read/write tracks into root files. 
 *  
 *  Author: Benjamin Schwenker, Göttingen University 
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

class KalmanAlignmentInputProvider {

 public:
  
  // Create TBTrack 
  TBTrack MakeTBTrack( AlignEvent& event, TBDetector& detector );
  
  // Fill AlignEvent 
  void FillEvent(TBTrack& RecoTrack, AlignEvent& event);
};

} // Namespace

#endif
