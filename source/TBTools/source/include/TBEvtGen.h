#ifndef TBEvtGen_H
#define TBEvtGen_H 1

#include "TBTrack.h"

#include <TRandom.h>
#include <TRandom3.h>

namespace depfet { 


/** Class TBEvtGen
 *  
 *  The TBEvtGen simulates the trigger of the eudet 
 *  system. 
 *
 *  The class decides whether a beam particle gets 
 *  recorded on disk or not.
 *   
 *  Author: Benjamin Schwenker, Göttingen University 
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

class TBEvtGen {

 public: 
   
  // Simulation time (sec)
  double time; 
  // Last trigger time (sec)
  double lasttrg;
  // Beam Intensity  (1/sec)
  double intensity; 
  // Integration time (sec)
  double readouttime; 
  // Trigger mode
  int trgmode; 
  
 public:
   
  // Constructor
  TBEvtGen( );

  // Destructor 
  ~TBEvtGen( );
         
  /** ACCEPT
   *  
   * A simulated track is accepted for detector
   * simulation if it generates a trigger or 
   * passes m26 sensors within integration time. 
   */
  bool  ACCEPT(TBTrack& TruthTrack);

 private: 
  
  /** HasCoincidence 
   *
   *  Flags coincident scinti signals to TLU
   */
  bool  HasCoincidence(TBTrack& TruthTrack);

  /** HasHit
   *
   *  Flags hit at layer ipl
   */
  bool  HasHit(TBTrack& TruthTrack, int ipl);
  
  /** Set seed for internal random number generator
   *
   */
  void SetSeed( int seed );
  
  // Internal random number generator (Mersenne Twister is safe choice)
  TRandom3 * myRng;
   
};


} // Namespace

#endif 

