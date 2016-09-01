#ifndef SeedGenerator_H
#define SeedGenerator_H 1

// Include DEPFETTrackTools header files
#include "TBDetector.h"
#include "TBHit.h"
#include "TBTrackState.h"

namespace depfet { 


/** Class SeedGenerator 
 * 
 *  The SeedGenarator builds an initial reference track state 
 *  from measured hits. This reference state should be used 
 *  to initialize and linearize a more advanced track fitter. 
 *  
 *  The SeedGenerator is very simple:  
 *
 *  - For a single hit, the reference state points along 
 *  the z-axis and crosses the given hit. 
 *  
 *  - For a hit pair, the reference state crosses both 
 *  hits on two different detectors.
 *
 *  - In case no specific data is available, the reference 
 *  state will follow the z axis (-> beam axis).  
 *    
 *  Author: Benjamin Schwenker, Göttingen University 
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */
  	
class SeedGenerator {
    
 public:
   
  // Constructor
  SeedGenerator( double acharge, double amom );

  // Create a seed track state 
  TBTrackState CreateSeedTrack(TBHit Hit, TBDetector& Detector);   
    
  // Create a seed track state 
  TBTrackState CreateSeedTrack(TBHit FirstHit, TBHit SecondHit, TBDetector& Detector);   
  
  // Create a seed track state  
  TBTrackState CreateSeedTrack();

 private:

  double charge; 
  double mom; 
  
};


} // Namespace

#endif 

