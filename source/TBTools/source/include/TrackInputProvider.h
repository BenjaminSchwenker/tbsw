#ifndef _TrackInputProvider_h
#define _TrackInputProvider_h
  	
// Include LCIO classes
#include <IMPL/TrackImpl.h>

// DEPFETTrackTools includes
#include "TBTrack.h"
#include "TBDetector.h"
 
/** Class TrackInputProvider
 *   
 *  Take a TBTrack and convert it to an lcio::Track, or vice versa.
 *  The lcio track is a persistency class used to read/write tracks into 
 *  lcio files. 
 *    
 *  Note: After creating a TBTrack from an lcio::TrackImpl, a refit is required. 
 *  
 *  @Author B. Schwenker, University of Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

namespace depfet {
   	
class TrackInputProvider {
    
 public:
  
  /** Build TBTrack 
   * 
   * The lcio::Track class is used to store reco tracks persistently in 
   * LCIO files.  
   *   
   */
  TBTrack MakeTBTrack( lcio::Track * lciotrk, TBDetector& detector );
  	
  /** Build lcio::Track 
   * 
   * The lcio::Track class is used to store reco tracks persistently in 
   * LCIO files.  
   *   
   */
  lcio::TrackImpl * MakeLCIOTrack( TBTrack& trk );
  
};

} // Namespace
  
#endif
