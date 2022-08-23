// ////////////////////////////////////////////////////////////////////////// //
//                                                                            //
//    EventEmbedder - Marlin Processor                                          //
// ////////////////////////////////////////////////////////////////////////// //


#ifndef EventEmbedder_H
#define EventEmbedder_H 1


// Include basic C
#include <vector>

// Include LCIO classes
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

namespace depfet {

//! EventEmbedder Processor
/*! This is a simple processor to add user defined collections from used defined 
 *  events from a different file to the current event.  
 *
 *  Author: B.Schwenker, Universität Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

class EventEmbedder : public marlin::Processor {

 public:

//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new EventEmbedder ; }
   
//!Constructor - set processor description and register processor parameters
   EventEmbedder();
   
//!Method called at the beginning of data processing - used for initialization
   virtual void init();
   
//!Method called for each run - used for run header processing
   virtual void processRunHeader(LCRunHeader * run);
   
//!Method called for each event - used for event data processing
   virtual void processEvent(LCEvent * evt);
   
//!Method called after each event - used for data checking
   virtual void check(LCEvent * evt);
   
//!Method called after all data processing
   virtual void end();
   
protected:
     
//!Method printing processor parameters
   void printProcessorParams() const;

//!Called by the processEvent() 
   void copy_digits(LCEvent * evt , LCCollectionVec * digitCollection );
         
   std::string _DigitCollectionName;       
   std::string _triggerFileName;    
   std::string _otherFileName;       
 
 private:

   //! internally used as storage for input decoding
   UTIL::BitField64 _inputDecodeHelper;
   CellIDEncodeConstructHelper _outputEncoderHelper; 
   
   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
   
   // List of external triggers
   std::vector<int> _triggers;   
   
}; // Class

} // Namespace

#endif 



