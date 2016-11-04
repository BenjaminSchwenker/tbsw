/////////////////////////////////////////////////////////  //
//                                                         //
//    TelUnpackerizer - Marlin Processor                   //
/////////////////////////////////////////////////////////  //

#ifndef TelUnpacke_H
#define TelUnpacker_H 1

#include "TBDetector.h"

// Include LCIO classes
#include <lcio.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerRawDataImpl.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include basic C
#include <vector>
#include <string>

namespace depfet {


class TelUnpacker : public marlin::Processor {

 public:

//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new TelUnpacker ; }

//!Constructor - set processor description and register processor parameters
   TelUnpacker();

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
      
// PROCESSOR PARAMETERS
   
//! Output collection name
   std::string _outputCollectionName;
    
//! Input data collection name  
   std::string _inputCollectionName;

//! Number of floats per digit 
   int m_modulus;

//! Use digits from these detectors
   std::vector<int >  _filterIDs;
 	       
 private: 
    
   // Handle to detector data sheets 
   TBDetector _detector;    
    
   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
   
}; // Class

} // Namespace

#endif 

