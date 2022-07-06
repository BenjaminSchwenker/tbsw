/////////////////////////////////////////////////////////  //
//                                                         //
//    HitsFilterProcessor - Marlin Processor               //
/////////////////////////////////////////////////////////  //

#ifndef HitsFilterProcessor_H
#define HitsFilterProcessor_H 1



// Include LCIO classes
#include <lcio.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>


// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include basic C
#include <vector>
#include <string>

namespace depfet {


class HitsFilterProcessor : public marlin::Processor {

 public:

   //!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new HitsFilterProcessor ; }

   //!Constructor - set processor description and register processor parameters
   HitsFilterProcessor();

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

   //! Only use hits from these these sensorIDs
   std::vector<int >  _filterIDs;
 	       
 private: 
   CellIDEncodeConstructHelper _outputEncoderHelper; 
    
   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
   int _hitElements;
   
}; // Class

} // Namespace

#endif 

