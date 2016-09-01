// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    CogHitMaker - Marlin Processor - Compute center of gravity hits from clusters         //
// ///////////////////////////////////////////////////////////////////////////////////////  //


#ifndef CogHitMaker_H
#define CogHitMaker_H 1

// Include DEPFETTrackTools 
#include "TBDetector.h"

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include <marlin/Exceptions.h>

// Include basic C
#include <string>



namespace depfet {

 //! CogHitMaker processor
 /*!  This processor uses the center of gravity algorithm to compute hits from 
   *  pixel clusters. The processor takes a cluster collection as input and 
   *  produces a hit collection as output. 
   *  
   *  
   *  Author: Benjamin Schwenker, Göttingen University 
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */

class CogHitMaker : public marlin::Processor {
  
 public:
   
//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new CogHitMaker ; }
   
//!Constructor - set processor description and register processor parameters
   CogHitMaker();
   
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
   
//!Method printing processor parameters
   void printProcessorParams() const;
   
protected:
      
//! Input cluster collection name
   std::string  _clusterCollectionName; 
   
//! Output hit collection name
   std::string  _hitCollectionName;
   
//! Cluster quality selection - use only good clusters
   int _clusterQualitySelect;
   
 private:
   
   // Handle to detector data 
   TBDetector _detector;    
    
   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
  
}; // Class

} // Namespace

#endif 

