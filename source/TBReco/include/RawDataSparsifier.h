// //////////////////////////////////////////////////////////////////////////  //
//                                                                             //
//    RawDataSparsifier - Marlin Processor                                     //
// //////////////////////////////////////////////////////////////////////////  //


#ifndef RawDataSparsifier_H
#define RawDataSparsifier_H 1

// Include basic C
#include <vector>


// Include LCIO classes

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"




namespace depfet {

//! RawDataSparsifier processor
  /*! This processor performs zero supression on full frame data. 
   * 
   *  The input data is organized in a collection of type TrackerData  
   *  named 'data'. Each element represents data from ons sensor for 
   *  one trigger.
   *  
   *  The output of the processor is a collection of sparsified or 
   *  zero suppressed digits.
   *  
   *  Author: Benjamin Schwenker, Göttingen University 
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de> 
   *  
   */

class RawDataSparsifier : public marlin::Processor {

 public:

//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new RawDataSparsifier ; }

//!Constructor - set processor description and register processor parameters
   RawDataSparsifier();

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

//! Full frame data collection name
   std::string _DataCollectionName;   

//! Noise collection name.
   std::string _NoiseCollectionName;

//! Status collection name. 
   std::string _StatusCollectionName;   

//! Sparsified data collection name
   std::string _SparseDataCollectionName;

//! Signal threshold   
   std::vector<float> _ThresholdVec;
   
//! Use SNR cut in zero suppression
   bool _useSNRCut;

 private:

   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
   
}; // Class

} // Namespace

#endif 



