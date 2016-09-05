// ///////////////////////////////////////////////////////////////////////////////////////     //
//                                                                                             //
//    DEPFETHotPixelMasker - Marlin Processor                                                  //
// ///////////////////////////////////////////////////////////////////////////////////////     //

#ifndef DEPFETHotPixelMasker_H
#define DEPFETHotPixelMasker_H 1

#include "TBDetector.h"

// Include LCIO classes
#include <EVENT/LCParameters.h>
#include <IMPL/LCCollectionVec.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include <marlin/ProcessorMgr.h>
#include <marlin/Exceptions.h>

// Include basic C
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>

namespace depfet {

 //! DEPFETHotPixelMasker Processor 
  /*! Remove hits from masked pixels 
   *  
   *  Author: B.Schwenker, Universität Göttingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */


class DEPFETHotPixelMasker : public marlin::Processor {

 public:

//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new DEPFETHotPixelMasker ; }

//!Constructor - set processor description and register processor parameters
   DEPFETHotPixelMasker();

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



  

//! Method printing processor parameters
   void printProcessorParams() const;

// Processor Parameters 

//! Input hit collection names.
/*!
 */
   std::vector< std::string >  _inputHitCollectionNameVec;   


//! Input status data collection name
   std::string _statusCollectionName;

 

 private:
  
  // Handle to detector data sheets 
  TBDetector _detector;     

  
    
  double _timeCPU; //!< CPU time
  int    _nRun ;   //!< Run number
  int    _nEvt ;   //!< Event number
  
 
     
}; // Class

} // Namespace

#endif 



