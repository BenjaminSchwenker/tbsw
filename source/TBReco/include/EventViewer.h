// ////////////////////////////////////////////////////////////////////////// //
//                                                                            //
//    EventViewer - Marlin Processor                                          //
// ////////////////////////////////////////////////////////////////////////// //


#ifndef EventViewer_H
#define EventViewer_H 1


// Include basic C
#include <vector>

// Include LCIO classes
#include <IMPL/LCCollectionVec.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include ROOT classes
#include <TFile.h>
#include <TH2D.h>

#include "TBDetector.h"


namespace depfet {

//! EventViewer Processor
/*! This is a simple processor to select interesting events in pixel sensors
 *  and dump event data - full frame or zero suppressed - into a TH2 histogram. 
 *  The processor outputs a TFile containing all histograms for browsing.  
 * 
 *  The main use case of this processor is to find and discplay interesting 
 *  events with large clusters from delta electrons or ...
 * 
 *  At the moment, following formats are supported
 * 
 *  A) FULL RAW FRAME:  NxM matrix of raw  signals 
 *  B) FULL CORR FRAME: NXM matrix of corrected signals 
 *  C) ZERO SUPPRESSED: Zero suppressed pixel signals  
 *  
 *  Author: B.Schwenker, Universität Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

class EventViewer : public marlin::Processor {

 public:

//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new EventViewer ; }
   
//!Constructor - set processor description and register processor parameters
   EventViewer();
   
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
   
    
//!Calculates external trigger value 
   bool checkExternalTrigger( LCEvent * evt );
   
//!Method dumping data event frames to 2D histos
   void dumpDataEvent( LCEvent * evt );
   
//!Method dumping raw data event frames to 2D histos
   void dumpRawDataEvent(  LCEvent * evt );
   
//!Method dumping MIMOSA26 data event frames to 2D histos
   void dumpZeroSuppEvent(  LCEvent * evt );
   
//!Method printing processor parameters
   void printProcessorParams() const;
   
//! Input collection names
   std::string _RawDataCollectionName;              
   std::string _DataCollectionName;    
   
// Processor Parameters 
   
   std::string _triggerFileName;    //!< dump all events listed here 
   std::string _rootFileName;       //!< store 2D histos here 
   std::string _selectDataType;     //!< select data type: DATA or RAWDATA or ZEROSUPP
   std::string _selectTrigger; 	    //!< "EXTERNAL" or "ALL"
   int _maxTriggers;                //!< maximum number triggers 
   bool _useStatusMap; 	            //!< if yes: use status map for masking	 	 
   int _displaySensorID;            //!< select modules to be displayed 
     
 private:
   
   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
   int    _triggerCounter;
      
   // List of external triggers
   std::vector<int> _triggers;   
   
   // ROOT_OUTPUT variables
   TFile * _rootFile;

   //! Handle to detector data 
   TBDetector _detector;   

   
}; // Class

} // Namespace

#endif 



