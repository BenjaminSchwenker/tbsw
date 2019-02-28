// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    RawHitDQM - Marlin Processor                                                    //
// ///////////////////////////////////////////////////////////////////////////////////////  //


#ifndef RawHitDQM_H
#define RawHitDQM_H 1

// Include DEPFETTrackTools 
#include "TBDetector.h"

// Include LCIO classes
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include <marlin/Exceptions.h>

// Include ROOT classes
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

// Include basic C
#include <string>
#include <vector>
#include <map>


namespace depfet {

 //! RawHitDQM 
 /*!  This processor produces DQM histos from one or more hit collections.
   *   
   *  
   *  Author: Benjamin Schwenker, Göttingen University 
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */

class RawHitDQM : public marlin::Processor {
  
 public:
   
//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new RawHitDQM ; }
   
//!Constructor - set processor description and register processor parameters
   RawHitDQM();
   
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
      
//! Histogram booking
   void bookHistos();

//! Input hit collection names.
   std::vector< std::string >  _inputHitCollectionNameVec; 
       
//! ROOT output file name  
   std::string _rootFileName;  
     
 private:
   
   // Handle to detector data 
   TBDetector _detector;    
   
   // Handle to root file
   TFile * _rootFile;
    
   std::map< std::string, TH1D *> _histoMap;
   std::map< std::string, TH2D *> _histoMap2D; 
   
   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number

}; // Class

} // Namespace

#endif 

