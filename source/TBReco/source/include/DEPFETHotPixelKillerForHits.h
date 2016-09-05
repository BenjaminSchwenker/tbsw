// ///////////////////////////////////////////////////////////////////////////////////////     //
//                                                                                             //
//    DEPFETHotPixelKillerForHits - Marlin Processor                                                  //
// ///////////////////////////////////////////////////////////////////////////////////////     //

#ifndef DEPFETHotPixelKillerForHits_H
#define DEPFETHotPixelKillerForHits_H 1

#include "TBDetector.h"

// Include LCIO classes
#include <EVENT/LCParameters.h>
#include <IMPL/LCCollectionVec.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include <marlin/ProcessorMgr.h>
#include <marlin/Exceptions.h>

// Include ROOT classes
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>

// Include basic C
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>

namespace depfet {

 //! DEPFETHotPixelKillerForHits Processor 
  /*! The processor generates a hot pixel mask in lcio format. This mask can be 
   *  during in the clustering processor to filter hot pixels. Moreover, the 
   *  HotPixelKiller processor generates a root file for easy visualisation of the 
   *  hot pixel mask and the noise occupancy of pixels. 
   * 
   *  The HotPixelKillerForHits processor needs an input collection containing hit data.
   *  The user can specify a maximum hit occupancy for good pixels. The user can 
   *  specify an offline zero suppression threshold. This threshold can be used in 
   *  case the signal amplitude is still available.   
   *   
   *  The processor also checks if zero suppressed pixels are not duplicated and have 
   *  a valid address.  
   *  
   *  In addition, the processor allows some minitoring of pixel signals using 
   *  root ntuples.  
   *  
   *  Author: B.Schwenker, Universität Göttingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */


class DEPFETHotPixelKillerForHits : public marlin::Processor {

 public:

//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new DEPFETHotPixelKillerForHits ; }

//!Constructor - set processor description and register processor parameters
   DEPFETHotPixelKillerForHits();

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


//! Read  event 
/*! Read data and accumulate hits 
 */
   void readEventData(LCEvent * evt);
  
//! Compute hot pixel mask 
/*! This method is called at the end of the processor to compute the 
 *  hot pixel mask
 */
   void computeMask();

//! This method is used to fill root occupancy ntuple   
   void fillOccupancyTuple();

//! This method is used to initialize algorithms  
   void initializeAlgorithms(); 

//! Method printing processor parameters
   void printProcessorParams() const;

// Processor Parameters 

//! Input hit collection names.
/*! A vector containing all the collection names to be used in hit 
 *  correlations.
 */
   std::vector< std::string >  _inputHitCollectionNameVec;   

//! Output noise DB file name
   std::string _noiseDBFileName; 

//! Output root file name containing all ntuples 
   std::string _rootFileName;       

//! Maximum pixel occupancy
/*! This is the maximum allowed frequency (hit occupancy) at which a 
 *  good pixel should fire. If a pixel fires too often, then it means 
 *  that there is something weird with this and it should be better masked out.
 */
   float _maxOccupancy;
    
 private:
  
  // Handle to detector data sheets 
  TBDetector _detector;     
   
  // Count pixel hits (-> mean firing frequency)
  std::vector < FloatVec > _hitCounter;
  
  // Status flags for all pixels in all detectors 
  std::vector < ShortVec > _status;
  
  double _timeCPU; //!< CPU time
  int    _nRun ;   //!< Run number
  int    _nEvt ;   //!< Event number
  
  //! ROOT_Output
  TFile * _rootFile;
  TTree * _rootOccTree; 
  TTree * _rootEventTree;  
  
  int _rootEventNumber; 
  int _rootDetectorID;
  int _rootPlaneNumber;  
  int _rootCol;                
  int _rootRow;      
  int _rootStatus;      
  double _rootHitFrequency; 
  int _rootNHits; 
  int _rootNGoodHits;
     
}; // Class

} // Namespace

#endif 



