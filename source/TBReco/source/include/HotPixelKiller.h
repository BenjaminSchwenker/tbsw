// ///////////////////////////////////////////////////////////////////////////////////////     //
//                                                                                             //
//  HotPixelKiller   - Marlin Processor                                                        //
// ///////////////////////////////////////////////////////////////////////////////////////     //

#ifndef HotPixelKiller_H
#define HotPixelKiller_H 1

// TBTools includes 
#include "TBDetector.h"

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
#include <map>
#include <vector>


namespace depfet {

 //! HotPixelKiller Processor 
  /*! The processor generates a hot pixel mask in lcio format. This mask can be used
   *  during the clustering processor to filter hot pixels. Moreover, the 
   *  HotPixelKiller processor produces a couple of DQM histos in a root file.
   * 
   *  The HotPixelKiller needs an input collection containing zero suppressed data.
   *  The user can specify a maximum hit occupancy for "normal" pixels. The user can 
   *  specify an offline zero suppression threshold. 
   *   
   *  The processor also checks if zero suppressed pixels are not duplicated and have 
   *  a valid pair of cellIDs.  
   *  
   *  Author: B.Schwenker, Universität Göttingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */


class HotPixelKiller : public marlin::Processor {

 public:

//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new HotPixelKiller ; }

//!Constructor - set processor description and register processor parameters
   HotPixelKiller();

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

//! Histogram booking
   void bookHistos();

//! Accumulate hits
/*! Accumulate hits to prepare computing of mask later
 */
   void accumulateHits(LCEvent * evt);
  
//! Compute hot pixel mask 
/*! This method is called at the end of the processor to compute the 
 *  hot pixel mask
 */
   void computeMask();

//! This method is used to initialize algorithms  
   void initializeAlgorithms(LCEvent * evt); 

//! Method printing processor parameters
   void printProcessorParams() const;

// Processor Parameters 

//! Input rawdata collection name
   std::string _zeroSuppressedDataCollectionName;   

//! Output atatus collection name 
   std::string _statusDataCollectionName;  

//! Output noise DB file name
   std::string _noiseDBFileName; 

//! Output root file name containing all ntuples 
   std::string _rootFileName;  

//! Offline ZS threshold
/*! Digits having a signal below this threshold are ignored. The same 
 *  threshold should be used later in the clusterizer. 
 */
   float _offlineZSCut;

//! Maximum pixel occupancy
/*! This is the maximum allowed hit rate or occupancy for nornal pixels.
 *  If a pixel fires more frequently, it is masked as HOT. 
 */
   float _maxOccupancy;
    
//! EventsForMask
/*! This is the number of events used for computing the hotpixel mask. 
 */
   int _eventsForMask; 

 private:
  
  // Handle to detector data sheets 
  TBDetector _detector;     
   
  // Modules to be processed 
  std::vector<int> _planeNumbers;  
  int _noOfDetector; 
  
  // Count  hits for all pixels on all  detectors 
  std::vector < FloatVec > _hitCounter;
  
  // Status mask for all pixels on all detectors 
  std::vector < ShortVec > _status;
  
  double _timeCPU; //!< CPU time
  int    _nRun ;   //!< Run number
  int    _nEvt ;   //!< Event number

  //! ROOT_Output
  TFile * _rootFile;
  TTree * _rootEventTree;  

  int _rootEventNumber; 
  int _rootDetectorID;  
  int _rootNHits; 
  int _rootNGoodHits;
   
  std::map< std::string, TH1D *> _histoMap;
  std::map< std::string, TH2D *> _histoMap2D;  

     
}; // Class

} // Namespace

#endif 


