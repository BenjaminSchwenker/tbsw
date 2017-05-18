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


// Include basic C
#include <string>
#include <map>
#include <vector>



namespace depfet {

  //! HotPixelKiller Processor 
  /*! The processor generates a hot pixel mask. This mask is used for clustering 
   *  to exclude hot pixels from entering into clusters. 
   * 
   *  The HotPixelKiller needs an input collection containing zero suppressed data or 
   *  digits. The user can specify a maximum hit rate for "normal" pixels. For The user
   *  can specify an offline zero suppression threshold. 
   *   
   *  The processor also checks if zero suppressed pixels are duplicated and have 
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
    
    //! Method printing processor parameters
    void printProcessorParams() const;
    
    // Processor Parameters 
    
    //! Input rawdata collection name
    std::string _zeroSuppressedDataCollectionName;   
    
    //! Output noise DB file name
    std::string _noiseDBFileName; 
    
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
    
   private:
    
    // Handle to detector data sheets 
    TBDetector _detector;     
     
    // Count  hits for all pixels 
    std::map < int, FloatVec > _hitCounterMap;
    
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
       
  }; // Class

} // Namespace

#endif 


