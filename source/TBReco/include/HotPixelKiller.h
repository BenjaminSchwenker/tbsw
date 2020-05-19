// ///////////////////////////////////////////////////////////////////////////////////////     //
//                                                                                             //
//  HotPixelKiller   - Marlin Processor                                                        //
// ///////////////////////////////////////////////////////////////////////////////////////     //

#ifndef HotPixelKiller_H
#define HotPixelKiller_H 1


// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include <marlin/ProcessorMgr.h>
#include <marlin/Exceptions.h>

#include <UTIL/CellIDDecoder.h>

// Include basic C
#include <string>
#include <map>
#include <vector>



namespace depfet {

  //! HotPixelKiller Processor 
  /*! The processor generates a hot/dead pixel mask. This mask is used for clustering 
   *  to exclude hot/dead pixels from entering into clusters. 
   * 
   *  The HotPixelKiller needs an input collection containing zero suppressed data or 
   *  digits. The user can specify a maximum pixel occupancy (hits/events) for "good"
   *  pixels. For The user can specify an offline zero suppression threshold.
   * 
   *  The pixel occupancy is shaped by the smooth beam profile. In order to remove the smooth
   *  beam profile we can normalize and mask only pixels having far less or far more hits 
   *  than the median hits in 5x5 neighbor pixels. A normal working pixel has normed occupancy
   *  close to 1.  
   *   
   *  The processor also checks if zero suppressed pixels are duplicated and have 
   *  a valid pair of u/v cellIDs.  
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
    
    //! Method for computing the median of std::vector<double>
    double get_median(std::vector<double>& v) const;
    
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
    /*! This is the maximum allowed pixel occupancy (hits/events) for good pixels.
     *  If a pixel fires more frequently, it is masked as HOT. 
     */
    float _maxOccupancy;
    
    //! Minimum pixel occupancy
    /*! This is the minimum allowed pixel occupancy (hits/events) for good pixels.
     *  If a pixel fires less frequently, it is masked as DEAD. 
     */
    float _minOccupancy;

    //! Maximum normed pixel occupancy
    /*! This is the maximum allowed normed pixel occupancy for good pixels.
     *  If a pixel fires N times more frequently than median of its neihbors, it is masked as HOT. 
     */
    float _maxNormedOccupancy;
    
    //! Minimum normed pixel occupancy
    /*! This is the minimum allowed normed pixel occupancy for good pixels.
     *  If a pixel fires N times less frequently than median of its neighbors, it is masked as DEAD. 
     */
    float _minNormedOccupancy;

    //! Masking on normalized occupancy
    /*! Normalize the occupancy before applying min/max cuts for pixel masking.  
     *  A Normal working pixel has a normalized occupancy close to 1. 
     */
    bool _maskNormalized;                 

   private:
    
    //! internally used as storage for input decoding
    UTIL::BitField64 _inputDecodeHelper;
     
    // Count  hits for all pixels 
    std::map < int, FloatVec > _hitCounterMap;
    
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
       
  }; // Class

} // Namespace

#endif 


