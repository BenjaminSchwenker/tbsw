/////////////////////////////////////////////////////////  //
//                                                         //
//    NIUnpacker - Marlin Processor                        //
/////////////////////////////////////////////////////////  //

#ifndef NIUnpacker_H
#define NIUnpacker_H 1

// Include LCIO classes
#include <lcio.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerRawDataImpl.h>

// Include Marlin classes
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include basic C
#include <vector>
#include <string>

namespace eudaqinput {
  
  
  class NIUnpacker : public marlin::Processor {
    typedef std::vector<unsigned char> datavect;
    typedef std::vector<unsigned char>::const_iterator datait;
    
   public:
   
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new NIUnpacker ; }
    
    //!Constructor - set processor description and register processor parameters
    NIUnpacker();
    
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
    
    //! Method to unpack source (raw data) -> result (digits)
    bool UnpackRawCollection(lcio::LCCollectionVec * result, lcio::LCCollectionVec * source);

    //! Method to decode raw frame 
    void DecodeFrame(lcio::TrackerDataImpl* zsFrame, size_t len, datait it, int frame, unsigned pivotpixel);
         
    //!Method printing processor parameters
    void printProcessorParams() const;
      
    // PROCESSOR PARAMETERS
     
    //! Output collection name
    std::string _outputCollectionName;
    
    //! Input data collection name  
    std::string _inputCollectionName;
    	        
   private: 
     
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
   
  }; // Class

} // Namespace

#endif 
