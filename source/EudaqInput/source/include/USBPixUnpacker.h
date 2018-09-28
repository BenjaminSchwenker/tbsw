/////////////////////////////////////////////////////////  //
//                                                         //
//    USBPixUnpacker - Marlin Processor                    //
/////////////////////////////////////////////////////////  //

#ifndef USBPixUnpacker_H
#define USBPixUnpacker_H 1

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

#include <bitset>
#include <algorithm>

namespace eudaqinput {
  
  
  class USBPixUnpacker : public marlin::Processor {
    typedef std::vector<unsigned char> datavect;
    typedef std::vector<unsigned char>::const_iterator datait;
    
   public:
    
    struct APIXPix {
	  int x, y, tot, lv1, channel;
	  APIXPix(int x, int y, int tot, int lv1, int channel): 
        x(x), y(y), tot(tot), lv1(lv1), channel(channel) {};
    };
    
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new USBPixUnpacker ; }
    
    //!Constructor - set processor description and register processor parameters
    USBPixUnpacker();
    
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
    
    template<typename T>
    T selectBits(T val, int offset, int length) {
	  return (val >> offset) & ((1ull << length) - 1);
    }
    
    //! Get a vector with all channels (sensors) which are contained in data 
    std::vector<int> getChannels(std::vector<unsigned char> const & data);
    
    //! Method to unpack source (raw data) -> result (digits)
    bool UnpackRawCollection(lcio::LCCollectionVec * result, lcio::LCCollectionVec * source);
    
    //! Method to decode raw frame 
    std::vector<APIXPix> decodeFEI4Data(std::vector<unsigned char> const & data);
         
    //!Method printing processor parameters
    void printProcessorParams() const;
      
    // PROCESSOR PARAMETERS
     
    //! Output collection name
    std::string _outputCollectionName;
    
    //! Input data collection name  
    std::string _inputCollectionName;
    
    //! Unpack data from this channel
    int  _filterChannel;
    	        
   private: 
     
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
    
    
   
  }; // Class

} // Namespace

#endif 
