
/////////////////////////////////////////////////////////  //
//                                                         //
//    DEPFETUnpacker - Marlin Processor             //
/////////////////////////////////////////////////////////  //

#ifndef DEPFETUnpacker_H
#define DEPFETUnpacker_H 1

// Include LCIO classes
#include <lcio.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/LCGenericObjectImpl.h>

// Include Marlin classes
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include basic C
#include <vector>
#include <string>

#include "depfetInterpreter.h"

namespace eudaqinput {
  
  
  class DEPFETUnpacker : public marlin::Processor {
    typedef std::vector<unsigned char> datavect;
    typedef std::vector<unsigned char>::const_iterator datait;
    
   public:
   
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new DEPFETUnpacker ; }
    
    //!Constructor - set processor description and register processor parameters
    DEPFETUnpacker();
    
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
      
    // PROCESSOR PARAMETERS
     
    //! Output collection name
    std::string _outputCollectionName;
    
    //! Input data collection name  
    std::string _inputCollectionName;

    //! Reset the sensorID to this number
    int _resetSensorID; 
    	        
   private: 
    DepfetInterpreter interpreter;
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
    int _modID;
    int _dheID;
    int _dhpID;
    Mapping mapping;
    std::string _mappingString;
    bool _swapAxes;

   
  }; // Class

} // Namespace

#endif 
