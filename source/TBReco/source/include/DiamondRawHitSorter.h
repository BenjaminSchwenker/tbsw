/////////////////////////////////////////////////////////  //
//                                                         //
//    DiamondRawHitSorter - Marlin Processor                   //
/////////////////////////////////////////////////////////  //

#ifndef DiamondRawHitSorter_H
#define DiamondRawHitSorter_H 1


// Include LCIO classes
#include <lcio.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include basic C
#include <vector>
#include <string>

namespace depfet {


class DiamondRawHitSorter : public marlin::Processor {

 public:

   //!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new DiamondRawHitSorter ; }

   //!Constructor - set processor description and register processor parameters
   DiamondRawHitSorter();

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
   
   int getPixelType(int x, int y);
   std::vector<int> getPixelCoordinates(int x, int y);
   std::vector<int> getPixelCoordinates(int x, int y, int type);
   std::vector<int> getRectCoord(int x, int y);
   std::vector<int> getHexCoord(int x, int y);
      
   // PROCESSOR PARAMETERS
   
   //! Output collection name
   std::string _outputCollectionName;
    
   //! Input data collection name  
   std::string _inputCollectionName;
   
   //! Select type of diamond pixel 
   int m_type;
          
 private: 
    
   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
   
   int layoutStartX = 65;
   int layoutStartY = 0;
   int layoutEndX = 79;
   int layoutEndY = 80;
   int metalStartX = 68;
   int metalEndX = 79;
   int metalStartY = 0;
   int metalEndY = 64;
   int hexStartY = 31;
    		
   int rectPeriod = 5;
   int hexPeriod = 4;
   
}; // Class

} // Namespace

#endif 

