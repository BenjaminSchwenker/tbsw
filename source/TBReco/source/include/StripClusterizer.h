/////////////////////////////////////////////////////////  //
//                                                         //
//    StripClusterizer - Marlin Processor                  //
/////////////////////////////////////////////////////////  //

#ifndef StripClusterizer_H
#define StripClusterizer_H 1

#include "TBDetector.h"

// Include LCIO classes
#include <lcio.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include basic C
#include <vector>
#include <string>

namespace depfet {

typedef std::vector<FloatVec> ClusterCandVec;

class StripClusterizer : public marlin::Processor {

 public:

//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new StripClusterizer ; }

//!Constructor - set processor description and register processor parameters
   StripClusterizer();

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
   
// This method is called to initialize the status information,
// namely linking status data to strip data. 
   void initializeStatus( LCEvent * event );
   
// Called by the processEvent() once for the clustering
// of sparsified digits. 
   void clusterize( LCEvent * evt , LCCollectionVec * clusterCollection );  
   
// This method is called inside the clusterize() method in order to 
// determine if a strip at cell should be added 
// to cluster.
   bool areNeighbours( FloatVec &cluster, int cell, int isV, int m_accept ); 
   
// This method is called inside the clusterize() method in order to 
// determine if strip at cell is already part 
// of cluster.
   bool isDuplicated( FloatVec &cluster, int cell, int isV) ;
   
// Checks if any other pixel group (apart from base group) neighbours cell. 
// If so, merge with base group.   
   void checkForMerge( int cell, int isV, 
           ClusterCandVec::iterator baseGroup,
           ClusterCandVec::iterator lastGroup); 
           
//!Method printing processor parameters
   void printProcessorParams() const;
      
// PROCESSOR PARAMETERS
   
//! Input sparsified data collection name
   std::string _sparseDataCollectionName;
    
//! Input status data collection name
   std::string _statusCollectionName;
   
//! Output hit cluster collection name
   std::string _clusterCollectionName;
   
//! Output original pixel data collection name 
   std::string _dummyCollectionName;
   
//! Minimimum signal for zero suppression
   float _sparseZSCut;
   
//! Minimimum signal for seed pixel in clusters
   float _sparseSeedCut;
   
//! Minimum signal for total cluster charge
   float _sparseClusterCut;
   
//! Accept strip clusters with gaps 
   int m_acceptGaps; 

//! Number of samples per hit
   int m_samples; 
   	       
 private: 
    
   // Handle to detector data sheets 
   TBDetector _detector;    
   
   //! Status ready switch
   bool _isStatusReady;
    
   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
   
}; // Class

} // Namespace

#endif 

