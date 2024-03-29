/////////////////////////////////////////////////////////  //
//                                                         //
//    PixelClusterizer - Marlin Processor                  //
/////////////////////////////////////////////////////////  //

#ifndef PixelClusterizer_H
#define PixelClusterizer_H 1



// Include LCIO classes
#include <lcio.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include basic C
#include <vector>
#include <string>


// Include ROOT classes
#include <TH2F.h>

namespace depfet {

  typedef std::vector<FloatVec> Pix_GroupVector;
  
  class PixelClusterizer : public marlin::Processor {
   
   public:
    
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new PixelClusterizer ; }

    //!Constructor - set processor description and register processor parameters
    PixelClusterizer();
    
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
      
    // Called by the processEvent() once for the clustering
    // of sparsified pixels. 
    void clusterize( LCEvent * evt , LCCollectionVec * clusterCollection );  
    
    // This method is called inside the clusterize() method in order to 
    // determine if the pixel cell at address (iU,iV) should be added 
    // to the pixel group passed as first argument.  
    bool areNeighbours( FloatVec &group, int iU, int iV, int planeNumber, float time ); 
       
    // This method is called inside the clusterize() method in order to 
    // determine if the pixel cell with address iU/iV is already part 
    // of pixel group passed as first argument.  
    bool isDuplicated( FloatVec &group, int iU, int iV ) ;
    
    // Checks if any other pixel group (apart from base group) neighbours (iU,iV). 
    // If so, merge with base group.   
    void checkForMerge( int iU, int iV, int planeNumber, float time,
           Pix_GroupVector::iterator baseGroup,
           Pix_GroupVector::iterator lastGroup); 
           
    //!Method printing processor parameters
    void printProcessorParams() const;
      
    // PROCESSOR PARAMETERS
    
    //! Input sparsified data collection name
    std::string _sparseDataCollectionName;
    
    //! Output hit cluster collection name
    std::string _clusterCollectionName;
    
    //! Output original pixel data collection name 
    std::string _dummyCollectionName;
    
    //! Name of noiseDB file with masks
    std::string  _noiseDBFileName;
    
    //! Minimimum signal for zero suppression
    float _sparseZSCut;
    
    //! Minimimum signal for seed pixel in clusters
    float _sparseSeedCut;
    
    //! Minimum signal for total cluster charge
    float _sparseClusterCut;

    //! Max absolute time difference between neighboring digits
    float _absoluteNeighborTimeCut;
    
    //! m_acceptDiagonalClusters: 
    //!   = 0: Add pixels which have a side in common with a pixel cell 
    //!        in the list
    //!   = 1: A common corner suffices
    //!    = 2: Max distance is a missing pixel in a iV or a iU
    //!    = 3: Max distance is a missing diagonal pixel 
   int m_acceptDiagonalClusters; 
   	       
  private: 
    
   
   
   //! internally used as storage for input decoding
   UTIL::BitField64 _inputDecodeHelper;
   CellIDEncodeConstructHelper _orginalOutputEncoderHelper;
   CellIDEncodeConstructHelper _clusterOutputEncoderHelper;
   // Pixel mask to filter brocken (hot) channels  
   // Key is sensorID 
   std::map< int, TH2F *> _DB_Map_Mask;
    
   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
   int _hitElements; //!< Integers per hit
   
}; // Class

} // Namespace

#endif 

