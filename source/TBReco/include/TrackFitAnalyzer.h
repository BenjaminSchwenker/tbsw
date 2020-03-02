// ////////////////////////////////////////////////////////////////////////// //
//                                                                            //
//    TrackFitAnalyzer - Marlin Processor                                          //
// ////////////////////////////////////////////////////////////////////////// //

#ifndef TrackFitAnalyzer_H
#define TrackFitAnalyzer_H 1


// Include ROOT classes
#include <TFile.h>
#include <TTree.h>


// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include basic C
#include <string>
#include <vector>
#include <map>

namespace depfet {

  //! TrackFitAnalyzer processor 
  /*! The Processor produces a tree with local track variables 
   *  and - if matched to the track in the trackfinding  - local 
   *  hit variables for all sensors in the telescope. The tree can 
   *  be used to analyze the resolution, alignment, in pixel charge 
   *  collection etc. It should - however - not be used to infer
   *  the hit efficiency.   
   *    
   *  Author: Benjamin Schwenker, Göttingen University 
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
  
  
  class TrackFitAnalyzer : public marlin::Processor {
  
   public:
     
    //!Method that returns a new instance of this processor
    Processor*  newProcessor() { return new TrackFitAnalyzer ; }
    
    //!Constructor - set processor description and register processor parameters
    TrackFitAnalyzer();
    
    //!Method called at the beginning of data processing - used for initialization
    void init();
    
    //!Method called for each run - used for run header processing
    void processRunHeader(LCRunHeader * run);
    
    //!Method called for each event - used for event data processing
    void processEvent(LCEvent * evt);
    
    //!Method called after each event - used for data checking
    void check(LCEvent * evt);
    
    //!Method called after all data processing
    void end();
    
    //!Method printing processor parameters
    void printProcessorParams() const;
      
   protected:
    
    //! Histogram booking
    void bookHistos();
    
    //! Processor Parameters 
    
    //! Input track collection name
    std::string _inputTrackCollectionName;
      
    //! ROOT output file name  
    std::string _rootFileName;
    
    //! Reference plane number, a track is required to have a hit on the reference plane 
    int _iref; 
    
    //! Ignore data from these sensors
    std::vector<int >  _ignoreIDVec;
         
   private:
   
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
     
    // Handle to root file
    TFile * _rootFile;

     /** One entry per track point */
     TTree * _rootHitTree;
        
     int _rootEventNumber;             // Event number from lcio file
     int _rootRunNumber;               // Run number from lcio file 
     int _rootTrackWithRefHit;         // Track has a hit on the reference plane
     double _rootTrackChi2 ;           // Chi2 value from fit of track
     int _rootTrackNDF;                // Number of degrees of freedom of track fit
     int _rootTrackNHits;              // Number of telescope hits used for track fitting 
     
     int _rootSensorID;                // SensorID from lcio file (this is typically NOT the plane number!!)
     int _rootPlaneNumber;             // Plane number 
     double _rootFitMomentum;          // Estimated track momentum from fit in GeV
     double _rootFitU;                 // Estimated track intersection u coordimate in mm, in local uvw coordinates        
     double _rootFitV;                 // Estimated track intersection v coordimate in mm, in local uvw coordinates                  
     double _rootFitdUdW;              // Estimated track slope du/dw in radians, in local uvw coordinates       
     double _rootFitdVdW;              // Estimated track slope dv/dw in radians, in local uvw coordinates    
     double _rootFitErrorU;            // Estimated 1x sigma uncertainty for track intersection u coordinate
     double _rootFitErrorV;            // Estimated 1x sigma uncertainty for track intersection v coordinate 
     int _rootFitCellU;                // Estimated track intersection u coordinate lies on this u cell      
     int _rootFitCellV;                // Estimated track intersection v coordinate lies on this v cell         
     double _rootFitCellUCenter;       // Central coordinate of cell 'FitCellU' in mm        
     double _rootFitCellVCenter;       // Central coordinate of cell 'FitCellV' in mm 
     int _rootTrackPixelType;          // PixelType of pixel cell intersected by track 
   
     int _rootHasHit;             // Track can be matched to a hit (== 0) 
     int _rootHitQuality;              // GoodCluster == 0, BadCluster != 0
     double _rootHitLocalChi2;         // Chi2 value from hit-track residual 
     double _rootPullResidualU;        // Standardized residual in u direction, should have mean = 0 and rms = 1
     double _rootPullResidualV;        // Standardized residual in v direction, should have mean = 0 and rms = 1       
     double _rootHitU;                 // Hit coordinate u reconstructed from cluster in mm, in local uvw coordinates       
     double _rootHitV;                 // Hit coordinate v reconstructed from cluster in mm, in local uvw coordinates     
     double _rootHitClusterCharge;     // Sum over all charges in the cluster 
     double _rootHitSeedCharge;        // Highest charge in cluster
     int _rootHitSeedPixelType;        // PixelType of seed pixel cell 
     int _rootHitSize;                 // Number of hit cells (pixels/strips) in cluster
     int _rootHitSizeU;                // Number of hit cells along u direction in cluster
     int _rootHitSizeV;                // Number of hit cells along v direction in cluster
     int _rootHitCellU;                // Hit u coordinate lies on this u cell
     int _rootHitCellV;                // Hit v coordinate lies on this v cell
     int _rootHitSeedCellU;            // Cluster seed pixel uCell  
     int _rootHitSeedCellV;            // Cluster seed pixel vCell 
     
  }; // Class

} // Namespace

#endif 

