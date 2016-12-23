// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    MCTreeProducer - Marlin Processor                                                     //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef MCTreeProducer_H
#define MCTreeProducer_H 1

// DEPFETTrackTools includes
#include "TBDetector.h"

// Include basic C
#include <vector>
#include <string>
#include <map>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include ROOT classes
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

namespace depfet {

  //! MCTreeProducer  
  /*! 
   *  The task of this processor is to study the relation between SimTrackerHits 
   *  and clusters.  
   *   
   *  The processor matches SimTrackerHits to reco hits and fills ntuples 
   *  for a final resolution analysis in root. The output is a root file 
   *  with a tuples for final user analysis. 
   *   
   *  Author: B.Schwenker, Universität Göttingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
   
  
  class MCTreeProducer : public marlin::Processor {
   
   public:
   
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new MCTreeProducer ; }
    
    //!Constructor - set processor description and register processor parameters
    MCTreeProducer();
    
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
    
    //!Method printing processor parameters
    void printProcessorParams() const;
   
    // Processor Parameters
    
    //! Input Track collection name
    std::string _simHitColName;
    
    //! Input DUT TrackerHit collection name
    std::string _hitColName;   
     
    //! Alignment DB file name 
    std::string _alignmentDBFileName;
    
    //! ROOT output file name  
    std::string _rootFileName;  
    
    //! DUT plane number, counting sensors in gear file along the beam line starting at zero 
    int _idut; 
      
    //! Max residual for hit-track matching in mm
    double _maxResidualU; 
    double _maxResidualV; 
    
   private:
    
    // ROOT_OUTPUT 
    TFile * _rootFile;
     
    /** One entry per cluster on the DUT */
    TTree * _rootClusterTree;
    
    // Variables in hit tree       
    int _rootEventNumber;             // Event number from lcio file
    int _rootRunNumber;               // Run number from lcio file 
    int _rootSensorID;                // SensorID from lcio file (this is typically NOT the plane number!!)
    int _rootClusterID;               // Cluster ID 
    int _rootClusterIDU;              // Cluster ID for u cluster 
    int _rootClusterIDV;              // Cluster ID for v cluster  
    double _rootClusterPosU;          // Cluster U position [mm]     
    double _rootClusterPosV;          // Cluster V position [mm]   
    double _rootClusterSigmaU;        // Sigma for cluster U position [mm]
    double _rootClusterSigmaV;        // Sigma for cluster V position [mm]   
    double _rootClusterCharge;        // Sum over all charges in the cluster 
    double _rootSeedCharge;           // Highest charge in cluster
    int _rootClusterSize;             // Number of hit cells (pixels/strips) in cluster
    int _rootClusterSizeU;            // Number of hit cells along u direction in cluster
    int _rootClusterSizeV;            // Number of hit cells along v direction in cluster
    int _rootClusterStartCellIdU;     // U Id of start cell
    int _rootClusterStartCellIdV;     // V Id of start cell  
    double _rootClusterStartCellPosU; // U position of start cell 
    double _rootClusterStartCellPosV; // V position of start cell 
    double _rootTrackMomentum;        // Track momentum [GeV/c]  
    double _rootTrackCharge;          // Track charge [e]
    double _rootTrackPosU;            // Track U position [mm] 
    double _rootTrackPosV;            // Track V position [mm]                   
    double _rootTrackdUdW;            // Track slope [rad]     
    double _rootTrackdVdW;            // Track slope [rad]     
    
    // Handle to detector data 
    TBDetector  _detector;    
     
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
     
  }; // Class

} // Namespace

#endif 



