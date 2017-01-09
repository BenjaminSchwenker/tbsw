// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    GoeClusterCalibratorFromMC - Marlin Processor                                         //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef GoeClusterCalibratorFromMC_H
#define GoeClusterCalibratorFromMC_H 1

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
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TMath.h>


namespace depfet {

  //! GoeClusterCalibratorFromMC  
  /*! 
   *  The task of this processor is to create a clusterDB for all clusters
   *  contained in the input cluster collection. The clusterDB will be used
   *  later to compute a 2D hit measurement together with a 2x2 covariance
   *  matrix for any given cluster.  
   *   
   *  Author: B.Schwenker, Universität Göttingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
   
  class GoeClusterCalibratorFromMC : public marlin::Processor {
   
   public:
   
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new GoeClusterCalibratorFromMC ; }
    
    //!Constructor - set processor description and register processor parameters
    GoeClusterCalibratorFromMC();
    
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
    
    //! Input TrackerHit collection name
    std::string _hitColName;   
       
    //! ROOT output file name  
    std::string _clusterDBFileName;   
    
    //! Minimum number of clusters occurances 
    int _minClusters; 
      
    //! Max residual for hit-track matching in mm
    double _maxResidualU; 
    double _maxResidualV; 
    
   private:
    
    // ROOT_OUTPUT 
    TFile * _rootFile;
    
    std::map<std::string, int>  _clusterMap;  
    
    std::map< std::string, TH1F *> _histoMapU;
    std::map< std::string, TH1F *> _histoMapV;
    std::map< std::string, TH2F *> _histoMapUV; 
     
    /** One entry per cluster */
    TTree * _rootClusterTree;
    
    // Variables in cluster tree      
    std::string _rootClusterID;       // Cluster ID 
    double _rootTrackMomentum;        // Track momentum [GeV/c]  
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



