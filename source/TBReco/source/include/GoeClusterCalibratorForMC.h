// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    GoeClusterCalibratorForMC - Marlin Processor                                         //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef GoeClusterCalibratorForMC_H
#define GoeClusterCalibratorForMC_H 1

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
#include <TH1F.h>
#include <TH2F.h>

namespace depfet {

  //! GoeClusterCalibratorForMC  
  /*! 
   *  The task of this processor is to create a clusterDB for all clusters
   *  contained in the input cluster collection. After calibration, the 
   *  clusterDB will be used by the GoeHitMaker processor to compute a 2D 
   *  position measurement (hit) in local sensor coordinates together with
   *  a 2x2 covariance matrix for all cluster registered in the clusterDB.
   * 
   *  This version of the ClusterCalibrator requires a collection of simHits 
   *  for creating the clusterDB. The final clusterDB can be used both for 
   *  simulated detector clusters and real detector clusters. In the later 
   *  case, the accuracy of the clusterDB depends on the accuracy of the 
   *  detector simulation (digitizer). 
   *    
   *  Author: B.Schwenker, Universität Göttingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
   
  class GoeClusterCalibratorForMC : public marlin::Processor {
   
   public:
   
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new GoeClusterCalibratorForMC ; }
    
    //!Constructor - set processor description and register processor parameters
    GoeClusterCalibratorForMC();
    
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
   
    // Processor Parameters
    
    //! Input Track collection name
    std::string _simHitColName;
    
    //! Input TrackerHit collection name
    std::string _hitColName;   
       
    //! ROOT output file name  
    std::string _clusterDBFileName;   
    
    //! Minimum number of clusters occurances 
    int _minClusters; 
      
    //! Position of steps for software ADC 
    std::vector<int> _swADCSteps; 

    //! Max residual for hit-track matching in mm
    double _maxResidualU; 
    double _maxResidualV; 

    //! Ignore clusters from these sensorIDs 
    std::vector<int>  _ignoreIDVec;
    
   private:
  
    // Intermediate histos to compute calibrated measurements 
    // Key is cluster label
    std::map<std::string, int>   _sensorMap;  
    std::map<std::string, TH1F *>  _clusterUMap;
    std::map<std::string, TH1F *>  _clusterVMap;
    std::map<std::string, TH2F *>  _clusterUVMap;
    
    std::map< std::string, TH1F *> _histoMap;
    
    // Handle to detector data 
    TBDetector  _detector;    
     
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
     
  }; // Class

} // Namespace

#endif 



