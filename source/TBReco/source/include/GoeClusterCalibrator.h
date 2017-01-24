// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    GoeClusterCalibrator - Marlin Processor                                         //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef GoeClusterCalibrator_H
#define GoeClusterCalibrator_H 1

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

  //! GoeClusterCalibrator  
  /*! 
   *  The task of this processor is to create a clusterDB for all clusters
   *  contained in the input track collection. After the calibration is 
   *  finished, the clusterDB is used by the GoeHitMaker to compute a 2D  
   *  position measurement (hit) together with a 2x2 covariance matrix 
   *  for all cluster shapes registered in the clusterDB .  
   *   
   *  This version of the ClusterCalibrator requires a collection of reco
   *  tracks in a fully aligned telescope for creating the clusterDB. A full 
   *  calibration schema involves bootstraping and iterations:
   *
   *  1) In the first iteration, clustering must be carried out with an 
   *     heuristic method like the center of gravity (bootstraping).  
   * 
   *     Center of gravity hits must also be used to as input to tracking 
   *     and track based alignment processors. 
   * 
   *  2) After track finding, a Kalman filter/smoother is used to estimate 
   *     local track states for all clusters in the track. 
   * 
   *  3) The cluster calibrator uses a sample of 'cluster:track-estimator' pairs
   *     to compute the clusterDB.
   *     
   *  4) Now, the process can be iterated. This time using the clusterDB to 
   *     compute hits instead of the center of gravity. 
   *   
   *  In many cases it will not be necessary to repeat the alignment processor 
   *  using the hits from the clusterDB.    
   *  
   *  Author: B.Schwenker, Universität Göttingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
   
  class GoeClusterCalibrator : public marlin::Processor {
    
   public:
    
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new GoeClusterCalibrator ; }
    
    //!Constructor - set processor description and register processor parameters
    GoeClusterCalibrator();
    
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
    
    //! Input track collection name
    std::string _inputTrackCollectionName;
    
    //! AlignmentDB file name 
    std::string _alignmentDBFileName;
       
    //! ROOT output file name  
    std::string _clusterDBFileName;   
    
    //! Minimum number of clusters occurances 
    int _minClusters; 
    
    //! Minimum variance of clusters covariance matrix
    float _minVarianceU; 
    float _minVarianceV;
        
   private:
    
    // Intermediate histos to compute averaged covariance matrix
    // Key is sensorID 
    std::map< int, TH1F *> _trackVarUMap;
    std::map< int, TH1F *> _trackVarVMap;
    std::map< int, TH1F *> _trackCovUVMap;
    
    // Intermediate histos to compute calibrated measurements 
    // Outer key is sensorID, inner key is clusterID 
    std::map<int, std::map<std::string, int> >  _sensorMap;  
    std::map< int, std::map<std::string, TH1F *> > _clusterUMap;
    std::map< int, std::map<std::string, TH1F *> > _clusterVMap;
    std::map< int, std::map<std::string, TH2F *> > _clusterUVMap;
    
    std::map< std::string, TH1F *> _histoMap;
          
    // Handle to detector data 
    TBDetector  _detector;    
     
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
     
  }; // Class

} // Namespace

#endif 



