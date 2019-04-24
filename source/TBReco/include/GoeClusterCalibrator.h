// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    GoeClusterCalibrator - Marlin Processor                                         //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef GoeClusterCalibrator_H
#define GoeClusterCalibrator_H 1

// Include basic C
#include <vector>
#include <string>
#include <map>
#include <set>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include ROOT classes
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>


namespace depfet {

  //! GoeClusterCalibrator  
  /*! 
   *  The task of this processor is to create a clusterDB for clusters
   *  contained in the input cluster collection. After calibration, the 
   *  clusterDB will be used by the GoeHitMaker processor to compute a 2D 
   *  position measurement (hit) in local sensor coordinates together with
   *  a 2x2 covariance matrix for all cluster registered in the clusterDB.
   *   
   *  This version of the ClusterCalibrator requires a collection of reco
   *  tracks in a fully aligned telescope for creating the clusterDB.
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
    
    //! Name of clusterDB output file   
    std::string _clusterDBFileName;   
       
    //! Minimum number of clusters occurances 
    int _minClusters; 
    
    //! Periodicity for vCells used for clusterDB
    int _vCellPeriod;
     
    //! Periodicity for uCells used for clusterDB
    int _uCellPeriod; 
    
    //! Max number of eta bins  
    int _maxEtaBins; 
    
    //! Minimum variance of clusters covariance matrix
    float _minVarianceU; 
    float _minVarianceV;
    
    //! Ignore clusters from these sensorIDs 
    std::vector<int >  _ignoreIDVec;
       
    //! Name of temporary file for collector output
    std::string _collectorOutputFileName;   
 
   private:
    
    /** Collector output file */ 
    TFile * _rootCollectorOutputFile;
    /** Histogram for track covariance matrix element UU */
    TH1F * _trackVarUHisto;
    /** Histogram for track covariance matrix element VV */
    TH1F * _trackVarVHisto;
    /** Histogram for track covariance matrix element UV */
    TH1F * _trackCovUVHisto;
    /** Histograms for track incident angle DuDw */
    TH1F * _trackDuDwHisto;
    /**Histograms for track incident angle DvDw */
    TH1F * _trackDvDwHisto;
    /**Container for histos*/
    std::map< std::string, TH1F *> _histoMap;
    /** Name of cluster tree */
    TTree * m_rootTree; 
    /** Name of cluster type */
    std::string m_typeName;
    /** Eta value of cluster for sector (+,+)*/
    float m_clusterEtaPP;
    /** Eta value of cluster for sector (-,+)*/
    float m_clusterEtaNP;
    /** Eta value of cluster for sector (+,-)*/
    float m_clusterEtaPN;
    /** Eta value of cluster for sector (-,-)*/
    float m_clusterEtaNN; 
    /** Position offset u of cluster */
    float m_positionOffsetU;
    /** Position offset v of cluster */
    float m_positionOffsetV;
            
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
    
    std::set<int> _setOfPlaneNumbers; //!< Set of plane numbers to be corrected  
 
  }; // Class

} // Namespace

#endif 



