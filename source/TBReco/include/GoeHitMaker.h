// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    GoeHitMaker - Marlin Processor - Compute pixel hits from clusters                     //
// ///////////////////////////////////////////////////////////////////////////////////////  //


#ifndef GoeHitMaker_H
#define GoeHitMaker_H 1

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include <marlin/Exceptions.h>

// Include basic C
#include <string>

// Include ROOT classes
#include <TH1F.h>



namespace depfet {

  //! GoeHitMaker processor
  /*! This processor computes 2D hits from pixel clusters. The 
   *  processor uses a clusterDB file with pre-computed offsets 
   *  and covariance matrix entries. 
   *  
   *  The needed clusterDB is either produced by the GoeClusterCalibrator
   *  or by the GoeClusterCalibratorForMC processors. 
   *  
   *  Clusters not found in the clusterDB are currently ignored. They are
   *  not available for tracking or DUT efficiency studies. The fraction 
   *  of ignored clusters is called clusterDB coverage inefficiency and is
   *  computed for all sensors. Coverage efficiencies are finally written 
   *  to the log file.  
   * 
   *  Author: Benjamin Schwenker, Göttingen University 
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
   
  class GoeHitMaker : public marlin::Processor {
  
   public:
   
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new GoeHitMaker ; }
    
    //!Constructor - set processor description and register processor parameters
    GoeHitMaker();
    
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
    
    //!Method printing processor parameters
    void printProcessorParams() const;
    
   protected:
     
    //! Input cluster collection name
    std::string  _clusterCollectionName; 
     
    //! Output hit collection name
    std::string  _hitCollectionName;
     
    //! Name of clusterDB file 
    std::string  _clusterDBFileName;

    //! Sigma U correction factors
    //! One value for sizeU = 1,2,3,... 
    std::vector<float> _sigmaUCorrections;
    
    //! Sigma V correction factors 
    //! One value for sizeV = 1,2,3,... 
    std::vector<float> _sigmaVCorrections; 
    
    //! Use Center of Gravity in case position is 
    //! not available from clusterDB. 
    bool _useCoGFallback; 

   private:
   
    //!Method searching for clusterID id on sensor sensorID in clusterDB. Returns success.  
    bool searchDB(int sensorID, std::string id, double& u, double& v, double& sig2_u, double& sig2_v, double& cov_uv);
     
    

    // Store cluster calibration
    TH1F * m_DB_Weight;
    TH1F * m_DB_U; 
    TH1F * m_DB_V; 
    TH1F * m_DB_Sigma2_U;
    TH1F * m_DB_Sigma2_V; 
    TH1F * m_DB_Cov_UV;
    TH1F * m_DB_Types;  
    
    // Map for eta bin edges stored by type name 
    std::map<std::string, std::vector<double> > m_etaBinEdgesMap;
    
    //! Periodicity for vCells used for clusterDB
    int _vCellPeriod {1};
    
    //! Periodicity for uCells used for clusterDB
    int _uCellPeriod {1}; 
    
    //! Incident angle into sensor DuDw
    double _thetaU {0}; 
    
    //! Incident angle into sensor DvDw
    double _thetaV {0};
    
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
    
    std::map< int, int> _countAllMap;   //!< Number of clusters per sensor
    std::map< int, int> _countCalMap;   //!< Number of calibrated clusters per sensor
    
  }; // Class

} // Namespace

#endif 
