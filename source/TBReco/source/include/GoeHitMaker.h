// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    GoeHitMaker - Marlin Processor - Compute pixel hits from clusters                     //
// ///////////////////////////////////////////////////////////////////////////////////////  //


#ifndef GoeHitMaker_H
#define GoeHitMaker_H 1

// Include TBTools 
#include "TBDetector.h"

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
    
    //!Method searching for clusterID id on sensor sensorID in clusterDB. Returns success.  
    bool searchDB(int sensorID, std::string id, double& u, double& v, double& sig2_u, double& sig2_v, double& cov_uv);
      
    //! Input cluster collection name
    std::string  _clusterCollectionName; 
     
    //! Output hit collection name
    std::string  _hitCollectionName;
     
    //! Name of clusterDB file 
    std::string  _clusterDBFileName;
    
    //! Clusters having more cells in u/v are considered bad
    int _maxSizeU;
    int _maxSizeV;
    
   private:
    
    // Handle to detector data 
    TBDetector _detector;  

    // Store pre coomputed cluster measurements 
    // Key is sensorID 
    std::map< int, TH1F *> _DB_Map_ID;
    std::map< int, TH1F *> _DB_Map_U; 
    std::map< int, TH1F *> _DB_Map_V; 
    std::map< int, TH1F *> _DB_Map_Sigma2_U;
    std::map< int, TH1F *> _DB_Map_Sigma2_V; 
    std::map< int, TH1F *> _DB_Map_Cov_UV;
      
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
    
    std::map< int, int> _countAllMap;   //!< Number of clusters per sensor
    std::map< int, int> _countBadMap;   //!< Number of bad clusters per sensor
    std::map< int, int> _countCalMap;   //!< Number of calibrated clusters per sensor
    
  }; // Class

} // Namespace

#endif 

