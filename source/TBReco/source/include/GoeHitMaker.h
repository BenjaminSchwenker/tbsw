// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    GoeHitMaker - Marlin Processor - Compute pixel hits from clusters                //
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
#include <TFile.h>
#include <TH1F.h>



namespace depfet {

  //! GoeHitMaker processor
  /*! This processor computes 2D hits from pixel clusters. The 
   *  processor uses a clusterDB file with pre-computed offsets 
   *  and covariance matrix entries. 
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
    
   private:
    
    // Handle to detector data 
    TBDetector _detector;    
    
    //! ClusterDB 
    TFile * _clusterDBFile;
    TH1F *_DB_U;         
    TH1F *_DB_V;       
    TH1F *_DB_Sigma_U;   
    TH1F *_DB_Sigma_V;  
    TH1F *_DB_Cov_UV;   
    
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
    
  }; // Class

} // Namespace

#endif 

