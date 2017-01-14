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
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>

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
       
    //! ROOT output file name  
    std::string _rootFileName;  
    
    //! DUT plane number, counting sensors in gear file along the beam line starting at zero 
    int _idut; 

    //! Minimum number of clusters occurances 
    int _minClusters; 
      
    //! Max residual for hit-track matching in mm
    double _maxResidualU; 
    double _maxResidualV; 
    
    //! Create clusterDB from hitmaker flag 
    bool _createDBFromHitMaker;

    //! Histogram binning parameters
    std::vector<float >  _binningU;
    std::vector<float >  _binningV;
    std::vector<float >  _binningTu;
    std::vector<float >  _binningTv;
    std::vector<float >  _binningMom;
    
   private:
    
    // ROOT_OUTPUT 
    TFile * _rootFile;
    
    std::map<std::string, int>  _clusterMap;  
    
    std::map< std::string, TH1F *> _histoMapU;
    std::map< std::string, TH1F *> _histoMapV;
    std::map< std::string, TH1F *> _histoMapTu;
    std::map< std::string, TH1F *> _histoMapTv;
    std::map< std::string, TH1F *> _histoMapMom;

    
    std::map< std::string, TH2F *> _histoMapU_V; 
    std::map< std::string, TH2F *> _histoMapU_Tu; 
    std::map< std::string, TH2F *> _histoMapU_Tv; 
    std::map< std::string, TH2F *> _histoMapU_Mom; 
    std::map< std::string, TH2F *> _histoMapV_Tu; 
    std::map< std::string, TH2F *> _histoMapV_Tv; 
    std::map< std::string, TH2F *> _histoMapV_Mom;  
    std::map< std::string, TH2F *> _histoMapTv_Tu; 
    std::map< std::string, TH2F *> _histoMapTv_Mom; 
    std::map< std::string, TH2F *> _histoMapTu_Mom; 
   
     
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



