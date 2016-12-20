// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    TrackFitValidation - Marlin Processor                                                     //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef TrackFitValidation_H
#define TrackFitValidation_H 1

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
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

namespace depfet {

  //! TrackFitValidation  
  /*! 
   *  The task of this processor is to study the relation between SimTrackerHits 
   *  and estimated track states.  
   *   
   *  The processor matches SimTrackerHits to estimated local track states and 
   *  histograms parameter pulls. 
   *   
   *  Author: B.Schwenker, Universität Göttingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
   
  
  class TrackFitValidation : public marlin::Processor {
   
   public:
   
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new TrackFitValidation ; }
    
    //!Constructor - set processor description and register processor parameters
    TrackFitValidation();
    
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
    
    //! Input SimTrackerHit collection name
    std::string _simHitColName;
    
    //! Input Track collection name
    std::string _trackColName; 
     
    //! Alignment DB file name 
    std::string _alignmentDBFileName;
    
    //! ROOT output file name  
    std::string _rootFileName;  
      
    //! Max residual for hit-track matching in mm
    double _maxResidualU; 
    double _maxResidualV; 
    
   private:
    
    // ROOT_OUTPUT 
    TFile * _rootFile; 
    std::map< std::string, TH1D* > _histoMap;
    std::map< std::string, TH2D* > _histoMap2D;
    
    // Handle to detector data 
    TBDetector  _detector;    
     
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
     
  }; // Class

} // Namespace

#endif 



