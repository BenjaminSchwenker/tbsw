// ////////////////////////////////////////////////////////////////////////// //
//                                                                            //
//    TrackFitDQM - Marlin Processor                                          //
// ////////////////////////////////////////////////////////////////////////// //

#ifndef TrackFitDQM_H
#define TrackFitDQM_H 1


// TBTools includes   
#include "TBDetector.h"

// Include ROOT classes
#include <TFile.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TH2D.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include basic C
#include <string>
#include <vector>
#include <map>

namespace depfet {

  //! TrackFitDQM processor 
  /*! The Processor produces data driven DQM histos to test the  
   *  quality of track fitting and detector alignment. All histos
   *  are written into a root file. 
   *    
   *  Author: Benjamin Schwenker, Göttingen University 
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */
  
  
  class TrackFitDQM : public marlin::Processor {
  
   public:
     
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new TrackFitDQM ; }
    
    //!Constructor - set processor description and register processor parameters
    TrackFitDQM();
    
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
    
    //! Histogram booking
    void bookHistos();
    
    //! Processor Parameters 
    
    //! Input track collection name
    std::string _inputTrackCollectionName;
      
    //! AlignmentDB file name 
    std::string _alignmentDBFileName;
    
    //! ROOT output file name  
    std::string _rootFileName;
         
   private:
   
    double _timeCPU; //!< CPU time
    int    _nRun ;   //!< Run number
    int    _nEvt ;   //!< Event number
    
    // Handle to aligned detector data 
    TBDetector  _detector;     
    
    // Handle to root file
    TFile * _rootFile;
    
    std::map< std::string, TH1D *> _histoMap;
    std::map< std::string, TH2D *> _histoMap2D;  
    std::map< std::string, TProfile *> _profileMap; 
     
  }; // Class

} // Namespace

#endif 

