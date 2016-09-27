// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    StripDUTAnalyzer - Marlin Processor                                                    //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef StripDUTAnalyzer_H
#define StripDUTAnalyzer_H 1


// DEPFETTrackTools includes
#include "TBDetector.h"


// Include basic C
#include <vector>
#include <string>
#include <map>

// Include LCIO classes
#include "lcio.h"
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include ROOT classes
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

namespace depfet {

//! StripDUTAnalyzer Processor 
/*! 
 *  The task of this processor is to analyze the intrinsic performance of a
 *  device under test (DUT) using a reference telescope. Instead of making 
 *  histograms directly, the processor outputs a TFile with flat  
 *  TTree objects. This allows simple TTree Draw() commands to make final plots. 
 *  wiht root scripts. The main advantage of this approach is more user flexibility 
 *  on final analysis cuts and histogram layout.  
 *  
 *  Relevant DUT benchmark variables for study are the hit detection efficiency,
 *  the occupancy, signal spectra, spatial resolution and in pixel (strip) studies 
 *  of charge collection efficiency and charge sharing.  
 *  
 *  
 *  Author: B.Schwenker, Universität Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */
   
  
class StripDUTAnalyzer : public marlin::Processor {
   
 public:
   
//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new StripDUTAnalyzer ; }
    
//!Constructor - set processor description and register processor parameters
   StripDUTAnalyzer();
   
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
   
//! A function to read DUT status map 
/*! Vector of status values from DUT  
 */
   bool getDUTStatus(LCEvent * evt, ShortVec & statusVec) ;
   
// Processor Parameters
   
//! Input Track collection name
   std::string _trackColName;
   
//! Input DUT TrackerHit collection name
   std::string _hitColName;   
   
//! Input DUT status data collection name
   std::string _statusColName;
   
//! Alignment DB file name 
   std::string _alignmentDBFileName;
   
//! ROOT output file name  
   std::string _rootFileName;  
   
//! DUT plane number 
   int _idut; 
   

      
//! Max hit2track distance for hit-track matching  
   double _hitdistmax; 
     
// ROOT_OUTPUT 
   TFile * _rootFile;
   TTree * _rootHitTree;
   TTree * _rootTrackTree;
   TTree * _rootEventTree;    
    
   // Common hit/track/event tree  
   int _rootEventNumber;
   int _rootRunNumber;  
   int _rootDetectorID;
   int _rootNTelTracks;  
   int _rootNDUTHits; 
   
   // Variables in hit tree          
   int _rootHitQuality;       // GoodCluster == 0
   double _rootHitU;          // in mm        
   double _rootHitV;          // in mm  
   double _rootHitSigmaU; 
   double _rootHitSigmaV; 
   double _rootHitClusterChargeU; 
   double _rootHitSeedChargeU;
   double _rootHitClusterChargeV; 
   double _rootHitSeedChargeV;
   int _rootHitSize;  
   int _rootHitSizeU;     
   int _rootHitSizeV;    
   int _rootHitCellU;
   int _rootHitCellV; 
   int _rootHitHasTrack;             // Matched to track == 0     
   double _rootHitFitMom; 
   double _rootHitFitSigmaMom;              
   double _rootHitFitU;              // Track impact point [mm] - in local frame       
   double _rootHitFitV;              
   double _rootHitFitdUdW;           // Track direction tangent [rad] - in local frame       
   double _rootHitFitdVdW;      
   double _rootHitFitSigmaU;          // rms error u
   double _rootHitFitSigmaV;          // rms error v 
   double _rootHitPullResidualU; 
   double _rootHitPullResidualV;             
   int _rootHitFitCellU;               // DUT  - readout column         
   int _rootHitFitCellV;               // DUT  - readout row   
   double _rootHitFitCellCenterU;        // Hit  center in u[mm]        
   double _rootHitFitCellCenterV;        // Hit  center in v[mm] 
   double _rootHitTrackChi2 ; 
   int _rootHitTrackNDF; 
   double _rootHitLocalChi2; 
   
 
   // Variables in track tree  
   int _rootTrackHasHit; 
   double _rootTrackFitMom;  
   int _rootTrackNDF;          // Track degrees of freedof (ndof) 
   double _rootTrackChi2;      // Track chisq
   double _rootTrackFitU ; 
   double _rootTrackFitV ; 
   double _rootTrackFitdUdW; 
   double _rootTrackFitdVdW; 
   int _rootTrackFitCellU; 
   int _rootTrackFitCellV;  
   double _rootTrackFitCellCenterU ; 
   double _rootTrackFitCellCenterV ; 
   int _rootTrackCellQualityU;
   int _rootTrackCellQualityV;        
   
   // Variables in event tree
   int _rootEventNDUTTracks; 
   int _rootEventNMatched; 
   
 private:
   
   // Handle to detector data 
   TBDetector  _detector;    
   
   // Few counter to show the final summary
   
   //! Number of event w/o input hit
   /*! It could happen that no input hit have passed the selection
    *  criteria posed by previous processors. Or, no particle passed 
    *  the DUT sensor in this event. 
    */
   int _noOfEventWOInputHit;

   //! Number of event w/o input track
   /*! It could happen that no particle track passed the selection
    *  criteria posed by previous processors. Or, no particle passed 
    *  the telescope sensors. 
    */
   int _noOfEventWOInputTrack;
   
   //! Total number of acccepted tracks
   /*! This is the total number of tracks passing DUTTool track selection
    *  cuts.
    */
   int _noOfTracks;       

   //! Total number of DUT hits 
   /*! This is the total number of DUT hits the processor was able to
    *  reconstruct.
    */
   int _noOfHits;     

   //! Total number of matched tracks
   /*! This is the total number of tracks the processor was able to
    *  match with DUT hits.
    */
   int _noOfMatchedTracks;     

   //! Number of hits with hit no track in telescope
   /*! A DUT hit without associated track may 
    *  indicate a lost track (desync or tracking)
    */
   int _noOfHitsWOTrack;      

   //! Number of hits with track in telescope
   /*! A track indicates synchronisation of telescope, 
    *  helps debugging
    */
   int _noOfHitsWTrack;    

   //! Last matched event
   /*! Event number of last matched hit track 
    *  pair. 
    */
   int _iLastMatchedEvent;    

   //! First matched event
   /*! Event number of first matched hit track 
    *  pair. 
    */
   int _iFirstMatchedEvent;     
   
   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
   
}; // Class

} // Namespace

#endif 



