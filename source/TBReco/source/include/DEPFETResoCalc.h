// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    DEPFETResoCalc - Marlin Processor - Calculate position resolution functions           //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef DEPFETResoCalc_H
#define DEPFETResoCalc_H 1



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
#include <IMPL/TrackerHitImpl.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include ROOT classes
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

namespace depfet {

//! DEPFETResoCalc  
/*! 
 *  The task of this processor is to study the relation between true hit position 
 *  and the digitized detector response simply called a hit. In particular, the 
 *  includes a sampling of the measurement position errors, the landau signal and
 *  hit efficiency or fake rates. 
 *   
 *  This processor requires two inputs: 
 *   
 *  - Collection of reconstructed hits: it means simulated particle crossings 
 *    after full digitization and clustering. Each simulated hit has an estimate 
 *    for hit coordinates and a covariance matrix. 
 *  
 *  - Collection of truth hits: The truth hit has the coordinates of the true 
 *    particle crossing. Each truth hit is calculated from Geant4 steps of the 
 *    particle crossing the detectors sensitive volume.   
 *   
 *  The processor matches simulated hits to truth hits and fills ntuples 
 *  for a final resolution analysis in root. The output is a root file 
 *  with a tuples for final user analysis. 
 *   
 *  Author: B.Schwenker, Universität Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */
   
  
class DEPFETResoCalc : public marlin::Processor {
   
 public:
   
//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new DEPFETResoCalc ; }
    
//!Constructor - set processor description and register processor parameters
   DEPFETResoCalc();
   
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
   
//! Input truth hit collection  
   std::string _truthHitColName;
   
//! Input reco hit collection  
   std::string _recoHitColName;   
   
//! ROOT output file name  
   std::string _rootFileName;  
   
//! DUT plane number 
   int _idut; 
   
//! Max distance for matching 
   double _hitdistmax;
   
//! Telescope error along u axis 
   double _telerrorU; 

//! Telescope error along v axis 
   double _telerrorV; 

//! Run number to identify mc sample
   int _runnumber;  
   

// ROOT_OUTPUT 
   TFile * _rootFile;
   TTree * _rootHitTree;
   TTree * _rootTrackTree;
   TTree * _rootEventTree;     
   
   // Common hit/track/event tree  
   int _rootEventNumber;
   int _rootRunNumber;  
   int _rootDetectorID;
   int _rootStartGate; 
   int _rootNTelTracks;  
   int _rootNDUTHits; 
   
   // Variables in hit tree          
   int _rootHitQuality;       // GoodCluster == 0
   double _rootHitU;         // in mm        
   double _rootHitV;         // in mm  
   double _rootHitCharge; 
   double _rootHitSeedCharge;
   int _rootHitSize;  
   int _rootHitSizeCol;     
   int _rootHitSizeRow;    
   int _rootHitCol;
   int _rootHitRow; 
   int _rootHitHasTrack;             // Matched to track == 0         
   double _rootHitFitU;              // Track impact point [mm] - in local frame       
   double _rootHitFitV;              
   double _rootHitFitdUdW;           // Track direction tangent [rad] - in local frame       
   double _rootHitFitdVdW;  
   double _rootHitFitErrorU;         // rms error u
   double _rootHitFitErrorV;         // rms error v     
   double _rootHitPullResidualU; 
   double _rootHitPullResidualV;  
   int _rootHitFitCol;               // DUT pixel - readout column         
   int _rootHitFitRow;               // DUT pixel - readout row   
   double _rootHitFitPixU;           // Hit pixel, center in u[mm]        
   double _rootHitFitPixV;           // Hit pixel, center in v[mm] 
   double _rootHitTrackChi2; 
   int _rootHitTrackNDF; 
   double _rootHitLocalChi2; 
   
    
   // Variables in track tree  
   int _rootTrackHasHit; 
   int _rootTrackNDF;          // Track degrees of freedof (ndof) 
   double _rootTrackChi2;      // Track chisq
   double _rootTrackFitU ; 
   double _rootTrackFitV ; 
   double _rootTrackFitdUdW; 
   double _rootTrackFitdVdW; 
   int _rootTrackFitCol; 
   int _rootTrackFitRow;  
   double _rootTrackFitPixU ; 
   double _rootTrackFitPixV ; 
   double _rootTrack1x1Charge;     // Pixle charge at Xing point 
   double _rootTrack3x3Charge;     // Cluster charge at Xing point 
   int _rootTrack1x1Quality;       // Good Xing = 0   
   int _rootTrack3x3Quality;       // Good Xing = 0  
   double _rootTrackSeedCharge;   
   
   // Variables in event tree
   int _rootEventNDUTTracks; 
   int _rootEventNMatched; 
   
 private:
   
   // Handle to detector data 
   TBDetector  _detector;    
     
   // Few counter to show the final summary
   
   //! Total number of acccepted tracks
   /*! This is the total number of tracks passing DUTTool track selection
    *  cuts.
    */
   int _noOfTruthHits;       
   
   //! Total number of reco hits 
   /*! This is the total number of reco hits the processor was able to
    *  reconstruct.
    */
   int _noOfRecoHits;     
   
   //! Total number of matched hits
   /*! This is the total number of hits the processor was able to
    *  match with truth hits.
    */
   int _noOfMatchedHits;     
   
   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
   
}; // Class

} // Namespace

#endif 



