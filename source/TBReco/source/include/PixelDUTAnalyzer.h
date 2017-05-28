// ///////////////////////////////////////////////////////////////////////////////////////  //
//                                                                                          //
//    PixelDUTAnalyzer - Marlin Processor                                                    //
// ///////////////////////////////////////////////////////////////////////////////////////  //

#ifndef PixelDUTAnalyzer_H
#define PixelDUTAnalyzer_H 1


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

//! PixelDUTAnalyzer Processor 
/*! 
 *  The task of this processor is to analyze the intrinsic performance of a
 *  device under test (DUT) using a reference telescope. Instead of making 
 *  histograms directly, the processor outputs a TFile with flat 
 *  TTree objects allowing simple TTree Draw() to make final plots. The 
 *  main advantage of this approach is more user flexibility on final analysis 
 *  cuts and histogram layout.  
 * 
 *  The DUT can be a planar pixel sensor, single sided strip sensor, or double 
 *  sided strip sensor. For single sided strips, please deactivate the cut on
 *  the maximum allowed residual on the long side of the strip. 
 *  
 *  Note: It is assumed that the trackfinder that input tracks do not include 
 *  hits on the selected DUT plane. Please set steering parameters in the 
 *  trackfinder processor accordingly.  
 *  
 *  Relevant DUT benchmark variables that can be studied with this processor are 
 *  the hit detection efficiency, fake rate, signal spectra, spatial resolution
 *  and in pixel/strip studies of charge collection efficiency and charge sharing.  
 *  
 *  
 *  Author: B.Schwenker, Universität Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */
   
  
class PixelDUTAnalyzer : public marlin::Processor {
   
 public:
   
//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new PixelDUTAnalyzer ; }
    
//!Constructor - set processor description and register processor parameters
   PixelDUTAnalyzer();
   
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
   std::string _trackColName;
   
//! Input DUT TrackerHit collection name
   std::string _hitColName;   
   
//! Alignment DB file name 
   std::string _alignmentDBFileName;
   
//! ROOT output file name  
   std::string _rootFileName;  
   
//! DUT plane number, counting sensors in gear file along the beam line starting at zero 
   int _idut; 
      
//! Max residual for hit-track matching in 
   double _maxResidualU;  // in mm, local DUT uvw coordinate frame 
   double _maxResidualV;  // in mm, local DUT uvw coordinate frame 
     
// ROOT_OUTPUT 
   TFile * _rootFile;

   /** One entry per cluster on the DUT */
   TTree * _rootHitTree;
   /** One entry per track in the reference telescope */
   TTree * _rootTrackTree;
   /** One entry per event */
   TTree * _rootEventTree;    
    
   // Variables in hit tree       
   int _rootEventNumber;             // Event number from lcio file
   int _rootRunNumber;               // Run number from lcio file 
   int _rootDEPFETGoodEvent;         // DEPFET good event flag
   int _rootDEPFETStartGate;         // DEPFET startgate
   int _rootSensorID;                // SensorID from lcio file (this is typically NOT the plane number!!)
   int _rootNTelTracks;              // Number of tracks in reference telescope in same event as hit
   int _rootNDUTHits;                // Number of DUT hits in the same event as hit   
   int _rootHitQuality;              // GoodCluster == 0, BadCluster != 0
   double _rootHitU;                 // Hit coordinate u reconstructed from DUT cluster in mm, in local DUT uvw coordinates       
   double _rootHitV;                 // Hit coordinate v reconstructed from DUT cluster in mm, in local DUT uvw coordinates     
   double _rootHitClusterCharge;     // Sum over all charges in the cluster 
   double _rootHitSeedCharge;        // Highest charge in cluster
   int _rootHitSize;                 // Number of hit cells (pixels/strips) in cluster
   int _rootHitSizeU;                // Number of hit cells along u direction in cluster
   int _rootHitSizeV;                // Number of hit cells along v direction in cluster
   int _rootHitCellU;                // Hit u coordinate lies on this u cell
   int _rootHitCellV;                // Hit v coordinate lies on this v cell
   int _rootHitHasTrack;             // Hit can be matched to track (== 0)     
   double _rootHitFitMomentum;       // Estimated track momentum from fit, only filled in case HasTrack==0            
   double _rootHitFitU;              // Estimated track intersection u coordimate in mm, in local DUT uvw coordinates        
   double _rootHitFitV;              // Estimated track intersection v coordimate in mm, in local DUT uvw coordinates                  
   double _rootHitFitdUdW;           // Estimated track slope du/dw in radians, in local DUT uvw coordinates       
   double _rootHitFitdVdW;           // Estimated track slope dv/dw in radians, in local DUT uvw coordinates    
   double _rootHitFitErrorU;         // Estimated 1x sigma uncertainty for track intersection u coordinate
   double _rootHitFitErrorV;         // Estimated 1x sigma uncertainty for track intersection v coordinate
   double _rootHitPullResidualU;     // Standardized residual in u direction, should have mean = 0 and rms = 1
   double _rootHitPullResidualV;     // Standardized residual in v direction, should have mean = 0 and rms = 1        
   int _rootHitFitCellU;             // Estimated track intersection u coordinate lies on this u cell      
   int _rootHitFitCellV;             // Estimated track intersection v coordinate lies on this v cell         
   double _rootHitFitCellUCenter;    // Central coordinate of cell 'FitCellU' in mm        
   double _rootHitFitCellVCenter;    // Central coordinate of cell 'FitCellV' in mm 
   double _rootHitTrackChi2 ;        // Chi2 value from fit of reference track
   double _rootHitLocalChi2;         // Chi2 value from hit-track residual on device under test 
   int _rootHitTrackNDF;             // Number of degrees of freedom of track fit
   int _rootHitTrackNHits;           // Number of telescope hits used for track fitting 
   
   
   // Variables in track tree  
   int _rootTrackHasHit;             // Track can be matched to a DUT hit (== 0) 
   double _rootTrackFitMomentum;     // Estimated track momentum from fit    
   int _rootTrackNDF;                // Number of degrees of freedom of track fit
   double _rootTrackChi2;            // Chi2 value from fit of reference track
   double _rootTrackLocalChi2;       // Chi2 value from hit-track residual on device under test 
   double _rootTrackFitU ;           // Estimated track intersection u coordimate in mm, in local DUT uvw coordinates 
   double _rootTrackFitV ;           // Estimated track intersection v coordimate in mm, in local DUT uvw coordinates 
   double _rootTrackFitdUdW;         // Estimated track slope du/dw in radians, in local DUT uvw coordinates     
   double _rootTrackFitdVdW;         // Estimated track slope dv/dw in radians, in local DUT uvw coordinates     
   int _rootTrackFitCellU;           // Estimated track intersection u coordinate lies on this u cell  
   int _rootTrackFitCellV;           // Estimated track intersection v coordinate lies on this v cell   
   int _rootTrackNHits;              // Number of telescope hits used for track fitting 
   double _rootTrackFitCellUCenter;  // Central coordinate of cell 'FitCellU' in mm 
   double _rootTrackFitCellVCenter;  // Central coordinate of cell 'FitCellV' in mm 
   double _rootTrackSeedCharge;      // Highest charge in cluster, only filled if cluster matched
   
   
   
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



