// //////////////////////////////////////////  //
//                                             //
//    FastTracker - Marlin Processor           //
// //////////////////////////////////////////  //


#ifndef FastTracker_H
#define FastTracker_H 1

// Include TBTools header files
#include "TBDetector.h"
#include "TBTrack.h"
#include "HitFactory.h"

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include basic C
#include <vector>

namespace depfet {

//! Check if tracks trk1 and trk2 have common hits
   bool check_incompatible( TBTrack& trk1, TBTrack& trk2 );
   
//! Is track trk1 better than trk2
   bool compare_tracks( TBTrack& trk1, TBTrack& trk2 );

//! Mark hits in track as used
   void mark_hits( TBTrack& trk, std::vector<std::vector<int>>&  usedIDs );

//! FastTracker Processor
/*! This processor performs track finding and track fitting in beam test 
 *  experiments. The processor allows a track finding in all sensors 
 *  specified in the gear file. Multiple hit collections can be given as 
 *  as input. The processor outputs a single track collection which can 
 *  be used in downstream processors. 
 *
 *  Track finder starts by registering all hits into a HitStore. The HitStore 
 *  sorts hits first by sensorID and assigns two indices to each hit. The HitID is a
 *  unique identifier for all hits in an event and will be used to find overlaps 
 *  between tracks. The SectorID assigns each hit to a small sector from a network of 
 *  sectors covering the full sensor. The sectors are computed from the residual cuts 
 *  supplied as steering parameters. 
 *  
 *  The second step is to build seed tracks from  pairs of sensors along the beam 
 *  line. The user is free to select any pair of sensors. In order to reduce combinatorial 
 *  background, the user can select a cut on the maximal angle between the direction 
 *  of the seed track and the telescope z axis, which is typically also the beam axis.  
 *  For tracks inside a magnetic field, the user has to provide an initial gues for the 
 *  particle momentum and charge. All seed tracks are followed through all sensor along 
 *  the beam line and compatible hits are picked up. The user can impose cuts 
 *  on the minimum number of hits and misses. Track candidates passing all cuts are stored
 *  in a TrackCollector. Finally, a opitmal set of non overlapping tracks is selected from
 *  the collector and stored in the final track  collection.
 *  
 *  Author: Benjamin Schwenker, Göttingen University 
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de> 
 */


class FastTracker : public marlin::Processor {
   
 public:
   
//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new FastTracker ; }
   
//!Constructor - set processor description and register processor parameters
   FastTracker();
   
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
       
//! Method printing processor parameters
   void printProcessorParams() const;
   
//! Called by the processEvent() to find track candidates  
   void findTracks( std::list<TBTrack>& TrackCollector , HitFactory& HitStore , int firstplane , int secondplane, int idir);  

//! Called by the processEvent() to add tracks to trackcollector using hits in hitstore
   void findTracks( std::list<TBTrack>& TrackCollector , HitFactory& HitStore , int seedplane); 

//! Method for building a track candidate from a seed. The finished candidate is added to the track collector.
   void buildTrackCand(TBTrack& trk, HitFactory& HitStore, std::list<TBTrack>& TrackCollector, int idir);

// Processor Parameters
   
//! Input hit collection names
   std::vector< std::string >  _hitCollectionNameVec;  

//! Output collection name for not used hits 
   std::string _notUsedhitCollectionName;
   
//! Output track collection name
   std::string _trackCollectionName;
   
//! Alignment DB file name 
   std::string _alignmentDBFileName;
      
//! Det Selection: Ignore hits from bad detectors
   std::vector<int >  _BadDetectorIDs;

//! Max residual for adding a hit to a track candidate
   std::vector<float >  _maxResidualU;
   std::vector<float >  _maxResidualV;   

//! Hit Selection: Discard bad hits
   int _hitQualitySelect;
   
//! Quality criteria for seed tracks 
   int _minHits; 
   int _maxGap;  
   float _maxSlope; 
   float _maxTrkChi2; 
   float _outlierChi2Cut; 

//! Perform double hit seeding using these plane 
   int _firstPass_firstPlane;
   int _firstPass_secondPlane;
   int _secondPass_firstPlane;
   int _secondPass_secondPlane;

//! Perform single hit seeding using these planes 
   std::vector<int>  _singleHitSeedingPlanes; 
   
//! Parameters for beam particles 
   std::vector<float> _momentum_list; 
   double _mass;
   double _charge;
      
 private: 
             
   // Handle to detector data sheets 
   TBDetector _detector;       
     
   // Total number of pixel sensors 
   int _nTelPlanes;
    
   // Total number of active sensors 
   int _nActivePlanes;
    
   // Active flag for sensors 
   std::vector<bool> _isActive;
   
   // Few counters to show the final summary
   
   //! Number of event with enough firing planes
   /*! 
    */
   int _noOfEventMinHits;


   //! Number of event w/o input hit
   /*! It could happen that no input hit has passed the selection
    *  criteria posed by previous processors.
    */
   int _noOfEventWOInputHit;
   
   //! Number of events with out tracks
   /*! This is the number of events in which the processor wasn't
    *  able to reconstruct any tracks
    */
   int _noOfEventWOTrack;
   
   //! Total number of reconstructed tracks
   /*! This is the total number of tracks the processor was able to
    * reconstruct.
    */
   int _noOfTracks;    

   
   int _noOfCandTracks;    
   int _noOfAmbiguousHits;

   //! Total number of failed final fits
   int _noOfFailedFits;    
   
   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
   
}; // Class

} // Namespace

#endif 
