// FastTracker implementation file
//                       
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


// TBTools includes 
#include "FastTracker.h"
#include "SeedGenerator.h"
#include "TBKalmanB.h"
#include "TBHit.h"
#include "TrackInputProvider.h"

// ROOT includes
#include <TMath.h>

// C++ includes 
#include <list>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>

// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;
using namespace CLHEP;

namespace depfet {

//
// Instantiate this object
//
FastTracker aFastTracker ;

//
// Constructor
//
FastTracker::FastTracker() : Processor("FastTracker")
{ 

// Processor description
   _description = "FastTracker: Performs track finding and fitting for beam tests";
   
//   
// First of all we need to register the input/output collections
   
   std::vector< string > inputHitCollectionNameVecExample;
   inputHitCollectionNameVecExample.push_back( "hits" );
   
   registerInputCollections (LCIO::TRACKERHIT, "InputHitCollectionNameVec",
                            "Input hit collection names",
                            _hitCollectionNameVec, inputHitCollectionNameVecExample );
   
   registerOutputCollection(LCIO::TRACK,"OutputTrackCollectionName",
                           "Collection name for fitted tracks",
                           _trackCollectionName, string ("fittracks"));
    
   registerOutputCollection (LCIO::TRACKERHIT, "HitCollectionName",
                            "Name of not used hit collection",
                            _notUsedhitCollectionName, 
                            string("unusedhits"));
   
// 
// Next, initialize the processor paramters
   
   registerProcessorParameter ("AlignmentDBFileName",
                             "This is the name of the LCIO file with the alignment constants (add .slcio)",
                             _alignmentDBFileName, static_cast< string > ( "alignmentDB.slcio" ) );   
   
   std::vector<float> initMomentumList;
   initMomentumList.push_back(4.0);
   registerProcessorParameter ("ParticleMomentum", "List of particle momenta [GeV]",
                              _momentum_list, initMomentumList );
   
   registerProcessorParameter ("ParticleMass", "Particle mass [GeV]",
                              _mass,  static_cast < double > (0.139));
   
   registerProcessorParameter ("ParticleCharge", "Particle charge [e]",
                              _charge,  static_cast < double > (+1));
   
   registerProcessorParameter ("MaxTrackChi2",
                              "Maximum track chisq in track fit",
                              _maxTrkChi2,  static_cast < float > (1000.) );
   
   std::vector<float> initMaxResidual;
   initMaxResidual.push_back(1);
   registerProcessorParameter ("MaxResidualU", "Maximum hit residual [mm]. Only used for pre-seleciton of hits",
                              _maxResidualU, initMaxResidual );
   
   registerProcessorParameter ("MaxResidualV", "Maximum hit residual [mm]. Only used for pre-seleciton of hits",
                              _maxResidualV, initMaxResidual );
   
   registerProcessorParameter("OutlierChi2Cut",  
                             "Maximum chi2 increment for any hit added to a track candidate",                       
                             _outlierChi2Cut, static_cast <float> (50) ); 
    
   registerProcessorParameter ("MinimumHits",
                              "Minimum number of hits in track",
                              _minHits,  static_cast < int > (2));
   
   std::vector<int> initSingleHitSeeding;
   registerProcessorParameter ("SingleHitSeeding", "Seed tracks using hits from these planes",
                              _singleHitSeedingPlanes, initSingleHitSeeding );
   
   registerProcessorParameter ("ForwardPass_FirstPlane",
                              "First plane for seeding the forward Kalman filter pass. Put -1 to deactivate pass",
                              _firstPass_firstPlane,  static_cast < int > (-1));
   
   registerProcessorParameter ("ForwardPass_SecondPlane",
                              "Second plane for seeding the forward Kalman filter pass. Put -1 to deactivate pass",
                              _firstPass_secondPlane,  static_cast < int > (-1));

   registerProcessorParameter ("BackwardPass_FirstPlane",
                              "First plane for seeding the backward Kalman filter pass. Put -1 to deactivate pass",
                              _secondPass_firstPlane,  static_cast < int > (-1)); 

   registerProcessorParameter ("BackwardPass_SecondPlane",
                              "Second plane for seeding the backward Kalman filter pass. Put -1 to deactivate pass",
                              _secondPass_secondPlane,  static_cast < int > (-1)); 
     
   registerProcessorParameter ("MaximumGap",
                              "Maximum number of consecutive missing hits in track",
                              _maxGap,  static_cast < int > (0));
   
   registerProcessorParameter ("MaximumSlope",
                              "Maximum slope of candiate tracks in Z=0 plane (rad). Set cut to zero to deactivate.",
                              _maxSlope,  static_cast < float > (0.1));
   
   registerProcessorParameter("HitQualitySelection",
                             "To use only GoodQuality write 0 here",
                             _hitQualitySelect, static_cast<int> ( 0 ));

   std::vector<int> initBadDetectorIDs;
   registerProcessorParameter ("ExcludeDetector",
                              "Enter plane numbers of bad detectors",
                              _BadDetectorIDs, initBadDetectorIDs);
   
   
   
}

//
// Method called at the beginning of data processing
//
void FastTracker::init() {
  
  // Print set parameters
  printProcessorParams();
   
  // CPU time start
  _timeCPU = clock()/1000;
    
  // Initialize variables
  _nRun = 0 ;
  _nEvt = 0 ;
   
  // Initialize all the counters
  _noOfEventWOInputHit = 0;
  _noOfEventWOTrack = 0;
  _noOfTracks  = 0;
  _noOfEventMinHits = 0;
  _noOfFailedFits = 0;  
  _noOfCandTracks = 0;   
  _noOfAmbiguousHits = 0;
  
  // Read detector constants from gear file
  _detector.ReadGearConfiguration();    
            
  // Read alignment data base file 
  _detector.ReadAlignmentDB( _alignmentDBFileName );      
        
  _nTelPlanes =  _detector.GetNSensors();
  _nActivePlanes = _nTelPlanes; 
  _isActive.resize(_nTelPlanes, true);
  
  if( _BadDetectorIDs.size() ) {
    streamlog_out ( MESSAGE3 ) <<  _BadDetectorIDs.size() << " detectors are excluded!!" << endl;
  }
    
  for(int iBad=0; iBad<(int)_BadDetectorIDs.size(); iBad++) {
    int ipl = _BadDetectorIDs[iBad];
    _isActive[ipl] = false; 
    _nActivePlanes--;
  }
   
  // Print out geometry information
  streamlog_out ( MESSAGE3 ) <<  "Telescope configuration with " << _nTelPlanes << " planes" << endl;
   
  for(int ipl=0; ipl < _nTelPlanes; ipl++) {
    stringstream ss ;
      
    if(_isActive[ipl]) {
        ss << "Active  plane" ;
    } else {
        ss << "Passive plane" ;
    }
    ss << "  ID = " << _detector.GetDet(ipl).GetDAQID()
       << "  at Z [mm] = " << setprecision(3) <<  _detector.GetDet(ipl).GetNominal().GetPosition()[2]; 
      
    streamlog_out( MESSAGE3 ) <<  ss.str() << endl;
  }
    
  // Check consistency of track selection
  if ( _minHits < 0 ) _minHits = 0;
  if ( _maxGap < 0 ) _maxGap = 0;  
  
  // Deactivate hit removal in final fit
  if (_outlierChi2Cut <= 0) {
    _outlierChi2Cut = numeric_limits< double >::max();
  }
  
  // Basic sanity check for particle hypthesis 
  if (_mass <= 0) {
    streamlog_out ( MESSAGE3 ) <<  "Bad steering file: negative particle mass " << _mass << " set" << endl;  
    _mass = 0.000511; 
  }
  
  if (_charge == 0) {
    streamlog_out ( MESSAGE3 ) <<  "Bad steering file: zero charge set." << endl;  
    _charge = -1; 
  }
  
  if ( (int)_momentum_list.size() == 0 ) {
    streamlog_out ( MESSAGE3 ) <<  "Bad steering file: zero charge set." << endl;  
    _momentum_list.push_back(4.0); 
  }
  
  if ( (int)_maxResidualU.size() == 0 ) {
    _maxResidualU.resize(_nTelPlanes, 1.0);
    streamlog_out ( MESSAGE3 ) <<  "Bad steering file: Parameter maxResidualU not set. Using default value " <<  _maxResidualU[0] << "mm." << endl;  
  } else if ( (int)_maxResidualU.size() != _nTelPlanes ) {
    double max = _maxResidualU[0];
    _maxResidualU.resize(_nTelPlanes, max);  
    streamlog_out ( MESSAGE3 ) <<  "Bad steering file: Paremeter maxResidualU not set for all planes. Appending with value " 
                               << _maxResidualU[0] << "mm." << endl; 
  }
  
  if ( (int)_maxResidualV.size() == 0 ) {
    _maxResidualV.resize(_nTelPlanes, 1.0);  
    streamlog_out ( MESSAGE3 ) <<  "Bad steering file: Parameter maxResidualU not set. Using default value " <<  _maxResidualV[0] << "mm." << endl;  
  } else if ( (int)_maxResidualV.size() != _nTelPlanes ) { 
    double max = _maxResidualV[0];
    _maxResidualV.resize(_nTelPlanes, max);  
    streamlog_out ( MESSAGE3 ) <<  "Bad steering file: Paremeter maxResidualV not set for all planes. Appending with value " 
                               << _maxResidualV[0] << "mm." << endl;
  }
  
  streamlog_out ( MESSAGE3 ) <<  "Use residual cuts for " << _nTelPlanes << " planes" << endl;
   
  for(int ipl=0; ipl < _nTelPlanes; ipl++) {
    stringstream ss ;
      
    if(_isActive[ipl]) {
        ss << "Active  plane" ;
    } else {
        ss << "Passive plane" ;
    }
    ss << "  ID = " << _detector.GetDet(ipl).GetDAQID()
       << "  maxU/mm= " << setprecision(3) <<  _maxResidualU[ipl]
       << "  maxV/mm= " << setprecision(3) <<  _maxResidualV[ipl]; 
      
    streamlog_out( MESSAGE3 ) <<  ss.str() << endl;
  }
       
}  

//
// Method called for each run
//
void FastTracker::processRunHeader(LCRunHeader * run)
{
// Print run number
   streamlog_out(MESSAGE3) << "Processing run: "
                           << (run->getRunNumber())
                           << std::endl << std::endl;
   
   _nRun++ ;
}

//
// Method called for each event
//
void FastTracker::processEvent(LCEvent * evt)
{
   
   // Print event number
   if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                   << (evt->getEventNumber())
                                                                   << std::endl << std::endl;
    
   
   streamlog_out(MESSAGE2) << "Events processed: " << (evt->getEventNumber())
                                                   << std::endl << std::endl;
   
   _nEvt ++ ;
    
   // Read input hit collections
   // ===========================
   // 
   // Loop over all hits and copy hits to HitStore 
   
   // Manage hits in this events 
   HitFactory HitStore(_detector);   
   
   for ( size_t iCol = 0 ; iCol < _hitCollectionNameVec.size(); ++iCol ) {
     
     try {
       
       LCCollectionVec * col = dynamic_cast < LCCollectionVec * > (evt->getCollection( _hitCollectionNameVec.at( iCol ) ) );
       
       for ( int iExt = 0 ; iExt < (int) col->size() ; ++iExt ) {
                 
         // Built a TBHit
         TrackerHitImpl * lciohit = dynamic_cast< TrackerHitImpl* > ( col->getElementAt( iExt ) );
         TBHit RecoHit ( lciohit  );        
         
         // We have to find plane number of the hit 
         int daqid = RecoHit.GetDAQID();
         int ipl = _detector.GetPlaneNumber(daqid);
            
         if (ipl == -99 ) continue; 

         // Ignore hit, iff plane not declared as active (i.e. plane excluded by user)
         if(! _isActive[ipl]) 
         { 
           streamlog_out ( MESSAGE2 ) << " ignore sensor!! " << endl;
           continue ;
         }
          
         // Ignore hit, iff hit quality is bad and quality filter active  
         if ( !_hitQualitySelect && RecoHit.GetQuality() != _hitQualitySelect )   
         { 
           streamlog_out ( MESSAGE2 ) << " bad quality hit on plane " << ipl << endl;
           continue ;
         }   
         
         streamlog_out(MESSAGE2) << "Use hit at plane " << ipl << ": "
                                 << RecoHit.GetCoord() 
                                 << RecoHit.GetCov()  
                                 << endl;
           
         // Add hit to store
         HitStore.AddRecoHit(RecoHit);  
                    
       } // End for hit loop      
          
     } catch (lcio::DataNotAvailableException& e) {
       streamlog_out ( MESSAGE2 ) << "Not able to get collection "
                                  << _hitCollectionNameVec.at( iCol )
                                  << "\nfrom event " << evt->getEventNumber()
                                  << " in run " << evt->getRunNumber()  << endl;
       
     }  
   }   
   
   if( HitStore.GetNHits() == 0 )  ++_noOfEventWOInputHit; 
   
   streamlog_out ( MESSAGE2 ) << "Total of " << HitStore.GetNHits() << " good hits." << endl;
   
   
   // Count number of firing planes 
   int firingPlanes=0;     
   for(int ipl=0;ipl<_nTelPlanes;ipl++)  { 
     if (_isActive[ipl]  && HitStore.GetNHits(ipl) > 0) 
       firingPlanes++;
   }
   
   if (firingPlanes>= _minHits) _noOfEventMinHits++;
   
   
   // Combinatorial Track Finder   
   //=========================================================
   // Track candidates are seeded from pairs of hits from two 
   // different sensor planes. A hit pair defines a reference 
   // trajectory (or seed track) for searching compatible hits.
   // 
   // All seed tracks are followed through the entire telescope. 
   // At each sensor plane, we add the best compatible hit to 
   // the seed track. Compatible hits lie in a rectangular field
   // of length MaxResidual around the seed track. In case multiple  
   // hits are compatible, the track finder selects always the 
   // closest hit. Hits are presorted to speed up the process of 
   // building candidate tracks from seed tracks. 
   // 
   // After all hits are assigned to a seed track, a Kalman Filter 
   // based track fit is performed. If the fit quality is good, the 
   // track candidate is stored in the TrackCollector. 
   
   // Container for candidate tracks 
   std::list<TBTrack> TrackCollector; 
   
   // First track finder pass      
   //=========================================================
   // Seed tracks are constructed from hits in two planes. 
   // Track seeds are followed in the beam direction.
   
   // Seed forward Kalman filter using hits on two sensor planes on upstream side
   if ( _firstPass_firstPlane != _firstPass_secondPlane) {
     
     findTracks(TrackCollector , HitStore , _firstPass_firstPlane , _firstPass_secondPlane,+1); 
     
     streamlog_out ( MESSAGE2 ) << "Total of " << TrackCollector.size() << " forward candidate tracks found" << endl;
   }
   
   // Seed backward Kalman filter using hits in two sensor planes on downstream side
   if ( _secondPass_firstPlane != _secondPass_secondPlane) {
     
     findTracks(TrackCollector , HitStore ,  _secondPass_firstPlane  , _secondPass_secondPlane,-1); 
     
     streamlog_out ( MESSAGE2 ) << "Total of " << TrackCollector.size() << " candidate tracks found" << endl;
   }
   
   // Single hit finding     
   //=========================================================
   // Seed forward Kalman filter using hits from a single seed plane and propagate along beam axis 
   // along the z axis. 
   for ( int seedplane : _singleHitSeedingPlanes ) {
     
     findTracks(TrackCollector , HitStore ,  seedplane); 
     
     streamlog_out ( MESSAGE2 ) << "Total of " << TrackCollector.size() << " candidate tracks found" << endl;
   } 
   
   // Final Track Selection  
   // ========================================================
   // So far, the TrackCollector may hold candidate tracks with 
   // one or more identical hits. In other words, the assignment 
   // of at least one hit to a track candidate is ambigiuous.
   // The ambiguity will be resolved by selecting the candidate 
   // with more hits. If candidates have the same number of hits, 
   // the track chi2 is used as a second criterion. Clearly, this 
   // makes sense if there are at least three active layers. 
   //  
   // After track selection, the TrackCollector holds all final
   // tracks.    
    
   
   if ( _nActivePlanes>=3 ) {
     
     // Ok, let's sort tracks hypotheses   
     TrackCollector.sort(compare_tracks);   
     
     // Loop over all tracks    
     for(list<TBTrack>::iterator ctrack=TrackCollector.begin(); ctrack!=TrackCollector.end(); ++ctrack) 
     {
       // We keep the current track candidate, but we delate all other candidates with overlap. 
       list<TBTrack>::iterator otrack = ctrack;
       ++otrack; 
       while ( otrack != TrackCollector.end() ) {   
         // Delete incompatible track
         if (  check_incompatible(*ctrack,*otrack)  ) { 
           otrack = TrackCollector.erase(otrack);
           streamlog_out ( MESSAGE2 ) << "   erase track" << endl; 
         } else {
           ++otrack;
         } 
       } 
     }
   }
    
   
   streamlog_out ( MESSAGE2 ) << "Total of " << TrackCollector.size() << " tracks found" << endl;
   
   // Store final tracks to LCIO
   // =========================================================== 
   // Accepted tracks are stored in collection of type LCIO::Tracks.
   // LCIO Tracks are linked to original measured hits contributing 
   // to the track.    
      
   // Create output track collection
   LCCollectionVec *fittrackvec = new LCCollectionVec(LCIO::TRACK);
   
   TrackInputProvider TrackLCIOWriter; 
   
   // Set flag for storing track hits in track collection
   LCFlagImpl flag(fittrackvec->getFlag());
   flag.setBit( LCIO::TRBIT_HITS );
   fittrackvec->setFlag(flag.getFlag());
   
   // Count stored tracks
   int nStoredTracks=0;

   // Remember the hit ids of used hits
   vector<vector<int>> usedIDs;
   usedIDs.resize(_nTelPlanes, vector<int>(0, 0));

   
   for(list<TBTrack>::iterator ctrack=TrackCollector.begin(); ctrack!=TrackCollector.end(); ++ctrack) 
   {
      
     // Print the track candidate 
     streamlog_out(MESSAGE2) << "Printing final track: " << std::endl;  
     for(int ipl=0;ipl<_nTelPlanes;ipl++) {
       if ( (*ctrack).GetTE(ipl).HasHit() ) {
         streamlog_out(MESSAGE2) << " at plane " << ipl << ": "
                                 << (*ctrack).GetTE(ipl).GetHit().GetCoord() 
                                 << endl;
       } 
     } 
     
     // Mark all hits used in this track
     mark_hits( *ctrack, usedIDs ); 
     
     // Convert TBTrack to LCIO::Track  
     TrackImpl* lciotrack = TrackLCIOWriter.MakeLCIOTrack( *ctrack );
       
     // Add lcio track to lcio collection 
     fittrackvec->addElement(lciotrack);
      
     // Increment the total track counter      
     _noOfTracks++;
     // Incremnent track in this event
     nStoredTracks++;
   }
   
   if ( nStoredTracks == 0 ) 
   {
     ++_noOfEventWOTrack;
   }
     
   evt->addCollection(fittrackvec,_trackCollectionName); 
   
   // Create collection of unused hits 
   // ============================================================= 
   // Unused hits are input hits not being part of any track in the 
   // final track selection. Putting them in a seperate collection 
   // may be helpfull for estimating track finding efficiency. 
   
   // Create collection for unused hits 
   LCCollectionVec * notUsedHitCollection = new LCCollectionVec(LCIO::TRACKERHIT) ;
   
   if (true) {  
     for(int ipl=0;ipl<_nTelPlanes;ipl++)  { 
       for (int ihit = 0; ihit < HitStore.GetNHits(ipl); ihit++ ) {  
         // Get hit    
         TBHit& hit = HitStore.GetRecoHitFromID(ihit, ipl); 
         // Check hit was not used
         if ( std::find(usedIDs[ipl].begin(), usedIDs[ipl].end(), hit.GetUniqueID() ) == usedIDs[ipl].end() )
         {
           // Print the track candidate 
           streamlog_out(MESSAGE2) << "Printing at plane " << ipl << "  not used hit: " << hit.GetCoord() 
                                   << endl;
           
           notUsedHitCollection->push_back( hit.MakeLCIOHit() );
         }
       }
     
     }
   }
    
   // Store collection of not used hits 
   evt->addCollection( notUsedHitCollection, _notUsedhitCollectionName );   

   return;
}


//
// Method called after each event to check the data processed
//
void FastTracker::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void FastTracker::end()
{ 
  
  // Print the summer
  streamlog_out(MESSAGE3) << "Total number of processed events:     " << setw(10) << setiosflags(ios::right) << _nEvt << resetiosflags(ios::right) << endl
                          << "Total number of event w/o input hit:  " << setw(10) << setiosflags(ios::right) << _noOfEventWOInputHit 
                          << resetiosflags(ios::right) << endl
                          << "Total number of event w/o track:      " << setw(10) << setiosflags(ios::right) << _noOfEventWOTrack
                          << resetiosflags(ios::right) << endl
                          << "Total number of reconstructed tracks: " << setw(10) << setiosflags(ios::right) << _noOfTracks << resetiosflags(ios::right)
                          << resetiosflags(ios::right) << endl
                          << "Number of events with " << _minHits << " firing planes: " << setw(9) << setiosflags(ios::right) << _noOfEventMinHits << resetiosflags(ios::right) << endl
                          << "Number of failed final fits: " << setw(10) << setiosflags(ios::right) <<  _noOfFailedFits << resetiosflags(ios::right) << endl
                          << "Number of ambiguous hits " << _noOfAmbiguousHits << " for number of cand tracks " << _noOfCandTracks << resetiosflags(ios::right) << endl                     
                          << endl; 

   
  // CPU time end
  _timeCPU = clock()/1000 - _timeCPU;
  
  streamlog_out(MESSAGE3) << std::endl
                          << " "
                          << "Time per event: "
                          << std::setiosflags(std::ios::fixed | std::ios::internal )
                          << std::setprecision(3)
                          << _timeCPU/_nEvt
                          << " ms"
                          << std::endl
                          << std::setprecision(3)
                          << std::endl
                          << " "
                          << "Processor succesfully finished!"
                          << std::endl;
   
  
  
  
}

//
// Called by the processEvent() to add tracks to trackcollector using hits in hitstore
//
void FastTracker::findTracks( std::list<TBTrack>& TrackCollector , HitFactory& HitStore , int firstplane , int secondplane, int idir) 
{
   // Check plane numbers are valid
   if (firstplane < 0 || secondplane < 0) return; 
   if (firstplane >= _nTelPlanes || secondplane >= _nTelPlanes) return; 
   
   streamlog_out ( MESSAGE2 ) << "First active plane is " << firstplane  << " and has " 
                              << HitStore.GetNHits(firstplane) << " good hits." << endl;
   streamlog_out ( MESSAGE2 ) << "Second active plane is " << secondplane << " and has " 
                              << HitStore.GetNHits(secondplane) << " good hits." << endl;
   
   // Scan momentum hypotheses 
   for (int imom = 0; imom < (int) _momentum_list.size() ; imom++ ) {
     
     double my_momentum = _momentum_list[imom];
     double my_charge = _charge; 
     
     // This case refers to the antiparticle
     if ( my_momentum  < 0 ) {
       my_momentum*=-1; 
       my_charge *=-1; 
     } 
     
     // Configue seed track generator  
     SeedGenerator TrackSeeder(my_charge,my_momentum);
     
     for (int ihit = 0; ihit < HitStore.GetNHits(firstplane); ihit++ ) {
       for (int jhit = 0; jhit < HitStore.GetNHits(secondplane); jhit++ ) {
       
         // Init new candidate track  
         TBTrack trk(_detector);
         trk.SetMass( _mass );
         trk.SetCharge( my_charge );
         trk.SetMomentum( my_momentum ); 
         
         // Compute seed track 
         TBHit& firsthit = HitStore.GetRecoHitFromID(ihit, firstplane);  
         TBHit& secondhit = HitStore.GetRecoHitFromID(jhit, secondplane);
         TBTrackState Seed = TrackSeeder.CreateSeedTrack(firsthit, secondhit, _detector);   
         
         // Skip all candidate tracks with a very 
         // large angle relative to Z axis. 
         if ( _maxSlope > 0 ) {         
           double dxdz = Seed.GetPars()[0][0];
           double dydz = Seed.GetPars()[1][0]; 
           if ( ( std::abs(dxdz) > _maxSlope ) || ( std::abs(dydz) > _maxSlope )  ) 
           {
             streamlog_out ( MESSAGE2 ) << "Slope too large. Skipping track candidate!! " << endl; 
             continue; 
           }
         }  

         trk.SetReferenceState(Seed);
         
         buildTrackCand(trk, HitStore, TrackCollector,idir);
         
       } // End track seeding
     } // End track seeding  
   } // End momentum scan 
}


//
// Called by the processEvent() to add tracks to trackcollector using hits in hitstore
//
void FastTracker::findTracks( std::list<TBTrack>& TrackCollector , HitFactory& HitStore , int seedplane) 
{
   // Check plane numbers are valid
   if (seedplane < 0 || seedplane >= _nTelPlanes) return; 
     
   // Compute direction of filter: forward (+1) or backward (-1)
   int idir = 1;
   if (seedplane > _nTelPlanes/2 ) idir = -1;
   
   streamlog_out ( MESSAGE2 ) << "Seed plane is " << seedplane  << " and has " 
                              << HitStore.GetNHits(seedplane) << " good hits. Filter runs in direction " << idir << endl;
    
   // Loop over different momentum hypothesis 
   
   for (int imom = 0; imom < (int) _momentum_list.size() ; imom++ ) {
     
     double my_momentum = _momentum_list[imom];
     double my_charge = _charge; 
     
     // This case refers to the antiparticle
     if ( my_momentum  < 0 ) {
       my_momentum*=-1; 
       my_charge *=-1; 
     } 
     
     // Configue seed track generator  
     SeedGenerator TrackSeeder(my_charge,my_momentum);
     
     for (int ihit = 0; ihit < HitStore.GetNHits(seedplane); ihit++ ) {
       
       // Init new candidate track  
       TBTrack trk(_detector);
       trk.SetMass( _mass );
       trk.SetCharge( my_charge );
       trk.SetMomentum( my_momentum ); 
         
       // Compute seed track 
       TBHit& seedhit = HitStore.GetRecoHitFromID(ihit, seedplane);  
       TBTrackState Seed = TrackSeeder.CreateSeedTrack(seedhit, _detector); 
       
       // Skip all candidate tracks with a very 
       // large angle relative to Z axis. 
       if ( _maxSlope > 0 ) {         
         double dxdz = Seed.GetPars()[0][0];
         double dydz = Seed.GetPars()[1][0]; 
         if ( ( std::abs(dxdz) > _maxSlope ) || ( std::abs(dydz) > _maxSlope )  ) 
         {
           streamlog_out ( MESSAGE2 ) << "Slope too large. Skipping track candidate!! " << endl; 
           continue; 
         }
       }         
       trk.SetReferenceState(Seed);
       
       buildTrackCand(trk, HitStore, TrackCollector, idir);
     }
   }      
}

void FastTracker::buildTrackCand(TBTrack& trk, HitFactory& HitStore, std::list<TBTrack>& TrackCollector, int idir)
{
  
  // Configure Kalman track fitter
  TBKalmanB TrackFitter(_detector);
   
  // Extrapolate seed to all planes
  bool exerr = TrackFitter.ExtrapolateSeed(trk);
  if ( exerr ) { // just skip the track 
    _noOfFailedFits++;
    streamlog_out ( MESSAGE2 ) << "Extrapolation of track seed into telescope failed!!" << endl; 
    return;
  }  
         
  // Add hits to candidate track  
  //============================
       
  // Counts consecutive missing hits
  int ngap = 0 ;    
         
  // Count number of added hits 
  int nhits = 0; 
               
  // Count hit ambiguities per track  
  int nAmbiguousHits = 0;   
         
  // This is the initial track state vector
  HepMatrix x0(5,1,0);
         
  // This is the initial covariance matrix 
  HepSymMatrix C0(5,0);
  for (int i = 0; i < 2; i++) {
    C0[i][i] = 1E-2;
    C0[i+2][i+2] = 1E1; 
  } 
  C0[4][4] = 1;  
          
  double finderChi2 = 0; 
  
  int istart, istop;   
  if (idir > 0) {
    istart = 0; 
    istop = _nTelPlanes; 
  } else {
    istart = _nTelPlanes-1; 
    istop = -1; 
  } 
  
  // Follow track along beam direction  
  for(int ipl=istart; ipl!=istop; ipl+=idir) { 
            
    // This is the reference state on the current track element
    HepMatrix& xref = trk.GetTE(ipl).GetState().GetPars();
                     
    // Skip passive sensors    
    if( _isActive[ipl] && trk.GetTE(ipl).IsCrossed() ) { 
                  
      // This is the filtered local track state using all hits 
      // encountered so far 
      HepMatrix x = xref + x0;
             
      // Get extrapolated intersection coordinates
      double u = x[2][0]; 
      double v = x[3][0]; 
             
      // Fast preselection of hit candidates compatible to   
      // predicted intersection coordinates. 
      vector<int> HitIdVec = HitStore.GetCompatibleHitIds(ipl, u, v, _maxResidualU[ipl], _maxResidualV[ipl]);
             
      // Now, we select the best hit candidate 
      int ncandhits = HitIdVec.size();
      int besthitid=-1;
      double bestdist = numeric_limits< double >::max();
             
      for (int icand = 0; icand < ncandhits; ++icand ) 
      {           
        // Get reco hit at plane ipl 
        int hitid = HitIdVec[icand];
        TBHit & RecoHit = HitStore.GetRecoHitFromID(hitid, ipl);   
        double uhit = RecoHit.GetCoord()[0][0];              
        double vhit = RecoHit.GetCoord()[1][0]; 
       
        // Discard hits with too large residuals
        if ( std::abs(u - uhit) >= _maxResidualU[ipl] && _maxResidualU[ipl] > 0) continue; 
        if ( std::abs(v - vhit) >= _maxResidualV[ipl] && _maxResidualV[ipl] > 0) continue; 
               
        // This hit could(?) belong to the track 
        nAmbiguousHits++;
               
        // Remember hit with smallest residual 
        double hitdist = 0; 
        if ( _maxResidualU[ipl] > 0 )  hitdist += std::abs(u - uhit); 
        if ( _maxResidualV[ipl] > 0 )  hitdist += std::abs(v - vhit); 
       
        if( hitdist < bestdist )
        {
          bestdist = hitdist;
          besthitid=hitid;
        }
      } 
              
      if ( besthitid!=-1 )  {
        // Add closest hit to candidate track
        // This is a simple greedy selection and may be wrong in 
        // there are hit ambiguities or the reference track is 
        // too bad.
               
        TBHit& BestHit = HitStore.GetRecoHitFromID(besthitid, ipl);
               
        if ( TrackFitter.GetPredictedChi2(x, C0, BestHit) < _outlierChi2Cut ) {
          double hitchi2 = TrackFitter.FilterHit(BestHit, xref, x0, C0);
          BestHit.SetUniqueID(besthitid);             
          trk.GetTE(ipl).SetHit(BestHit);               
          // Some bookkeeping   
          finderChi2 += hitchi2; 
          nhits++;  
          ngap = 0;
        }          
      } else {
        // No matching hit found on this sensor
        ngap++;    
        if( ngap > _maxGap ) {
          streamlog_out(MESSAGE1) << "Too many missing hits. Skip seed track! " << endl;  
          return;       
        }          
      }
    } // End is active
           
    // Extrapolate filtered state to next track element 
    int inext = ipl+idir;
    if (inext!=istop)  {
      HepMatrix& nxref = trk.GetTE(inext).GetState().GetPars();
      exerr = TrackFitter.PropagateState(trk.GetTE(ipl), trk.GetTE(inext), xref, nxref, x0, C0); 
      if ( exerr ) { // just skip the track 
        _noOfFailedFits++;
        streamlog_out ( MESSAGE2 ) << "Extrapolation of track seed into telescope failed!!" << endl; 
        return;
      }  
    }
  
  } // End sensor loop  
                                   
  // Reject track candidate if number of hits too small
  if ( nhits<_minHits ) {
    streamlog_out ( MESSAGE1 ) << "Number of hits too small. Skipping candidate track!" << endl;
    return;
  }
         
  // Reject hit if total chisq gets too large
  if( finderChi2>_maxTrkChi2  ||  std::isnan(finderChi2)  || finderChi2 < 0 ) { 
    streamlog_out ( MESSAGE1 ) << "Bad chisq. Skipping track candidate!" << endl;
    return;      
  }
                        
  // update the running counters  
  _noOfCandTracks++;    
  _noOfAmbiguousHits += nAmbiguousHits;
              
  // Ok, we keep this track for final selection
  trk.SetChiSqu(finderChi2);
  TrackFitter.SetNdof(trk);      
  TrackCollector.push_back( trk ); 
  
  return;
}

//
// Method printing processor parameters
//
void FastTracker::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "FastTracker Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}

// 
// Check if tracks trk1 and trk2 have common hits
// 
bool check_incompatible( TBTrack& trk1, TBTrack& trk2 )
{
  // Get all track elements (TEs) in trk1
  std::vector<TBTrackElement>& TEVec1 = trk1.GetTEs();
  int nTE1 = (int) TEVec1.size();  
  
  // Get all track elements (TEs) in trk2
  std::vector<TBTrackElement>& TEVec2 = trk2.GetTEs();
  int nTE2 = (int) TEVec2.size();
  
  if (nTE1!=nTE2) {  
    return false; 
  }  
  
  // Loop over track elements 
    
  for(int iTE=0;iTE<nTE1;++iTE) {
     
    // Check both tracks have hits
    if ( TEVec1[iTE].HasHit() && TEVec2[iTE].HasHit() ) {
      
      // Get unique hitId for comparison
      int hitId1 = TEVec1[iTE].GetHit().GetUniqueID();          
      int hitId2 = TEVec2[iTE].GetHit().GetUniqueID();       
      
      // Check if reco(!) hits are same  
      if ( hitId1==hitId2 && hitId1!=-1 && hitId2!=-1 ) {
        return true; 
      }  
    } 	 
  }   
  return false; 
}
   
// 
// Is track trk1 better than trk2
// 
bool compare_tracks ( TBTrack& trk1, TBTrack& trk2 )
{
  // More hits are better 
  if ( trk1.GetNumHits() > trk2.GetNumHits() ) {
     return true;
  } else if ( trk1.GetNumHits() < trk2.GetNumHits() ) {
     return false; 
  } else { 
     // Smaller chisq is better
     if (trk1.GetChiSqu() <= trk2.GetChiSqu() ) {
        return true; 
     } else {
        return false; 
     }
  }
}


// 
// mark hits in track as used
// 
void mark_hits ( TBTrack& trk, vector<vector<int>>&  usedIDs )
{ 
  // Loop over track elements 
  for(TBTrackElement& te : trk.GetTEs() ) {     
    // Check both tracks have hits
    if ( te.HasHit() ) {
      // Get unique hitId 
      int hitId = te.GetHit().GetUniqueID();
      int ipl = te.GetDet().GetPlaneNumber();
      if (hitId >= 0) usedIDs[ipl].push_back(hitId);   
    } 	 
  }   
  return; 
}


} // Namespace

