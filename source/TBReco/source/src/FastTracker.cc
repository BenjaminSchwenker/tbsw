// FastTracker implementation file
//                       
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


// TBTools includes 
#include "FastTracker.h"
#include "SeedGenerator.h"
#include "GenericTrackFitter.h"
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
   
   registerProcessorParameter("MaxResidual",  
                             "Maximum hit residual [mm]",                       
                             _maxResidual, static_cast <float> (0.5) );
    
   registerProcessorParameter("OutlierChi2Cut",  
                             "Chi2 cut for removal of bad hits",                       
                             _outlierChi2Cut, static_cast <float> (50) ); 

   registerProcessorParameter("OutlierIterations",  
                             "Fitting iteration for removal of bad hits",                       
                             _outlierIterations, static_cast <int> (0) ); 
   
   registerProcessorParameter ("MinimumHits",
                              "Minimum number of hits in track",
                              _minHits,  static_cast < int > (2));

   registerProcessorParameter ("PassOne_FirstPlane",
                              "Build track seeds from first plane and second plane",
                              _firstPass_firstPlane,  static_cast < int > (-1));
  
   registerProcessorParameter ("PassOne_SecondPlane",
                              "Build track seeds from first plane and second plane",
                              _firstPass_secondPlane,  static_cast < int > (-1));

   registerProcessorParameter ("PassTwo_FirstPlane",
                              "Build track seeds from first plane and second plane",
                              _secondPass_firstPlane,  static_cast < int > (-1)); 

   registerProcessorParameter ("PassTwo_SecondPlane",
                              "Build track seeds from first plane and second plane",
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
    _outlierIterations = 0; // no iteration, one pass 
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
           //cout << " ignore sensor!! " << endl;
           continue ;
         }
          
         // Ignore hit, iff hit quality is bad and quality filter active  
         if ( !_hitQualitySelect && RecoHit.GetQuality() != _hitQualitySelect )   
         { 
           //cout << " bad quality hit on plane " << ipl << endl;
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
   // different subdetectors. A hit pair defines a reference 
   // trajectory (or seed track) for searching compatible hits
   // on subdetectors. Tracks are seeded from both ends of 
   // the telescope to increase the track finding efficiency.
   // 
   // All seed tracks are followed through the entire detector. 
   // At each subdetector, we add the best compatible hit to 
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
   
   // We need at least two active planes to seed a straight track
   if ( _firstPass_firstPlane != _firstPass_secondPlane) {
     
     findTracks(TrackCollector , HitStore , _firstPass_firstPlane , _firstPass_secondPlane); 
    
     streamlog_out ( MESSAGE2 ) << "Total of " << TrackCollector.size() << " forward candidate tracks found" << endl;
   }
   
   // Second track finder pass      
   //=========================================================
   // Seed tracks are constructed from hits in two planes. 
   // Track seeds are followed in the beam direction.
   
  
   // We need at least three active planes to find new straight track; otherwise, we would just repeat forward pass;)
   if ( _secondPass_firstPlane != _secondPass_secondPlane) {
     
     findTracks(TrackCollector , HitStore ,  _secondPass_firstPlane  , _secondPass_secondPlane ); 
     
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
         if ( check_incompatible(*ctrack,*otrack) ) { 
           otrack = TrackCollector.erase(otrack);
           streamlog_out ( MESSAGE1 ) << "   erase track" << endl; 
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
   
   int nStoredTracks=0;
   
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
                          << "Number of events with " << _minHits << " firing planes: " << setw(9) << setiosflags(ios::right) << _noOfEventMinHits << resetiosflags(ios::right)
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
void FastTracker::findTracks( std::list<TBTrack>& TrackCollector , HitFactory& HitStore , int firstplane , int secondplane) 
{
   
   streamlog_out ( MESSAGE2 ) << "First active plane is " << firstplane  << " and has " 
                              << HitStore.GetNHits(firstplane) << " good hits." << endl;
   streamlog_out ( MESSAGE2 ) << "Second active plane is " << secondplane << " and has " 
                              << HitStore.GetNHits(secondplane) << " good hits." << endl;
   
   // Check plane numbers are valid
   if (firstplane < 0 || secondplane < 0) return; 
 
   // Configure Kalman track fitter
   GenericTrackFitter TrackFitter(_detector);
   TrackFitter.SetNumIterations(_outlierIterations+1);
   TrackFitter.SetOutlierCut(_outlierChi2Cut); 

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
         trk.SetReferenceState(Seed);
         
         // Option: Skip all candidate tracks with a very 
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
       
         // Extrapolate seed to all planes
         bool exerr = TrackFitter.ExtrapolateSeed(trk);
         if ( exerr ) { // just skip the track 
           continue;
         }  
         
         // Add hits to candidate track  
         //============================
       
         // Counts consecutive missing hits
         int ngap = 0 ;    
         
         // Count number of added hits 
         int nhits = 0; 
       
         // Flag success of track building
         bool isValid = true;   
                              
         // Follow track along beam direction    
         for(int ipl=0;ipl<_nTelPlanes;ipl++)  { 
                
           // Skip passive sensors    
           if( _isActive[ipl] && trk.GetTE(ipl).IsCrossed() ) { 
                  
             // Get extrapolated intersection coordinates
             double u = trk.GetTE(ipl).GetState().GetPars()[2][0]; 
             double v = trk.GetTE(ipl).GetState().GetPars()[3][0]; 
             
             // Fast preselection of hit candidates compatible to   
             // predicted intersection coordinates. 
             vector<int> HitIdVec = HitStore.GetCompatibleHitIds(ipl, u, v, _maxResidual);
             
             // Now, we select the best hit candidate 
             int ncandhits = HitIdVec.size();
             int besthitid=-1;
             double bestdist = numeric_limits< double >::max();
             
             for (int icand = 0; icand < ncandhits; ++icand ) 
             {    
             
               // Get reco hit at plane ipl 
               int hitid = HitIdVec[icand];
               TBHit & RecoHit = HitStore.GetRecoHitFromID(hitid, ipl);   
               
               // Calculate hit2track distance
               double uhit = RecoHit.GetCoord()[0][0];              
               double vhit = RecoHit.GetCoord()[1][0]; 
               double hitdist = std::abs( u - uhit ) + std::abs( v - vhit ) ;
                     
               if( hitdist < bestdist )
               {
                 bestdist = hitdist;
                 besthitid=hitid;
               }
             } 
              
             // Check iff good hit found 
             if ( besthitid!=-1 && bestdist< 2*_maxResidual )  {
              
               // Add hit to candidate track 
               TBHit& BestHit = HitStore.GetRecoHitFromID(besthitid, ipl);
               BestHit.SetUniqueID(besthitid);            
               trk.GetTE(ipl).SetHit(BestHit);
                 
               // Some bookkeeping    
               nhits++;  
               ngap = 0;
             
             } else {
           
               // Ok, missing hit for that seed track 
               ngap++;    
               if( ngap > _maxGap ) {
                 streamlog_out(MESSAGE1) << "Too many missing hits. Skip seed track! " << endl;  
                 isValid = false;
                 break;       
               }
              
             }
           } // End is active 
          
         } // End subdetector loop 
                
         // Check if track candidate is valid, otherwise proceed to 
         // next. 
         // =============================
       
         if( !isValid ) {
           continue;
         }  
       
         if ( nhits<_minHits ) {
           continue;
         }
         
         // Fit track candidate 
         bool trkerr = TrackFitter.Fit(trk); 
         if( trkerr ) {
           streamlog_out ( MESSAGE1 ) << "Fit failed. Skipping candidate track!" << endl;
           continue;
         }
         
         // Reject track candidate if number of hits too small
         // Can change in outlier rejection of fitter
         if(  trk.GetNumHits() < _minHits  ) { 
           streamlog_out ( MESSAGE1 ) << "Number of hits too small. Skipping candidate track!" << endl;
           continue;      
         }   
         
         // Reject hit if total chisq gets too large
         if(  trk.GetChiSqu() >= _maxTrkChi2  ) { 
           streamlog_out ( MESSAGE1 ) << "Track chisq too big. Skipping candidate track!" << endl;
           continue;      
         }
              
         // Ok, we keep this track for final selection      
         TrackCollector.push_back( trk ); 
       } 
     } // End track seeding  
   } // End momentum scan 
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
    //cout << "ERR: (CTRACK) At least one track is not complete." << endl;    
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

} // Namespace

