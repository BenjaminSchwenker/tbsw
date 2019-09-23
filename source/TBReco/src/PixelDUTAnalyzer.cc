// PixelDUTAnalyzer Processor  
// 
// See PixelDUTAnalyzer.h for full documentation of processor. 
// 
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


#include "PixelDUTAnalyzer.h"

// TBTools includes
#include "TBDetector.h"
#include "TBTrack.h"
#include "TBHit.h"
#include "PixelCluster.h"
#include "GenericTrackFitter.h"
#include "TrackInputProvider.h"
#include "Det.h"
#include "Utilities.h"

// Include basic C
#include <iostream>
#include <limits>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <string>
#include <algorithm>

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCGenericObjectImpl.h>

// Used namespaces
using namespace std; 
using namespace lcio ;
using namespace marlin ;

namespace depfet {


//
// Instantiate this object
//
PixelDUTAnalyzer aPixelDUTAnalyzer ;

//
// Constructor
//
PixelDUTAnalyzer::PixelDUTAnalyzer() : Processor("PixelDUTAnalyzer"),_inputDecodeHelper("")
{

// Processor description
  _description = "PixelDUTAnalyzer: DUT Analysis Processor" ;
   
//
// Input collections 
   
  registerInputCollection( LCIO::TRACK,
                           "TrackCollection" ,
                           "Name of telescope track collection"  ,
                           _trackColName ,
                           std::string("tracks") ) ;
   
  registerInputCollection( LCIO::TRACKERHIT,
                           "HitCollection" ,
                           "Name of DUT hit collection"  ,
                           _hitColName ,
                           std::string("hit") ) ;
    
  registerInputCollection( LCIO::TRACKERDATA,
                           "DigitCollection" ,
                           "Name of unpacked DUT digit collection"  ,
                           _digitColName ,
                           std::string("zsdata_dut") ) ;
  
  // Processor parameters
  
  registerProcessorParameter ("MetaInfoCollection",
                             "Name of optional unpacker meta info collection for DUT. Set empty string if it does not exist",
                             _metaDataColName  , static_cast< string > ( "" ) );  
  
  registerProcessorParameter ("DUTPlane",
                              "Plane number of DUT along the beam line",
                              _idut,  static_cast < int > (3));
  
  registerProcessorParameter ("ReferencePlane",
                              "Plane number of teference plane. Use only tracks having a hit on the reference plane to measure DUT efficiency. Put -1 to deactivate.",
                              _iref,  static_cast < int > (-1));
        
  registerProcessorParameter ("MaxResidualU",
                              "Maximum u residual for matching DUT hits to telescope track [mm]. Put -1 to deactivate cut.",
                              _maxResidualU,  static_cast < double > (0.2));
  
  registerProcessorParameter ("MaxResidualV",
                              "Maximum v residual for matching DUT hits to telescope track [mm]. Put -1 to deactivate cut.",
                              _maxResidualV,  static_cast < double > (0.2));
  
  registerProcessorParameter( "RootFileName",
                               "Output root file name",
                               _rootFileName, std::string("histos.root"));
   
   std::vector<std::string> initTestPixels;
   registerProcessorParameter ("DUTTestPixels", 
                               "List of pixels in format 'row1:col1 row2:col2 ...'. A flag 'hasTestPixels' is written to output tree if these pixels are found in an event",
                               _testPixels, initTestPixels); 
   
   

}

//
// Method called at the beginning of data processing
//
void PixelDUTAnalyzer::init() {
   
   // Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
   _timeCPU = clock()/1000 ;
   
   // Initialize all the counters
   _noOfEventWOInputHit = 0;
   _noOfEventWOInputTrack = 0;  
   _noOfTracks = 0;
   _noOfHits = 0;   
   _noOfMatchedTracks = 0;      
   _noOfHitsWOTrack = 0; 
   _noOfHitsWTrack = 0; 
   _iLastMatchedEvent = 0; 
   _iFirstMatchedEvent = -1; 

   if (_iref < 0 && _iref >= TBDetector::GetInstance().GetNSensors()) {
     streamlog_out ( MESSAGE3 )  << "Steering parameter reference plane has invalid value and will be ignored."
                                 << endl << endl;  
     _iref = -1;
   }
   
   // Print set parameters
   printProcessorParams();
   
   // Print out geometry information  
   streamlog_out ( MESSAGE3 )  << "D.U.T. plane  ID = " << TBDetector::Get(_idut).GetSensorID()
                               << "  at position = " << _idut 
                               << endl << endl;
    
      
   bookHistos();
   
}

//
// Method called for each run
//
void PixelDUTAnalyzer::processRunHeader(LCRunHeader * run)
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
void PixelDUTAnalyzer::processEvent(LCEvent * evt)
{
  
  // Print event number
  if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Event number is: "
                                                                << (evt->getEventNumber())
                                                                << std::endl << std::endl;
  
  streamlog_out(MESSAGE2) << "Event number is: " << (evt->getEventNumber())
                                                   << std::endl << std::endl;
  
  _nEvt ++ ;
  
  // Configure Kalman track fitter
  GenericTrackFitter TrackFitter(TBDetector::GetInstance());
  TrackFitter.SetNumIterations(1); 
     
  // Load DUT module    
  const Det & dut = TBDetector::Get(_idut);   
  
  // 
  //  Get DUT unpacker meta info, if exists
  //   
  bool isGoodEvent = true;    
  try {
    LCCollectionVec* metainfo = dynamic_cast < LCCollectionVec * > (evt->getCollection( _metaDataColName )) ;
    LCGenericObjectImpl* metaobj = dynamic_cast<LCGenericObjectImpl* > ( metainfo->getElementAt(0) );
    isGoodEvent = metaobj->getIntVal(0);            
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out(MESSAGE2) << " Meta info not available "
                            << endl << endl;
  }     
  
  streamlog_out(MESSAGE2) << "Unpacker flags isGoodEvent= " << isGoodEvent << endl; 
  
  //
  // Get digit collection 
  //
  int nDUTDigits = 0;   
  bool hasTestPixels = true;
  
  try {
    LCCollection* digitcol = evt->getCollection( _digitColName ) ;
    CellIDDecoder<TrackerDataImpl> DigitDecoder(digitcol, &_inputDecodeHelper);
    // Search for digits from DUT 
    for (int iDet = 0; iDet < digitcol->getNumberOfElements(); iDet++) {    
      TrackerDataImpl * digits = dynamic_cast<TrackerDataImpl* > ( digitcol->getElementAt(iDet) );
      int sensorID = DigitDecoder( digits ) ["sensorID"];
      if ( sensorID ==  dut.GetSensorID()  ) { 
        const FloatVec &pixVector = digits->getChargeValues();
         
        // Read the number of pixels on DUT in the event
        nDUTDigits = pixVector.size()/3;
               
        std::set<char> delims{':'};
        for (auto pixelstr: _testPixels) {
          auto tokens = splitpath(pixelstr, delims);   
          int test_col = std::stoi(tokens.at(1));
          int test_row = std::stoi(tokens.at(0));
          bool flag = false; 
          
          for (int iPix = 0; iPix < nDUTDigits; iPix++) {   
            int col = static_cast<int> (pixVector[iPix * 3]);
            int row = static_cast<int> (pixVector[iPix * 3 + 1]);
            if (col == test_col and row == test_row) {
              flag = true;
              break;
            } 
          }
          if (not flag) {
            hasTestPixels = false;
            break;
          }
        }
      }
    }
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out(MESSAGE2) << "Not able to get collection "
                            << _digitColName
                            << " from event " << evt->getEventNumber()
                            << " in run " << evt->getRunNumber()  << endl << endl;   
  }     
  
  streamlog_out(MESSAGE2) << "Total of " << nDUTDigits << " DUT digits in collection " << _digitColName << endl; 
     
  //
  // Get telescope track collection
  //
  
  LCCollection* trackcol = nullptr;
  int nTrack = 0;   
  try {
    trackcol = evt->getCollection( _trackColName ) ;
    nTrack = trackcol->getNumberOfElements();
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out(MESSAGE2) << "Not able to get collection "
                            << _trackColName
                            << " from event " << evt->getEventNumber()
                            << " in run " << evt->getRunNumber()  << endl << endl;   
  }  
   
  streamlog_out(MESSAGE2) << "Total of " << nTrack << " tracks in collection " << _trackColName << endl; 
  
  // 
  // Get DUT hit collection 
  // 
  
  LCCollection* hitcol = 0;
  int nHit = 0; 
  try {
    hitcol = evt->getCollection( _hitColName ) ;
    nHit = hitcol->getNumberOfElements();
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out(MESSAGE2) << "Not able to get collection "
                            << _hitColName
                            << " from event " << evt->getEventNumber()
                            << " in run " << evt->getRunNumber() << endl << endl;
  } 
 
  streamlog_out(MESSAGE2) << "Total of " << nHit << " hit(s) in collection " << _hitColName << endl;
  
  // Read telescope tracks and DUT hits 
  // ----------------------------------
  

  TrackStore.clear();
  HitStore.clear();
  // Get tracks from the lcio file 
  TrackInputProvider TrackLCIOReader;  
  int nTrackWithRefHit = 0;     

  for(int itrk=0; itrk< nTrack ; itrk++)
  {  
    // Retrieve track from LCIO 
    Track * lciotrk = dynamic_cast<Track*> (trackcol->getElementAt(itrk));
      
    // Convert LCIO -> TB track  
    TBTrack trk = TrackLCIOReader.MakeTBTrack( lciotrk, TBDetector::GetInstance() );  
    
    // Check that track has no hit on the DUT to avoid bias of residuals and efficiency
    if ( trk.GetTE(_idut).HasHit()  ) {
      streamlog_out ( MESSAGE3 ) << "Track has already hit on dut plane. Danger to bias final results" << endl;
    } 
    
    // Count tracks having hit on reference (timing) plane
    if (_iref >= 0) {
      if ( trk.GetTE(_iref).HasHit()  ) nTrackWithRefHit += 1;   
    }
    
    // Refit track in current alignment. Restrict smoothing to DUT plane.
    bool trkerr = TrackFitter.Fit(trk, _idut);
    if ( trkerr ) {
      streamlog_out ( MESSAGE2 ) << "Fit failed. Skipping track!" << endl;
      continue;
    } 
    
    // Compute intersection point of track with DUT (at w=0 plane)
    double u_trk = trk.GetTE(_idut).GetState().GetPars()[2];
    double v_trk = trk.GetTE(_idut).GetState().GetPars()[3];
    
    // Check that track intersects with the DUT sensitive volume 
    if (  dut.isPointOutOfSensor( u_trk, v_trk ) ) {
      streamlog_out(MESSAGE2) << "Ignore track which does not intersect the DUT."<< std::endl; 
      continue; 
    }
        
    TrackStore.push_back(std::move(trk));
  } // End track loop
    
  // Get hits from the lcio file   
  for(int ihit=0; ihit< nHit ; ihit++)
  {
    // Built a TBHit
    TrackerHitImpl * lciohit = dynamic_cast<TrackerHitImpl*>( hitcol->getElementAt(ihit) ) ;
    TBHit RecoHit ( lciohit  );        
         
    // We have to find plane number of the hit 
    int sensorid = RecoHit.GetSensorID();      
    int ipl = TBDetector::GetInstance().GetPlaneNumber(sensorid);  
      
    // Store all hits on the DUT module  
    if( dut.GetPlaneNumber() == ipl )
    {         
        streamlog_out(MESSAGE2) << " DUT hit at plane " << ipl 
                                << "   U [mm]= " << RecoHit.GetCoord()[0]
                                << "   V [mm]= " << RecoHit.GetCoord()[1] 
                                << endl;
      
        HitStore.push_back( std::move(RecoHit) );
    }
  } 
  
  streamlog_out(MESSAGE2) << "Total of " << TrackStore.size() << " good track(s)" << endl; 
  _noOfTracks += TrackStore.size();   
  
  
  streamlog_out(MESSAGE2) << "Total of " << HitStore.size() << " DUT hit(s)" << endl; 
  _noOfHits += HitStore.size(); 
  if ( HitStore.size() == 0)  ++_noOfEventWOInputHit;  
  
  // 
  // Match DUT hits with telescope track 
  int nMatch=0;	
  
  // Record for each hit a matched track 
  vector<int> hit2track(HitStore.size(), -1);
   
  // Record for each track a matched hit 
  vector<int> track2hit(TrackStore.size(), -1); 
  
  // Continue matching tracks and hits until all tracks are matched 
  // or no hit is close enough to a track!!

  { 
  double distmin=numeric_limits<double >::max();
  int besttrk=-1;
  int besthit=-1;   

  do{
    besttrk=-1;
    besthit=-1;
    
    distmin=numeric_limits<double >::max();
    
    // Find hit/track pair with minimum chi2 distance.  
    for(int itrk=0;itrk<(int)TrackStore.size(); itrk++)
    {
            
      // If matched, skip track 
      if (track2hit[itrk] >= 0) continue;  
      
      for(int ihit=0; ihit< (int)HitStore.size() ; ihit++)
      {
         
        // If matched, skip hit 
        if (hit2track[ihit] >= 0) continue;   
        
        TBHit& RecoHit = HitStore[ihit];
        double uhit = RecoHit.GetCoord()[0];
        double vhit = RecoHit.GetCoord()[1];

        TBTrack& trk = TrackStore[itrk]; 
        double utrk = trk.GetTE(_idut).GetState().GetPars()[2];
        double vtrk = trk.GetTE(_idut).GetState().GetPars()[3];

        // Skip all DUT hits with too large residuum 
        if ( std::abs(utrk-uhit) >= _maxResidualU && _maxResidualU > 0 ) continue;  
        if ( std::abs(vtrk-vhit) >= _maxResidualV && _maxResidualV > 0 ) continue; 
        
        // Finally, we will use a simple 2D distance to select best matching hit
        double hitdist = 0; 
        if ( _maxResidualU > 0 )  hitdist += std::abs(utrk-uhit); 
        if ( _maxResidualV > 0 )  hitdist += std::abs(vtrk-vhit); 
                      
        if( hitdist<distmin )
        {
          distmin=hitdist;
          besthit=ihit;
          besttrk=itrk;
        }
      }
    }
     
    streamlog_out(MESSAGE2) << "In matching loop: best hit " << besthit << " to best track " << besttrk << endl; 
    streamlog_out(MESSAGE2) << "  distmin: " <<  distmin  << endl; 
    
    // Check if a match was found
    if( besttrk>-1 &&  besthit>-1   )
    {   
      streamlog_out(MESSAGE2) << "  match found!!!"   << endl;
      nMatch++;
      hit2track[besthit] = besttrk;
      track2hit[besttrk] = besthit;   
      
      _iLastMatchedEvent = evt->getEventNumber(); 
      if (_iFirstMatchedEvent < 0)
        _iFirstMatchedEvent = evt->getEventNumber();
      
    } 
  
  } // End of do loop of matching DUT hits to fitted positions
  while( besttrk>-1 &&  besthit>-1);
  
  }
  _noOfMatchedTracks += nMatch; 
  
  streamlog_out(MESSAGE2) << nMatch << " DUT hits matched to tracks " << endl;
  streamlog_out(MESSAGE2) << HitStore.size() << " DUT hits not matched to any track " << endl;
  streamlog_out(MESSAGE2) << TrackStore.size() << " Tracks not matched to any DUT hit "<< endl;
  
  // Fill event tree
  _rootRunNumber = evt->getRunNumber();  
  _rootEventNumber = evt->getEventNumber();  
  _rootSensorID = dut.GetSensorID();        
  _rootNTelTracks = nTrack; 
  _rootNTelTracksWithRefHit = nTrackWithRefHit; 
  _rootNDUTDigits = nDUTDigits;
  _rootDUTHasTestPixels = hasTestPixels;
  _rootDUTGoodEvent = isGoodEvent; 
  
  _rootFile->cd("");
  _rootEventTree->Fill();  
  
  // Fill hit tree 
  
  streamlog_out(MESSAGE2) << "Start fill hit tree" << endl; 
  
  for(int ihit=0;ihit<(int)HitStore.size(); ++ihit)
  {
      
    TBHit& hit = HitStore[ihit];
    
    // Get uv hit coordinate of hit in mm 
    _rootHitU = hit.GetCoord()[0];         
    _rootHitV = hit.GetCoord()[1];   
    
    // Compute pixel cell which covers the hit coordinates 
    _rootHitCellU = dut.GetUCellFromCoord( _rootHitU, _rootHitV );  
    _rootHitCellV = dut.GetVCellFromCoord( _rootHitU, _rootHitV ); 
    
    // Cluster shape variables   
    PixelCluster Cluster = hit.GetCluster();
    _rootHitSeedCellU = Cluster.getUSeed();     
    _rootHitSeedCellV = Cluster.getVSeed();    
    _rootHitSeedPixelType = dut.GetPixelType(_rootHitSeedCellV, _rootHitSeedCellU);   
       
    _rootHitQuality = 0; 
    _rootHitClusterCharge = Cluster.getCharge() ; 
    _rootHitSeedCharge = Cluster.getSeedCharge() ; 
    _rootHitSize = Cluster.getSize();  
    _rootHitSizeU = Cluster.getUSize();     
    _rootHitSizeV = Cluster.getVSize();     
    
    // Add variables for matched track 
    if ( hit2track[ihit] >= 0 ) {  
     
      TBTrack& trk = TrackStore[hit2track[ihit]];
      _rootHitHasTrack = 0;  // matched   
       
      // Check track has a hit on reference (timing) plane
      if ( _iref >= 0 ) {
        if ( trk.GetTE(_iref).HasHit() ) {
          streamlog_out ( MESSAGE2 ) << "Track has hit on reference plane." << endl;
          _rootHitHasTrackWithRefHit = 0;
        }  
      } 
      
      auto p = trk.GetTE(_idut).GetState().GetPars();
      auto C = trk.GetTE(_idut).GetState().GetCov();  
      

      _rootHitLocalChi2 = TrackFitter.GetPredictedChi2(p, C, hit);
      _rootHitFitMomentum = trk.GetMomentum();   
           
      // Get predicted hit coordinates 
      double pu = p[2];
      double pv = p[3];
        
      // Get readout channels  
      int fitcol = dut.GetUCellFromCoord( pu, pv );     
      int fitrow = dut.GetVCellFromCoord( pu, pv );           
       
      _rootHitFitdUdW = p[0];     
      _rootHitFitdVdW = p[1];    
      _rootHitFitU = p[2];           
      _rootHitFitV = p[3];          
          
      _rootHitFitErrorU  = TMath::Sqrt(C(2,2));    
      _rootHitFitErrorV  = TMath::Sqrt(C(3,3));  
      _rootHitPullResidualU = (hit.GetCoord()[0] - p[2]) / TMath::Sqrt( C(2,2) + hit.GetCov()(0,0) ) ;   
      _rootHitPullResidualV = (hit.GetCoord()[1] - p[3]) / TMath::Sqrt( C(3,3) + hit.GetCov()(1,1) ) ;  
                                 
      _rootHitFitCellU = fitcol;      
      _rootHitFitCellV = fitrow;    
      _rootHitFitCellUCenter = dut.GetPixelCenterCoordU( fitrow, fitcol ); 
      _rootHitFitCellVCenter = dut.GetPixelCenterCoordV( fitrow, fitcol );                                        
      _rootHitTrackChi2 = trk.GetChiSqu(); 
      _rootHitTrackNDF = trk.GetNDF();
      _rootHitTrackNHits = trk.GetNumHits(); 
      
     
      
    } else {

      _rootHitHasTrack = -1; // no match with track
      // These are dummy values, always query hasTrack 
      _rootHitHasTrackWithRefHit = -1; 
      _rootHitLocalChi2 = -1;
      _rootHitFitMomentum = -1;            
      _rootHitFitU = -1;           
      _rootHitFitV = -1;           
      _rootHitFitdUdW = -1;      
      _rootHitFitdVdW = -1;        
      _rootHitFitErrorU = -1;  
      _rootHitFitErrorV = -1;                          
      _rootHitFitCellU = -1;      
      _rootHitFitCellV = -1;      
      _rootHitFitCellUCenter = -1;  
      _rootHitFitCellVCenter = -1;                                       
      _rootHitTrackChi2 = -1;  
      _rootHitTrackNDF = -1; 
      _rootHitTrackNHits = -1;  
      _rootHitPullResidualU = -1;   
      _rootHitPullResidualV = -1;  
    }
    
    // Fill tree with set variables 
    _rootFile->cd("");
    _rootHitTree->Fill();
  }

  // Fill track tree 
  streamlog_out(MESSAGE2) << "Start fill track tree" << endl; 
   
  for(int itrk=0;itrk<(int)TrackStore.size(); ++itrk)
  {
    TBTrack& trk = TrackStore[itrk];
    
    // Check track has a hit on reference (timing) plane
    _rootTrackWithRefHit = -1;
    if ( _iref >= 0 ) {
      if ( trk.GetTE(_iref).HasHit() ) {
        streamlog_out ( MESSAGE2 ) << "Track has hit on reference plane." << endl;
        _rootTrackWithRefHit = 0;
      }  
    } 
    
    auto p = trk.GetTE(_idut).GetState().GetPars();
    auto C = trk.GetTE(_idut).GetState().GetCov();  
           
    // Get predicted hit coordinates 
    double pu = p[2];
    double pv = p[3];
        
    // Get readout channels  
    int fitcellu = dut.GetUCellFromCoord( pu, pv );     
    int fitcellv = dut.GetVCellFromCoord( pu, pv );       
    
    _rootTrackPixelType = dut.GetPixelType(fitcellv, fitcellu);  
    _rootTrackFitMomentum = trk.GetMomentum();      
    _rootTrackFitdUdW = p[0];     
    _rootTrackFitdVdW = p[1];    
    _rootTrackFitU = p[2];           
    _rootTrackFitV = p[3];    
    _rootTrackFitCellU = fitcellu;      
    _rootTrackFitCellV = fitcellv;    
    _rootTrackFitCellUCenter = dut.GetPixelCenterCoordU( fitcellv, fitcellu ); 
    _rootTrackFitCellVCenter = dut.GetPixelCenterCoordV( fitcellv, fitcellu );                                        
    _rootTrackChi2 = trk.GetChiSqu(); 
    _rootTrackNDF = trk.GetNDF();  
    _rootTrackNHits = trk.GetNumHits();    
           
    if ( track2hit[itrk] >= 0  ) {
      TBHit& hit = HitStore[ track2hit[itrk] ];  
      PixelCluster Cluster = hit.GetCluster();
      _rootTrackLocalChi2 = TrackFitter.GetPredictedChi2(p, C, hit);
      _rootTrackHasHit = 0;  // match
      _rootTrackSeedCharge = Cluster.getSeedCharge() ; 
      _rootTrackHitU = hit.GetCoord()[0];
      _rootTrackHitV = hit.GetCoord()[1];
    } else {
      _rootTrackLocalChi2 = -1;
      _rootTrackHasHit = -1; // no match
      _rootTrackSeedCharge = -1;
      _rootTrackHitU = -1;
      _rootTrackHitV = -1;
    }
       
    _rootFile->cd("");
    _rootTrackTree->Fill(); 
  }
  
  
  if ( nTrack == 0 )  _noOfHitsWOTrack += (int)HitStore.size(); 
  else  _noOfHitsWTrack += (int)HitStore.size();  
   
  return;
}


//
// Method called after each event to check the data processed
//
void PixelDUTAnalyzer::check( LCEvent * )
{
}

//
// Method called after all data processing
//
void PixelDUTAnalyzer::end()
{

   // Print the summer
   streamlog_out(MESSAGE3) << endl << endl 
                           << "Total number of processed events:     " << setiosflags(ios::right) << _nEvt 
                           << resetiosflags(ios::right) << endl
                           << "Total number of events wo DUT hits:   " << setiosflags(ios::right) << _noOfEventWOInputHit 
                           << resetiosflags(ios::right) << endl
                           << "Total number of events wo Tel track:  " << setiosflags(ios::right) << _noOfEventWOInputTrack 
                           << resetiosflags(ios::right) << endl
                           << "Total number of accepted Tel tracks:  " << setiosflags(ios::right) << _noOfTracks 
                           << resetiosflags(ios::right) << endl
                           << "Total number of DUT hits:             " << setiosflags(ios::right) << _noOfHits
                           << resetiosflags(ios::right) << endl
                           << "Total number of matched tracks:       " << setiosflags(ios::right) << _noOfMatchedTracks 
                           << resetiosflags(ios::right) << endl
                           << "Total number of hits w/o track:       " << setiosflags(ios::right) <<  _noOfHitsWOTrack
                           << resetiosflags(ios::right) << endl 
                           << "Total number of hits with track:      " << setiosflags(ios::right) <<  _noOfHitsWTrack 
                           << resetiosflags(ios::right) << endl  
                           << "First match hit/track in event:       "  << setiosflags(ios::right) <<  _iFirstMatchedEvent 
                           << resetiosflags(ios::right) << endl  
                           << "Last match hit/track in event:        "  << setiosflags(ios::right) <<  _iLastMatchedEvent 
                           << resetiosflags(ios::right) << endl  
                           << endl << endl; 
     
   
   // CPU time end
   _timeCPU = clock()/1000 - _timeCPU;

   // Print message
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

    
   // Close root  file
   _rootFile->Write();
   _rootFile->Close();
   delete _rootFile;

}


//
// Method printing processor parameters
//
void PixelDUTAnalyzer::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "PixelDUTAnalyzer Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}


 // ROOT_OUTPUT
void PixelDUTAnalyzer::bookHistos()
{   
   
   _rootFile = new TFile( _rootFileName.c_str(),"recreate");
   _rootFile->cd("");
   
   // 
   // Hit Tree  
   _rootHitTree = new TTree("Hit","Hit info");
   _rootHitTree->Branch("iRun"            ,&_rootRunNumber        ,"iRun/I");
   _rootHitTree->Branch("iEvt"            ,&_rootEventNumber      ,"iEvt/I");
   _rootHitTree->Branch("sensorID"        ,&_rootSensorID       ,"sensorID/I");
   _rootHitTree->Branch("nTelTracks"      ,&_rootNTelTracks       ,"nTelTracks/I"); 
   _rootHitTree->Branch("nDutDigits"      ,&_rootNDUTDigits          ,"nDutDigits/I");
   _rootHitTree->Branch("hasTestPixels"   ,&_rootDUTHasTestPixels, "hasTestPixels/O");
   _rootHitTree->Branch("isGoodEvent"     ,&_rootDUTGoodEvent,  "isGoodEvent/O");
   _rootHitTree->Branch("clusterQuality"  ,&_rootHitQuality   ,"clusterQuality/I");
   _rootHitTree->Branch("u_hit"           ,&_rootHitU             ,"u_hit/D");
   _rootHitTree->Branch("v_hit"           ,&_rootHitV             ,"v_hit/D");     
   _rootHitTree->Branch("clusterCharge"   ,&_rootHitClusterCharge ,"clusterCharge/D");
   _rootHitTree->Branch("seedCharge"      ,&_rootHitSeedCharge       ,"seedCharge/D");
   _rootHitTree->Branch("sizeU"           ,&_rootHitSizeU        ,"sizeU/I");
   _rootHitTree->Branch("sizeV"           ,&_rootHitSizeV        ,"sizeV/I");
   _rootHitTree->Branch("size"            ,&_rootHitSize         ,"size/I");
   _rootHitTree->Branch("hasTrack"        ,&_rootHitHasTrack         ,"hasTrack/I");  
   _rootHitTree->Branch("hasTrackWithRefHit", &_rootHitHasTrackWithRefHit,"hasTrackWithRefHit/I"); 
   _rootHitTree->Branch("u_fit"           ,&_rootHitFitU             ,"u_fit/D");
   _rootHitTree->Branch("v_fit"           ,&_rootHitFitV             ,"v_fit/D"); 
   _rootHitTree->Branch("dudw_fit"        ,&_rootHitFitdUdW          ,"dudw_fit/D");
   _rootHitTree->Branch("dvdw_fit"        ,&_rootHitFitdVdW          ,"dvdw_fit/D");    
   _rootHitTree->Branch("u_fiterr"        ,&_rootHitFitErrorU        ,"u_fiterr/D");
   _rootHitTree->Branch("v_fiterr"        ,&_rootHitFitErrorV        ,"v_fiterr/D");   
   _rootHitTree->Branch("pull_resu"       ,&_rootHitPullResidualU    ,"pull_resu/D");
   _rootHitTree->Branch("pull_resv"       ,&_rootHitPullResidualV    ,"pull_resv/D");  
   _rootHitTree->Branch("cellU_fit"       ,&_rootHitFitCellU           ,"cellU_fit/I");
   _rootHitTree->Branch("cellV_fit"       ,&_rootHitFitCellV           ,"cellV_fit/I");
   _rootHitTree->Branch("cellU_hit"       ,&_rootHitCellU              ,"cellU_hit/I");
   _rootHitTree->Branch("cellV_hit"       ,&_rootHitCellV              ,"cellV_hit/I");
   _rootHitTree->Branch("cellU_seed"      ,&_rootHitSeedCellU          ,"cellU_seed/I");
   _rootHitTree->Branch("cellV_seed"      ,&_rootHitSeedCellV          ,"cellV_seed/I");
   _rootHitTree->Branch("cellUCenter_fit" ,&_rootHitFitCellUCenter  ,"cellUCenter_fit/D");
   _rootHitTree->Branch("cellVCenter_fit" ,&_rootHitFitCellVCenter  ,"cellVCenter_fit/D");                                      
   _rootHitTree->Branch("trackChi2"       ,&_rootHitTrackChi2      ,"trackChi2/D");
   _rootHitTree->Branch("trackNdof"       ,&_rootHitTrackNDF       ,"trackNdof/I");
   _rootHitTree->Branch("trackNHits"      ,&_rootHitTrackNHits     ,"trackNHits/I");  
   _rootHitTree->Branch("momentum"        ,&_rootHitFitMomentum      ,"momentum/D");    
   _rootHitTree->Branch("localChi2"       ,&_rootHitLocalChi2       ,"localChi2/D"); 
   _rootHitTree->Branch("pixeltype"       ,&_rootHitSeedPixelType     ,"pixeltype/I");  
    
   // 
   // Track Tree 
   _rootTrackTree = new TTree("Track","Track info");
   _rootTrackTree->Branch("iRun"            ,&_rootRunNumber      ,"iRun/I");
   _rootTrackTree->Branch("iEvt"            ,&_rootEventNumber    ,"iEvt/I");
   _rootTrackTree->Branch("sensorID"        ,&_rootSensorID        ,"sensorID/I");
   _rootTrackTree->Branch("nTelTracks"      ,&_rootNTelTracks     ,"nTelTracks/I"); 
   _rootTrackTree->Branch("nTelTracksWithRefHit" ,&_rootNTelTracksWithRefHit ,"nTelTracksWithRefHit/I");
   _rootTrackTree->Branch("nDutDigits"      ,&_rootNDUTDigits       ,"nDutDigits/I");
   _rootTrackTree->Branch("isGoodEvent"     ,&_rootDUTGoodEvent,  "isGoodEvent/O");
   _rootTrackTree->Branch("hasTestPixels"   ,&_rootDUTHasTestPixels, "hasTestPixels/O");
   _rootTrackTree->Branch("hasHit"          ,&_rootTrackHasHit         ,"hasHit/I");
   _rootTrackTree->Branch("u_hit"           ,&_rootTrackHitU             ,"u_hit/D");
   _rootTrackTree->Branch("v_hit"           ,&_rootTrackHitV             ,"v_hit/D");
   _rootTrackTree->Branch("hasRefHit"       ,&_rootTrackWithRefHit     ,"hasRefHit/I");
   _rootTrackTree->Branch("momentum"        ,&_rootTrackFitMomentum    ,"momentum/D");                                                           
   _rootTrackTree->Branch("u_fit"           ,&_rootTrackFitU           ,"u_fit/D");
   _rootTrackTree->Branch("v_fit"           ,&_rootTrackFitV           ,"v_fit/D");
   _rootTrackTree->Branch("dudw_fit"        ,&_rootTrackFitdUdW        ,"dudw_fit/D");
   _rootTrackTree->Branch("dvdw_fit"        ,&_rootTrackFitdVdW        ,"dvdw_fit/D");
   _rootTrackTree->Branch("cellU_fit"       ,&_rootTrackFitCellU       ,"cellU_fit/I");
   _rootTrackTree->Branch("cellV_fit"       ,&_rootTrackFitCellV       ,"cellV_fit/I");
   _rootTrackTree->Branch("cellUCenter_fit" ,&_rootTrackFitCellUCenter ,"cellUCenter_fit/D");
   _rootTrackTree->Branch("cellVCenter_fit" ,&_rootTrackFitCellVCenter ,"cellVCenter_fit/D");
   _rootTrackTree->Branch("trackChi2"       ,&_rootTrackChi2           ,"trackChi2/D");
   _rootTrackTree->Branch("trackNdof"       ,&_rootTrackNDF            ,"trackNdof/I");
   _rootTrackTree->Branch("trackNHits"      ,&_rootTrackNHits          ,"trackNHits/I");  
   _rootTrackTree->Branch("seedCharge"      ,&_rootTrackSeedCharge     ,"seedCharge/D");  
   _rootTrackTree->Branch("localChi2"       ,&_rootTrackLocalChi2      ,"localChi2/D"); 
   _rootTrackTree->Branch("pixeltype"      ,&_rootTrackPixelType     ,"pixeltype/I");     

   // 
   // Event Summay Tree 
   _rootEventTree = new TTree("Event","Event info");
   _rootEventTree->Branch("iRun"            ,&_rootRunNumber      ,"iRun/I");
   _rootEventTree->Branch("iEvt"            ,&_rootEventNumber    ,"iEvt/I");
   _rootEventTree->Branch("sensorID"        ,&_rootSensorID       ,"sensorID/I");   
   _rootEventTree->Branch("nTelTracks"      ,&_rootNTelTracks     ,"nTelTracks/I"); 
   _rootEventTree->Branch("nTelTracksWithRefHit" ,&_rootNTelTracksWithRefHit ,"nTelTracksWithRefHit/I");  
   _rootEventTree->Branch("nDutDigits"      ,&_rootNDUTDigits       ,"nDutDigits/I");
   _rootEventTree->Branch("hasTestPixels"   ,&_rootDUTHasTestPixels, "hasTestPixels/O");
   _rootEventTree->Branch("isGoodEvent"     ,&_rootDUTGoodEvent,  "isGoodEvent/O");
   
}

} // Namespace



