// DUTTreeProducer Processor  
// 
// See DUTTreeProducer.h for full documentation of processor. 
// 
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


#include "DUTTreeProducer.h"

// DEPFETTrackTools includes
#include "TBTrack.h"
#include "TBHit.h"
#include "PixelCluster.h"
#include "GenericTrackFitter.h"
#include "TrackInputProvider.h"
#include "Det.h"
#include "Utilities.h"
#include "MatrixDecoder.h"
#include "DEPFET.h" 

// Include basic C
#include <iostream>
#include <limits>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <string>
#include <algorithm>

// Include LCIO classes
#include <UTIL/CellIDDecoder.h>
#include <IMPL/LCFlagImpl.h>

// Include CLHEP classes
#include <CLHEP/Matrix/Vector.h>

// Used namespaces
using namespace std; 
using namespace CLHEP; 
using namespace lcio ;
using namespace marlin ;

namespace depfet {


//
// Instantiate this object
//
DUTTreeProducer aDUTTreeProducer ;

//
// Constructor
//
DUTTreeProducer::DUTTreeProducer() : Processor("DUTTreeProducer")
{

// Processor description
  _description = "DUTTreeProducer: DUT Analysis Processor" ;
   
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
  
  
   
  registerInputCollection (LCIO::TRACKERRAWDATA, "StatusCollection",
                           "Name of DUT status collection",
                           _statusColName, string("status")); 
  
  
   
 
     
  
  // Processor parameters:
  
  registerProcessorParameter ("AlignmentDBFileName",
                             "This is the name of the LCIO file with the alignment constants (add .slcio)",
                             _alignmentDBFileName, static_cast< string > ( "alignmentDB.slcio" ) );     
  
  registerProcessorParameter ("DUTPlane",
                              "Plane number of DUT along the beam line",
                              _idut,  static_cast < int > (3));
       
  registerProcessorParameter ("HitDistMax",
                              "Maximum hit2track distance for hit matching [mm]",
                              _hitdistmax,  static_cast < double > (2.0));
                               
  
  
  registerProcessorParameter ("SingleTrackEvents",
                              "Use only events with one track in telescope",
                              _singleTrackEvents,  static_cast < bool > (false));
  
  registerProcessorParameter( "RootFileName",
                               "Output root file name",
                               _rootFileName, std::string("dut_histos.root"));
   
}

//
// Method called at the beginning of data processing
//
void DUTTreeProducer::init() {
   
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
   
   // Print set parameters
   printProcessorParams();
   
   // Read detector constants from gear file
   _detector.ReadGearConfiguration();    
            
   // Read alignment data base file 
   _detector.ReadAlignmentDB( _alignmentDBFileName );    
   
   // Load DUT module    
   Det & dut = _detector.GetDet(_idut); 
          
   // Print out geometry information  
   streamlog_out ( MESSAGE3 )  << "D.U.T. plane  ID = " << dut.GetDAQID()
                               << "  at position = " << _idut 
                               << endl << endl;
    
      
   bookHistos();
   
}

//
// Method called for each run
//
void DUTTreeProducer::processRunHeader(LCRunHeader * run)
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
void DUTTreeProducer::processEvent(LCEvent * evt)
{
  
  // Print event number
  if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Event number is: "
                                                                << (evt->getEventNumber())
                                                                << std::endl << std::endl;
  
  streamlog_out(MESSAGE2) << "Event number is: " << (evt->getEventNumber())
                                                   << std::endl << std::endl;
  
  _nEvt ++ ;
  
  // Configure Kalman track fitter
  GenericTrackFitter TrackFitter(_detector);
  TrackFitter.SetNumIterations(2); 
     
  // Load DUT module    
  Det & dut = _detector.GetDet(_idut);   
  
  
    
  // Decoding of DUT matrix
  MatrixDecoder matrixDecoder(dut.GetNColumns(), dut.GetNRows()); 
  
  
  
  ShortVec statusVec; 
  bool isDUTStatusOk = getDUTStatus( evt, statusVec );
       
  //
  // Get telescope track collection
  //
  
  LCCollection* trackcol = NULL;
  bool isTrackok = true;   
  try {
    trackcol = evt->getCollection( _trackColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out(MESSAGE2) << "Not able to get collection "
                            << _trackColName
                            << " from event " << evt->getEventNumber()
                            << " in run " << evt->getRunNumber()  << endl << endl;   
    
    isTrackok = false; 
  }  
   
  int nTrack = 0;
  if(isTrackok) nTrack = trackcol->getNumberOfElements();
  if ( nTrack == 0)  ++_noOfEventWOInputTrack;   
  
  streamlog_out(MESSAGE2) << "Total of " << nTrack << " tracks in collection " << _trackColName << endl; 
  
  // 
  // Get DUT hit collection 
  // 
  
  LCCollection* hitcol = NULL;
  bool isDUTok=true;
  try {
    hitcol = evt->getCollection( _hitColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out(MESSAGE2) << "Not able to get collection "
                            << _hitColName
                            << " from event " << evt->getEventNumber()
                            << " in run " << evt->getRunNumber() << endl << endl;
    isDUTok=false;
  } 
  
  int nHit = 0;
  if(isDUTok) nHit = hitcol->getNumberOfElements();
    
  streamlog_out(MESSAGE2) << "Total of " << nHit << " hit(s) in collection " << _hitColName << endl;
  
  
  // Load telescope tracks and DUT hits 
  // ----------------------------------
  
  std::vector< TBTrack > TrackStore;
  std::vector<TBHit> HitStore;    
  
  
  if ( _singleTrackEvents  &&  nTrack > 1 ) {
    
    streamlog_out(MESSAGE2) << "Track Finder found more than one track in telescope."
                            <<" SingleTrackEvent requested, discard event!!" 
                            << endl; 
  } else {
    
    // Get tracks from the lcio file 
    TrackInputProvider TrackLCIOReader;  
    
    for(int itrk=0; itrk< nTrack ; itrk++)
    {
    
      // Retrieve track from LCIO 
      Track * lciotrk = dynamic_cast<Track*> (trackcol->getElementAt(itrk));
      
      // Convert LCIO -> TB track  
      TBTrack trk = TrackLCIOReader.MakeTBTrack( lciotrk, _detector );  
      
      // Refit track in nominal alignment
      bool trkerr = TrackFitter.Fit(trk);
      if ( trkerr ) {
        streamlog_out ( MESSAGE1 ) << "Fit failed. Skipping track!" << endl;
        continue;
      } 
      
      
      TrackStore.push_back(trk);                         
       
              
    } // End track loop
    
    // Get hits from the lcio file 
      
    for(int ihit=0; ihit< nHit ; ihit++)
    {
      // Built a TBHit
      TrackerHitImpl * lciohit = dynamic_cast<TrackerHitImpl*>( hitcol->getElementAt(ihit) ) ;
      TBHit RecoHit ( lciohit  );        
         
      // We have to find plane number of the hit 
      int daqid = RecoHit.GetDAQID();      
      int ipl = _detector.GetPlaneNumber(daqid);  
      
      // Store all hits on the DUT module  
      if( dut.GetPlaneNumber() == ipl )
      {         
        streamlog_out(MESSAGE2) << " DUT hit at plane " << ipl 
                                << "   U [mm]= " << RecoHit.GetCoord()[0][0]
                                << "   V [mm]= " << RecoHit.GetCoord()[1][0] 
                                << endl;
      
        HitStore.push_back( RecoHit );  
      }
    } 
    
     
  } // endif singletrack events
  
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
  double distmin=numeric_limits<double >::max();
   
  do{
    int besttrk=-1;
    int besthit=-1;
    
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
        double uhit = RecoHit.GetCoord()[0][0];
        double vhit = RecoHit.GetCoord()[1][0];

        TBTrack& trk = TrackStore[itrk]; 
        double utrk = trk.GetTE(_idut).GetState().GetPars()[2][0];
        double vtrk = trk.GetTE(_idut).GetState().GetPars()[3][0];
        
        double hitdist = std::abs(utrk-uhit) + std::abs(vtrk-vhit);
                      
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
    
    // Check if best match is good enough
    if( distmin < _hitdistmax  )
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
  
  while( distmin < _hitdistmax );
  
  _noOfMatchedTracks += nMatch; 
  
  streamlog_out(MESSAGE2) << nMatch << " DUT hits matched to tracks " << endl;
  streamlog_out(MESSAGE2) << HitStore.size() << " DUT hits not matched to any track " << endl;
  streamlog_out(MESSAGE2) << TrackStore.size() << " Tracks not matched to any DUT hit "<< endl;
  

  // Fill hit tree 

  streamlog_out(MESSAGE2) << "Start fill hit tree" << endl; 
  
  for(int ihit=0;ihit<(int)HitStore.size(); ++ihit)
  {
    
    _rootRunNumber = evt->getRunNumber();  
    _rootEventNumber = evt->getEventNumber();  
    _rootDetectorID = dut.GetDAQID();       
    _rootNTelTracks = nTrack; 
    _rootNDUTHits = (int)HitStore.size(); 
    
    
    TBHit& hit = HitStore[ihit];
    
    _rootHitU = hit.GetCoord()[0][0];         
    _rootHitV = hit.GetCoord()[1][0];   
    
    _rootHitCol = dut.GetColumnFromCoord( _rootHitU, _rootHitV );  
    _rootHitRow = dut.GetRowFromCoord( _rootHitU, _rootHitV );  

    // Cluster shape variables   
    
    PixelCluster Cluster = hit.GetCluster();
       
    _rootHitQuality = 0; 
    _rootHitCharge = Cluster.getCharge() ; 
    _rootHitSeedCharge = Cluster.getSeedCharge() ; 
    _rootHitSize = Cluster.getSize();  
    _rootHitSizeCol = Cluster.getUSize();     
    _rootHitSizeRow = Cluster.getVSize();     
   
    

    // Add variables for matched track 
    if ( hit2track[ihit] >= 0 ) {  
     
      _rootHitHasTrack = 0;  // matched         
       
      TBTrack& trk = TrackStore[hit2track[ihit]];      
      HepMatrix p = trk.GetTE(_idut).GetState().GetPars();
      HepSymMatrix C = trk.GetTE(_idut).GetState().GetCov();  
      double hitchi2 = TrackFitter.GetPredictedChi2(p, C, hit);

      _rootHitFitMomentum = trk.GetMomentum();   
           
      // Get predicted hit coordinates 
      double pu = p[2][0];
      double pv = p[3][0];
        
      // Get readout channels  
      int fitcol = dut.GetColumnFromCoord( pu, pv );     
      int fitrow = dut.GetRowFromCoord( pu, pv );           
       
      _rootHitFitdUdW = p[0][0];     
      _rootHitFitdVdW = p[1][0];    
      _rootHitFitU = p[2][0];           
      _rootHitFitV = p[3][0];          
          
      _rootHitFitErrorU  = TMath::Sqrt(C[2][2]);    
      _rootHitFitErrorV  = TMath::Sqrt(C[3][3]);  
      _rootHitPullResidualU = (hit.GetCoord()[0][0] - p[2][0]) / TMath::Sqrt( C[2][2] + hit.GetCov()[0][0] ) ;   
      _rootHitPullResidualV = (hit.GetCoord()[1][0] - p[3][0]) / TMath::Sqrt( C[3][3] + hit.GetCov()[1][1] ) ;  
                                 
      _rootHitFitCol = fitcol;      
      _rootHitFitRow = fitrow;    
      _rootHitFitPixU = dut.GetPixelCenterCoordU( fitrow, fitcol ); 
      _rootHitFitPixV = dut.GetPixelCenterCoordV( fitrow, fitcol );                                        
      _rootHitTrackChi2 = trk.GetChiSqu(); 
      _rootHitTrackNDF = trk.GetNDF();
      _rootHitLocalChi2 = hitchi2;   
      
     
      
    } else {

      _rootHitHasTrack = -1; // no match
      _rootHitFitMomentum = -1;       
      // These are dummy values, always query hasTrack      
      _rootHitFitU = -1;           
      _rootHitFitV = -1;           
      _rootHitFitdUdW = -1;      
      _rootHitFitdVdW = -1;        
      _rootHitFitErrorU = -1;  
      _rootHitFitErrorV = -1;                          
      _rootHitFitCol = -1;      
      _rootHitFitRow = -1;      
      _rootHitFitPixU = -1;  
      _rootHitFitPixV = -1;                                       
      _rootHitTrackChi2 = -1;  
      _rootHitTrackNDF = -1; 
      _rootHitLocalChi2 = -1;  
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
    
    _rootRunNumber = evt->getRunNumber();  
    _rootEventNumber = evt->getEventNumber();  
    _rootDetectorID = dut.GetDAQID();   
    _rootNTelTracks = nTrack; 
    _rootNDUTHits = (int)HitStore.size(); 
   
       
    TBTrack& trk = TrackStore[itrk];

    _rootTrackFitMomentum = trk.GetMomentum();   
    
    HepMatrix p = trk.GetTE(_idut).GetState().GetPars();
    HepSymMatrix C = trk.GetTE(_idut).GetState().GetCov();  
           
    // Get predicted hit coordinates 
    double pu = p[2][0];
    double pv = p[3][0];
        
    // Get readout channels  
    int fitcol = dut.GetColumnFromCoord( pu, pv );     
    int fitrow = dut.GetRowFromCoord( pu, pv );           
    int Xpixel = matrixDecoder.getIndexFromXY(fitcol, fitrow);
        
    _rootTrackFitdUdW = p[0][0];     
    _rootTrackFitdVdW = p[1][0];    
    _rootTrackFitU = p[2][0];           
    _rootTrackFitV = p[3][0];    
    _rootTrackFitCol = fitcol;      
    _rootTrackFitRow = fitrow;    
    _rootTrackFitPixU = dut.GetPixelCenterCoordU( fitrow, fitcol ); 
    _rootTrackFitPixV = dut.GetPixelCenterCoordV( fitrow, fitcol );                                        
    _rootTrackChi2 = trk.GetChiSqu(); 
    _rootTrackNDF = trk.GetNDF();   
    _rootTrack1x1Quality = 0;    
    _rootTrack3x3Quality = 0;  
            
    if (  isDUTStatusOk  ) {
      
      if ( ( fitcol >= 0 )  &&  ( fitcol < dut.GetNColumns() ) &&
               ( fitrow >= 0 )  &&  ( fitrow < dut.GetNRows() ) ) {
         
        _rootTrack1x1Quality = statusVec[Xpixel];
        
      }
               
      for (int iRow = fitrow - (3 / 2); iRow <= fitrow + (3 / 2); iRow++) {
        for (int iCol =  fitcol - (3 / 2); iCol <= fitcol + (3 / 2); iCol++) {
          // Always check we are still within the sensor boundaries!!!
          if ( ( iCol >= 0 )  &&  ( iCol < dut.GetNColumns() ) &&
               ( iRow >= 0 )  &&  ( iRow < dut.GetNRows() ) ) {
                   
              int index = matrixDecoder.getIndexFromXY(iCol, iRow);
              
              bool isGood = true; 
              if (isDUTStatusOk) {
                isGood = ( statusVec[index] == DEPFET::GOODPIXEL );
              }
                
              if ( !isGood  ) {
                // Bad pixel in cluster degrades Xing quality
                _rootTrack3x3Quality++;  
              }     
          }
        }
      }
            
    }      
    
    if ( track2hit[itrk] >= 0  ) {
      
      _rootTrackHasHit = 0;  // match   
      
      TBHit& hit = HitStore[track2hit[itrk]];  
      PixelCluster Cluster = hit.GetCluster();
      _rootTrackSeedCharge = Cluster.getSeedCharge() ; 

        
    } else {
      _rootTrackHasHit = -1; // no match
      _rootTrackSeedCharge = -1;
    }
     
    
    
    
    _rootFile->cd("");
    _rootTrackTree->Fill(); 
  }
  
  
  if ( nTrack == 0 )  _noOfHitsWOTrack += (int)HitStore.size(); 
  else  _noOfHitsWTrack += (int)HitStore.size();  
  
  // Fill event tree
  _rootNDUTHits = (int)HitStore.size(); 
  _rootEventNDUTTracks = (int) TrackStore.size(); 
  _rootNTelTracks = nTrack; 
  _rootEventNMatched = nMatch; 
  
  
  _rootFile->cd("");
  _rootEventTree->Fill();  
   
  return;
}


//
// Method called after each event to check the data processed
//
void DUTTreeProducer::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void DUTTreeProducer::end()
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
     
    
   // Load DUT module    
   Det & dut = _detector.GetDet(_idut); 

   // DUT fake rate 
   double dutNPixels = dut.GetNColumns()*dut.GetNRows(); 
   double dutFakeRate  = static_cast<double> (_noOfHits - _noOfMatchedTracks) / static_cast<double> ( _nEvt ) / dutNPixels;
   
   // Telescope track effi 
   double telTrackEffi = static_cast<double> (_noOfMatchedTracks) / static_cast<double> (_noOfHits);
   
   // Telescope track fake rate
   double telTrackFakeRate = static_cast<double> ( _noOfTracks - _noOfMatchedTracks ) / static_cast<double> ( _nEvt );
   
   // Print perfomance summer
   streamlog_out(MESSAGE3) << std::setprecision(3)
                           << "DUT Fake Rate is ~ " << setiosflags(ios::right) << dutFakeRate
                           << resetiosflags(ios::right) << endl
                           << std::setprecision(3) 
                           << "Telescope Track Efficiency is ~ "  << setiosflags(ios::right) << telTrackEffi
                           << resetiosflags(ios::right) << endl
                           << std::setprecision(0)
                           << std::setprecision(3) 
                           << "Telescope Track Fake Rate is ~ "  << setiosflags(ios::right) << telTrackFakeRate
                           << resetiosflags(ios::right) << endl
                           << std::setprecision(0)
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

//! A function to read DUT status map 
/*! Vector of status values for DUT pixels 
 */
bool DUTTreeProducer::getDUTStatus(LCEvent * evt, ShortVec & statusVec) {
  
  // Load DUT module    
  Det & dut = _detector.GetDet(_idut); 
  
  bool isOk = false; 
  
  try {  
     
    LCCollectionVec * statusCollection = dynamic_cast < LCCollectionVec * > (evt->getCollection(_statusColName)); 
    CellIDDecoder<TrackerRawDataImpl> Decoder(statusCollection);    
     
    for ( size_t iDet = 0; iDet < statusCollection->size() ; iDet++) {
        
      // Get pixel status matrix 
      TrackerRawDataImpl * status = dynamic_cast<TrackerRawDataImpl*> (statusCollection->getElementAt(iDet));    
           
      // DAQ ID for pixel detector	
      int sensorID =  Decoder(status)["sensorID"];
         
      if (sensorID == dut.GetDAQID()) { 
        statusVec = status->getADCValues(); 
        isOk = true; 
        break;
      }
    }  
    
  } catch(DataNotAvailableException &e){} 
  
  return isOk; 
}
   

   




//
// Method printing processor parameters
//
void DUTTreeProducer::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "DUTTreeProducer Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}


 // ROOT_OUTPUT
void DUTTreeProducer::bookHistos()
{   
   
   _rootFile = new TFile( _rootFileName.c_str(),"recreate");
   _rootFile->cd("");
   
   // 
   // Hit Tree  
   _rootHitTree = new TTree("Hit","Hit info");
   _rootHitTree->Branch("iRun"            ,&_rootRunNumber        ,"iRun/I");
   _rootHitTree->Branch("iEvt"            ,&_rootEventNumber      ,"iEvt/I");
   _rootHitTree->Branch("det"             ,&_rootDetectorID       ,"det/I");
   _rootHitTree->Branch("nTelTracks"      ,&_rootNTelTracks       ,"nTelTracks/I"); 
   _rootHitTree->Branch("nduthits"        ,&_rootNDUTHits          ,"nduthits/I");



   _rootHitTree->Branch("clusterQuality"  ,&_rootHitQuality   ,"clusterQuality/I");
   _rootHitTree->Branch("u_hit"           ,&_rootHitU             ,"u_hit/D");
   _rootHitTree->Branch("v_hit"           ,&_rootHitV             ,"v_hit/D");     
   _rootHitTree->Branch("clusterCharge"   ,&_rootHitCharge    ,"clusterCharge/D");
   _rootHitTree->Branch("seedCharge"      ,&_rootHitSeedCharge       ,"seedCharge/D");
   _rootHitTree->Branch("sizeCol"         ,&_rootHitSizeCol   ,"sizeCol/I");
   _rootHitTree->Branch("sizeRow"         ,&_rootHitSizeRow   ,"sizeRow/I");
   _rootHitTree->Branch("size"            ,&_rootHitSize      ,"size/I");
   _rootHitTree->Branch("hasTrack"        ,&_rootHitHasTrack         ,"hasTrack/I");   
   _rootHitTree->Branch("u_fit"           ,&_rootHitFitU             ,"u_fit/D");
   _rootHitTree->Branch("v_fit"           ,&_rootHitFitV             ,"v_fit/D"); 
   _rootHitTree->Branch("dudw_fit"        ,&_rootHitFitdUdW          ,"dudw_fit/D");
   _rootHitTree->Branch("dvdw_fit"        ,&_rootHitFitdVdW          ,"dvdw_fit/D");    
   _rootHitTree->Branch("u_fiterr"        ,&_rootHitFitErrorU        ,"u_fiterr/D");
   _rootHitTree->Branch("v_fiterr"        ,&_rootHitFitErrorV        ,"v_fiterr/D");   
   _rootHitTree->Branch("pull_resu"       ,&_rootHitPullResidualU    ,"pull_resu/D");
   _rootHitTree->Branch("pull_resv"       ,&_rootHitPullResidualV    ,"pull_resv/D");  
   _rootHitTree->Branch("col_fit"         ,&_rootHitFitCol           ,"col_fit/I");
   _rootHitTree->Branch("row_fit"         ,&_rootHitFitRow           ,"row_fit/I");
   _rootHitTree->Branch("col_hit"         ,&_rootHitCol           ,"col_hit/I");
   _rootHitTree->Branch("row_hit"         ,&_rootHitRow           ,"row_hit/I");
   _rootHitTree->Branch("u_pixel"         ,&_rootHitFitPixU          ,"u_pixel/D");
   _rootHitTree->Branch("v_pixel"         ,&_rootHitFitPixV          ,"v_pixel/D");                                      
   _rootHitTree->Branch("chi2"            ,&_rootHitTrackChi2      ,"chi2/D");
   _rootHitTree->Branch("ndof"            ,&_rootHitTrackNDF       ,"ndof/I");
   _rootHitTree->Branch("chi2pred"        ,&_rootHitLocalChi2        ,"chi2pred/D");  
   _rootHitTree->Branch("momentum"        ,&_rootHitFitMomentum      ,"momentum/D");    

  
   
    
   // 
   // Track Tree 
   _rootTrackTree = new TTree("Track","Track info");
   _rootTrackTree->Branch("iRun"            ,&_rootRunNumber      ,"iRun/I");
   _rootTrackTree->Branch("iEvt"            ,&_rootEventNumber    ,"iEvt/I");
   _rootTrackTree->Branch("det"             ,&_rootDetectorID     ,"det/I");
   _rootTrackTree->Branch("nTelTracks"      ,&_rootNTelTracks     ,"nTelTracks/I"); 
   _rootTrackTree->Branch("nduthits"        ,&_rootNDUTHits       ,"nduthits/I");
   
   _rootTrackTree->Branch("hasHit"          ,&_rootTrackHasHit         ,"hasHit/I");
   _rootTrackTree->Branch("momentum"        ,&_rootTrackFitMomentum    ,"momentum/D");    
   _rootTrackTree->Branch("seedCharge"      ,&_rootTrackSeedCharge     ,"seedCharge/D");                                                        
   _rootTrackTree->Branch("1x1Quality"      ,&_rootTrack1x1Quality     ,"1x1Quality/I");
   _rootTrackTree->Branch("3x3Quality"      ,&_rootTrack3x3Quality     ,"3x3Quality/I");
   _rootTrackTree->Branch("u_fit"           ,&_rootTrackFitU           ,"u_fit/D");
   _rootTrackTree->Branch("v_fit"           ,&_rootTrackFitV           ,"v_fit/D");
   _rootTrackTree->Branch("dudw_fit"        ,&_rootTrackFitdUdW        ,"dudw_fit/D");
   _rootTrackTree->Branch("dvdw_fit"        ,&_rootTrackFitdVdW        ,"dvdw_fit/D");
   _rootTrackTree->Branch("col_fit"         ,&_rootTrackFitCol         ,"col_fit/I");
   _rootTrackTree->Branch("row_fit"         ,&_rootTrackFitRow         ,"row_fit/I");
   _rootTrackTree->Branch("u_pixel"         ,&_rootTrackFitPixU        ,"u_pixel/D");
   _rootTrackTree->Branch("v_pixel"         ,&_rootTrackFitPixV        ,"v_pixel/D");
   _rootTrackTree->Branch("chi2"            ,&_rootTrackChi2           ,"chi2/D");
   _rootTrackTree->Branch("ndof"            ,&_rootTrackNDF            ,"ndof/I");
   
   // 
   // Event Summay Tree 
   _rootEventTree = new TTree("Event","Event info");
   _rootEventTree->Branch("iRun"            ,&_rootRunNumber      ,"iRun/I");
   _rootEventTree->Branch("iEvt"            ,&_rootEventNumber    ,"iEvt/I");
   _rootEventTree->Branch("det"             ,&_rootDetectorID     ,"det/I");    
   _rootEventTree->Branch("nTelTracks"      ,&_rootNTelTracks     ,"nTelTracks/I"); 
   _rootEventTree->Branch("nduthits"           ,&_rootNDUTHits          ,"nduthits/I");
   _rootEventTree->Branch("nDUTTracks"      ,&_rootEventNDUTTracks     ,"nDUTTracks/I");
   _rootEventTree->Branch("matched"         ,&_rootEventNMatched       ,"matched/I");
   
}

} // Namespace



