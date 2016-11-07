// X0ImageProduce
// 
// See X0ImageProducer.h for full documentation of processor. 
// 
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


#include "X0ImageProducer.h"

// DEPFETTrackTools includes
#include "TBTrack.h"
#include "TBHit.h"
#include "TBKalmanMSC.h"
#include "TrackInputProvider.h"
#include "Utilities.h"

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
X0ImageProducer aX0ImageProducer ;

//
// Constructor
//
X0ImageProducer::X0ImageProducer() : Processor("X0ImageProducer")
{

// Processor description
  _description = "X0ImageProducer: X/X0 measurement for EUDET/AIDA telescope data" ;
   
//
// Input collections 
   
  registerInputCollection( LCIO::TRACK,
                           "DownStreamTrackCollection" ,
                           "Name of downstream track collection"  ,
                           _downStreamTrackColName ,
                           std::string("down_tracks") ) ;  

  registerInputCollection( LCIO::TRACK,
                           "UpStreamTrackCollection" ,
                           "Name of upstream track collection"  ,
                           _upStreamTrackColName ,
                           std::string("up_tracks") ) ;  
  
  // Processor parameters:
  
  registerProcessorParameter ("AlignmentDBFileName",
                             "This is the name of the LCIO file with the alignment constants (add .slcio)",
                             _alignmentDBFileName, static_cast< string > ( "alignmentDB.slcio" ) );     
  
  registerProcessorParameter ("DUTPlane",
                              "Plane number of DUT along the beam line",
                              _idut,  static_cast < int > (3));
                               
  registerProcessorParameter( "RootFileName",
                               "Output root file name",
                               _rootFileName, std::string("X0.root"));

  registerProcessorParameter ("MaxDist",
                              "Maximum distance between up and downstream tracks at dut plane [mm]",
                              _maxDist,  static_cast < double > (0.1));

}

//
// Method called at the beginning of data processing
//
void X0ImageProducer::init() {
   
   // Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
   _timeCPU = clock()/1000 ;
   
   
        
   // Print set parameters
   printProcessorParams();
   
   // Read detector constants from gear file
   _detector.ReadGearConfiguration();    
            
   // Read alignment data base file 
   _detector.ReadAlignmentDB( _alignmentDBFileName );    
   
   // Load DUT module    
   Det & dut = _detector.GetDet(_idut); 
          
   // Print out geometry information  
   streamlog_out ( MESSAGE3 )  << "Scatter DUT plane  ID = " << dut.GetDAQID()
                               << "  at position = " << _idut 
                               << endl << endl;
    
      
   bookHistos();
   
}

//
// Method called for each run
//
void X0ImageProducer::processRunHeader(LCRunHeader * run)
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
void X0ImageProducer::processEvent(LCEvent * evt)
{
  
  // Print event number
  if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                << (evt->getEventNumber())
                                                                << std::endl << std::endl;
  
  streamlog_out(MESSAGE2) << "Events processed: " << (evt->getEventNumber())
                                                   << std::endl << std::endl;
  
  _nEvt ++ ; 
 
       
  //
  // Get telescope track collection
  //
  
  LCCollection* downtrackcol = 0;   
  try {
    downtrackcol = evt->getCollection( _downStreamTrackColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out(MESSAGE2) << "Not able to get collection "
                            << _downStreamTrackColName
                            << " from event " << evt->getEventNumber()
                            << " in run " << evt->getRunNumber()  << endl << endl;   
     
    throw SkipEventException(this);
  }  
  
  LCCollection* uptrackcol = 0;
  try {
    uptrackcol = evt->getCollection( _upStreamTrackColName ) ;
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out(MESSAGE2) << "Not able to get collection "
                            << _upStreamTrackColName
                            << " from event " << evt->getEventNumber()
                            << " in run " << evt->getRunNumber()  << endl << endl;   
    
    throw SkipEventException(this);
  } 
  
  // Configure Kalman track fitter
  TBKalmanMSC TrackFitterMSC(_detector);
  TrackInputProvider TrackLCIOReader;  
  
  // Store tracks 
  std::vector< TBTrack > downTrackStore;    
    
  // Loop over tracks in input track collection
    
  int nDownTracks = downtrackcol->getNumberOfElements();  
  for(int itrk=0; itrk< nDownTracks ; itrk++) {
    
    // Retrieve track from LCIO 
    Track * lciotrk = dynamic_cast<Track*> (downtrackcol->getElementAt(itrk));
      
    // Convert LCIO -> TB track  
    TBTrack trk = TrackLCIOReader.MakeTBTrack( lciotrk, _detector ); 
    
    // Refit track in nominal alignment
    bool trkerr = TrackFitterMSC.ProcessTrack(trk, 1, 0);
    if ( trkerr ) {
      streamlog_out ( MESSAGE1 ) << "Fit failed. Skipping track!" << endl;
      continue;
    } 
          
    downTrackStore.push_back(trk);   
  }

  // Store tracks 
  std::vector< TBTrack > upTrackStore;    
    
  // Loop over tracks in input track collection
    
  int nUpTracks = uptrackcol->getNumberOfElements();  
  for(int itrk=0; itrk< nUpTracks ; itrk++) {
     
    // Retrieve track from LCIO 
    Track * lciotrk = dynamic_cast<Track*> (uptrackcol->getElementAt(itrk));
      
    // Convert LCIO -> TB track  
    TBTrack trk = TrackLCIOReader.MakeTBTrack( lciotrk, _detector ); 

    // Refit track in nominal alignment
    bool trkerr = TrackFitterMSC.ProcessTrack(trk, -1, 0);  
    if ( trkerr ) {
      streamlog_out ( MESSAGE1 ) << "Fit failed. Skipping track!" << endl;
      continue;
    }  
    
    upTrackStore.push_back(trk);   
  }
  
  
  // 
  // Match upstream and downstream tracks at DUT plane 
  int nMatch=0;	
  vector< vector<int> > up2down(upTrackStore.size() );
  vector< vector<int> > down2up(downTrackStore.size() );
   
  // Continue matching tracks and hits until all tracks are matched 
  // or no hit is close enough to a track!! 
  double distmin=numeric_limits<double >::max();
   
  Det dut = _detector.GetDet(_idut);

  do{
    int bestup=-1;
    int bestdown=-1;
    
    distmin=numeric_limits<double >::max();
      
    for(int iup=0;iup<(int)upTrackStore.size(); iup++)
    {
            
      // If matched, skip track 
      if ( up2down[iup].size() >= 0 ) continue;    
      
      for(int idown=0; idown< (int)downTrackStore.size() ; idown++)
      {

        // If matched, skip track 
        if (down2up[idown].size() >= 0) continue; 
        
        TBTrack& uptrack = upTrackStore[iup];
        TBTrack& downtrack = downTrackStore[idown];
        
        // In and OutStates of the reconstructed Track at the current detector
        TBTrackState& inState=uptrack.GetTE(_idut).GetState();
        TBTrackState& OutState=downtrack.GetTE(_idut).GetState(); 
        
        double u_in = InState.GetPars()[2][0];
        double v_in = InState.GetPars()[3][0];
        double u_out = OutState.GetPars()[2][0];
        double v_out = OutState.GetPars()[3][0];
        
        double hitdist = std::abs(u_in-u_out) + std::abs(v_in-v_out);
                      
        if( hitdist<distmin )
        {
          distmin=hitdist;
          bestup=iup;
          bestdown=idown;
        }
      }
    }
    
    streamlog_out(MESSAGE2) << "In matching loop: best up " << bestup << " to best down " << bestdown << endl; 
    streamlog_out(MESSAGE2) << "  distmin: " <<  distmin  << endl; 
    
    // Check if best match is good enough
    if( distmin < _maxDist  )
    {   

      streamlog_out(MESSAGE2) << "  match found!!!"   << endl;
      nMatch++;
      up2down[bestup].push_back( bestdown );
      down2up[bestdown].push_back( bestup );     
    } 
  
  } // End matching loop
  while( distmin < _maxDist );
  
  // Fill event tree
  _rootRunNumber = evt->getRunNumber();  
  _rootEventNumber = evt->getEventNumber();  
  _rootnUpTracks = upTrackStore.size();  
  _rootnDownTracks = downTrackStore.size(); 
  _rootNMatched = nMatch; 
  _rootFile->cd("");
  _rootEventTree->Fill();    

 

  for(int iup=0;iup<(int)upTrackStore.size(); iup++)
  {
  
    // Check upstream track is matched 
    if ( up2down[iup].size() != 1 ) continue; 
    
    TBTrack& uptrack = upTrackStore[iup];
    TBTrack& downtrack = downTrackStore[ up2down[iup][0] ];

    // comboChi2 is chi2 combination of track in upstream and downstream telescope arm
    double comboChi2 = uptrack.GetChiSqu()+downtrack.GetChiSqu(); 
    
    
    //MSC Analysis for the reconstructed angles
    //Here we use the In and Out State and the GetScatterKinks function of the TBKalmanMSC Class
     
    // In and OutStates of the reconstructed Track at the current detector
    TBTrackState& InState=uptrack.GetTE(_idut).GetState();
    TBTrackState& OutState=downtrack.GetTE(_idut).GetState(); 
    
    //Angles and angle errors
    HepMatrix theta(2,1,0);
    HepSymMatrix Cov(2,0);
    theta = TrackFitterMSC.GetScatterKinks(dut, InState, OutState); 
    Cov = TrackFitterMSC.GetScatterKinkCov(dut, InState, OutState);
    
    // Get the track parameters of the fitted track on the current sensor
    // The u and v positions are needed for a position-resolved measurement
    HepMatrix p_in = InState.GetPars();
    HepMatrix p_out = OutState.GetPars();
     	
    // Fill root variables 
    _rootDaqID = dut.GetDAQID(); 
    _rootPlaneID = _idut;
    _root_momentum = uptrack.GetMomentum(); 
    _rootTrackHits = uptrack.GetNumHits() + downtrack.GetNumHits(); 
    _rootTrackProbUp = TMath::Prob(uptrack.GetChiSqu(),uptrack.GetNDF());
    _rootTrackProbDown = TMath::Prob(downtrack.GetChiSqu(),downtrack.GetNDF());
    _rootTrackProbCombo = TMath::Prob( comboChi2 ,uptrack.GetNDF()+downtrack.GetNDF());
    
    _root_u_in = p_in[2][0]; 
    _root_v_in = p_in[3][0];
    _root_u_out = p_out[2][0]; 
    _root_v_out = p_out[3][0];
    _root_u = 0.5*(p_in[2][0] + p_out[2][0]); 
    _root_v = 0.5*(p_in[3][0] + p_out[3][0]);  
    _root_dudw = 0.5*(p_in[0][0] + p_out[0][0]);    
    _root_dvdw = 0.5*(p_in[1][0] + p_out[1][0]);   
   
    _root_angle1 = theta[0][0];
    _root_angle2 = theta[1][0];
    _root_angle1_err = TMath::Sqrt(Cov[0][0]);
    _root_angle2_err = TMath::Sqrt(Cov[1][1]);
    
    _rootMscTree->Fill();       
    
  } // end track loop 		
     
   
  
  return;
}




//
// Method called after each event to check the data processed
//
void X0ImageProducer::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void X0ImageProducer::end()
{

   // Print the summer
   streamlog_out(MESSAGE3) << endl << endl 
                           << "Total number of processed events:     " << setiosflags(ios::right) << _nEvt 
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

    
   // Close root  files
   _rootFile->Write();
 
}

//
// Method printing processor parameters
//
void X0ImageProducer::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "X0ImageProducer Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}


 // ROOT_OUTPUT
void X0ImageProducer::bookHistos() {   
   
  _rootFile = new TFile( _rootFileName.c_str(),"recreate");
   
  // 
  // Event Summary Tree 
  _rootEventTree = new TTree("Event","Event info");
  _rootEventTree->Branch("iRun"            ,&_rootRunNumber      ,"iRun/I");
  _rootEventTree->Branch("iEvt"            ,&_rootEventNumber    ,"iEvt/I");
  _rootEventTree->Branch("nDownTracks"      ,&_rootnDownTracks     ,"nDownTracks/I");
  _rootEventTree->Branch("nUpTracks"      ,&_rootnUpTracks     ,"nUpTracks/I"); 
  _rootEventTree->Branch("nMatched"      ,&_rootNMatched     ,"nMatched/I");

  
    
  // 
  // Track Tree of the whole track 
  _rootMscTree = new TTree("MSCTree","Multiple scattering data");
  _rootMscTree->Branch("iRun"            ,&_rootRunNumber       ,"iRun/I");
  _rootMscTree->Branch("iEvt"            ,&_rootEventNumber     ,"iEvt/I"); 
  _rootMscTree->Branch("daqid"           ,&_rootDaqID           ,"daqid/I"); 
  _rootMscTree->Branch("ipl"             ,&_rootPlaneID         ,"ipl/I"); 
  _rootMscTree->Branch("nhits"           ,&_rootTrackHits       ,"nhits/I"); 
  _rootMscTree->Branch("prob_up"         ,&_rootTrackProbUp     ,"prob_up/D"); 
  _rootMscTree->Branch("prob_down"       ,&_rootTrackProbDown   ,"prob_down/D"); 
  _rootMscTree->Branch("prob_combo"      ,&_rootTrackProbCombo  ,"prob_combo/D"); 

  _rootMscTree->Branch("u"            ,&_root_u             ,"u/D");
  _rootMscTree->Branch("v"            ,&_root_v             ,"v/D");
  _rootMscTree->Branch("dudw"         ,&_root_dudw          ,"dudw/D");
  _rootMscTree->Branch("dvdw"         ,&_root_dvdw          ,"dvdw/D"); 


  _rootMscTree->Branch("u_in"         ,&_root_u_in          ,"u_in/D");
  _rootMscTree->Branch("v_in"         ,&_root_v_in          ,"v_in/D");
  _rootMscTree->Branch("u_out"        ,&_root_u_out         ,"u_out/D");
  _rootMscTree->Branch("v_out"        ,&_root_v_out         ,"v_out/D");

  _rootMscTree->Branch("theta1_val"      ,&_root_angle1      ,"theta1_val/D"); 
  _rootMscTree->Branch("theta2_val"      ,&_root_angle2      ,"theta2_val/D");
  _rootMscTree->Branch("theta1_err"      ,&_root_angle1_err  ,"theta1_err/D");
  _rootMscTree->Branch("theta2_err"      ,&_root_angle2_err  ,"theta2_err/D");
  _rootMscTree->Branch("momentum"        ,&_root_momentum    ,"momentum/D");
  

}


} // Namespace


