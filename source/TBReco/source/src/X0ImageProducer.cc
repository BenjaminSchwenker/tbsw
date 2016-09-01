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
#include "GenericTrackFitter.h"
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
                           "TrackCollection" ,
                           "Name of track collection"  ,
                           _trackColName ,
                           std::string("tracks") ) ;    
  
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
 
}

//
// Method called at the beginning of data processing
//
void X0ImageProducer::init() {
   
   // Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
   _timeCPU = clock()/1000 ;
   
   // Initialize all the counters
   _noOfTracks = 0;
        
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
  
  // Configure Kalman track fitter
  GenericTrackFitter TrackFitter(_detector);
  TBKalmanMSC TrackFitterMSC(_detector);
  TrackFitter.SetNumIterations(2);
   
  
  // Store tracks 
  std::vector< TBTrack > TrackStore;
      
  int nTrack = 0;
  if(isTrackok) nTrack = trackcol->getNumberOfElements();
  if ( nTrack == 0)  ++_noOfEventWOInputTrack;  
    
  // Loop over tracks in input track collection and 
  // fill root trees.
    
  TrackInputProvider TrackLCIOReader;  
    
  for(int itrk=0; itrk< nTrack ; itrk++) {
     
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
      
  }
  
  // Fill event tree
  _rootRunNumber = evt->getRunNumber();  
  _rootEventNumber = evt->getEventNumber();  
  _rootNTelTracks = TrackStore.size();  
  _rootFile->cd("");
  _rootEventTree->Fill();    
  
  // Analyze tracks   
  int nsensor = _detector.GetNSensors();   

  for(int itrk=0;itrk<(int)TrackStore.size(); ++itrk) {
     
    TBTrack& trk = TrackStore[itrk];
     
    TBTrack ForwardTrack(trk);
    TBTrack BackwardTrack(trk); 
    
    //Forward and backward direction
    int forwarddir = 1;
    int backwarddir =-1;
    //No bias: Hit on central plane won't be used
    bool forwardbias = 0;
    bool backwardbias = 0;
    
    // Try to run a forward Kalman Filter (without bias) on ForwardTrack
    //
    // Only the hits on sensor planes 0,1,2 (first telescope arm) are used
    for( int ipl = 3; ipl < nsensor; ipl++ )
    {
	ForwardTrack.GetTE(ipl).RemoveHit();
    }
    // Fit the track with the hits that are left (in this case the hits in the upstream telescope arm)
    TrackFitterMSC.ProcessTrack(ForwardTrack, forwarddir, forwardbias);

    // Try to run a backward Kalman Filter (without bias) on BackwardTrack
    //
    // Only the hits on sensors planes nsensor-3, nsensor-2, nsensor-1 (second telescope arm will be used)
    for( int ipl = 0; ipl < nsensor-3; ipl++ )
    {
	BackwardTrack.GetTE(ipl).RemoveHit();
    }
     
    // Fit the track with the hits that are left (in this case the hits in the upstream telescope arm)
    TrackFitterMSC.ProcessTrack(BackwardTrack, backwarddir, backwardbias);   
    
    // comboChi2 is chi2 combination of track in upstream and downstream telescope arm
    double comboChi2 = ForwardTrack.GetChiSqu()+BackwardTrack.GetChiSqu(); 
    
    
    //MSC Analysis for the reconstructed angles
    //Here we use the In and Out State and the GetScatterKinks function of the TBKalmanMSC Class
    
    Det dut = _detector.GetDet(_idut);
     
    // In and OutStates of the reconstructed Track at the current detector
    TBTrackState& InState=ForwardTrack.GetTE(_idut).GetState();
    TBTrackState& OutState=BackwardTrack.GetTE(_idut).GetState(); 
    
    //Angles and angle errors
    HepMatrix theta(2,1,0);
    HepSymMatrix Cov(2,0);
    theta = TrackFitterMSC.GetScatterKinks(dut, InState, OutState); 
    Cov = TrackFitterMSC.GetScatterKinkCov(dut, InState, OutState);
    
    // Get the track parameters of the fitted track on the current sensor
    // The u and v positions are needed for a position-resolved measurement
    HepMatrix p = trk.GetTE(_idut).GetState().GetPars();
    HepMatrix p_in = InState.GetPars();
    HepMatrix p_out = OutState.GetPars();
     	
    // Fill root variables 
    _rootDaqID = dut.GetDAQID(); 
    _rootPlaneID = _idut;
    _root_momentum = trk.GetMomentum(); 
    _rootTrackHits = trk.GetNumHits();
    _rootTrackChi2 = trk.GetChiSqu(); 
    _rootTrackProb = TMath::Prob(trk.GetChiSqu(),trk.GetNDF());
    _rootTrackProbUp = TMath::Prob(ForwardTrack.GetChiSqu(),ForwardTrack.GetNDF());
    _rootTrackProbDown = TMath::Prob(BackwardTrack.GetChiSqu(),BackwardTrack.GetNDF());
    _rootTrackProbCombo = TMath::Prob( comboChi2 ,ForwardTrack.GetNDF()+BackwardTrack.GetNDF());
    
    _root_x = trk.GetTE(0).GetState().GetPars()[2][0];   
    _root_y = trk.GetTE(0).GetState().GetPars()[3][0];   
    _root_dxdz = trk.GetTE(0).GetState().GetPars()[0][0];  
    _root_dydz = trk.GetTE(0).GetState().GetPars()[1][0];
    _root_u = p[2][0]; 
    _root_v = p[3][0]; 
    _root_u_in = p_in[2][0]; 
    _root_v_in = p_in[3][0];
    _root_u_out = p_out[2][0]; 
    _root_v_out = p_out[3][0];
    _root_dudw = p[0][0]; 
    _root_dvdw = p[1][0]; 
   
    _root_angle1 = theta[0][0];
    _root_angle2 = theta[1][0];
    _root_angle1_err = TMath::Sqrt(Cov[0][0]);
    _root_angle2_err = TMath::Sqrt(Cov[1][1]);
    
    _rootMscTree->Fill();       
    
  } // end track loop 		
     
  streamlog_out(MESSAGE2) << "Total of " << nTrack << " tracks in collection " << _trackColName << endl;      
  
  streamlog_out(MESSAGE2) << "Total of " << TrackStore.size() << " good tracks" << endl; 
  _noOfTracks += TrackStore.size();   
  
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
                           << "Total number of accepted Tel tracks:  " << setiosflags(ios::right) << _noOfTracks 
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
  _rootEventTree->Branch("nTelTracks"      ,&_rootNTelTracks     ,"nTelTracks/I");
    
  // 
  // Track Tree of the whole track 
  _rootMscTree = new TTree("MSCTree","Multiple scattering data");
  _rootMscTree->Branch("iRun"            ,&_rootRunNumber       ,"iRun/I");
  _rootMscTree->Branch("iEvt"            ,&_rootEventNumber     ,"iEvt/I"); 
  _rootMscTree->Branch("daqid"           ,&_rootDaqID           ,"daqid/I"); 
  _rootMscTree->Branch("ipl"             ,&_rootPlaneID         ,"ipl/I"); 
  _rootMscTree->Branch("nhits"           ,&_rootTrackHits       ,"nhits/I"); 
  _rootMscTree->Branch("chi2"            ,&_rootTrackChi2       ,"chi2/D"); 
  _rootMscTree->Branch("prob"            ,&_rootTrackProb       ,"prob/D"); 
  _rootMscTree->Branch("prob_up"         ,&_rootTrackProbUp     ,"prob_up/D"); 
  _rootMscTree->Branch("prob_down"       ,&_rootTrackProbDown   ,"prob_down/D"); 
  _rootMscTree->Branch("prob_combo"      ,&_rootTrackProbCombo  ,"prob_combo/D"); 

  _rootMscTree->Branch("x"            ,&_root_x             ,"x/D");
  _rootMscTree->Branch("y"            ,&_root_y             ,"y/D");
  _rootMscTree->Branch("dxdz"         ,&_root_dxdz          ,"dxdz/D");
  _rootMscTree->Branch("dydz"         ,&_root_dydz          ,"dydz/D");

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


