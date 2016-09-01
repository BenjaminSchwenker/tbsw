// KalmanAligner Processor 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// Local includes 
#include "KalmanAligner.h"

// DEPFETTrackTools includes
#include "TBTrack.h"
#include "TrackInputProvider.h"
#include "Utilities.h"
#include "ThreeDModel.h"
#include "KalmanAlignmentInputProvider.h"
#include "KalmanAlignmentAlgorithm2.h"

// C++ includes
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <string>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <Exceptions.h>

// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;
using namespace CLHEP;

namespace depfet {

//
// Instantiate this object
//
KalmanAligner aKalmanAligner ;

//
// Constructor
//
KalmanAligner::KalmanAligner() : Processor("KalmanAligner")
{
   
// Processor description
  _description = "KalmanAligner: Tracking telescope alignment using the Kalman Alignment Algorithm (KAA) with annealing";
   
//
// Input collections  
  registerInputCollection(LCIO::TRACK,"TrackCollectionName",
                          "Track collection for alignment",
                          _trackCollectionName,std::string("aligntracks"));
   
// 
// Processor parameters
  
  registerProcessorParameter ("AlignmentDBFileName",
                             "This is the name of the LCIO file with the alignment constants (add .slcio)",
                             _alignmentDBFileName, static_cast< string > ( "alignmentDB.slcio" ) ); 
  
  registerProcessorParameter ("AlignConfigFileName",
                             "Name of the alignment config file",
                             _alignConfigFileName, std::string("cfg/align.cfg"));
  
  registerProcessorParameter( "RootFileName",
                             "Root file containing alignment validation plots",
                             _rootFileName, std::string("KalmanAlign-TB2011.root"));
  
  registerProcessorParameter ("UpdateAlignment",
                              "Update lcio alignmentDB using alignment results (true/false)?",
                              _updateAlignment, static_cast <bool> (false) ); 
  
  
  
 
                                
}

//
// Method called at the beginning of data processing
//
void KalmanAligner::init() {
  
  // Initialize variables
  _nRun = 0 ;
  _nEvt = 0 ;
  
  // Print set parameters
  printProcessorParams();
  
  // CPU time start
  _timeCPU = clock()/1000;
  
  // Read detector constants from gear file
  _detector.ReadGearConfiguration();    
  
  // Read alignment data base file 
  _detector.ReadAlignmentDB( _alignmentDBFileName );       
    
  //////////////////////////////////////////////////////////////////////
  // Alignment Data I/O 
  
  _nKAATracks = 0;
      
  // Open buffer file for storing alignment data 
  alignment_data = new TFile(_rootFileName.c_str(), "RECREATE");
  if (alignment_data == 0 || alignment_data->IsOpen() != kTRUE) {
    streamlog_out ( ERROR4) << "Could not open alignment data file!!" << endl;
    exit(-1);
  }
     
  // Create alignment "event" which records all information from one track
  myEvent = new AlignEvent;
  // Create a tree where the alignment info is being stored
  AlignTree = new TTree("AlignTree", "Alignment data", 0);
  AlignTree->Branch("AlignEvent", & myEvent, 64000, 0);
  // Auto save every MB
  AlignTree->SetAutoSave(1000000);
  
}

//
// Method called for each run
//
void KalmanAligner::processRunHeader(LCRunHeader * run)
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
void KalmanAligner::processEvent(LCEvent * evt)
{
    
  //////////////////////////////////////////////////////////////////////  
  // Process next event
  ++_nEvt;
   
  if ( _nEvt % 1000 == 0 ) {
    streamlog_out( MESSAGE3 ) << "Processing event "
                              << evt->getEventNumber() << " in run "
                              << evt->getRunNumber() << endl; 
                              
    streamlog_out( MESSAGE3 ) << "Currently having " << _nKAATracks << " tracks" << endl;
  }
    
  
  LCCollection* collection;
  try {
      collection = evt->getCollection(_trackCollectionName);
  } catch (DataNotAvailableException& e) {
      throw SkipEventException(this);
  }
  
 
  
  // Main loop over all tracks
  
  int nTracks = collection->getNumberOfElements(); 
  TrackInputProvider TrackLCIOReader;  

  for (int itrk = 0; itrk < nTracks; itrk++) {
    
    // Retrieve track from LCIO 
    Track * lciotrk = dynamic_cast<Track*> (collection->getElementAt(itrk));
    
    // Convert LCIO -> TB track  
    TBTrack rectrack = TrackLCIOReader.MakeTBTrack( lciotrk, _detector );  
    
    // Added one track
    ++_nKAATracks;
    
    // Fill alignment container   
    alignment_data->cd("");
    
    KalmanAlignmentInputProvider kaip;
    kaip.FillEvent(rectrack, *myEvent);
    
    AlignTree->Fill();
                               
  } // End loop over all tracks 
        
}

//
// Method called after each event to check the data processed
//
void KalmanAligner::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void KalmanAligner::end()
{
   
  ////////////////////////////////////////////////////////////
  // Try to fit alignment corrections from track residuals  
  
  TBDetector tmp_detector = _detector;

  cout << " Total of " << _nKAATracks << " tracks found" << endl;
  
  streamlog_out ( MESSAGE3 ) << endl;
  streamlog_out ( MESSAGE3 ) << "Starting alignment ..." << endl;
  
  KalmanAlignmentAlgorithm2 Aligner;
  AlignableDet reco_const = Aligner.Fit(tmp_detector, alignment_data, _alignConfigFileName );
  
  bool error_fim = Aligner.AlignDetector(tmp_detector, reco_const);
  if ( error_fim ) {
    streamlog_out ( MESSAGE3 ) << "Alignment failed!" << endl;
  } 
    
  ////////////////////////////////////////////////////////////
  // Print detector geometrie after alignment   
  
  streamlog_out ( MESSAGE3 ) << "Detector geometry after alignment procedure: " << endl;   
  tmp_detector.Print();    
  
  ////////////////////////////////////////////////////////////
  // Visualize alignment corrections 
  int nsensor = _detector.GetNSensors();   
  
  TH1D * hxshift = new TH1D("hxshift", "X Shift", nsensor, 0, nsensor); 
  hxshift->SetXTitle("plane number");
  hxshift->SetYTitle("x shift [mm]");
  
  TH1D * hyshift = new TH1D("hyshift", "Y Shift", nsensor, 0, nsensor); 
  hyshift->SetXTitle("plane number");
  hyshift->SetYTitle("y shift [mm]");
  
  TH1D * hzshift = new TH1D("hzshift", "Z Shift", nsensor, 0, nsensor);  
  hzshift->SetXTitle("plane number");
  hzshift->SetYTitle("z shift [mm]");
  
  TH1D * hxrot = new TH1D("hxrot", "X Rotation", nsensor, 0, nsensor); 
  hxrot->SetXTitle("plane number");
  hxrot->SetYTitle("x rot [rad]");
  
  TH1D * hyrot = new TH1D("hyrot", "Y Rotation", nsensor, 0, nsensor); 
  hyrot->SetXTitle("plane number");
  hyrot->SetYTitle("y rot [rad]");
  
  TH1D * hzrot = new TH1D("hzrot", "Z Rotation", nsensor, 0, nsensor);  
  hzrot->SetXTitle("plane number");
  hzrot->SetYTitle("z rot [rad]");         

  streamlog_out ( MESSAGE3 )  << endl << "Final geometry constants:" << endl << endl; 
    
  HepSymMatrix reco_cov = reco_const.alignmentCovariance;
  
  for ( int ipl = 0; ipl < nsensor; ipl++ ) {
    
    // Print final geometry constants 
    // ------------------------------  
    
    // This is the position vector of the sensor
    HepVector pos_f = tmp_detector.GetDet(ipl).GetNominal().GetPosition(); 
    
    // This is the rotation matrix of the sensor; it 
    // contains a discrete and a continuous factor. 
    HepMatrix Rot_f = tmp_detector.GetDet(ipl).GetNominal().GetRotation();

    // This is the discrete factor of sensor rotation. 
    HepMatrix DRot = tmp_detector.GetDet(ipl).GetDiscrete().GetRotation();
    
    // This is finally the continous factor of the rotation
    HepMatrix CRot_f = Rot_f*DRot.T(); 
    
    // Euler angles are defined wrt. the continous rotation 
    double alpha_f, beta_f, gamma_f; 
    GetAnglesKarimaki(CRot_f, alpha_f, beta_f, gamma_f); 
    
    // Fill root histos
    hxshift->SetBinContent(ipl+1,pos_f[0]);
    hyshift->SetBinContent(ipl+1,pos_f[1]);
    hzshift->SetBinContent(ipl+1,pos_f[2]);
    hxrot->SetBinContent(ipl+1,alpha_f);
    hyrot->SetBinContent(ipl+1,beta_f);
    hzrot->SetBinContent(ipl+1,gamma_f);
     
    // Print 
    streamlog_out ( MESSAGE3 ) << endl << "Sensor plane " << ipl << endl << endl;  
    streamlog_out ( MESSAGE3 ) << endl << "  final x [mm] " << pos_f[0] << " +/- " << std::sqrt(reco_cov[ipl*6+0][ipl*6+0]) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  final y [mm] " << pos_f[1] << " +/- " << std::sqrt(reco_cov[ipl*6+1][ipl*6+1]) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  final z [mm] " << pos_f[2] << " +/- " << std::sqrt(reco_cov[ipl*6+2][ipl*6+2]) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  final alpha [rad] " << alpha_f << " +/- " << std::sqrt(reco_cov[ipl*6+3][ipl*6+3]) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  final beta  [rad] " << beta_f << " +/- " << std::sqrt(reco_cov[ipl*6+4][ipl*6+4]) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  final gamma [rad] " << gamma_f << " +/- " << std::sqrt(reco_cov[ipl*6+5][ipl*6+5]) << endl; 

    // Print incremental alignment corrections 
    // ---------------------------------------  
    
    // Initial alignment 
    HepMatrix Rot_i = _detector.GetDet(ipl).GetNominal().GetRotation(); 
    HepVector pos_i = _detector.GetDet(ipl).GetNominal().GetPosition(); 

    // Diff rotation  
    HepMatrix diffRot = Rot_f*Rot_i.T();
    double dalpha, dbeta, dgamma; 
    GetAnglesKarimaki(diffRot, dalpha, dbeta, dgamma);
       
    // Diff shift   
    HepVector dr = pos_f - pos_i;
     
    // Print 
    streamlog_out ( MESSAGE3 ) << endl << "Sensor plane " << ipl << endl << endl;  
    streamlog_out ( MESSAGE3 ) << endl << "  correction dx [mm] " << dr[0] << " +/- " << std::sqrt(reco_cov[ipl*6+0][ipl*6+0]) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  correction dy [mm] " << dr[1] << " +/- " << std::sqrt(reco_cov[ipl*6+1][ipl*6+1]) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  correction dz [mm] " << dr[2] << " +/- " << std::sqrt(reco_cov[ipl*6+2][ipl*6+2]) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  correction dalpha [rad] " << dalpha << " +/- " << std::sqrt(reco_cov[ipl*6+3][ipl*6+3]) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  correction dbeta [rad] " << dbeta << " +/- " << std::sqrt(reco_cov[ipl*6+4][ipl*6+4]) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  correction dgamma [rad] " << dgamma << " +/- " << std::sqrt(reco_cov[ipl*6+5][ipl*6+5]) << endl; 
      
  }
    
  // Write validation histos 
  alignment_data->Write();
  alignment_data->Close();
  delete alignment_data;
     
  //////////////////////////////////////////////////////////////////////  
  // Create aligned detector 
      
  _detector = tmp_detector; 
   
  if ( _updateAlignment ) { 
    _detector.WriteAlignmentDB( ); 
  } else {
    streamlog_out ( MESSAGE3 ) << endl;
    streamlog_out ( MESSAGE3 ) << "NO UPDATE OF ALIGNMENT DB LCIO FILE" << endl; 
  }
         
  streamlog_out ( MESSAGE3 ) << endl;
  streamlog_out ( MESSAGE3 ) << "Successfully finished" << endl;
  
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

 
}


//
// Method printing processor parameters
//
void KalmanAligner::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "KalmanAligner Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}

} // Namespace
