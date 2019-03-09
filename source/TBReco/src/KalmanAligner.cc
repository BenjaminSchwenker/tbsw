// KalmanAligner Processor 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// Local includes 
#include "KalmanAligner.h"

// TBTools includes
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
                          _trackCollectionName,std::string("tracks"));
   
// 
// Processor parameters
  
  registerProcessorParameter ("AlignmentDBFileName",
                             "This is the name of the file with the alignment constants (add .root)",
                             _alignmentDBFileName, static_cast< string > ( "alignmentDB.root" ) ); 

  registerProcessorParameter ("LogLevel", "LogLever during alignment",
                              _logLevel,  static_cast < int > (2));
  
  registerProcessorParameter ("MaxTracks", "Maximum number of tracks passed to alignemnt",
                              _maxTracks,  static_cast < int > (70000)); 

  registerProcessorParameter ("UseBeamConstraint", "Use beam model to constrain track fitting",
                              _useBC,  static_cast < bool > (false));
  
  registerProcessorParameter ("pValueCut", "P-Value cut for tracks used during alignment",
                              _pValueCut,  static_cast < double > (0.0)); 

  registerProcessorParameter ("DeviationCut", "Reject alignment corrections exceeding DeviationCut*Sigma",
                              _deviationCut,  static_cast < double > (0.0)); 

  registerProcessorParameter ("AnnealingTracks",
                              "Number of tracks before the annealign is turned OFF",
                              _annealingTracks,  static_cast < int > (4000));

  registerProcessorParameter ("AnnealingFactor", "Scale factor for annealing schedule",
                              _annealingFactor,  static_cast < double > (10000.));
  
  std::vector<float> initErrorsShiftX;
  initErrorsShiftX.push_back(0.0);
  registerProcessorParameter("ErrorsShiftX", "Initial errors on alignment x shift [mm] for sensors ordered along beam line.",
                              _errorsShiftX, initErrorsShiftX );

  std::vector<float> initErrorsShiftY;
  initErrorsShiftY.push_back(0.0);
  registerProcessorParameter("ErrorsShiftY", "Initial errors on alignment y shift [mm] for sensors ordered along beam line.",
                              _errorsShiftY, initErrorsShiftY );
  
  std::vector<float> initErrorsShiftZ;
  initErrorsShiftZ.push_back(0.0);
  registerProcessorParameter("ErrorsShiftZ", "Initial errors on alignment z shift [mm] for sensors ordered along beam line.",
                              _errorsShiftZ, initErrorsShiftZ );

  std::vector<float> initErrorsAlpha;
  initErrorsAlpha.push_back(0.0);
  registerProcessorParameter("ErrorsAlpha", "Initial errors on alignment alpha [rad] for sensors ordered along beam line.",
                              _errorsAlpha, initErrorsAlpha );

  std::vector<float> initErrorsBeta;
  initErrorsBeta.push_back(0.0);
  registerProcessorParameter("ErrorsBeta", "Initial errors on alignment beta [rad] for sensors ordered along beam line.",
                              _errorsBeta, initErrorsBeta );

  std::vector<float> initErrorsGamma;
  initErrorsGamma.push_back(0.0);
  registerProcessorParameter("ErrorsGamma", "Initial errors on alignment gamma [rad] for sensors ordered along beam line.",
                              _errorsGamma, initErrorsGamma );

  registerProcessorParameter ("UpdateAlignment",
                              "Update alignment DB using offset corrections (true/false)?",
                              _updateAlignment, static_cast <bool> (true) ); 

  registerProcessorParameter ("NewAlignment",
                              "Start alignment from scratch (true/false)?",
                              _newAlignment, static_cast <bool> (false) ); 
                                
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
  if(!_newAlignment) _detector.ReadAlignmentDB( _alignmentDBFileName );
  // This is needed, because if the AlignmentDB is not read, the detector construct doesn't know the alignmentDB name
  else  _detector.SetAlignmentDBName( _alignmentDBFileName );     
  
  if ( (int)_errorsShiftX.size() != _detector.GetNSensors() ) {
    _errorsShiftX.resize(_detector.GetNSensors(), 0.0);
    streamlog_out ( MESSAGE3 ) <<  "Bad steering file: Parameter errorsShiftX has wrong size. Resize using default error 0.0." << endl;  
  } 
    
  if ( (int)_errorsShiftY.size() != _detector.GetNSensors() ) {
    _errorsShiftY.resize(_detector.GetNSensors(), 0.0);
    streamlog_out ( MESSAGE3 ) <<  "Bad steering file: Parameter errorsShiftY has wrong size. Resize using default error 0.0." << endl;  
  } 
  
  if ( (int)_errorsShiftZ.size() != _detector.GetNSensors() ) {
    _errorsShiftZ.resize(_detector.GetNSensors(), 0.0);
    streamlog_out ( MESSAGE3 ) <<  "Bad steering file: Parameter errorsShiftZ has wrong size. Resize using default error 0.0." << endl;  
  } 
  
  if ( (int)_errorsAlpha.size() != _detector.GetNSensors() ) {
    _errorsAlpha.resize(_detector.GetNSensors(), 0.0);
    streamlog_out ( MESSAGE3 ) <<  "Bad steering file: Parameter errorsAlpha has wrong size. Resize using default error 0.0." << endl;  
  } 
  
  if ( (int)_errorsBeta.size() != _detector.GetNSensors() ) {
    _errorsBeta.resize(_detector.GetNSensors(), 0.0);
    streamlog_out ( MESSAGE3 ) <<  "Bad steering file: Parameter errorsBeta has wrong size. Resize using default error 0.0." << endl;  
  } 
  
  if ( (int)_errorsGamma.size() != _detector.GetNSensors() ) {
    _errorsGamma.resize(_detector.GetNSensors(), 0.0);
    streamlog_out ( MESSAGE3 ) <<  "Bad steering file: Parameter errorsGamma has wrong size. Resize using default error 0.0." << endl;  
  } 
  
  //////////////////////////////////////////////////////////////////////
  // Alignment Data I/O 
  
  _nKAATracks = 0;
      
  // Open buffer file for storing alignment data 
  alignment_data = new TFile("tmpAlignTracks.root", "RECREATE");
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
void KalmanAligner::check( LCEvent * )
{
}

//
// Method called after all data processing
//
void KalmanAligner::end()
{
  
  ///////////////////////////////////////////////////////////
  // Construct the initial alignment state 
  
  int nSensors = _detector.GetNSensors();
  AlignableDet AlignState(nSensors);
  
  for (int iSensor=0; iSensor < nSensors; iSensor++){     

    SensorAlignmentParameters alignParams;
    alignParams << 0, 0, 0 , 0, 0, 0;

    SensorAlignmentCovariance alignCov = SensorAlignmentCovariance::Zero();   
    alignCov(0,0) = _errorsShiftX[iSensor]*_errorsShiftX[iSensor];
    alignCov(1,1) = _errorsShiftY[iSensor]*_errorsShiftY[iSensor];
    alignCov(2,2) = _errorsShiftZ[iSensor]*_errorsShiftZ[iSensor];
    alignCov(3,3) = _errorsAlpha[iSensor]*_errorsAlpha[iSensor]; 
    alignCov(4,4) = _errorsBeta[iSensor]*_errorsBeta[iSensor];
    alignCov(5,5) = _errorsGamma[iSensor]*_errorsGamma[iSensor];
   
    AlignState.SetAlignState(iSensor, alignParams);
    AlignState.SetAlignCovariance(iSensor, alignCov);     
  }
   
  ////////////////////////////////////////////////////////////
  // Try to fit alignment corrections from track residuals  
  
  TBDetector tmp_detector = _detector;
  
  streamlog_out ( MESSAGE3 ) << " Total of " << _nKAATracks << " tracks found" << endl;
  streamlog_out ( MESSAGE3 ) << endl;
  streamlog_out ( MESSAGE3 ) << "Starting alignment ..." << endl;

  KalmanAlignmentAlgorithm2 Aligner;
  AlignableDet reco_const = Aligner.Fit(tmp_detector, alignment_data, AlignState, _maxTracks, _annealingTracks, _annealingFactor,  _pValueCut, _deviationCut, _useBC, _logLevel );
  
  bool error_fim = Aligner.AlignDetector(tmp_detector, reco_const);
  if ( error_fim ) {
    streamlog_out ( MESSAGE3 ) << "Alignment failed!" << endl;
  } 

  // Close alignment_data file
  alignment_data->Close();
  delete alignment_data;
    
  ////////////////////////////////////////////////////////////
  // Print alignment results     
  
  streamlog_out ( MESSAGE3 )  << endl << "Alignment constants:" << endl << endl; 
    
  for ( int ipl = 0; ipl < _detector.GetNSensors(); ipl++ ) {
    
    // Print final geometry constants 
    // ------------------------------  
    
    // This is the position vector of the sensor
    auto pos_f = tmp_detector.GetDet(ipl).GetNominal().GetPosition(); 
    
    // This is the rotation matrix of the sensor; it 
    // contains a discrete and a continuous factor. 
    auto Rot_f = tmp_detector.GetDet(ipl).GetNominal().GetRotation();

    // This is the discrete factor of sensor rotation. 
    auto DRot = tmp_detector.GetDet(ipl).GetDiscrete().GetRotation();
    
    // This is finally the continous factor of the rotation
    auto CRot_f = Rot_f*DRot.transpose(); 
    
    // Euler angles are defined wrt. the continous rotation 
    double alpha_f, beta_f, gamma_f; 
    GetAnglesKarimaki(CRot_f, alpha_f, beta_f, gamma_f); 
    
    // Print 
    streamlog_out ( MESSAGE3 ) << endl << "Sensor plane " << ipl << endl << endl;  
    streamlog_out ( MESSAGE3 ) << endl << "  final x [mm] " << pos_f[0] << " +/- " << std::sqrt( reco_const.GetAlignCovariance(ipl)(0,0) )  << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  final y [mm] " << pos_f[1] << " +/- " << std::sqrt( reco_const.GetAlignCovariance(ipl)(1,1) ) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  final z [mm] " << pos_f[2] << " +/- " << std::sqrt( reco_const.GetAlignCovariance(ipl)(2,2) ) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  final alpha [rad] " << alpha_f << " +/- " << std::sqrt( reco_const.GetAlignCovariance(ipl)(3,3) ) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  final beta  [rad] " << beta_f << " +/- " << std::sqrt( reco_const.GetAlignCovariance(ipl)(4,4) ) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  final gamma [rad] " << gamma_f << " +/- " << std::sqrt( reco_const.GetAlignCovariance(ipl)(5,5) ) << endl; 

    // Print incremental alignment corrections 
    // ---------------------------------------  
    
    // Initial alignment 
    auto Rot_i = _detector.GetDet(ipl).GetNominal().GetRotation(); 
    auto pos_i = _detector.GetDet(ipl).GetNominal().GetPosition(); 

    // Diff rotation  
    auto diffRot = Rot_f*Rot_i.transpose();
    double dalpha, dbeta, dgamma; 
    GetAnglesKarimaki(diffRot, dalpha, dbeta, dgamma);
       
    // Diff shift   
    auto dr = pos_f - pos_i;
     
    // Print 
    streamlog_out ( MESSAGE3 ) << endl << "Sensor plane " << ipl << endl << endl;  
    streamlog_out ( MESSAGE3 ) << endl << "  correction dx [mm] " << dr[0] << " +/- " << std::sqrt( reco_const.GetAlignCovariance(ipl)(0,0) ) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  correction dy [mm] " << dr[1] << " +/- " << std::sqrt( reco_const.GetAlignCovariance(ipl)(1,1) ) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  correction dz [mm] " << dr[2] << " +/- " << std::sqrt( reco_const.GetAlignCovariance(ipl)(2,2) ) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  correction dalpha [rad] " << dalpha << " +/- " << std::sqrt( reco_const.GetAlignCovariance(ipl)(3,3) ) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  correction dbeta [rad] " << dbeta << " +/- " << std::sqrt( reco_const.GetAlignCovariance(ipl)(4,4) ) << endl; 
    streamlog_out ( MESSAGE3 ) << endl << "  correction dgamma [rad] " << dgamma << " +/- " << std::sqrt( reco_const.GetAlignCovariance(ipl)(5,5) ) << endl; 
      
  }
       
  //////////////////////////////////////////////////////////////////////  
  // Create aligned detector 
      
  _detector = tmp_detector; 
   
  if ( _updateAlignment ) { 
    _detector.WriteAlignmentDB( ); 
  } else {
    streamlog_out ( MESSAGE3 ) << endl;
    streamlog_out ( MESSAGE3 ) << "NO UPDATE OF ALIGNMENT DB" << endl; 
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
