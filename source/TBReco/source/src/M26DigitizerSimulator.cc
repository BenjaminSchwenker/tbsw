// M26DigitizerSimulator
// Validation tool for track fitting routines 
//                       
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// user includes
#include "M26DigitizerSimulator.h"


// DEPFETTrackTools includes 
#include "ParticleGun.h"
#include "GenericTrackFitter.h"
#include "TBHit.h"
#include "SeedGenerator.h"
#include "Utilities.h"
#include "KalmanAlignmentAlgorithm2.h"

// marlin includes
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"

// system includes
#include <string> 
#include <list>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace marlin;
using namespace CLHEP;

namespace depfet {

M26DigitizerSimulator::M26DigitizerSimulator ():DataSourceProcessor  ("M26DigitizerSimulator") {
  
  _description =
    "Validation tool M26 digitizer model";
   
  registerProcessorParameter ("NumberOfTracks",
                              "How many tracks for simulation",
                              _maxTracks,  static_cast < int > (1000)); 
    
  registerProcessorParameter ("ParticleMomentum", "Particle momentum [GeV]",
                              _momentum,  static_cast < double > (120.0));
   
  registerProcessorParameter ("ParticleMass", "Particle mass [GeV]",
                              _mass,  static_cast < double > (0.139));
   
  registerProcessorParameter ("ParticleCharge", "Particle charge [e]",
                              _charge,  static_cast < double > (+1));
  
  registerProcessorParameter ("GunPositionX", "X position of particle gun [mm]",
                              _GunXPosition,  static_cast < double > (0));
  
  registerProcessorParameter ("GunPositionY", "Y position of particle gun [mm]",
                              _GunYPosition,  static_cast < double > (0));

  registerProcessorParameter ("GunPositionZ", "Z position of particle gun [mm]",
                              _GunZPosition,  static_cast < double > (-5000));
  
  registerProcessorParameter ("GunRotationX", "X rotation of particle gun [rad]",
                              _GunRotX,  static_cast < double > (0)); 
  
  registerProcessorParameter ("GunRotationY", "Y rotation of particle gun [rad]",
                              _GunRotY,  static_cast < double > (0)); 

  registerProcessorParameter ("GunRotationZ", "Z rotation of particel gun [rad]",
                              _GunRotZ,  static_cast < double > (0)); 
  
  registerProcessorParameter ("GunSpotSizeX", "Smearing of X vertex position at beam collimator [mm]",
                              _GunSpotSizeX,  static_cast < double > (1)); 

  registerProcessorParameter ("GunSpotSizeY", "Smearing of Y vertex position at beam collimator [mm]",
                              _GunSpotSizeY,  static_cast < double > (1)); 
  
  registerProcessorParameter ("GunDivergenceX", "RMS track slope in XZ plane [rad]",
                              _GunDivergenceX,  static_cast < double > (0.0001)); 
  
  registerProcessorParameter ("GunDivergenceY", "RMS track slope in YZ plane [rad]",
                              _GunDivergenceY,  static_cast < double > (0.0001)); 
  
  registerProcessorParameter ("GunCorrelationX", "Beam correlation coefficient X",
                              _GunCorrelationX,  static_cast < double > (0.0)); 
  
  registerProcessorParameter ("GunCorralationY", "Beam correlation coefficient Y",
                              _GunCorrelationY,  static_cast < double > (0.0)); 
  
  registerProcessorParameter ("HitDigitizer", 
                              "Choose digitizer model for pixel detectors:\n"
                              "0: Digital Mimosa26 digitizer\n"
                              "1: Gaussian smearing digitizer\n"
                              "2: Box digitizer",
                              _digi_type,  static_cast < int > (1)); 
  
  registerProcessorParameter ("MultipleScattering", "Choose model of multiple scattering in the telescope:\n"
                              "0: Highland\n"
                              "1: Moliere",
                              _mscmodel,  static_cast < int > (0));  
  
  registerProcessorParameter( "RootFileName", "Output root file name",
                              _rootFileName, std::string("TrackingTest.root"));
  
}


M26DigitizerSimulator * M26DigitizerSimulator::newProcessor () {  
  return new M26DigitizerSimulator;
}


void M26DigitizerSimulator::init () {
  
  printParameters ();
  
  // Read detector constants from gear file
  _detector.ReadGearConfiguration();    
  
  // Print out geometry information
  streamlog_out ( MESSAGE3 ) <<  "Telescope configuration with " << _detector.GetNSensors() << " planes" << endl << endl;
   
  for(int ipl=0; ipl < _detector.GetNSensors(); ipl++) {
      
    streamlog_out( MESSAGE3 ) << "ID = " << _detector.GetDet(ipl).GetDAQID()
                              << " at Z [mm] = " <<  _detector.GetDet(ipl).GetNominal().GetPosition()[2]
                              << endl; 
                               
    // Print detector position  
    _detector.GetDet(ipl).GetNominal().PrintHepMatrix();  
  }
   
  // Book correlation histograms   
  bookHistos();   
}


void M26DigitizerSimulator::readDataSource (int Ntrig) {
       
   // Configure ParticleGun
   ParticleGun TBGun; 
   TBGun.positionX = _GunXPosition; 
   TBGun.positionY = _GunYPosition; 
   TBGun.positionZ = _GunZPosition; 
   TBGun.rotationX = _GunRotX; 
   TBGun.rotationY = _GunRotY; 
   TBGun.rotationZ = _GunRotZ; 
   TBGun.divergenceX = _GunDivergenceX; 
   TBGun.divergenceY = _GunDivergenceY;
   TBGun.correlationX = _GunCorrelationX;
   TBGun.correlationY = _GunCorrelationY;
   TBGun.spotsizeX = _GunSpotSizeX; 
   TBGun.spotsizeY = _GunSpotSizeY;
   TBGun.mscmodel = _mscmodel;
   TBGun.digi_type = _digi_type;
   TBGun.efficiency = 1;    
 
   // Configure Kalman track fitter
   GenericTrackFitter TrackFitter(_detector);
   TrackFitter.SetNumIterations(2);
   TrackFitter.SetOutlierCut(10000); 
      
   // Configure seed generator 
   SeedGenerator TrackSeeder(_charge, _momentum);
   
   int nsensor = _detector.GetNSensors(); 
   int nParameters = 6; 
       
   // Create a misaligned detector for simulation 
   // of tracks
   TBDetector misdetector = _detector;
         
   AlignableDet MisalignStoreMod(nsensor,nParameters);
   
     
   for ( int ipl = 0; ipl < nsensor; ipl++ ) {    
     // Shifts, mm
     MisalignStoreMod.alignmentParameters[ipl*nParameters + 0] = gRandom->Uniform(-0, 0);
     MisalignStoreMod.alignmentParameters[ipl*nParameters + 1] = gRandom->Uniform(-0, 0);
     MisalignStoreMod.alignmentParameters[ipl*nParameters + 2] = gRandom->Uniform(-0, 0);
     // Tilts, rad 
     MisalignStoreMod.alignmentParameters[ipl*nParameters + 3] = gRandom->Uniform(-0, 0);
     MisalignStoreMod.alignmentParameters[ipl*nParameters + 4] = gRandom->Uniform(-0, 0);
     MisalignStoreMod.alignmentParameters[ipl*nParameters + 5] = gRandom->Uniform(-0.0001, 0.0001);
   }
     
   // Use same mechnism for misalignment and alignment 
   KalmanAlignmentAlgorithm2 ModMisAligner; 
   bool miserr_mod = ModMisAligner.AlignDetector(misdetector, MisalignStoreMod);
   if ( miserr_mod ) {
     cout << "BIG ERROR: Misalignmemnt of detector failed. Skipping!!" << endl;
    
   } 
      
   // This is the nominal detector. Tracking 
   // is done in this (wrong) detector
   TBDetector nomdetector = misdetector;  
       
   
   for(int eventNumber = 0; eventNumber < _maxTracks; ++eventNumber) { 
       
     // Print event number 
     if (eventNumber%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                       << eventNumber
                                                       << std::endl << std::endl;
     
     streamlog_out(MESSAGE1) << "Events processed: " << eventNumber << endl;
     
     // Store simulated hits along the particle 
     // trajectory for track fitting 
     vector<TBHit> HitStore;  
      
     // Generate a beam particle 
     TBTrack TruthTrack = TBGun.SimulateTrajectory(misdetector,_mass, _charge, _momentum);
     
    
     // Simulate hits along particle trajectory
     TBGun.SimulateTrackHits(TruthTrack, HitStore);
     
     // Start rekonstruction of particle 
     // trajectory from hits. 
     
     // Require at least 3 hits 
     if ( HitStore.size() < 3 ) continue; 
           
     // Init reco track  
     TBTrack RecoTrack(misdetector);
     RecoTrack.SetMass( _mass );
     RecoTrack.SetCharge( _charge );
     RecoTrack.SetMomentum( _momentum ); 
      
     // Add hits to track 
     for (int ihit=0; ihit<(int)HitStore.size(); ++ihit) {
       int daqid = HitStore[ihit].GetDAQID();
       int ipl = _detector.GetPlaneNumber(daqid); 
       RecoTrack.GetTE(ipl).SetHit(HitStore[ihit]); 
     }
     
     // Compute seed track 
     TBTrackState Seed = TrackSeeder.CreateSeedTrack(HitStore[0], HitStore[1], misdetector);   
     RecoTrack.SetReferenceState(Seed);
      
     // Try to run track fitter
     bool trkerr = TrackFitter.Fit(RecoTrack);
     if ( trkerr ) {
       streamlog_out ( MESSAGE3 ) << "Fit with " << RecoTrack.GetNumHits() << " hits failed. Skipping track!" << endl;
       continue;
     }  
     
     // Fill hit tree
     
     for(int ipl = 0; ipl < _detector.GetNSensors(); ipl++)  { 
          
       if ( RecoTrack.GetTE(ipl).HasHit() ) {
         
         Det& dut = RecoTrack.GetTE(ipl).GetDet();
         TBHit& hit = RecoTrack.GetTE(ipl).GetHit();
         
         _rootDetectorID = dut.GetDAQID(); 
         _rootHitU = hit.GetCoord()[0][0];         
         _rootHitV = hit.GetCoord()[1][0];   
         _rootHitCol = dut.GetColumnFromCoord( _rootHitU, _rootHitV );  
         _rootHitRow = dut.GetRowFromCoord( _rootHitU, _rootHitV );  
         _rootHitQuality = 0; 
         _rootHitCharge = -1; 
         _rootHitSeedCharge = -1; 
         _rootHitSize = -1;  
         _rootHitSizeCol = -1;     
         _rootHitSizeRow = -1;      
         
         // matched track
         _rootHitHasTrack = 0;    	     
         _rootHitLocalChi2 = RecoTrack.GetTE(ipl).GetChiSqu();  
                  
         _rootHitTruthU = TruthTrack.GetTE(ipl).GetState().GetPars()[2][0];
         _rootHitTruthV = TruthTrack.GetTE(ipl).GetState().GetPars()[3][0];  

         HepMatrix p = RecoTrack.GetTE(ipl).GetState().GetPars();
         HepSymMatrix C = RecoTrack.GetTE(ipl).GetState().GetCov();  
           
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
         _rootHitTrackChi2 = RecoTrack.GetChiSqu(); 
         _rootHitTrackNDF = RecoTrack.GetNDF();
         
         // Fill tree with set variables 
         _rootFile->cd("");
         _rootHitTree->Fill();
        
       }
         
     }
      
   } // End loop over events  
}

void M26DigitizerSimulator::end () {
   streamlog_out(MESSAGE3) << std::endl
                           << "Processor succesfully finished!"
                           << std::endl;
    // ROOT Output
  _rootFile->Write();
  _rootFile->Close();
}


void M26DigitizerSimulator::bookHistos() 
{
  
  streamlog_out ( MESSAGE4 ) <<  "Booking histograms" << endl;
  
  // ROOT Output 
  _rootFile = new TFile(_rootFileName.c_str(),"recreate"); 
  _rootFile->cd("");

  _rootHitTree = new TTree("Hit","Hit info");    
 
  _rootHitTree->Branch("det"             ,&_rootDetectorID       ,"det/I");
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
  _rootHitTree->Branch("u_truth"         ,&_rootHitTruthU             ,"u_truth/D");
  _rootHitTree->Branch("v_truth"         ,&_rootHitTruthV             ,"v_truth/D"); 

    
}


} // Namespace


