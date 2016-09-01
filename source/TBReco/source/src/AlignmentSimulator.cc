// AlignmentSimulator
// Validation tool for test beam alignment 
//                       
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include <CLHEP/Random/RandGamma.h>

// User includes
#include "AlignmentSimulator.h"

// DEPFETTrackTools includes 
#include "ParticleGun.h"
#include "TBEvtGen.h"
#include "GenericTrackFitter.h"
#include "TBHit.h"
#include "SeedGenerator.h"
#include "Utilities.h"
#include "ThreeDModel.h"
#include "KalmanAlignmentInputProvider.h"  
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

AlignmentSimulator::AlignmentSimulator ():DataSourceProcessor  ("AlignmentSimulator") {
  
  _description =
    "Validation tool for test beam alignment codes";
   
  registerProcessorParameter ("AlignConfigFileName",
                             "Name of the alignment config file",
                             _alignConfigFileName, std::string("cfg/align.cfg"));
  
  registerProcessorParameter ("NumberOfTracks",
                              "How many tracks for simulation",
                              _nTracks,  static_cast < int > (10000));
     
  registerProcessorParameter ("NumberOfSetups",
                              "How many misaligned detector setups",
                              _nMisalign,  static_cast < int > (100));   
  
  registerProcessorParameter ("RandomShiftErrorX",
                              "Uniform smear from -error,..,+error [mm] of sensor x position ",
                              _PositionAlignmentErrorX,  static_cast < double > (0.1));   
  
  registerProcessorParameter ("RandomShiftErrorY",
                              "Uniform smear from -error,..,+error [mm] of sensor y position ",
                              _PositionAlignmentErrorY,  static_cast < double > (0.1));   
  
  registerProcessorParameter ("RandomShiftErrorZ",
                              "Uniform smear from -error,..,+error [mm] of sensor z position ",
                              _PositionAlignmentErrorZ,  static_cast < double > (0.1));   
  
  registerProcessorParameter ("RandomRotationErrorU",
                              "Uniform smear from -error,..,+error [rad] of rotation around detector U axis ",
                              _RotationAlignmentErrorU,  static_cast < double > (0.1));   
  
  registerProcessorParameter ("RandomRotationErrorV",
                              "Uniform smear from -error,..,+error [rad] of rotation around detector V axis ",
                              _RotationAlignmentErrorV,  static_cast < double > (0.1)); 
  
  registerProcessorParameter ("RandomRotationErrorW",
                              "Uniform smear from -error,..,+error [rad] of rotation around detector W axis ",
                              _RotationAlignmentErrorW,  static_cast < double > (0.1));     
   
  std::vector<int> initDontMisAlignPlanes;
  registerProcessorParameter ("DontMisAlign",
                              "Enter plane numbers kept in ideal positions, not misaligned",
                              _DontMisAlignPlanes, initDontMisAlignPlanes);

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

  registerProcessorParameter ("GunSpotSizeX", "Smearing of X vertex position at beam collimator [mm]",
                              _GunSpotSizeX,  static_cast < double > (1)); 

  registerProcessorParameter ("GunSpotSizeY", "Smearing of Y vertex position at beam collimator [mm]",
                              _GunSpotSizeY,  static_cast < double > (1)); 
  
  registerProcessorParameter ("GunRotationX", "X rotation of particle gun [rad]",
                              _GunRotX,  static_cast < double > (0)); 
  
  registerProcessorParameter ("GunRotationY", "Y rotation of particle gun [rad]",
                              _GunRotY,  static_cast < double > (0)); 
   
  registerProcessorParameter ("GunRotationZ", "Z rotation of particle gun [rad]",
                              _GunRotZ,  static_cast < double > (0)); 
  
  registerProcessorParameter ("GunDivergenceX", "RMS track slope in XZ plane [rad]",
                              _GunDivergenceX,  static_cast < double > (0.0001)); 
  
  registerProcessorParameter ("GunDivergenceY", "RMS track slope in YZ plane [rad]",
                              _GunDivergenceY,  static_cast < double > (0.0001)); 

  registerProcessorParameter ("GunCorrelationX", "Beam correlation coefficient X",
                              _GunCorrelationX,  static_cast < double > (0.0)); 
  
  registerProcessorParameter ("GunCorralationY", "Beam correlation coefficient Y",
                              _GunCorrelationY,  static_cast < double > (0.0)); 
   
  registerProcessorParameter ("GunIntensity", "Number of particles per second",
                              _beamIntensity,  static_cast < double > (10000)); 

  registerProcessorParameter ("MaterialBeforeTelescope", "Thickness of material before telescope [X/X0]",
                              _BetheHeitlerT0,  static_cast < double > (0.0)); 
  
  registerProcessorParameter ("TriggerMode", "0:front sct, 1:back sct, 2: all sct",
                              _trgmode,  static_cast < int > (0)); 
  
  registerProcessorParameter ("HitDigitizer", 
                              "Choose digitizer model for pixel detectors:\n"
                              "0: Digital Mimosa26 digitizer\n"
                              "1: Gaussian smearing digitizer\n"
                              "2: Box digitizer",
                              _digi_type,  static_cast < int > (1)); 
     
  registerProcessorParameter ("DetectionEfficiency", "Hit detection efficiency of sub detectors [%]",
                              _hitEfficiency,  static_cast < double > (100)); 
  
  registerProcessorParameter ("MultipleScattering", "Choose model of multiple scattering in the telescope:\n"
                              "0: Highland\n"
                              "1: Moliere",
                              _mscmodel,  static_cast < int > (0));  
  
  registerProcessorParameter( "RootFileName", "Output root file name",
                              _rootFileName, std::string("AlignerTest.root"));
  
}


AlignmentSimulator * AlignmentSimulator::newProcessor () {  
  return new AlignmentSimulator;
}


void AlignmentSimulator::init () {
   
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
    
  // Rescale to intervall 0,..,1
  _hitEfficiency /= 100; 
  
  // Read reference sensors   
  
  _isMisAligned.resize(_detector.GetNSensors(), true);  

  for(int iRef=0; iRef<(int)_DontMisAlignPlanes.size(); ++iRef) {
    int ipl = _DontMisAlignPlanes[iRef];
    streamlog_out ( MESSAGE3 ) << "Do not misalign plane:  " << ipl << endl; 
    _isMisAligned[ipl] = false; 
  }

  
   
  // Book histograms   
  bookHistos();   
}


void AlignmentSimulator::readDataSource (int Ntrig) {
        
   // The idea is to simulate nMisalign variations 
   // of the nominal detector. These variations are
   // called misaligned detectors. Then, we simulate
   // nTracks tracks for each misaligned detector.
   // Then we try to to estimate the corrections to 
   // go from the nominal detector to the misaligned 
   // detector from track data. This step is called 
   // alignmnet. Finally, we compute some quality 
   // indicators for alignment corrections.  
    
   streamlog_out(MESSAGE3) << "Simulate " << _nMisalign << " misaligned detector(s) " 
                           << " with " << _nTracks << " tracks."<< endl; 
    
   
   
   streamlog_out(MESSAGE3)  << "Position error X [mm] "  << _PositionAlignmentErrorX << endl
                            << "Position error Y [mm] "  << _PositionAlignmentErrorY << endl
                            << "Position error Z [mm] "  << _PositionAlignmentErrorZ << endl
                            << "Rotation error U axis [rad] " << _RotationAlignmentErrorU << endl 
                            << "Rotation error V axis [rad] " << _RotationAlignmentErrorV << endl
                            << "Rotation error W axis [rad] " << _RotationAlignmentErrorW << endl
                            << endl << endl;   
   
   //////////////////////////////////////////////////////////////////////
   // Particle gun settings   
   // 
   // Please note that the each gun object initializes an internal random 
   // number generator to sample initial track state, scatter angles and 
   // errors of position measurements. 
   // 
   // In order to have independent track samples for each misalignment, the 
   // gun must be initialized here, before the loop over misaligned setups.   
   
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
   TBGun.efficiency = _hitEfficiency;    
     
   // Configure Trigger Simulation
   TBEvtGen TriggerSim;
   // Particler per sec 
   TriggerSim.intensity = _beamIntensity; 
   // Select trigger mode
   TriggerSim.trgmode = _trgmode;  
   
   // Main loop over misaligned setups  
    
   for(int iMisalign = 0; iMisalign < _nMisalign; ++iMisalign) { 
           
     streamlog_out(MESSAGE3) << endl << endl
                             << "Misaligned processed: " << iMisalign << " !!!!!!!!!!!"
                             << endl << endl; 
     
     
           
     //////////////////////////////////////////////////////////////////////
     // Detector Geometry 
     
     // Number of sub detectors/pixel modules
     int nsensor = _detector.GetNSensors(); 
      
     // This is the nominal detector. Initially, tracking 
     // is done in this (wrong) detector
     TBDetector nomdetector = _detector; 
     
     // Create a misaligned detector for simulation 
     // of tracks
     TBDetector misdetector = nomdetector;
     
     // We assume that sensors are randomly misaligned 
     // in shifts of sensor origin and rotations of 
     // sensor axes.  
     // 
     // We model all misalignments to be independent and
     // uncorrelated. That is, we do not model a joint 
     // displacement of telescope arms for example.
     // 
     // The starting point for small misalignments of 
     // sensors is the detector position as read from 
     // the gear file. This is the same situation as in 
     // test beams. 
      
    
     
     // Store true misalignment parameters for each detector for 
     // later comparison 
     int nAlignablesMod = nsensor;
     int nParameters = 6;   
 
     AlignableDet MisalignStoreMod(nAlignablesMod,nParameters);
     
     for ( int ipl = 0; ipl < nsensor; ipl++ ) { 
       
       // Reference sensors will not be misaligned. They remain 
       // exactly in the nominal position as given in the gear file.
        
       // Reference sensors can help to eliminate so called weak 
       // modes from the alignment problem. The term weak mode
       // refers to alignment corrections that cannot be determined
       // from cluster residuals.  
       
       if ( _isMisAligned[ipl] ) {
         // Shifts, mm
         MisalignStoreMod.alignmentParameters[ipl*nParameters + 0] = gRandom->Uniform(-_PositionAlignmentErrorX, _PositionAlignmentErrorX);
         MisalignStoreMod.alignmentParameters[ipl*nParameters + 1] = gRandom->Uniform(-_PositionAlignmentErrorY, _PositionAlignmentErrorY);
         MisalignStoreMod.alignmentParameters[ipl*nParameters + 2] = gRandom->Uniform(-_PositionAlignmentErrorZ, _PositionAlignmentErrorZ);
         // Tilts, rad 
         MisalignStoreMod.alignmentParameters[ipl*nParameters + 3] = gRandom->Uniform(-_RotationAlignmentErrorU, _RotationAlignmentErrorU);
         MisalignStoreMod.alignmentParameters[ipl*nParameters + 4] = gRandom->Uniform(-_RotationAlignmentErrorV, _RotationAlignmentErrorV);
         MisalignStoreMod.alignmentParameters[ipl*nParameters + 5] = gRandom->Uniform(-_RotationAlignmentErrorW, _RotationAlignmentErrorW);
       } 
        
     }
     
     // Use same mechnism for misalignment and alignment 
     KalmanAlignmentAlgorithm2 ModMisAligner; 
     bool miserr_mod = ModMisAligner.AlignDetector(misdetector, MisalignStoreMod);
     if ( miserr_mod ) {
        cout << "BIG ERROR: Misalignmemnt of detector failed. Skipping!!" << endl;
        continue;
     } 
      
     // Now, compute total misalignment of sensors
     HepVector truemisalign(nParameters*nsensor);
     
     for ( int ipl = 0; ipl < nsensor; ipl++ ) {
            
       // Initial detector  
       HepMatrix Rot_i = nomdetector.GetDet(ipl).GetNominal().GetRotation(); 
       HepVector pos_i = nomdetector.GetDet(ipl).GetNominal().GetPosition(); 
       
       // Misaligned detector 
       HepMatrix Rot_f = misdetector.GetDet(ipl).GetNominal().GetRotation();
       HepVector pos_f = misdetector.GetDet(ipl).GetNominal().GetPosition(); 
       
       // Diff rotation  
       HepMatrix dRot = Rot_f*Rot_i.T();
       
       // Diff shift   
       HepVector dr = pos_f - pos_i;
       
       truemisalign[ipl*nParameters + 0] = dr[0]; 
       truemisalign[ipl*nParameters + 1] = dr[1]; 
       truemisalign[ipl*nParameters + 2] = dr[2]; 
       
       double dalpha, dbeta, dgamma; 
       GetAnglesKarimaki(dRot, dalpha, dbeta, dgamma);
       truemisalign[ipl*nParameters + 3] = dalpha; 
       truemisalign[ipl*nParameters + 4] = dbeta; 
       truemisalign[ipl*nParameters + 5] = dgamma; 
     }
     
     
     //////////////////////////////////////////////////////////////////////
     // Alignment Data I/O 
     
     // Create a tree to buffer tracks for alignment 
     _rootFile->cd("");
     myAlignTree = new TTree("AlignTree", "Alignment data", 0);
     
     // Auto save every MB
     myAlignTree->SetAutoSave(1000000); 
     myAlignTree->Branch("AlignEvent", & myEvent, 64000, 0);  
     
     // Create a tree to buffer hits for pre- alignment  
     _rootFile->cd("");
     myHitTree = new TTree("HitTree","Hit data",0);
     
     // Auto save every MB
     myHitTree->SetAutoSave(1000000); 
     myHitTree->Branch("det", &myPlane, "det/I"); 
     myHitTree->Branch("u", &myHitU, "u/D");
     myHitTree->Branch("v", &myHitV, "v/D");
     
     //////////////////////////////////////////////////////////////////////
     // Simulate Tracks in misaligned (true) detector   
     
     // Configure Kalman track fitter
     GenericTrackFitter TrackFitter(_detector);
     TrackFitter.SetNumIterations(2); 
     
     // Configure seed generator 
     SeedGenerator TrackSeeder(_charge, _momentum);
     
     for(int iTrack = 0; iTrack < _nTracks; ++iTrack) {   
          
       //////////////////////////////////////////////////////////////////////
       // Simulate Track 
       
       // Store simulated hits along the particle 
       // trajectory for track fitting 
       vector<TBHit> HitStore;  

       // For electrons in PCMAG, Bethe Heitler model used to 
       // simulate energy loss by Bremsstrahlung
       double true_momentum = _momentum; 
       
       if ( _BetheHeitlerT0 > 0  ) {
         double t0 = _BetheHeitlerT0;
         double c = t0/TMath::Log(2);
         double u = double(RandGamma::shoot(c,1));
         double z = TMath::Exp(-u);    
         true_momentum *= z;  
       } 
       
       // Generate a beam particle 
       TBTrack TruthTrack = TBGun.SimulateTrajectory(misdetector,_mass, _charge, true_momentum);
       
       if (  !TriggerSim.ACCEPT( TruthTrack )  ) {
         continue; 
       } 
       
       // Simulate hits along particle trajectory
       TBGun.SimulateTrackHits(TruthTrack, HitStore);
        
       // Sample hit data for pre- alignment of 
       // detector
        
      
       for (int ihit=0; ihit<(int)HitStore.size(); ++ihit) {
         // Fill hit tree        
         int daqid = HitStore[ihit].GetDAQID();
         myPlane = nomdetector.GetPlaneNumber(daqid);
         myHitU = HitStore[ihit].GetCoord()[0][0];        
         myHitV = HitStore[ihit].GetCoord()[1][0];  
         myHitTree->Fill();  
       }
       
       // Sample track data for final alignment of 
       // detector.
       
       // Important: Track fitting is done in the 
       // nominal(!) detector. 

       // Require at least 2 hits for seeding 
       if ( HitStore.size() < 6 ) continue;           
       
       // Init reco track  
       TBTrack RecoTrack(nomdetector);
       RecoTrack.SetMass( _mass );
       RecoTrack.SetCharge( _charge );
       RecoTrack.SetMomentum( _momentum ); 
       
       // Add hits to track 
       for (int ihit=0; ihit<(int)HitStore.size(); ++ihit) {
         int daqid = HitStore[ihit].GetDAQID();
         int ipl = nomdetector.GetPlaneNumber(daqid); 
         RecoTrack.GetTE(ipl).SetHit(HitStore[ihit]); 
       }
       
       // Compute seed track 
       TBTrackState Seed = TrackSeeder.CreateSeedTrack(HitStore[0], HitStore[1], nomdetector);   
       RecoTrack.SetReferenceState(Seed);        

       // Try to run a Kalman Filter
       bool trkerr = TrackFitter.Fit(RecoTrack);
       if ( trkerr ) {
         streamlog_out ( MESSAGE3 ) << "Fit with " << RecoTrack.GetNumHits() << " hits failed. Skipping track!" << endl;
         continue;
       } 
            
       _histoMap["trkchi2_before"]->Fill( RecoTrack.GetChiSqu() ); 
       _histoMap["trkchi2ndof_before"]->Fill( RecoTrack.GetChiSqu()/RecoTrack.GetNDF()); 
          
       KalmanAlignmentInputProvider kaip;
       kaip.FillEvent(RecoTrack, *myEvent);
       
       // Fill track into alignment tree 
       myAlignTree->Fill();
       
     } // End track loop
     
     // Run alignment code
     
     streamlog_out ( MESSAGE3 ) << endl;
     streamlog_out ( MESSAGE3 ) << "Starting alignment..." << endl;
            
     KalmanAlignmentAlgorithm2 Aligner;
    
     // Try to fit alignment constants 
     
     TBDetector tmp_detector = nomdetector;
          
     AlignableDet reco = Aligner.Fit(tmp_detector, _rootFile, _alignConfigFileName );
     
     bool error = Aligner.AlignDetector(tmp_detector, reco);
     if ( error ) {
       streamlog_out ( MESSAGE3 ) << "Alignment failed!" << endl;
     } 
          
     // Compute total movement of sensors
     
     HepVector reco_state(6*nsensor);
     HepSymMatrix reco_cov = reco.alignmentCovariance;
     
     for ( int ipl = 0; ipl < nsensor; ipl++ ) {
            
       // Initial alignment 
       HepMatrix Rot_i = nomdetector.GetDet(ipl).GetNominal().GetRotation(); 
       HepVector pos_i = nomdetector.GetDet(ipl).GetNominal().GetPosition(); 
       
       // Final alignment  
       HepMatrix Rot_f = tmp_detector.GetDet(ipl).GetNominal().GetRotation();
       HepVector pos_f = tmp_detector.GetDet(ipl).GetNominal().GetPosition(); 
       
       // Diff rotation  
       HepMatrix dRot = Rot_f*Rot_i.T();
       
       // Diff shift   
       HepVector dr = pos_f - pos_i;
       
       reco_state[ipl*6 + 0] = dr[0]; 
       reco_state[ipl*6 + 1] = dr[1]; 
       reco_state[ipl*6 + 2] = dr[2]; 
       
       double dalpha, dbeta, dgamma; 
       GetAnglesKarimaki(dRot, dalpha, dbeta, dgamma);
       reco_state[ipl*6 + 3] = dalpha; 
       reco_state[ipl*6 + 4] = dbeta; 
       reco_state[ipl*6 + 5] = dgamma; 

     }
     
     // Histogram final geometry constants 
     
     for ( int ipl = 0; ipl < nsensor; ipl++ ) { 
        
       // This is the final position vector of the sensor
       HepVector pos_f = tmp_detector.GetDet(ipl).GetNominal().GetPosition(); 
    
       // This is the final rotation matrix of the sensor; it 
       // contains a discrete and a continuous factor. 
       HepMatrix Rot_f = tmp_detector.GetDet(ipl).GetNominal().GetRotation();
       
       // This is the discrete factor of sensor rotation. 
       HepMatrix DRot = tmp_detector.GetDet(ipl).GetDiscrete().GetRotation();
       
       // This is the continous factor of the rotation
       HepMatrix CRot_f = Rot_f*DRot.T(); 
       
       // Euler angles are defined wrt. the continous rotation 
       double alpha_f, beta_f, gamma_f; 
       GetAnglesKarimaki(CRot_f, alpha_f, beta_f, gamma_f); 

       _histoMap["hxshift"]->SetBinContent(ipl+1,pos_f[0]);
       _histoMap["hyshift"]->SetBinContent(ipl+1,pos_f[1]);
       _histoMap["hzshift"]->SetBinContent(ipl+1,pos_f[2]);
       _histoMap["hxrot"]->SetBinContent(ipl+1,alpha_f);
       _histoMap["hyrot"]->SetBinContent(ipl+1,beta_f);
       _histoMap["hzrot"]->SetBinContent(ipl+1,gamma_f);
     }
     
     
     // Evaluate alignment quality 
      
     streamlog_out ( MESSAGE1 ) << endl << "Final fit " << endl << endl; 
     
     for ( int ipl = 0; ipl < nsensor; ipl++ ) { 
       
       streamlog_out ( MESSAGE1 ) << endl << "Sensor plane " << ipl << endl << endl; 
       
       streamlog_out ( MESSAGE1 ) << "  reco dx " << reco_state[ipl*6+0] << " +/- " << std::sqrt(reco_cov[ipl*6+0][ipl*6+0]) << endl; 
       streamlog_out ( MESSAGE1 ) << "  true dx " << truemisalign[ipl*6+0] << endl;
       
       streamlog_out ( MESSAGE1 ) << endl; 
       streamlog_out ( MESSAGE1 ) << "  reco dy " << reco_state[ipl*6+1] << " +/- " << std::sqrt(reco_cov[ipl*6+1][ipl*6+1]) << endl;
       streamlog_out ( MESSAGE1 ) << "  true dy " << truemisalign[ipl*6+1] << endl; 
       
       streamlog_out ( MESSAGE1 ) << endl; 
       streamlog_out ( MESSAGE1 ) << "  reco dz " << reco_state[ipl*6+2] << " +/- " << std::sqrt(reco_cov[ipl*6+2][ipl*6+2]) << endl;
       streamlog_out ( MESSAGE1 ) << "  true dz " << truemisalign[ipl*6+2] << endl;
       
       streamlog_out ( MESSAGE1 ) << endl; 
       streamlog_out ( MESSAGE1 ) << "  reco dalpha " << reco_state[ipl*6+3] << " +/- " << std::sqrt(reco_cov[ipl*6+3][ipl*6+3]) << endl;
       streamlog_out ( MESSAGE1 ) << "  true dalpha " << truemisalign[ipl*6+3] << endl; 
       
       streamlog_out ( MESSAGE1 ) << endl;  
       streamlog_out ( MESSAGE1 ) << "  reco dbeta " << reco_state[ipl*6+4] << " +/- " << std::sqrt(reco_cov[ipl*6+4][ipl*6+4]) << endl;
       streamlog_out ( MESSAGE1 ) << "  true dbeta " << truemisalign[ipl*6+4] << endl;
       
       streamlog_out ( MESSAGE1 ) << endl;  
       streamlog_out ( MESSAGE1 ) << "  reco dgamma " << reco_state[ipl*6+5] << " +/- " << std::sqrt(reco_cov[ipl*6+5][ipl*6+5]) << endl;
       streamlog_out ( MESSAGE1 ) << "  true dgamma " << truemisalign[ipl*6+5] << endl; 
     }
     
     //Count wrong signs 
     int nwrong = 0; 
      
     for(int ipl = 0; ipl < nsensor; ipl++)  { 
       int offset = ipl*6; 
       if (truemisalign[offset+3]*reco_state[offset+3] <= 0) {
         if ( fabs(truemisalign[offset+3]) > 0.0001 ) ++nwrong;
       } 
       // Wrong sign for dbeta 
       if (truemisalign[offset+4]*reco_state[offset+4] <= 0) {
         if ( fabs(truemisalign[offset+4]) > 0.0001 ) ++nwrong;
       } 
     }
     
     streamlog_out ( MESSAGE1 ) <<" number of wrong signs " << nwrong << endl;  
     
     fill_align_pulls(reco_state, reco_cov, truemisalign);
      
     // Histogram data 
     
     _rootFile->cd("");
     
     string histoName = "covariance_mis"+to_string( iMisalign );
     string histoTitle ="Alignment Covariance "+to_string( iMisalign ); 
     int nBins = reco_cov.num_row(); 
     _histoMap2D[ histoName  ] = new TH2D(histoName.c_str(), histoTitle.c_str(), nBins, -0.5, nBins-0.5, nBins , -0.5, nBins-0.5); 
     _histoMap2D[ histoName  ]->SetXTitle( "" );
     _histoMap2D[ histoName  ]->SetYTitle( "");
     _histoMap2D[ histoName  ]->SetStats( false );
     _histoMap2D[ histoName  ]->GetYaxis()->SetTitleOffset(1.5);
     
     for (int ybin = 0; ybin< nBins; ybin++) {
       for (int xbin = 0; xbin < nBins; xbin++) {
         double corr = reco_cov[xbin][ybin]/(TMath::Sqrt(reco_cov[xbin][xbin])*TMath::Sqrt(reco_cov[ybin][ybin])  ) ;  
         _histoMap2D[ histoName ]->Fill(xbin,ybin, corr  );
       } 
     } 
     
     histoName = "alignerror"+to_string( iMisalign );
     histoTitle ="Align Error"+to_string( iMisalign ); 
     _histoMap[ histoName  ] = new TH1D(histoName.c_str(), histoTitle.c_str(), nBins, 0, nBins); 
     
     for (int xbin = 0; xbin < nsensor*6; xbin++) { 
       _histoMap[ histoName ]->Fill(xbin, reco_state[xbin] - truemisalign[xbin]  );
     }

     // Clean up  
     myAlignTree->Delete("all");
     delete myAlignTree;
     myHitTree->Delete("all");
     
     //////////////////////////////////////////////////////////////////////  
     // Create aligned detector 
     
     nomdetector = tmp_detector; 
     
     //////////////////////////////////////////////////////////////////////
     // Simulate tracks in misaligned detector   
     
     for(int iTrack = 0; iTrack < _nTracks; ++iTrack) { 
             
       
       //////////////////////////////////////////////////////////////////////
       // Simulate Track 
       
       // Store simulated hits along the particle 
       // trajectory for track fitting 
       vector<TBHit> HitStore;  

       // For electrons in PCMAG, Bethe Heitler model used to 
       // simulate energy loss by Bremsstrahlung
       double true_momentum = _momentum; 
     
       if ( _BetheHeitlerT0 > 0  ) {
         double t0 = _BetheHeitlerT0;
         double c = t0/TMath::Log(2);
         double u = double(RandGamma::shoot(c,1));
         double z = TMath::Exp(-u);    
         true_momentum *= z;  
       } 
       
       // Generate a beam particle, in misdetector(!!)
       TBTrack TruthTrack = TBGun.SimulateTrajectory(misdetector,_mass, _charge, true_momentum);
       
       if (  !TriggerSim.ACCEPT( TruthTrack )  ) {
         continue; 
       } 
       
       // Simulate hits along particle trajectory
       TBGun.SimulateTrackHits(TruthTrack, HitStore);
        
       // Start rekonstruction of particle 
       // trajectory from hits. 
       
       // Require at least 6 hits 
       //if ( HitStore.size() < 3 ) continue; 
       if ( HitStore.size() < 6 ) continue;
           
       // Init reco track, in nomdetector (now aligned!!)
       TBTrack RecoTrack(nomdetector);
       RecoTrack.SetMass( _mass );
       RecoTrack.SetCharge( _charge );
       RecoTrack.SetMomentum( _momentum ); 
       
       // Add hits to track 
       for (int ihit=0; ihit<(int)HitStore.size(); ++ihit) {
         int daqid = HitStore[ihit].GetDAQID();
         int ipl = nomdetector.GetPlaneNumber(daqid); 
         RecoTrack.GetTE(ipl).SetHit(HitStore[ihit]); 
       }
       
       // Compute seed track 
       TBTrackState Seed = TrackSeeder.CreateSeedTrack(HitStore[0], HitStore[1], nomdetector);   
       RecoTrack.SetReferenceState(Seed); 
            
       // Try to run a Kalman Filter
       bool trkerr = TrackFitter.Fit(RecoTrack);
       if ( trkerr ) {
         streamlog_out ( MESSAGE3 ) << "Fit with " << RecoTrack.GetNumHits() << " hits failed. Skipping track!" << endl;
         continue;
       } 
                
       _histoMap["trkchi2_after"]->Fill( RecoTrack.GetChiSqu() ); 
       _histoMap["trkchi2ndof_after"]->Fill( RecoTrack.GetChiSqu()/RecoTrack.GetNDF()); 
       _histoMap["chi2prob"]->Fill(TMath::Prob(RecoTrack.GetChiSqu(), RecoTrack.GetNDF()));
       _histoMap["mom_reco"]->Fill(RecoTrack.GetMomentum());
       
       fill_tracking_histos(RecoTrack); 
       
     } // End track loop

  
     streamlog_out ( MESSAGE3 ) << "Mean trk chi2/ndof " << _histoMap["trkchi2ndof_after"]->GetMean() << endl; 
    
   }  // End misalign loop  
     
   
}

void AlignmentSimulator::end () {
   streamlog_out(MESSAGE3) << std::endl
                           << "Processor succesfully finished!"
                           << std::endl;
    // ROOT Output
  _rootFile->Write();
  _rootFile->Close();
}

/**  Compute align pulls, difference to truth dividided by align errors
 *
 * This function uses the known truth, i.e. the true align parameters, and
 * compares the reconstructed ones to them, dividing by the align errors. This
 * "pull" distribution should be, in case of correct parameters and error
 * estimates, have mean zero and RMS one. With gaussian input errors
 * (i.e. gaussian hit errors), it should follow a gaussian
 * distribution. However, hits are distributed according to uniform
 * distribution, therefore the shape is not gaussian (but mean should still be
 * zero and RMS one).
 */
void AlignmentSimulator::fill_align_pulls(HepVector& reco_state, HepSymMatrix& reco_cov, 
                                         HepVector& truth_state)
{
  int nsensors = _detector.GetNSensors(); 
  std::string histoName;
  
  // This is the absolute alignment error 
  HepVector diff = reco_state-truth_state;
  
  for(int ipl = 0; ipl < nsensors; ipl++)  { 
     
    int offset = ipl*6; 
     
    // Fill correlation histos
    
    histoName = "hsim_dx_det"+to_string( ipl );
    _histoMap2D[ histoName  ]->Fill(truth_state[offset+0], reco_state[offset+0]);  
    
    histoName = "hsim_dy_det"+to_string( ipl );
    _histoMap2D[ histoName  ]->Fill(truth_state[offset+1], reco_state[offset+1]);  
             
    histoName = "hsim_dz_det"+to_string( ipl );
    _histoMap2D[ histoName  ]->Fill(truth_state[offset+2], reco_state[offset+2]); 
    
    histoName = "hsim_dalpha_det"+to_string( ipl );
    _histoMap2D[ histoName  ]->Fill( truth_state[offset+3], reco_state[offset+3]);  
    
    histoName = "hsim_dbeta_det"+to_string( ipl );
    _histoMap2D[ histoName  ]->Fill( truth_state[offset+4], reco_state[offset+4]);   
    
    histoName = "hsim_dgamma_det"+to_string( ipl );
    _histoMap2D[ histoName  ]->Fill( truth_state[offset+5], reco_state[offset+5]);  
      
    // Fill align parameter pulls
             
    histoName = "hpull_dx_det"+to_string( ipl );
    _histoMap[ histoName ]->Fill(diff[offset+0]/TMath::Sqrt(reco_cov[offset+0][offset+0]));
    
    histoName = "hpull_dy_det"+to_string( ipl );
    _histoMap[ histoName ]->Fill(diff[offset+1]/TMath::Sqrt(reco_cov[offset+1][offset+1]));
    
    histoName = "hpull_dz_det"+to_string( ipl );
    _histoMap[ histoName ]->Fill(diff[offset+2]/TMath::Sqrt(reco_cov[offset+2][offset+2]));
    
    histoName = "hpull_dalpha_det"+to_string( ipl );
    _histoMap[ histoName ]->Fill(diff[offset+3]/TMath::Sqrt(reco_cov[offset+3][offset+3])); 
    
    histoName = "hpull_dbeta_det"+to_string( ipl );
    _histoMap[ histoName ]->Fill(diff[offset+4]/TMath::Sqrt(reco_cov[offset+4][offset+4]));
    
    histoName = "hpull_dgamma_det"+to_string( ipl );
    _histoMap[ histoName ]->Fill(diff[offset+5]/TMath::Sqrt(reco_cov[offset+5][offset+5])); 
    
    // Fill align parameter errors
    
    histoName = "hsigma_dx_det"+to_string( ipl );
    _histoMap[ histoName  ]->Fill( TMath::Sqrt(reco_cov[offset+0][offset+0])) ; 
    
    histoName = "hsigma_dy_det"+to_string( ipl );
    _histoMap[ histoName  ]->Fill( TMath::Sqrt(reco_cov[offset+1][offset+1])) ;  
    
    histoName = "hsigma_dz_det"+to_string( ipl );
    _histoMap[ histoName  ]->Fill( TMath::Sqrt(reco_cov[offset+2][offset+2])) ; 
    
    histoName = "hsigma_dalpha_det"+to_string( ipl );
    _histoMap[ histoName  ]->Fill( TMath::Sqrt(reco_cov[offset+3][offset+3])) ;
    
    histoName = "hsigma_dbeta_det"+to_string( ipl );
    _histoMap[ histoName  ]->Fill( TMath::Sqrt(reco_cov[offset+4][offset+4])) ; 
    
    histoName = "hsigma_dgamma_det"+to_string( ipl );
    _histoMap[ histoName  ]->Fill( TMath::Sqrt(reco_cov[offset+5][offset+5])) ;
    
  }
}


/** Compute tracking parameter resolutions
 */   
void AlignmentSimulator::fill_tracking_histos(TBTrack & rectrk)
{
   
  
  int nsensors = _detector.GetNSensors(); 
  std::string histoName;
  
  double charge =  rectrk.GetCharge();
   
  for(int ipl = 0; ipl < nsensors; ipl++)  { 
    
    if ( rectrk.GetTE(ipl).IsCrossed() ) {
      
      // Kalman smoother results
      
      HepMatrix Kal_State = rectrk.GetTE(ipl).GetState().GetPars();
      HepSymMatrix Kal_Cov = rectrk.GetTE(ipl).GetState().GetCov();   

      double mom = std::abs(charge/Kal_State[4][0]); 
      
      // Local track parameters 
        
      histoName = "htrk_u_det"+to_string( ipl );
      _histoMap[ histoName ]->Fill(Kal_State[2][0]);  
      
      histoName = "htrk_v_det"+to_string( ipl );
      _histoMap[ histoName ]->Fill(Kal_State[3][0]);
      
      histoName = "htrk_tu_det"+to_string( ipl );
      _histoMap[ histoName ]->Fill(Kal_State[0][0]);
      
      histoName = "htrk_tv_det"+to_string( ipl );
      _histoMap[ histoName ]->Fill(Kal_State[1][0]); 

      histoName = "htrk_mom_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( mom ); 
      
      // Fill track parameter errors
      
      histoName = "hsigma_track_tu_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TMath::Sqrt(Kal_Cov[0][0]) ) ; 
      
      histoName = "hsigma_track_tv_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TMath::Sqrt(Kal_Cov[1][1]) ) ;  
      
      histoName = "hsigma_track_u_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TMath::Sqrt( Kal_Cov[2][2]) ); 
      
      histoName = "hsigma_track_v_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TMath::Sqrt( Kal_Cov[3][3]) );

      histoName = "hsigma_track_mom_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( mom*mom*TMath::Sqrt( Kal_Cov[4][4]) );  
      
      // Hit variables  
      
      if ( rectrk.GetTE(ipl).HasHit() ) {
       
        // Get hit residual
        HepMatrix r = rectrk.GetTE(ipl).GetHit().GetCoord();
        r -= rectrk.GetTE(ipl).GetState().GetXCoord();
        
        // Fill hit residuals 
        histoName = "hresidualU_det"+to_string( ipl );
        _histoMap[ histoName ]-> Fill(r[0][0]); 
        
        histoName = "hresidualV_det"+to_string( ipl );
        _histoMap[ histoName ]->Fill(r[1][0]);
        
      }
    }     
  }
    
}

void AlignmentSimulator::bookHistos() 
{
  
  streamlog_out ( MESSAGE4 ) <<  "Booking histograms" << endl;
  
  // ROOT Output 
  _rootFile = new TFile(_rootFileName.c_str(),"recreate"); 
    
  // Create alignment "event" which records all information from one track
  myEvent = new AlignEvent;
        
  int nsensors = _detector.GetNSensors(); 
   
  // Track chi2 histograms 
  _histoMap["trkchi2_before"] = new TH1D("htrkchi2_before", "", 100, 0, 10*_detector.GetNSensors() ); 
  _histoMap["trkchi2_before"]->SetYTitle("tracks"); 
  _histoMap["trkchi2_before"]->SetXTitle("track #chi^{2} before alignment");
  
  _histoMap["trkchi2ndof_before"] = new TH1D("htrkchi2ndof_before", "", 100, 0, 10); 
  _histoMap["trkchi2ndof_before"]->SetYTitle("tracks");
  _histoMap["trkchi2ndof_before"]->SetXTitle("track #chi^{2}/ndof before align"); 
   
  _histoMap["trkchi2_after"] = new TH1D("htrkchi2_after", "", 100, 0, 10*_detector.GetNSensors() );
  _histoMap["trkchi2_after"]->SetYTitle("tracks"); 
  _histoMap["trkchi2_after"]->SetXTitle("track #chi^{2} after align");
    
  _histoMap["trkchi2ndof_after"] = new TH1D("htrkchi2ndof_after", "", 100, 0, 10); 
  _histoMap["trkchi2ndof_after"]->SetYTitle("tracks");
  _histoMap["trkchi2ndof_after"]->SetXTitle("track #chi^{2}/ndof after align"); 
  
  _histoMap["chi2prob"] = new TH1D("hchi2prob", "", 100, 0, 1);
  _histoMap["chi2prob"]->SetXTitle("track p-value after align"); 
  _histoMap["chi2prob"]->SetYTitle("tracks");
  _histoMap["chi2prob"]->SetMinimum(0.);  	

  _histoMap["hxshift"] = new TH1D("hxshift", "", nsensors, 0, nsensors); 
  _histoMap["hxshift"]->SetXTitle("plane number");
  _histoMap["hxshift"]->SetYTitle("x shift [mm]");
  
  _histoMap["hyshift"] = new TH1D("hyshift", "", nsensors, 0, nsensors);
  _histoMap["hyshift"]->SetXTitle("plane number");
  _histoMap["hyshift"]->SetYTitle("y shift [mm]");
   
  _histoMap["hzshift"] = new TH1D("hzshift", "", nsensors, 0, nsensors); 
  _histoMap["hzshift"]->SetXTitle("plane number");
  _histoMap["hzshift"]->SetYTitle("z shift [mm]");
   
  _histoMap["hxrot"] = new TH1D("hxrot", "", nsensors, 0, nsensors); 
  _histoMap["hxrot"]->SetXTitle("plane number");
  _histoMap["hxrot"]->SetYTitle("x rotation [rad]");
  
  _histoMap["hyrot"] = new TH1D("hyrot", "", nsensors, 0, nsensors); 
  _histoMap["hyrot"]->SetXTitle("plane number");
  _histoMap["hyrot"]->SetYTitle("y rotation [rad]");
  
  _histoMap["hzrot"] = new TH1D("hzrot", "", nsensors, 0, nsensors);  
  _histoMap["hzrot"]->SetXTitle("plane number");
  _histoMap["hzrot"]->SetXTitle("z rotation [rad]");
  
  _histoMap["mom_reco"] = new TH1D("hmom_reco", "", 200, 0, 5*_momentum); 
  _histoMap["mom_reco"]->SetYTitle("tracks"); 
  _histoMap["mom_reco"]->SetXTitle("reco initial momentum [GeV]");
   
  // Create subdirs for detectors
  std::string dirName; 
  for (int ipl=0 ; ipl < _detector.GetNSensors(); ipl++) {
    std::string dirName = "Det"+to_string( ipl );
    _rootFile->mkdir(dirName.c_str());     
  }      
  
  // Detector histograms (pulls, residuals ...)
  for (int ipl=0 ; ipl < _detector.GetNSensors(); ipl++) {
    
    Det & adet = _detector.GetDet(ipl);
    
    dirName = "/Det"+to_string(ipl)+"/";
    _rootFile->cd(dirName.c_str());
    
    std::string histoName;
    double max; 
    double min; 
    
    // True offset parameter vs. fitted   
    histoName = "hsim_dx_det"+to_string( ipl );
    _histoMap2D[ histoName  ] = new TH2D(histoName.c_str(), "", 100, -5*_PositionAlignmentErrorX, 5*_PositionAlignmentErrorX,
                                                                              100, -5*_PositionAlignmentErrorX, 5*_PositionAlignmentErrorX); 
    _histoMap2D[ histoName  ]->SetXTitle( "True offset #Deltax, mm" );
    _histoMap2D[ histoName  ]->SetYTitle( "Fitted offset #Deltax, mm");
    _histoMap2D[ histoName  ]->SetStats( false );
    _histoMap2D[ histoName  ]->GetYaxis()->SetTitleOffset(1.5);
    
    histoName = "hsim_dy_det"+to_string( ipl );
    
    _histoMap2D[ histoName  ] = new TH2D(histoName.c_str(), "", 100, -5*_PositionAlignmentErrorY, 5*_PositionAlignmentErrorY,
                                                                              100, -5*_PositionAlignmentErrorY, 5*_PositionAlignmentErrorY); 
    _histoMap2D[ histoName  ]->SetXTitle( "True offset #Deltay, mm" );
    _histoMap2D[ histoName  ]->SetYTitle( "Fitted offset #Deltay, mm");
    _histoMap2D[ histoName  ]->SetStats( false );
    _histoMap2D[ histoName  ]->GetYaxis()->SetTitleOffset(1.5);
    
    histoName = "hsim_dz_det"+to_string( ipl );
    _histoMap2D[ histoName  ] = new TH2D(histoName.c_str(), "", 100, -5*_PositionAlignmentErrorZ, 5*_PositionAlignmentErrorZ,
                                                                              100, -5*_PositionAlignmentErrorZ, 5*_PositionAlignmentErrorZ); 
    _histoMap2D[ histoName  ]->SetXTitle( "True offset #Deltaz, mm" );
    _histoMap2D[ histoName  ]->SetYTitle( "Fitted offset #Deltaz, mm");
    _histoMap2D[ histoName  ]->SetStats( false );
    _histoMap2D[ histoName  ]->GetYaxis()->SetTitleOffset(1.5); 
    
    histoName = "hsim_dalpha_det"+to_string( ipl );
    _histoMap2D[ histoName  ] = new TH2D(histoName.c_str(), "", 100, -5*_RotationAlignmentErrorU, 5*_RotationAlignmentErrorU,
                                                                              100, -5*_RotationAlignmentErrorU, 5*_RotationAlignmentErrorU); 
    _histoMap2D[ histoName  ]->SetXTitle( "True offset #Delta#alpha, rad" );
    _histoMap2D[ histoName  ]->SetYTitle( "Fitted offset #Delta#alpha, rad");
    _histoMap2D[ histoName  ]->SetStats( false );
    _histoMap2D[ histoName  ]->GetYaxis()->SetTitleOffset(1.5);
     
    histoName = "hsim_dbeta_det"+to_string( ipl );
    _histoMap2D[ histoName  ] = new TH2D(histoName.c_str(), "", 100, -5*_RotationAlignmentErrorV, 5*_RotationAlignmentErrorV,
                                                                              100, -5*_RotationAlignmentErrorV, 5*_RotationAlignmentErrorV); 
    _histoMap2D[ histoName  ]->SetXTitle( "True offset #Delta#beta, rad" );
    _histoMap2D[ histoName  ]->SetYTitle( "Fitted offset #Delta#beta, rad");
    _histoMap2D[ histoName  ]->SetStats( false );
    _histoMap2D[ histoName  ]->GetYaxis()->SetTitleOffset(1.5);
         
    histoName = "hsim_dgamma_det"+to_string( ipl );
    _histoMap2D[ histoName  ] = new TH2D(histoName.c_str(), "", 100, -5*_RotationAlignmentErrorW, 5*_RotationAlignmentErrorW,
                                                                              100, -5*_RotationAlignmentErrorW, 5*_RotationAlignmentErrorW); 
    _histoMap2D[ histoName  ]->SetXTitle( "True offset #Delta#gamma, rad" );
    _histoMap2D[ histoName  ]->SetYTitle( "Fitted offset #Delta#gamma, rad");
    _histoMap2D[ histoName  ]->SetStats( false );
    _histoMap2D[ histoName  ]->GetYaxis()->SetTitleOffset(1.5);
    
    // Local track parameters 
       
    double  uBox = 1.1 * 0.5 * adet.GetModuleBoxSizeU();
    double  vBox = 1.1 * 0.5 * adet.GetModuleBoxSizeV();
    
    histoName = "htrk_u_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 2000, -uBox, +uBox); 
    _histoMap[ histoName ]->SetXTitle("intersect u [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks");      
    
    histoName = "htrk_v_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 2000, -vBox, +vBox);
    _histoMap[ histoName ]->SetXTitle("intersect v [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks");      
    
    histoName = "htrk_tu_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100000, -3, 3);
    _histoMap[ histoName ]->SetXTitle("slope du/dw [rad]");
    _histoMap[ histoName ]->SetYTitle("tracks");      
    
    histoName = "htrk_tv_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100000, -3, 3);
    _histoMap[ histoName ]->SetXTitle("slope dv/dw [rad]");
    _histoMap[ histoName ]->SetYTitle("tracks");  
    
    histoName = "htrk_mom_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, 0, 5*_momentum); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 
    _histoMap[ histoName ]->SetXTitle("reco momentum [GeV]"); 
    
    // Local align parameter errors 
     
    histoName = "hsigma_dx_det"+to_string( ipl );
    max = 5*_PositionAlignmentErrorX;
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 100000, 0, max); 
    _histoMap[ histoName  ]->SetXTitle("sigma #Deltax [mm]"); 
    _histoMap[ histoName ]->SetYTitle("misaligned setups");     

    histoName = "hsigma_dy_det"+to_string( ipl );
    max = 5*_PositionAlignmentErrorY;
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 100000, 0, max);
    _histoMap[ histoName  ]->SetXTitle("sigma #Deltay [mm]"); 
    _histoMap[ histoName ]->SetYTitle("misaligned setups");     

    histoName = "hsigma_dz_det"+to_string( ipl );
    max = 5*_PositionAlignmentErrorZ;
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 100000, 0, max);
    _histoMap[ histoName  ]->SetXTitle("sigma #Deltaz [mm]"); 
    _histoMap[ histoName ]->SetYTitle("misaligned setups");     

    histoName = "hsigma_dalpha_det"+to_string( ipl );
    max = 5*_RotationAlignmentErrorU;
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 100000, 0, max); 
    _histoMap[ histoName  ]->SetXTitle("sigma #Delta#alpha [rad]"); 
    _histoMap[ histoName ]->SetYTitle("misaligned setups");     

    histoName = "hsigma_dbeta_det"+to_string( ipl );
    max = 5*_RotationAlignmentErrorV;
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 100000, 0, max);
    _histoMap[ histoName  ]->SetXTitle("sigma #Delta#beta [rad]"); 
    _histoMap[ histoName ]->SetYTitle("misaligned setups");     

    histoName = "hsigma_dgamma_det"+to_string( ipl );
    max = 5*_RotationAlignmentErrorW;
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 100000, 0, max);
    _histoMap[ histoName  ]->SetXTitle("sigma #Delta#gamma [rad]");
    _histoMap[ histoName ]->SetYTitle("misaligned setups");  
        
    // Local align parameter pulls
    
    histoName = "hpull_dx_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, -10., 10.);
    _histoMap[ histoName  ]->SetXTitle("pull #Deltax");
    _histoMap[ histoName ]->SetYTitle("misaligned setups"); 
    
    histoName = "hpull_dy_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, -10., 10.);
    _histoMap[ histoName  ]->SetXTitle("pull #Deltay");
    _histoMap[ histoName ]->SetYTitle("misaligned setups");     

    histoName = "hpull_dz_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, -10., 10.); 
    _histoMap[ histoName  ]->SetXTitle("pull #Deltaz");
    _histoMap[ histoName ]->SetYTitle("misaligned setups");     

    histoName = "hpull_dalpha_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, -10., 10.);  
    _histoMap[ histoName  ]->SetXTitle("pull #Delta#alpha");
    _histoMap[ histoName ]->SetYTitle("misaligned setups");      

    histoName = "hpull_dbeta_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, -10., 10.);  
    _histoMap[ histoName  ]->SetXTitle("pull #Delta#beta");
    _histoMap[ histoName ]->SetYTitle("misaligned setups"); 
    
    histoName = "hpull_dgamma_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, -10., 10.);  
    _histoMap[ histoName  ]->SetXTitle("pull #Delta#gamma");
    _histoMap[ histoName ]->SetYTitle("misaligned setups"); 
    
    // Local track parameter errors - after alignment 
      
    histoName = "hsigma_track_u_det"+to_string( ipl );
    max = 100*adet.GetResolutionU();
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 8000, 0, max); 
    _histoMap[ histoName  ]->SetXTitle("error u [mm]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks");
    
    histoName = "hsigma_track_v_det"+to_string( ipl );
    max = 100*adet.GetResolutionV();
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 8000, 0, max);
    _histoMap[ histoName  ]->SetXTitle("error v [mm]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
    
    histoName = "hsigma_track_tu_det"+to_string( ipl );  
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 4000, 0, 0.01);
    _histoMap[ histoName  ]->SetXTitle("error du/dw [rad]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks");  
    
    histoName = "hsigma_track_tv_det"+to_string( ipl );
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 4000, 0, 0.01);  
    _histoMap[ histoName  ]->SetXTitle("error dv/dw [rad]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks");  
     
    histoName = "hsigma_track_mom_det"+to_string( ipl );
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 4000, 0, _momentum);
    _histoMap[ histoName  ]->SetXTitle("sigma p [GeV]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
    
    // Track residuals - after alignment
    
    histoName = "hresidualU_det"+to_string( ipl );
    min = -10*adet.GetResolutionU();
    max = +10*adet.GetResolutionU();
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, min, max); 
    _histoMap[ histoName ]->SetXTitle("u residual [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 
    
    histoName = "hresidualV_det"+to_string( ipl );
    min = -10*adet.GetResolutionV();
    max = +10*adet.GetResolutionV();
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, min, max); 
    _histoMap[ histoName ]->SetXTitle("v residual [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 
    
  }
  
}


} // Namespace


