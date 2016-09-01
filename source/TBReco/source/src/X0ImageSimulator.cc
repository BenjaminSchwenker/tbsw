// X0ImageSimulator
// Validation tool for material estimation
//
//                       
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// user includes
#include "X0ImageSimulator.h"

// DEPFETTrackTools includes 
#include "ParticleGun.h"
#include "GenericTrackFitter.h" 
#include "TBKalmanMSC.h" 
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

#include <CLHEP/Random/RandGamma.h>

using namespace std;
using namespace marlin;
using namespace CLHEP;

namespace depfet {

X0ImageSimulator::X0ImageSimulator ():DataSourceProcessor  ("X0ImageSimulator") {
  
  _description =
    "Validation tool for material estimation";

  registerProcessorParameter ("AlignmentDBFileName",
                             "This is the name of the LCIO file with the alignment constants (add .slcio)",
                             _alignmentDBFileName, static_cast< string > ( "alignmentDB.slcio" ) );  
   
  registerProcessorParameter ("NumberOfTracks",
                              "How many tracks for simulation",
                              _maxTracks,  static_cast < int > (250000)); 
    
  registerProcessorParameter ("ParticleMomentum", "Particle momentum [GeV]",
                              _momentum,  static_cast < double > (120.0));

  registerProcessorParameter ("ParticleMomentumError", "Random Error on Particle momentum [GeV]",
                              _momentum_error,  static_cast < double > (0.1));
   
  registerProcessorParameter ("ParticleMass", "Particle mass [GeV]",
                              _mass,  static_cast < double > (0.139));
   
  registerProcessorParameter ("ParticleCharge", "Particle charge [e]",
                              _charge,  static_cast < double > (+1));
  
  registerProcessorParameter ("MaterialBeforeTelescope", "Thickness of material before telescope [X/X0]",
                              _BetheHeitlerT0,  static_cast < double > (0.0)); 
  
  registerProcessorParameter ("UseTrueEnergy", "Reconstruct track using true true particle momentum",
                              _useTrueEnergy,  static_cast < bool > (false)); 
  
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
     
  registerProcessorParameter ("DetectionEfficiency", "Hit detection efficiency of telescope detectors [%]",
                              _hitEfficiency,  static_cast < double > (100)); 
  
  registerProcessorParameter ("MSCModel", "Choose model for multiple scattering:\n"
			      "0: Highland\n"
                              "1: Moliere",
                              _mscmodel,  static_cast < int > (0)); 

  
  
  registerProcessorParameter( "RootFileName", "Output root file name",
                              _rootFileName, std::string("X0Test.root"));

  registerProcessorParameter ("DUTPlane",
                              "Plane number of DUT along the beam line",
                              _idut,  static_cast < int > (3));

  
  
}


X0ImageSimulator * X0ImageSimulator::newProcessor () {  
  return new X0ImageSimulator;
}

void X0ImageSimulator::init () {
  
  printParameters ();
  
  // Rescale efficiency to intervall 0,..,1
  _hitEfficiency /= 100; 
  
  // Read detector constants from gear file
  _detector.ReadGearConfiguration();    
  
  // Read alignment data base file 
  _detector.ReadAlignmentDB( _alignmentDBFileName );
  
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


void X0ImageSimulator::readDataSource (int Ntrig) {


   
   //////////////////////////////////////////////////////////////////////
   // Particle gun settings   
       
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
      
   // Configure Kalman track fitter
   GenericTrackFitter TrackFitter(_detector);
   TBKalmanMSC TrackFitterMSC(_detector);
   TrackFitter.SetNumIterations(3);
   TrackFitter.SetOutlierCut(10000); 
   
   // Configure seed generator 
   SeedGenerator TrackSeeder(_charge, _momentum);
   
   
   int nsensor = _detector.GetNSensors(); 
   int nParameters = 6; 
            
   AlignableDet MisalignStoreMod(nsensor,nParameters);
     
   for ( int ipl = 0; ipl < nsensor; ipl++ ) {    
     // Shifts, mm
     MisalignStoreMod.alignmentParameters[ipl*nParameters + 0] = gRandom->Uniform(-0, 0);
     MisalignStoreMod.alignmentParameters[ipl*nParameters + 1] = gRandom->Uniform(-0, 0);
     MisalignStoreMod.alignmentParameters[ipl*nParameters + 2] = gRandom->Uniform(-0, 0);
     // Tilts, rad 
     MisalignStoreMod.alignmentParameters[ipl*nParameters + 3] = gRandom->Uniform(-0, 0);
     MisalignStoreMod.alignmentParameters[ipl*nParameters + 4] = gRandom->Uniform(-0, 0);
     MisalignStoreMod.alignmentParameters[ipl*nParameters + 5] = gRandom->Uniform(-0.01, 0.01);
   }
     
   // Use same mechnism for misalignment and alignment 
   KalmanAlignmentAlgorithm2 ModMisAligner; 
   bool miserr_mod = ModMisAligner.AlignDetector(_detector, MisalignStoreMod);
   if ( miserr_mod ) {
     cout << "BIG ERROR: Misalignmemnt of detector failed. Skipping!!" << endl;
    
   } 
   
   for(int eventNumber = 0; eventNumber < _maxTracks; ++eventNumber) { 
       
     // Print event number 
     if (eventNumber%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                       << eventNumber
                                                       << std::endl << std::endl;
     
     streamlog_out(MESSAGE1) << "Events processed: " << eventNumber << endl;
   
     
     // Store simulated hits along the particle 
     // trajectory for track fitting 
     vector<TBHit> HitStore;  
     
     // For electrons in PCMAG, Bethe Heitler model used to 
     // simulate energy loss by Bremsstrahlung
     double true_momentum = _momentum+_momentum_error;
     
     if ( _BetheHeitlerT0 > 0  ) {
       double t0 = _BetheHeitlerT0;
       double c = t0/TMath::Log(2);
       double u = double(RandGamma::shoot(c,1));
       double z = TMath::Exp(-u);    
       true_momentum *= z;  
     } 

     // Generate a beam particle 
     TBTrack TruthTrack = TBGun.SimulateTrajectory(_detector,_mass, _charge, true_momentum);

     // Simulate hits along particle trajectory
     TBGun.SimulateTrackHits(TruthTrack, HitStore);
     
     // Start rekonstruction of particle 
     // trajectory from hits. 
     
     int nsensor = _detector.GetNSensors(); 

     // Require at least 3 hits 
     if ( HitStore.size() < 6 ) continue; 
           
     // Init reco track  
     TBTrack RecoTrack(_detector);
     RecoTrack.SetMass( _mass );
     RecoTrack.SetCharge( _charge );
     RecoTrack.SetMomentum( _momentum );
     
     if (_useTrueEnergy) RecoTrack.SetMomentum( true_momentum );  
      
     // Add hits to track 
     for (int ihit=0; ihit<(int)HitStore.size(); ++ihit) {
       int daqid = HitStore[ihit].GetDAQID();
       int ipl = _detector.GetPlaneNumber(daqid);
       RecoTrack.GetTE(ipl).SetHit(HitStore[ihit]); 
     }
     
     // Compute seed track 
     TBTrackState Seed = TrackSeeder.CreateSeedTrack(HitStore[0], HitStore[1], _detector);   
     RecoTrack.SetReferenceState(Seed);

     // Try to run a Kalman Filter on RecoTrack
     bool trkerr = TrackFitter.Fit(RecoTrack);
     if ( trkerr ) {
       streamlog_out ( MESSAGE3 ) << "Fit with " << RecoTrack.GetNumHits() << " hits failed. Skipping track!" << endl;
       continue;
     }  

     //Copy the RecoTrack two times for the forward backward Kalman Filter pair
     TBTrack ForwardTrack(RecoTrack);
     TBTrack BackwardTrack(RecoTrack);
     
     //MSC Analysis

     //Parameters needed in the TrackfitterMSC

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
     
     double comboChi2 = ForwardTrack.GetChiSqu()+BackwardTrack.GetChiSqu();

     // definition of DUT
     Det dut = _detector.GetDet(_idut);
 
     // In and OutStates of the reconstructed Track at the current detector
     TBTrackState& InState=ForwardTrack.GetTE(_idut).GetState();
     TBTrackState& OutState=BackwardTrack.GetTE(_idut).GetState(); 
     
     // Calculate the Kinks and Covariances from the two trackstates and the sensor position/orientation
     HepMatrix theta(2,1,0);
     HepSymMatrix Cov (2,0);
     theta = TrackFitterMSC.GetScatterKinks(dut, InState, OutState); 
     Cov = TrackFitterMSC.GetScatterKinkCov(dut, InState, OutState);
     
     // Get the track parameters of the fitted track on the current sensor
     // The u and v positions are needed for a position-resolved measurement
     HepMatrix p = RecoTrack.GetTE(_idut).GetState().GetPars();
     HepMatrix p_in = InState.GetPars();
     HepMatrix p_out = OutState.GetPars();
     
     // Fill root variables 
     
     _rootRunNumber = 0;  
     _rootEventNumber = eventNumber;  
     
     _rootDaqID = dut.GetDAQID(); 
     _rootPlaneID = _idut;
     _root_momentum = RecoTrack.GetMomentum(); 
     _rootTrackHits = RecoTrack.GetNumHits();
     _rootTrackChi2 = RecoTrack.GetChiSqu(); 
     _rootTrackProb = TMath::Prob(RecoTrack.GetChiSqu(),RecoTrack.GetNDF());
     _rootTrackProbUp = TMath::Prob(ForwardTrack.GetChiSqu(),ForwardTrack.GetNDF());
     _rootTrackProbDown = TMath::Prob(BackwardTrack.GetChiSqu(),BackwardTrack.GetNDF());
     _rootTrackProbCombo = TMath::Prob( comboChi2 ,ForwardTrack.GetNDF()+BackwardTrack.GetNDF());
     
     _root_x = RecoTrack.GetTE(0).GetState().GetPars()[2][0];   
     _root_y = RecoTrack.GetTE(0).GetState().GetPars()[3][0];   
     _root_dxdz = RecoTrack.GetTE(0).GetState().GetPars()[0][0];  
     _root_dydz = RecoTrack.GetTE(0).GetState().GetPars()[1][0];
     _root_u = p[2][0]; 
     _root_v = p[3][0]; 
     _root_dudw = p[0][0]; 
     _root_dvdw = p[1][0]; 
     
     _root_u_in = p_in[2][0]; 
     _root_v_in = p_in[3][0];
     _root_u_out = p_out[2][0]; 
     _root_v_out = p_out[3][0];
     
     _root_angle1 = theta[0][0];
     _root_angle2 = theta[1][0];
     _root_angle1_err = TMath::Sqrt(Cov[0][0]);
     _root_angle2_err = TMath::Sqrt(Cov[1][1]);
     		     
     if(TruthTrack.GetTE(_idut).HasHit()) {			    
       TBHit MSC = TruthTrack.GetTE(_idut).GetHit();
       _root_angle1_truth = MSC.GetCoord()[0][0];
       _root_angle2_truth = MSC.GetCoord()[1][0];
     } else {
       _root_angle1_truth = -99;
       _root_angle2_truth = -99;
     }
        
     _root_truth_momentum = TruthTrack.GetMomentum(); 

     _rootMscTree->Fill(); 
       				     
     fill_tracking_pulls(RecoTrack, TruthTrack);
     

  } // End loop over events
   
  
  
}

void X0ImageSimulator::end () {
   streamlog_out(MESSAGE3) << std::endl
                           << "Processor succesfully finished!"
                           << std::endl;
    // ROOT Output
  _rootFile->Write();
  _rootFile->Close();
  delete _rootFile; 
}


/**  Compute tracking pulls, difference to truth dividided by tracking errors
 *
 * This function uses the known truth, i.e. the true track parameters, and
 * compares the reconstructed ones to them, dividing by the track errors. This
 * "pull" distribution should be, in case of correct parameters and error
 * estimates, have mean zero and RMS one. With gaussian input errors
 * (i.e. gaussian hit errors), it should follow a gaussian
 * distribution. However, hits are distributed according to uniform
 * distribution, therefore the shape is not gaussian (but mean should still be
 * zero and RMS one).
 */
void X0ImageSimulator::fill_tracking_pulls(TBTrack & reco, TBTrack & truth)
{
  int nsensors = _detector.GetNSensors(); 
  std::string histoName;

  for(int ipl = 0; ipl < nsensors; ipl++)  { 
    
    if ( reco.GetTE(ipl).IsCrossed() ) {
      
      HepMatrix Sim_State = truth.GetTE(ipl).GetState().GetPars();
         
      // Fill truth track parameters
      
      histoName = "htrk_tu_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(Sim_State[0][0]);  
    
      histoName = "htrk_tv_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(Sim_State[1][0]);  
             
      histoName = "htrk_u_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(Sim_State[2][0]); 
    
      histoName = "htrk_v_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(Sim_State[3][0]);  
      
      // Track fitter results 
      HepMatrix Kal_State = reco.GetTE(ipl).GetState().GetPars();
      HepSymMatrix Kal_Cov = reco.GetTE(ipl).GetState().GetCov();   
      
      // Fill track parameter pulls
      
      HepMatrix diff = Kal_State-Sim_State;
           
      histoName = "hp1_det"+to_string( ipl );
      _histoMap[ histoName ]->Fill(diff[0][0]/TMath::Sqrt(Kal_Cov[0][0]));
      
      histoName = "hp2_det"+to_string( ipl );
      _histoMap[ histoName ]->Fill(diff[1][0]/TMath::Sqrt(Kal_Cov[1][1]));
      
      histoName = "hp3_det"+to_string( ipl );
      _histoMap[ histoName ]->Fill(diff[2][0]/TMath::Sqrt(Kal_Cov[2][2]));
      
      histoName = "hp4_det"+to_string( ipl );
      _histoMap[ histoName ]->Fill(diff[3][0]/TMath::Sqrt(Kal_Cov[3][3])); 
      
      // Fill track parameter errors
      
      histoName = "herr_tu_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TMath::Sqrt(Kal_Cov[0][0]) ) ; 
      
      histoName = "herr_tv_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TMath::Sqrt(Kal_Cov[1][1]) ) ;  
      
      histoName = "herr_u_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TMath::Sqrt( Kal_Cov[2][2]) ); 
      
      histoName = "herr_v_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TMath::Sqrt( Kal_Cov[3][3]) );
      
      // Hit variables  
      
      if ( reco.GetTE(ipl).HasHit() ) {
        
        // Get local chi2   	     
        double  hitchi2 = reco.GetTE(ipl).GetChiSqu();  
        
        // Get hit residual
        HepMatrix r = reco.GetTE(ipl).GetHit().GetCoord();
        r -= reco.GetTE(ipl).GetState().GetXCoord();
        
        // Fill hit residuals 
        histoName = "hresU_det"+to_string( ipl );
        _histoMap[ histoName ]-> Fill(r[0][0]); 
        
        histoName = "hresV_det"+to_string( ipl );
        _histoMap[ histoName ]->Fill(r[1][0]);
        
        histoName = "hchi2inc_det"+to_string( ipl );  
        _histoMap[histoName]->Fill(hitchi2);  
        
        histoName = "hchi2incprob_det"+to_string( ipl );
        _histoMap[histoName]->Fill(TMath::Prob(hitchi2, 2));
        
      }
    }
       
  }
}


void X0ImageSimulator::bookHistos() 
{
  
  streamlog_out ( MESSAGE4 ) <<  "Booking histograms" << endl;
  
  // ROOT Output 
  _rootFile = new TFile(_rootFileName.c_str(),"recreate"); 

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
  _rootMscTree->Branch("u_in"         ,&_root_u_in          ,"u_in/D");
  _rootMscTree->Branch("v_in"         ,&_root_v_in          ,"v_in/D");
  _rootMscTree->Branch("u_out"        ,&_root_u_out         ,"u_out/D");
  _rootMscTree->Branch("v_out"        ,&_root_v_out         ,"v_out/D");
  _rootMscTree->Branch("dudw"         ,&_root_dudw          ,"dudw/D");
  _rootMscTree->Branch("dvdw"         ,&_root_dvdw          ,"dvdw/D"); 
  _rootMscTree->Branch("theta1_val"      ,&_root_angle1      ,"theta1_val/D"); 
  _rootMscTree->Branch("theta2_val"      ,&_root_angle2      ,"theta2_val/D");
  _rootMscTree->Branch("theta1_err"      ,&_root_angle1_err  ,"theta1_err/D");
  _rootMscTree->Branch("theta2_err"      ,&_root_angle2_err  ,"theta2_err/D");
  _rootMscTree->Branch("momentum"        ,&_root_momentum    ,"momentum/D");
  _rootMscTree->Branch("truth_momentum"  ,&_root_truth_momentum    ,"truth_momentum/D");
  _rootMscTree->Branch("theta1_true", &_root_angle1_truth);
  _rootMscTree->Branch("theta2_true", &_root_angle2_truth);
  
     
  // Create subdirs for detectors 
  for (int ipl=0 ; ipl < _detector.GetNSensors(); ipl++) {
    std::string dirName = "Det"+to_string( ipl );
    _rootFile->mkdir(dirName.c_str());     
  }      
  
  // Detector histograms (pulls, residuals ...)
  for (int ipl=0 ; ipl < _detector.GetNSensors(); ipl++) {
                    
    Det & adet = _detector.GetDet(ipl);
    
    std::string dirName = "/Det"+to_string(ipl)+"/";
    _rootFile->cd(dirName.c_str());
    
    std::string histoName;
    std::string histoTitle;
    double min, max; 
    
    // Local track parameters 
    
    double  uBox = 1.1 * 0.5 * adet.GetModuleBoxSizeU();
    double  vBox = 1.1 * 0.5 * adet.GetModuleBoxSizeV();
    
    histoName = "htrk_u_det"+to_string( ipl );
    histoTitle ="Track intersection u"; 
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), histoTitle.c_str(), 10000, -uBox, +uBox); 
    _histoMap[ histoName  ]->SetXTitle("intersect u [mm]");  
    
    histoName = "htrk_v_det"+to_string( ipl );
    histoTitle ="Track intersection v"; 
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), histoTitle.c_str(), 10000, -vBox, +vBox);
    _histoMap[ histoName  ]->SetXTitle("intersect v [mm]");  
   
    histoName = "htrk_tu_det"+to_string( ipl );
    histoTitle ="Track slope du/dw";   
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), histoTitle.c_str(), 100000, -3, 3); 
    _histoMap[ histoName  ]->SetXTitle("slope du/dw [rad]"); 
    
    histoName = "htrk_tv_det"+to_string( ipl );
    histoTitle ="Track slope dv/dw"; 
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), histoTitle.c_str(), 100000, -3, 3); 
    _histoMap[ histoName  ]->SetXTitle("slope dv/dw [rad]"); 
       
    // Local track parameter errors 
     
    histoName = "herr_u_det"+to_string( ipl );
    histoTitle ="RMS error intersect u"; 
    max = 100*adet.GetResolutionU();
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), histoTitle.c_str(), 8000, 0, max); 
    _histoMap[ histoName  ]->SetXTitle("error u [mm]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
    
    histoName = "herr_v_det"+to_string( ipl );
    histoTitle ="RMS error intersect v"; 
    max = 100*adet.GetResolutionV();
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), histoTitle.c_str(), 8000, 0, max);
    _histoMap[ histoName  ]->SetXTitle("error v [mm]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
    
    histoName = "herr_tu_det"+to_string( ipl );
    histoTitle ="RMS error slope du/dw";   
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), histoTitle.c_str(), 4000, 0, 0.01); 
    _histoMap[ histoName  ]->SetXTitle("error du/dw [rad]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
    
    histoName = "herr_tv_det"+to_string( ipl );
    histoTitle ="RMS error slope dv/dw"; 
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), histoTitle.c_str(), 4000, 0, 0.01);
    _histoMap[ histoName  ]->SetXTitle("error dv/dw [rad]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks");  
       
    // Local track parameter pulls
    
    histoName = "hp1_det"+to_string( ipl );
    histoTitle ="Pull slope du/dw"; 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), histoTitle.c_str(), 100, -4., 4.);
    
    histoName = "hp2_det"+to_string( ipl );
    histoTitle ="Pull slope dv/dw"; 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), histoTitle.c_str(), 100, -4., 4.);
    
    histoName = "hp3_det"+to_string( ipl );
    histoTitle ="Pull intersect u"; 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), histoTitle.c_str(), 100, -4., 4.); 
    
    histoName = "hp4_det"+to_string( ipl );
    histoTitle ="Pull intersect v"; 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), histoTitle.c_str(), 100, -4., 4.);  
    
    histoName = "hresU_det"+to_string( ipl );
    histoTitle ="U Residuals"; 
    min = -10*adet.GetResolutionU();
    max = +10*adet.GetResolutionU();
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), histoTitle.c_str(), 100, min, max);
    _histoMap[ histoName ]->SetXTitle("u residual [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 
    
    histoName = "hresV_det"+to_string( ipl );
    histoTitle ="V Residuals"; 
    min = -10*adet.GetResolutionV();
    max = +10*adet.GetResolutionV();
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), histoTitle.c_str(), 100, min, max); 
    _histoMap[ histoName ]->SetXTitle("v residual [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 
     
    histoName = "hchi2inc_det"+to_string( ipl );    
    _histoMap[histoName] = new TH1D(histoName.c_str(), "#chi^{2} increment", 100, 0, 20 );    
    _histoMap[histoName]->SetXTitle("#chi^{2} increment"); 
    _histoMap[histoName]->SetYTitle("tracks");  
    
    histoName = "hchi2incprob_det"+to_string( ipl );
    _histoMap[histoName] = new TH1D(histoName.c_str(), "#chi^{2} increment probability", 100, 0, 1 ); 
    _histoMap[histoName]->SetMinimum(0.);
    _histoMap[histoName]->SetXTitle("p-value"); 
    _histoMap[histoName]->SetYTitle("tracks");    
  }
  
} 

} // Namespace


