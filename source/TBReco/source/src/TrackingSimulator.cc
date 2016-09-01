// TrackingSimulator
// Validation tool for track fitting routines 
//                       
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// user includes
#include "TrackingSimulator.h"
#include "TBEvtGen.h"

// DEPFETTrackTools includes 
#include "ParticleGun.h"
#include "GenericTrackFitter.h"
#include "TBHit.h"
#include "SeedGenerator.h"
#include "Utilities.h"

// marlin includes
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"

#include <CLHEP/Random/RandGamma.h>

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

TrackingSimulator::TrackingSimulator ():DataSourceProcessor  ("TrackingSimulator") {
  
  _description =
    "Validation tool for telescope track fitting algorithms";
   
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
  
  registerProcessorParameter ("GunIntensity", "Number of particles per second",
                              _beamIntensity,  static_cast < double > (10000)); 
  
  registerProcessorParameter ("TriggerMode", "0:front sct, 1:back sct, 2: all sct",
                              _trgmode,  static_cast < int > (0)); 
   
  registerProcessorParameter ("MaterialBeforeTelescope", "Thickness of material before telescope [X/X0]",
                              _BetheHeitlerT0,  static_cast < double > (0.0)); 

  registerProcessorParameter ("UseTrueEnergy", "Reconstruct track using true true particle momentum",
                              _useTrueEnergy,  static_cast < bool > (false)); 

  registerProcessorParameter ("NumberOfFitIterations", "Iterate the track fit several times",
                              _numit,  static_cast < int > (1)); 
  
  registerProcessorParameter ("MinHits", "Minimum number of hits in telescope for fitting",
                              _minHits,  static_cast < int > (6)); 
   
  registerProcessorParameter ("HitDigitizer", 
                              "Choose digitizer model for pixel detectors:\n"
                              "0: Digital Mimosa26 digitizer\n"
                              "1: Gaussian smearing digitizer\n"
                              "2: Box digitizer",
                              _digi_type,  static_cast < int > (1)); 
     
  registerProcessorParameter ("DetectionEfficiency", "Hit detection efficiency of pixel detectors [%]",
                              _hitEfficiency,  static_cast < double > (100)); 
  
  registerProcessorParameter ("MultipleScattering", "Choose model of multiple scattering in the telescope:\n"
                              "0: Highland\n"
                              "1: Moliere",
                              _mscmodel,  static_cast < int > (0));  
  
  registerProcessorParameter( "RootFileName", "Output root file name",
                              _rootFileName, std::string("TrackingTest.root"));
  
}


TrackingSimulator * TrackingSimulator::newProcessor () {  
  return new TrackingSimulator;
}


void TrackingSimulator::init () {
  
  printParameters ();
  
  // Rescale efficiency to intervall 0,..,1
  _hitEfficiency /= 100; 
  
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


void TrackingSimulator::readDataSource (int Ntrig) {
       
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
   TriggerSim.intensity = _beamIntensity; 
   TriggerSim.trgmode = _trgmode;  
 
   // Configure Kalman track fitter
   GenericTrackFitter  TrackFitter(_detector); 
   TrackFitter.SetNumIterations(_numit);
   TrackFitter.SetOutlierCut(10000); 
      
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
     double true_momentum = _momentum; 
     
     if ( _BetheHeitlerT0 > 0  ) {
       double t0 = _BetheHeitlerT0;
       double c = t0/TMath::Log(2);
       double u = double(RandGamma::shoot(c,1));
       double z = TMath::Exp(-u);    
       true_momentum *= z;  
     } 
     
     // Generate a beam particle 
     TBTrack TruthTrack = TBGun.SimulateTrajectory(_detector,_mass, _charge, true_momentum);
      
     if (  !TriggerSim.ACCEPT( TruthTrack )  ) {
       continue; 
     } 
     
     // Simulate hits along particle trajectory
     TBGun.SimulateTrackHits(TruthTrack, HitStore);
     
     // Start rekonstruction of particle 
     // trajectory from hits. 
     
     // Require at least 2 simulated hits for seeding 
     if ( HitStore.size() < 2 ) continue; 
       
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
     
     // Configure seed generator 
     SeedGenerator TrackSeeder(_charge, _momentum);

     // Compute seed track 
     TBTrackState Seed = TrackSeeder.CreateSeedTrack(HitStore[0], HitStore[1], _detector);

     // Use truth state at plane 0 as seed
     //TBTrackState StateAtZ0;
     //StateAtZ0.Pars = TruthTrack.GetTE(0).GetState().GetPars();
     //StateAtZ0.Plane = TruthTrack.GetTE(0).GetDet().GetNominal();
     
     RecoTrack.SetReferenceState(Seed);
     
     // Try to run track fitter
     bool trkerr = TrackFitter.Fit(RecoTrack);
     if ( trkerr ) {
       streamlog_out ( MESSAGE3 ) << "Fit with " << RecoTrack.GetNumHits() << " hits failed. Skipping track!" << endl;
       continue;
     }   
     
     // Require at least 6 hits 
     if ( RecoTrack.GetNumHits() < _minHits ) continue; 
     
     fill_tracking_pulls(RecoTrack, TruthTrack);
     fill_track_chi2(RecoTrack);
 
     
  } // End loop over events  
   
  
   
  
   
}

void TrackingSimulator::end () {
   streamlog_out(MESSAGE3) << std::endl
                           << "Processor succesfully finished!"
                           << std::endl;
    // ROOT Output
  _rootFile->Write();
  _rootFile->Close();
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
void TrackingSimulator::fill_tracking_pulls(TBTrack & reco, TBTrack & truth)
{
  
  _histoMap["mom_truth"]->Fill(truth.GetMomentum());    
  _histoMap["mom_reco"]->Fill(reco.GetMomentum());
  
  int nsensors = _detector.GetNSensors(); 
  std::string histoName;
  
  double charge =  reco.GetCharge();

  for(int ipl = 0; ipl < nsensors; ipl++)  { 
    
    if ( truth.GetTE(ipl).IsCrossed() ) {
      
      // Truth track state 
      HepMatrix Sim_State = truth.GetTE(ipl).GetState().GetPars();  
      
      histoName = "htrk_dir_truth_det"+to_string( ipl ); 
      _histoMap2D[ histoName  ]->Fill(Sim_State[0][0],Sim_State[1][0]);  
      
      histoName = "htrk_tu_truth_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(Sim_State[0][0]);  
      
      histoName = "htrk_tv_truth_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(Sim_State[1][0]);  
      
      histoName = "htrk_u_truth_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(Sim_State[2][0]);  
      
      histoName = "htrk_v_truth_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(Sim_State[3][0]);  

      histoName = "htrk_mom_truth_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(std::abs(charge/Sim_State[4][0]));  
    }

    if ( reco.GetTE(ipl).IsCrossed() && truth.GetTE(ipl).IsCrossed() ) {
      
      // Truth track state 
      HepMatrix Sim_State = truth.GetTE(ipl).GetState().GetPars();
      
      // Track fitter results 
      HepMatrix Kal_State = reco.GetTE(ipl).GetState().GetPars();
      HepSymMatrix Kal_Cov = reco.GetTE(ipl).GetState().GetCov();   

      double mom = std::abs(charge/Kal_State[4][0]); 

      // Momentum calibration
      if (ipl == 0) _histoMap2D[ "eCal1"  ]->Fill( truth.GetMomentum(), Kal_State[0][0] );
      if (ipl == 0) _histoMap2D[ "eCal2"  ]->Fill( truth.GetMomentum(), Kal_State[1][0] );         

      // Fill reco track parameters
             
      histoName = "htrk_tu_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(Kal_State[0][0]);  
      
      histoName = "htrk_tv_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(Kal_State[1][0]); 
 
      histoName = "htrk_u_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(Kal_State[2][0]); 
      
      histoName = "htrk_v_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(Kal_State[3][0]); 
      
      histoName = "htrk_mom_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( mom ); 
 
      // Fill track parameter pulls
      
      HepMatrix diff = Kal_State-Sim_State;
       
      int ierr; 
      HepMatrix jchisq = diff.T()*Kal_Cov.inverse(ierr)*diff;
         
      histoName = "hJ_det"+to_string( ipl );
      _histoMap[histoName]->Fill(jchisq[0][0]);  
        
      histoName = "hJp_det"+to_string( ipl );
      _histoMap[histoName]->Fill(TMath::Prob(jchisq[0][0], 5));
           
      histoName = "hp1_det"+to_string( ipl );
      _histoMap[ histoName ]->Fill(diff[0][0]/TMath::Sqrt(Kal_Cov[0][0]));
      
      histoName = "hp2_det"+to_string( ipl );
      _histoMap[ histoName ]->Fill(diff[1][0]/TMath::Sqrt(Kal_Cov[1][1]));
      
      histoName = "hp3_det"+to_string( ipl );
      _histoMap[ histoName ]->Fill(diff[2][0]/TMath::Sqrt(Kal_Cov[2][2]));
      
      histoName = "hp4_det"+to_string( ipl );
      _histoMap[ histoName ]->Fill(diff[3][0]/TMath::Sqrt(Kal_Cov[3][3])); 

      histoName = "hp5_det"+to_string( ipl );
      _histoMap[ histoName ]->Fill(diff[4][0]/TMath::Sqrt(Kal_Cov[4][4])); 
      
      // Fill track parameter errors
      
      histoName = "hsigma_tu_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TMath::Sqrt(Kal_Cov[0][0]) ) ; 
      
      histoName = "hsigma_tv_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TMath::Sqrt(Kal_Cov[1][1]) ) ;  
      
      histoName = "hsigma_u_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TMath::Sqrt( Kal_Cov[2][2]) ); 
      
      histoName = "hsigma_v_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TMath::Sqrt( Kal_Cov[3][3]) );
      
      histoName = "hsigma_mom_det"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( mom*mom*TMath::Sqrt( Kal_Cov[4][4]) ); 
      
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

/**  Make a histogram of track chi2, chi2/ndof and chi2-Probability. 
 */
void TrackingSimulator::fill_track_chi2(TBTrack& trk)
{ 
  _histoMap["nhits"]->Fill(trk.GetNumHits());    
  _histoMap["ndf"]->Fill(trk.GetNDF());   
  _histoMap["chi2"]->Fill(trk.GetChiSqu());
  _histoMap["chi2ndof"]->Fill(trk.GetChiSqu()/trk.GetNDF());
  _histoMap["chi2prob"]->Fill(TMath::Prob(trk.GetChiSqu(), trk.GetNDF()));
}

void TrackingSimulator::bookHistos() 
{
  
  streamlog_out ( MESSAGE4 ) <<  "Booking histograms" << endl;
  
  // ROOT Output 
  _rootFile = new TFile(_rootFileName.c_str(),"recreate"); 
      
  // Track chi2 histograms 
  _histoMap["chi2"] = new TH1D("hchi2", "", 30*_detector.GetNSensors(), 0, 10*_detector.GetNSensors() ); 
  _histoMap["chi2"]->SetYTitle("tracks"); 
  _histoMap["chi2"]->SetXTitle("#chi^{2}");
  
  _histoMap["nhits"] = new TH1D("hnhits", "", _detector.GetNSensors()+2, 0, _detector.GetNSensors()+2); 
  _histoMap["nhits"]->SetYTitle("tracks"); 
  _histoMap["nhits"]->SetXTitle("hits"); 

  _histoMap["ndf"] = new TH1D("hndf", "", 2*_detector.GetNSensors(), 0, 2*_detector.GetNSensors()); 
  _histoMap["ndf"]->SetYTitle("tracks"); 
  _histoMap["ndf"]->SetXTitle("degrees of freedom"); 
  
  _histoMap["mom_truth"] = new TH1D("hmom_truth", "", 200, 0, 1.2*_momentum); 
  _histoMap["mom_truth"]->SetYTitle("tracks"); 
  _histoMap["mom_truth"]->SetXTitle("true initial momentum [GeV]"); 
  
  _histoMap["mom_reco"] = new TH1D("hmom_reco", "", 200, 0, 1.2*_momentum); 
  _histoMap["mom_reco"]->SetYTitle("tracks"); 
  _histoMap["mom_reco"]->SetXTitle("reco initial momentum [GeV]"); 
  
  _histoMap["chi2ndof"] = new TH1D("hchi2ndof", "", 100, 0, 10); 
  _histoMap["chi2ndof"]->SetYTitle("tracks"); 
  _histoMap["chi2ndof"]->SetXTitle("#chi^{2}/ndof");  
   
  _histoMap["chi2prob"] = new TH1D("hchi2prob", "", 100, 0, 1);
  _histoMap["chi2prob"]->SetXTitle("p-value"); 
  _histoMap["chi2prob"]->SetYTitle("tracks");
  _histoMap["chi2prob"]->SetMinimum(0.);  
  
  // Average beam directions  
  int nEBins = 60;
  double eMin = 0; 
  double eMax = 6; 
  int NSBins = 400;   
  double sMin = -0.04; 
  double sMax = +0.04;   
  
  _histoMap2D[ "eCal1"  ]  = new TH2D("eCal1","calibration1", nEBins,eMin,eMax,NSBins,sMin,sMax);
  _histoMap2D[ "eCal1"  ]->SetXTitle("true momentum at first tel.plane [GeV]"); 
  _histoMap2D[ "eCal1"  ]->SetYTitle("reco du/dw at first tel. plane [rad]");
  _histoMap2D[ "eCal1"  ]->SetStats( false );

  _histoMap2D[ "eCal2"  ]  = new TH2D("eCal2","calibration2", nEBins,eMin,eMax,NSBins,sMin,sMax);
  _histoMap2D[ "eCal2"  ]->SetXTitle("true momentum at first tel.plane [GeV]"); 
  _histoMap2D[ "eCal2"  ]->SetYTitle("reco dv/dw at first tel. plane [rad]");
  _histoMap2D[ "eCal2"  ]->SetStats( false );
   
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
    
    double min, max; 
    std::string histoName;

    // Local track parameters 
    
    double  uBox = 1.1 * 0.5 * adet.GetModuleBoxSizeU();
    double  vBox = 1.1 * 0.5 * adet.GetModuleBoxSizeV();
    
    histoName = "htrk_dir_truth_det"+to_string( ipl ); 
    _histoMap2D[ histoName  ] = new TH2D(histoName.c_str(), "",1000, -0.03, 0.03,1000, -0.03, 0.03); 
    _histoMap2D[ histoName  ]->SetXTitle( "du/dw [rad]" );
    _histoMap2D[ histoName  ]->SetYTitle( "dv/dw [rad]");
    _histoMap2D[ histoName  ]->SetStats( false );
    _histoMap2D[ histoName  ]->GetYaxis()->SetTitleOffset(1.5); 
    
    histoName = "htrk_u_det"+to_string( ipl ); 
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 10000, -uBox, +uBox); 
    _histoMap[ histoName  ]->SetXTitle("intersect u [mm]");  
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
    
    
    histoName = "htrk_v_det"+to_string( ipl );
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 10000, -vBox, +vBox);
    _histoMap[ histoName  ]->SetXTitle("intersect v [mm]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
   
    histoName = "htrk_tu_det"+to_string( ipl );
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 100000, -0.3, 0.3); 
    _histoMap[ histoName  ]->SetXTitle("slope du/dw [rad]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks");    

    histoName = "htrk_tv_det"+to_string( ipl ); 
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 100000, -0.3, 0.3); 
    _histoMap[ histoName  ]->SetXTitle("slope dv/dw [rad]");
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
    
    histoName = "htrk_mom_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, 0, 1.2*_momentum); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 
    _histoMap[ histoName ]->SetXTitle("reco momentum [GeV]"); 
    
    histoName = "htrk_tu_truth_det"+to_string( ipl ); 
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 100000, -0.3, 0.3); 
    _histoMap[ histoName  ]->SetXTitle("true slope du/dw [rad]"); 
    _histoMap[ histoName ]->SetYTitle("tracks");     

    histoName = "htrk_tv_truth_det"+to_string( ipl );
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 100000, -0.3, 0.3); 
    _histoMap[ histoName  ]->SetXTitle("true slope dv/dw [rad]"); 
    _histoMap[ histoName ]->SetYTitle("tracks");     

    histoName = "htrk_u_truth_det"+to_string( ipl ); 
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 10000, -uBox, +uBox); 
    _histoMap[ histoName  ]->SetXTitle("true intersect u [mm]");  
    _histoMap[ histoName ]->SetYTitle("tracks");         

    histoName = "htrk_v_truth_det"+to_string( ipl );
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 10000, -vBox, +vBox);
    _histoMap[ histoName  ]->SetXTitle("true intersect v [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks");         

    histoName = "htrk_mom_truth_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, 0, 1.2*_momentum); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 
    _histoMap[ histoName ]->SetXTitle("true momentum [GeV]"); 
    
    // Local track parameter errors 
     
    histoName = "hsigma_u_det"+to_string( ipl ); 
    max = 100*adet.GetResolutionU();
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 8000, 0, max); 
    _histoMap[ histoName  ]->SetXTitle("sigma u [mm]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
    
    histoName = "hsigma_v_det"+to_string( ipl );
    max = 100*adet.GetResolutionV();
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 8000, 0, max);
    _histoMap[ histoName  ]->SetXTitle("sigma v [mm]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
    
    histoName = "hsigma_tu_det"+to_string( ipl );
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 4000, 0, 0.01); 
    _histoMap[ histoName  ]->SetXTitle("sigma du/dw [rad]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
    
    histoName = "hsigma_tv_det"+to_string( ipl );
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 4000, 0, 0.01);
    _histoMap[ histoName  ]->SetXTitle("sigma dv/dw [rad]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks");  

    histoName = "hsigma_mom_det"+to_string( ipl );
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 4000, 0, _momentum);
    _histoMap[ histoName  ]->SetXTitle("sigma p [GeV]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
       
    // Local track parameter pulls
    
    histoName = "hp1_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, -4., 4.);
    _histoMap[ histoName  ]->SetXTitle("pull slope du/dw"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
    
    histoName = "hp2_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, -4., 4.);
    _histoMap[ histoName  ]->SetXTitle("pull slope dv/dw"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
    
    histoName = "hp3_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, -4., 4.); 
    _histoMap[ histoName  ]->SetXTitle("pull intersect u"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
    
    histoName = "hp4_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, -4., 4.);  
    _histoMap[ histoName  ]->SetXTitle("pull intersect v"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
     
    histoName = "hp5_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, -4., 4.);  
    _histoMap[ histoName  ]->SetXTitle("pull q/p"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
    
    histoName = "hJ_det"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, 0, 40);
    _histoMap[ histoName ]->SetXTitle("local state #chi^{2}");
    _histoMap[ histoName ]->SetYTitle("tracks"); 
    
    histoName = "hJp_det"+to_string( ipl );
    _histoMap[histoName] = new TH1D(histoName.c_str(), "", 100, 0, 1 ); 
    _histoMap[histoName]->SetMinimum(0.);
    _histoMap[histoName]->SetXTitle("local state p-value"); 
    _histoMap[histoName]->SetYTitle("tracks");    
    
    histoName = "hresU_det"+to_string( ipl );
    min = -10*adet.GetResolutionU();
    max = +10*adet.GetResolutionU();
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, min, max);
    _histoMap[ histoName ]->SetXTitle("u residual [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 
    
    histoName = "hresV_det"+to_string( ipl );
    min = -10*adet.GetResolutionV();
    max = +10*adet.GetResolutionV();
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, min, max); 
    _histoMap[ histoName ]->SetXTitle("v residual [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 
     
    histoName = "hchi2inc_det"+to_string( ipl );    
    _histoMap[histoName] = new TH1D(histoName.c_str(), "", 100, 0, 20 );    
    _histoMap[histoName]->SetXTitle("#chi^{2} increment"); 
    _histoMap[histoName]->SetYTitle("tracks");  
    
    histoName = "hchi2incprob_det"+to_string( ipl );
    _histoMap[histoName] = new TH1D(histoName.c_str(), "", 100, 0, 1 ); 
    _histoMap[histoName]->SetMinimum(0.);
    _histoMap[histoName]->SetXTitle("increment p-value"); 
    _histoMap[histoName]->SetYTitle("tracks");    
  }
  
}


} // Namespace


