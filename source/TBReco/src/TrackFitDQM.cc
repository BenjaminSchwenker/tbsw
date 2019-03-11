// TrackFitDQM implementation file
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// Local includes 
#include "TrackFitDQM.h"

// TBTools includes
#include "TBDetector.h"
#include "TBTrack.h"
#include "TrackInputProvider.h"
#include "GenericTrackFitter.h"
#include "PixelCluster.h"
#include "ThreeDModel.h"

// C++ includes
#include <iomanip>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <Exceptions.h>

#include <TMath.h>

// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;

namespace depfet {

//
// Instantiate this object
//
TrackFitDQM aTrackFitDQM ;

//
// Constructor
//
TrackFitDQM::TrackFitDQM() : Processor("TrackFitDQM")
{
   
// Processor description
  _description = "TrackFitDQM: DQM plots for monitoring the track fit quality";
   

//
// Input collections  
  registerInputCollection(LCIO::TRACK,"InputTrackCollectionName",
                          "Track input collection",
                          _inputTrackCollectionName,std::string("tracks"));
   
// 
// Processor parameters
  
  registerProcessorParameter( "RootFileName",
                              "Output root file name",
                              _rootFileName, 
                              std::string("Histos.root"));
                                 
}

//
// Method called at the beginning of data processing
//
void TrackFitDQM::init() {
  
  // Initialize variables
  _nRun = 0 ;
  _nEvt = 0 ;
  
  // Print set parameters
  printProcessorParams();
  
  // CPU time start
  _timeCPU = clock()/1000;
  
  // Book all needed histograms 
  bookHistos();
}

//
// Method called for each run
//
void TrackFitDQM::processRunHeader(LCRunHeader * run)
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
void TrackFitDQM::processEvent(LCEvent * evt)
{
    
  //////////////////////////////////////////////////////////////////////  
  // Process next event
  ++_nEvt;
   
  if ( _nEvt % 1000 == 0 ) {
    streamlog_out( MESSAGE3 ) << "Processing event "
                              << evt->getEventNumber() << " in run "
                              << evt->getRunNumber() << endl; 
                               
  }

  TrackInputProvider TrackIO; 
  
  GenericTrackFitter TrackFitter(TBDetector::GetInstance());
  TrackFitter.SetNumIterations(1); 
   
  LCCollection* inputCollection;
  try {
      inputCollection = evt->getCollection(_inputTrackCollectionName);
  } catch (DataNotAvailableException& e) {
      throw SkipEventException(this);
  }
  
  // Main loop over all tracks
  int nTracks = inputCollection->getNumberOfElements(); 
  _overviewHistoMap["hntracks"]->Fill(nTracks);
  
  for (int itrk = 0; itrk < nTracks; itrk++) {
    
    // Retrieve track from LCIO 
    Track * inputtrack = dynamic_cast<Track*> (inputCollection->getElementAt(itrk));
    
    // Convert LCIO -> TB track  
    TBTrack track = TrackIO.MakeTBTrack( inputtrack, TBDetector::GetInstance() );  
    
    // ReFit track 
    bool trkerr = TrackFitter.Fit(track);
    if ( trkerr ) {
      continue;
    }  
    
    // Fill DQM histos 
    // ===============

    //
    // Top level histograms 
    
    _overviewHistoMap["hnhits"]->Fill(track.GetNumHits());
    _overviewHistoMap["htrkchi2"]->Fill(track.GetChiSqu());
    _overviewHistoMap["htrkchi2ndof"]->Fill(track.GetChiSqu()/track.GetNDF());
    _overviewHistoMap["hchi2prob"]->Fill(TMath::Prob(track.GetChiSqu(), track.GetNDF()));
    _overviewHistoMap["hmom"]->Fill(track.GetMomentum());
    _overviewHistoMap["hcharge"]->Fill(track.GetCharge());

    //
    // Sensor level histograms 
  
    // Get number of sensors
    int nSens = TBDetector::GetInstance().GetNSensors(); 


    for (int ipl= 0; ipl< nSens; ++ipl) {  
        
      // Get sensor data 
      //------------------------
      TBTrackElement& TE = track.GetTE(ipl);  
           
      // Get local track parameters 
      double trk_tu = TE.GetState().GetPars()[0];  // rad
      double trk_tv = TE.GetState().GetPars()[1];  // rad
      double trk_u = TE.GetState().GetPars()[2];   // mm
      double trk_v = TE.GetState().GetPars()[3];   // mm
      double trk_qp = TE.GetState().GetPars()[4];   // 1/GeV
         
      double trk_charge = track.GetCharge();
      double trk_mom = std::abs(trk_charge/trk_qp); 
      auto& _histoMap=_perLayerHistoMap[ipl];
      auto& _histoMap2D=_perLayerHistoMap2D[ipl];
      auto& _profileMap=_perLayerProfileMap[ipl];

      // Fill track parameter errors
      

      _histoMap[ "hsigma2_tu_sensor"  ]->Fill( TE.GetState().GetCov()(0,0) ) ;
      _histoMap[ "hsigma2_tv_sensor"  ]->Fill( TE.GetState().GetCov()(1,1) ) ;

      _histoMap[ "hsigma2_u_sensor"  ]->Fill( TE.GetState().GetCov()(2,2) );
      _histoMap[ "hsigma2_v_sensor"  ]->Fill( TE.GetState().GetCov()(3,3) );

      _histoMap[ "hsigma2_qp_sensor"  ]->Fill( TE.GetState().GetCov()(4,4) );
      
      //
      // Fill beam profile histograms 
      
      _profileMap[ "htrk_dudw_vs_u_sensor" ]->Fill( trk_u, trk_tu );
      _profileMap[ "htrk_dvdw_vs_u_sensor" ]->Fill( trk_u, trk_tv );

      _profileMap[ "htrk_dudw_vs_v_sensor" ]->Fill( trk_v, trk_tu );
      _profileMap[ "htrk_dvdw_vs_v_sensor" ]->Fill( trk_v, trk_tv );

      _profileMap[ "htrk_mom_vs_u_sensor" ]->Fill( trk_u, trk_mom );
      _profileMap[ "htrk_mom_vs_v_sensor" ]->Fill( trk_v, trk_mom );

      _histoMap[ "htrk_u_sensor" ]->Fill(trk_u);
      _histoMap[ "htrk_v_sensor" ]->Fill(trk_v);

      _histoMap2D[ "hhitmap_sensor"  ]->Fill(trk_u,trk_v);

      _histoMap[ "htrk_tu_sensor"  ]->Fill(trk_tu);
      _histoMap[ "htrk_tv_sensor"  ]->Fill(trk_tv);

      _histoMap[ "htrk_mom_sensor"  ]->Fill(trk_mom);
	
       
      // Skip sensor w/o measurment
      if ( !TE.HasHit() ) continue;  
       
      // Get pixel residuals 
      double du = TE.GetHit().GetCoord()[0] - TE.GetState().GetPars()[2]; // mm 
      double dv = TE.GetHit().GetCoord()[1] - TE.GetState().GetPars()[3]; // mm
               
      double pull_u = du / TMath::Sqrt( TE.GetState().GetCov()(2,2) + TE.GetHit().GetCov()(0,0) ) ; 
      double pull_v = dv / TMath::Sqrt( TE.GetState().GetCov()(3,3) + TE.GetHit().GetCov()(1,1) ) ;  
      
      //PixelCluster Cluster = TE.GetHit().GetCluster(); 
      
      _histoMap[ "hresU_sensor" ]->Fill( du );
      _histoMap[ "hresV_sensor" ]->Fill( dv );
      
      _histoMap[ "hpull_resU_sensor" ]->Fill( pull_u );
      _histoMap[ "hpull_resV_sensor" ]->Fill( pull_v );
       
      // Fill alignment profile 
      
      _profileMap[ "hduvsu_sensor" ]->Fill( trk_u , du);
      _profileMap[ "hdvvsv_sensor" ]->Fill( trk_v , dv );
      _profileMap[ "hduvsv_sensor" ]->Fill( trk_v , du );
      _profileMap[ "hdvvsu_sensor" ]->Fill( trk_u , dv );
            

      _profileMap[ "hduvsthetau_sensor" ]->Fill( trk_tu , du );
      _profileMap[ "hdvvsthetav_sensor" ]->Fill( trk_tv , dv );
                  
    } // end sensor loop   
    
                       
  } // end track loop 
}

//
// Method called after each event to check the data processed
//
void TrackFitDQM::check( LCEvent * )
{
}

//
// Method called after all data processing
//
void TrackFitDQM::end()
{

  // Loop over all sensors
  for (int ipl=0 ; ipl < TBDetector::GetInstance().GetNSensors(); ipl++) {
    
    // Fill summary histos on telescope resolution 
    auto& _histoMap=_perLayerHistoMap[ipl];

    double sigma_u = TMath::Sqrt(_histoMap[ "hsigma2_u_sensor" ]->GetMean());
    double sigma_error_u = 0; 
    if (sigma_u > 0) sigma_error_u = 0.5*_histoMap[ "hsigma2_u_sensor" ]->GetMeanError()/sigma_u;

    double sigma_v = TMath::Sqrt(_histoMap[ "hsigma2_v_sensor" ]->GetMean());
    double sigma_error_v = 0; 
    if (sigma_v > 0) sigma_error_v = 0.5*_histoMap[ "hsigma2_v_sensor" ]->GetMeanError()/sigma_v;

    double sigma_tu = TMath::Sqrt(_histoMap[ "hsigma2_tu_sensor" ]->GetMean());
    double sigma_error_tu = 0; 
    if (sigma_tu > 0) sigma_error_tu = 0.5*_histoMap[ "hsigma2_tu_sensor" ]->GetMeanError()/sigma_tu;
    
    double sigma_tv = TMath::Sqrt(_histoMap[ "hsigma2_tv_sensor" ]->GetMean());
    double sigma_error_tv = 0; 
    if (sigma_tv > 0) sigma_error_tv = 0.5*_histoMap[ "hsigma2_tv_sensor" ]->GetMeanError()/sigma_tv;
    
    double pull_rms_u = _histoMap[ "hpull_resU_sensor" ]->GetRMS();
    double pull_rms_error_u = _histoMap[ "hpull_resU_sensor" ]->GetRMSError();

    double pull_rms_v = _histoMap[ "hpull_resV_sensor" ]->GetRMS();
    double pull_rms_error_v = _histoMap[ "hpull_resV_sensor" ]->GetRMSError();


    double res_rms_u = _histoMap[ "hresU_sensor" ]->GetRMS();
    double res_rms_error_u = _histoMap[ "hresU_sensor" ]->GetRMSError();
    
    double res_rms_v = _histoMap[ "hresV_sensor" ]->GetRMS();
    double res_rms_error_v = _histoMap[ "hresV_sensor" ]->GetRMSError();

    _overviewHistoMap["hfit_sigma_u"]->SetBinContent(ipl+1, sigma_u );
    _overviewHistoMap["hfit_sigma_u"]->SetBinError(ipl+1, sigma_error_u );
    
    _overviewHistoMap["hfit_sigma_v"]->SetBinContent(ipl+1, sigma_v );
    _overviewHistoMap["hfit_sigma_v"]->SetBinError(ipl+1, sigma_error_v );

    _overviewHistoMap["hfit_sigma_tu"]->SetBinContent(ipl+1, sigma_tu );
    _overviewHistoMap["hfit_sigma_tu"]->SetBinError(ipl+1, sigma_error_tu );

    _overviewHistoMap["hfit_sigma_tv"]->SetBinContent(ipl+1, sigma_tv );
    _overviewHistoMap["hfit_sigma_tv"]->SetBinError(ipl+1, sigma_error_tv );
    
    _overviewHistoMap["hfit_pull_rms_u"]->SetBinContent(ipl+1, pull_rms_u );
    _overviewHistoMap["hfit_pull_rms_u"]->SetBinError(ipl+1, pull_rms_error_u );
    
    _overviewHistoMap["hfit_pull_rms_v"]->SetBinContent(ipl+1, pull_rms_v );
    _overviewHistoMap["hfit_pull_rms_v"]->SetBinError(ipl+1, pull_rms_error_v );
    
    _overviewHistoMap["hfit_res_rms_u"]->SetBinContent(ipl+1, res_rms_u );
    _overviewHistoMap["hfit_res_rms_u"]->SetBinError(ipl+1, res_rms_error_u );
    
    _overviewHistoMap["hfit_res_rms_v"]->SetBinContent(ipl+1, res_rms_v );
    _overviewHistoMap["hfit_res_rms_v"]->SetBinError(ipl+1, res_rms_error_v );

    
    
    _rootFile->cd("");

    // Fill summary histos on telescope alignment
    // ------------------------------  
     
	// This is the position vector of the sensor in the aligned telescope geometry
	auto pos_f = TBDetector::Get(ipl).GetNominal().GetPosition(); 
		
	// This is the rotation matrix of the sensor in the aligned telescope geometry; it 
	// contains a discrete and a continuous factor. 
	auto Rot_f = TBDetector::Get(ipl).GetNominal().GetRotation();

	// This is the discrete factor of sensor rotation in the aligned telescope geometry. 
	auto DRot = TBDetector::Get(ipl).GetDiscrete().GetRotation();
		
	// This is finally the continous factor of the rotation in the aligned telescope geometry
	auto CRot_f = Rot_f*DRot.transpose(); 
		
	// Euler angles are defined wrt. the continous rotation in the aligned telescope geometry
	double alpha_f, beta_f, gamma_f; 
	GetAnglesKarimaki(CRot_f, alpha_f, beta_f, gamma_f); 
	//
	// Fill alignment histograms

    _overviewHistoMap["hxshift_diff"]->SetBinContent(ipl+1,pos_f[0]);
    _overviewHistoMap["hyshift_diff"]->SetBinContent(ipl+1,pos_f[1]);
    _overviewHistoMap["hzshift_diff"]->SetBinContent(ipl+1,pos_f[2]);
    _overviewHistoMap["hxrot_diff"]->SetBinContent(ipl+1,alpha_f);
    _overviewHistoMap["hyrot_diff"]->SetBinContent(ipl+1,beta_f);
    _overviewHistoMap["hzrot_diff"]->SetBinContent(ipl+1,gamma_f);
      
    
  }
  
  _rootFile->Write();
  _rootFile->Close();
   
  delete _rootFile;
   
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
void TrackFitDQM::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "TrackFitDQM Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}


void TrackFitDQM::bookHistos()
{   
  
  _rootFile = new TFile( _rootFileName.c_str(),"recreate");
  _rootFile->cd("");
  
  //
  // Top level histograms 
   
  _overviewHistoMap["hntracks"] = new TH1D("hntracks", "", 30, 0, 30);
  _overviewHistoMap["hntracks"]->SetYTitle("events");
  _overviewHistoMap["hntracks"]->SetXTitle("number of tracks per event");

  _overviewHistoMap["hnhits"] = new TH1D("hnhits", "", 14, 0, 14);
  _overviewHistoMap["hnhits"]->SetYTitle("tracks");
  _overviewHistoMap["hnhits"]->SetXTitle("number of hits per track");
  
  _overviewHistoMap["hmom"] = new TH1D("hmom", "", 100, 0, 0);
  _overviewHistoMap["hmom"]->SetYTitle("tracks");
  _overviewHistoMap["hmom"]->SetXTitle("track momentum [GeV]");
  
  _overviewHistoMap["hcharge"] = new TH1D("hcharge", "", 10, -2, 2);
  _overviewHistoMap["hcharge"]->SetYTitle("tracks");
  _overviewHistoMap["hcharge"]->SetXTitle("track charge [e]");
  
  _overviewHistoMap["htrkchi2"] = new TH1D("htrkchi2", "", 400, 0, 100 );
  _overviewHistoMap["htrkchi2"]->SetYTitle("tracks");
  _overviewHistoMap["htrkchi2"]->SetXTitle("track #chi^{2}");
   
  _overviewHistoMap["htrkchi2ndof"] = new TH1D("htrkchi2ndof", "", 200, 0, 10);
  _overviewHistoMap["htrkchi2ndof"]->SetYTitle("tracks");
  _overviewHistoMap["htrkchi2ndof"]->SetXTitle("track #chi^{2}/ndof");

  _overviewHistoMap["hchi2prob"] = new TH1D("hchi2prob", "", 100, 0, 1);
  _overviewHistoMap["hchi2prob"]->SetXTitle("track p-value");
  _overviewHistoMap["hchi2prob"]->SetYTitle("tracks");
  _overviewHistoMap["hchi2prob"]->SetMinimum(0.);
  
  // Get number of sensors
  int nSens = TBDetector::GetInstance().GetNSensors();
  
  _overviewHistoMap["hfit_sigma_u"] = new TH1D("hfit_sigma_u","",nSens,0,nSens);
  _overviewHistoMap["hfit_sigma_u"]->SetStats( false );
  _overviewHistoMap["hfit_sigma_u"]->SetXTitle("plane number");
  _overviewHistoMap["hfit_sigma_u"]->SetYTitle("track fit sigma u [mm]");
  
  _overviewHistoMap["hfit_sigma_v"] = new TH1D("hfit_sigma_v","",nSens,0,nSens);
  _overviewHistoMap["hfit_sigma_v"]->SetStats( false );
  _overviewHistoMap["hfit_sigma_v"]->SetXTitle("plane number");
  _overviewHistoMap["hfit_sigma_v"]->SetYTitle("track fit sigma v [mm]");

  _overviewHistoMap["hfit_sigma_tu"] = new TH1D("hfit_sigma_tu","",nSens,0,nSens);
  _overviewHistoMap["hfit_sigma_tu"]->SetStats( false );
  _overviewHistoMap["hfit_sigma_tu"]->SetXTitle("plane number");
  _overviewHistoMap["hfit_sigma_tu"]->SetYTitle("track fit sigma dudw [rad]");
  
  _overviewHistoMap["hfit_sigma_tv"] = new TH1D("hfit_sigma_tv","",nSens,0,nSens);
  _overviewHistoMap["hfit_sigma_tv"]->SetStats( false );
  _overviewHistoMap["hfit_sigma_tv"]->SetXTitle("plane number");
  _overviewHistoMap["hfit_sigma_tv"]->SetYTitle("track fit sigma dvdw [rad]");
   
  _overviewHistoMap["hfit_pull_rms_u"] = new TH1D("hfit_pull_rms_u","",nSens,0,nSens);
  _overviewHistoMap["hfit_pull_rms_u"]->SetStats( false );
  _overviewHistoMap["hfit_pull_rms_u"]->SetXTitle("plane number");
  _overviewHistoMap["hfit_pull_rms_u"]->SetYTitle("RMS pull u residual");
  
  _overviewHistoMap["hfit_pull_rms_v"] = new TH1D("hfit_pull_rms_v","",nSens,0,nSens);
  _overviewHistoMap["hfit_pull_rms_v"]->SetStats( false );
  _overviewHistoMap["hfit_pull_rms_v"]->SetXTitle("plane number");
  _overviewHistoMap["hfit_pull_rms_v"]->SetYTitle("RMS pull v residual");

  _overviewHistoMap["hfit_res_rms_u"] = new TH1D("hfit_res_rms_u","",nSens,0,nSens);
  _overviewHistoMap["hfit_res_rms_u"]->SetStats( false );
  _overviewHistoMap["hfit_res_rms_u"]->SetXTitle("plane number");
  _overviewHistoMap["hfit_res_rms_u"]->SetYTitle("RMS u residual");
  
  _overviewHistoMap["hfit_res_rms_v"] = new TH1D("hfit_res_rms_v","",nSens,0,nSens);
  _overviewHistoMap["hfit_res_rms_v"]->SetStats( false );
  _overviewHistoMap["hfit_res_rms_v"]->SetXTitle("plane number");
  _overviewHistoMap["hfit_res_rms_v"]->SetYTitle("RMS v residual");


  // Create subdir for alignment plots
  TDirectory *alignDir = _rootFile->mkdir("alignment");
  alignDir->cd();


  _overviewHistoMap["hxshift_diff"] = new TH1D("hxshift_diff", "", nSens,0,nSens);
  _overviewHistoMap["hxshift_diff"]->SetYTitle("x shift diff [mm]");
  _overviewHistoMap["hxshift_diff"]->SetXTitle("sensor");

  _overviewHistoMap["hyshift_diff"] = new TH1D("hyshift_diff", "", nSens,0,nSens);
  _overviewHistoMap["hyshift_diff"]->SetYTitle("y shift diff [mm]");
  _overviewHistoMap["hyshift_diff"]->SetXTitle("sensor");

  _overviewHistoMap["hzshift_diff"] = new TH1D("hzshift_diff", "", nSens,0,nSens);
  _overviewHistoMap["hzshift_diff"]->SetYTitle("z shift diff [mm]");
  _overviewHistoMap["hzshift_diff"]->SetXTitle("sensor");

  _overviewHistoMap["hxrot_diff"] = new TH1D("hxrot_diff", "", nSens,0,nSens);
  _overviewHistoMap["hxrot_diff"]->SetYTitle("x rot diff [rad]");
  _overviewHistoMap["hxrot_diff"]->SetXTitle("sensor");

  _overviewHistoMap["hyrot_diff"] = new TH1D("hyrot_diff", "", nSens,0,nSens);
  _overviewHistoMap["hyrot_diff"]->SetYTitle("y rot diff [rad]");
  _overviewHistoMap["hyrot_diff"]->SetXTitle("sensor");

  _overviewHistoMap["hzrot_diff"] = new TH1D("hzrot_diff", "", nSens,0,nSens);
  _overviewHistoMap["hzrot_diff"]->SetYTitle("z rot diff [rad]");
  _overviewHistoMap["hzrot_diff"]->SetXTitle("sensor");

  // Change current directory to root
  _rootFile->cd("");
   
  //
  // Sensor level histograms 
      
  // Loop over all sensors
  for (int ipl=0 ; ipl < nSens; ipl++) {
  
    // Create subdirs for sensors
    std::string dirName = "Sensor"+to_string( ipl );
    TDirectory *sensDir = _rootFile->mkdir(dirName.c_str());
    sensDir->cd();
    auto& _histoMap=_perLayerHistoMap[ipl];
    auto& _histoMap2D=_perLayerHistoMap2D[ipl];
    auto& _profileMap=_perLayerProfileMap[ipl];
    double max; 
    int nbins; 
    double safetyFactor = 1.1;
    std::string histoName;
      
    // Get handle to sensor data
    Det & Sensor = TBDetector::GetInstance().GetDet(ipl); 
    
    // Temporary histograms used to compute the mean variance
    // of track parameters. Histograms will not be added to the 
    // root output file.  
    
    histoName = "hsigma2_u_sensor"+to_string( ipl );
    _histoMap[ "hsigma2_u_sensor"  ] = new TH1D(histoName.c_str(), "", 1, 0, 1);
    _histoMap[ "hsigma2_u_sensor"  ]->StatOverflows();
    _histoMap[ "hsigma2_u_sensor"  ]->SetDirectory(0);
    
    histoName = "hsigma2_v_sensor"+to_string( ipl );  
    _histoMap[ "hsigma2_v_sensor"  ] = new TH1D(histoName.c_str(), "", 1, 0, 1);
    _histoMap[ "hsigma2_v_sensor"  ]->StatOverflows();
    _histoMap[ "hsigma2_v_sensor"  ]->SetDirectory(0);
    
    histoName = "hsigma2_tu_sensor"+to_string( ipl );
    _histoMap[ "hsigma2_tu_sensor"  ] = new TH1D(histoName.c_str(), "", 1, 0, 1);
    _histoMap[ "hsigma2_tu_sensor"  ]->StatOverflows();
    _histoMap[ "hsigma2_tu_sensor"  ]->SetDirectory(0);

    histoName = "hsigma2_tv_sensor"+to_string( ipl );
    _histoMap[ "hsigma2_tv_sensor"  ] = new TH1D(histoName.c_str(), "", 1, 0, 1);
    _histoMap[ "hsigma2_tv_sensor"  ]->StatOverflows();
    _histoMap[ "hsigma2_tv_sensor"  ]->SetDirectory(0);
    
    histoName = "hsigma2_qp_sensor"+to_string( ipl );
    _histoMap[ "hsigma2_qp_sensor"  ] = new TH1D(histoName.c_str(), "", 1, 0, 1);
    _histoMap[ "hsigma2_qp_sensor"  ]->StatOverflows();
    _histoMap[ "hsigma2_qp_sensor"  ]->SetDirectory(0);
    
    // Plot residuals U/V 
    
    histoName = "hresU_sensor"+to_string( ipl );
    max = 5*safetyFactor*( Sensor.GetSensitiveMaxU() - Sensor.GetSensitiveMinU()) /(Sensor.GetMaxUCell()+2); 
    _histoMap[ "hresU_sensor" ] = new TH1D(histoName.c_str(), "", 500, -max, +max);
    _histoMap[ "hresU_sensor" ]->SetXTitle("u residual [mm]");
    _histoMap[ "hresU_sensor" ]->SetYTitle("tracks");
    
    histoName = "hresV_sensor"+to_string( ipl );
    max = 5*safetyFactor*( Sensor.GetSensitiveMaxV() - Sensor.GetSensitiveMinV())/(Sensor.GetMaxVCell()+2); 
    _histoMap[ "hresV_sensor" ] = new TH1D(histoName.c_str(), "", 500, -max, +max);
    _histoMap[ "hresV_sensor" ]->SetXTitle("v residual [mm]");
    _histoMap[ "hresV_sensor" ]->SetYTitle("tracks");
    
    // Plot residuals pulls  

    histoName = "hpull_resU_sensor"+to_string( ipl );
    _histoMap[ "hpull_resU_sensor" ] = new TH1D(histoName.c_str(), "", 100, -4, +4);
    _histoMap[ "hpull_resU_sensor" ]->SetXTitle(" pull u residual");
    _histoMap[ "hpull_resU_sensor" ]->SetYTitle(" tracks");
    
    histoName = "hpull_resV_sensor"+to_string( ipl );
    _histoMap[ "hpull_resV_sensor" ] = new TH1D(histoName.c_str(), "", 100, -4, +4);
    _histoMap[ "hpull_resV_sensor" ]->SetXTitle(" pull v residual");
    _histoMap[ "hpull_resV_sensor" ]->SetYTitle(" tracks");
           
    // Plot residual profiles 
    
    histoName = "hduvsu_sensor"+to_string( ipl );
    double maxU = safetyFactor*Sensor.GetSensitiveMaxU(); 
    double minU = safetyFactor*Sensor.GetSensitiveMinU(); 
    nbins = 100; 
    _profileMap[ "hduvsu_sensor" ] = new TProfile(histoName.c_str(), "", nbins, minU, maxU);
    _profileMap[ "hduvsu_sensor" ]->SetXTitle("u [mm]");
    _profileMap[ "hduvsu_sensor" ]->SetYTitle("mean residual u [mm]");
      
    histoName = "hdvvsv_sensor"+to_string( ipl );
    double maxV = safetyFactor*Sensor.GetSensitiveMaxV(); 
    double minV = safetyFactor*Sensor.GetSensitiveMinV();  
    nbins = 100;  
    _profileMap[ "hdvvsv_sensor" ] = new TProfile(histoName.c_str(), "", nbins, minV, +maxV);
    _profileMap[ "hdvvsv_sensor" ]->SetXTitle("v [mm]");
    _profileMap[ "hdvvsv_sensor" ]->SetYTitle("mean residual v [mm]");
    
    histoName = "hduvsv_sensor"+to_string( ipl );
    nbins = 100;  
    _profileMap[ "hduvsv_sensor" ] = new TProfile(histoName.c_str(), "", nbins, minV, +maxV);
    _profileMap[ "hduvsv_sensor" ]->SetXTitle("v [mm]");
    _profileMap[ "hduvsv_sensor" ]->SetYTitle("mean residual u [mm]");
       
    histoName = "hdvvsu_sensor"+to_string( ipl );
    nbins = 100; 
    _profileMap[ "hdvvsu_sensor" ] = new TProfile(histoName.c_str(), "", nbins, minU, +maxU);
    _profileMap[ "hdvvsu_sensor" ]->SetXTitle("u [mm]");
    _profileMap[ "hdvvsu_sensor" ]->SetYTitle("mean residual v [mm]");
          
    histoName = "hduvsthetau_sensor"+to_string( ipl );
    max = 0;  // max track slope [rad]
    _profileMap[ "hduvsthetau_sensor" ] = new TProfile(histoName.c_str(), "", 100, -max, +max);
    _profileMap[ "hduvsthetau_sensor" ]->SetXTitle("du/dw [rad]");
    _profileMap[ "hduvsthetau_sensor" ]->SetYTitle("mean residual u [mm]");
    
    histoName = "hdvvsthetav_sensor"+to_string( ipl );
    max = 0;  // max track slope [rad]
    _profileMap[ "hdvvsthetav_sensor" ] = new TProfile(histoName.c_str(), "", 100, -max, +max);
    _profileMap[ "hdvvsthetav_sensor" ]->SetXTitle("dv/dw [rad]");
    _profileMap[ "hdvvsthetav_sensor" ]->SetYTitle("mean residual v [mm]");

    // Plot beam profiles
    TDirectory *beamDir = sensDir->mkdir("BeamProfile");
    beamDir->cd(); 
    
    histoName = "htrk_mom_sensor"+to_string( ipl );
    _histoMap[ "htrk_mom_sensor" ] = new TH1D(histoName.c_str(), "", 100, 0, 0);
    _histoMap[ "htrk_mom_sensor" ]->SetXTitle("track momentum [GeV]");
    _histoMap[ "htrk_mom_sensor" ]->SetYTitle("tracks");
    
    histoName = "htrk_tu_sensor"+to_string( ipl ); 
    _histoMap[ "htrk_tu_sensor" ] = new TH1D(histoName.c_str(), "", 100, 0, 0);
    _histoMap[ "htrk_tu_sensor" ]->SetXTitle("fit du/dw [rad]");
    _histoMap[ "htrk_tu_sensor" ]->SetYTitle("tracks");
    
    histoName = "htrk_tv_sensor"+to_string( ipl ); 
    _histoMap[ "htrk_tv_sensor" ] = new TH1D(histoName.c_str(), "", 100, 0, 0);
    _histoMap[ "htrk_tv_sensor" ]->SetXTitle("fit dv/dw [rad]");
    _histoMap[ "htrk_tv_sensor" ]->SetYTitle("tracks");
      
    
    nbins = 100;           
    histoName = "htrk_u_sensor"+to_string( ipl ); 
    _histoMap[ "htrk_u_sensor" ] = new TH1D(histoName.c_str(), "", nbins, minU, +maxU);
    _histoMap[ "htrk_u_sensor" ]->SetXTitle("fit u [mm]");
    _histoMap[ "htrk_u_sensor" ]->SetYTitle("tracks");

    histoName = "htrk_v_sensor"+to_string( ipl ); 
    _histoMap[ "htrk_v_sensor" ] = new TH1D(histoName.c_str(), "", nbins, minV, maxV);
    _histoMap[ "htrk_v_sensor" ]->SetXTitle("fit v [mm]");
    _histoMap[ "htrk_v_sensor" ]->SetYTitle("tracks");

    histoName = "hhitmap_sensor"+to_string( ipl );
    _histoMap2D["hhitmap_sensor"] = new TH2D(histoName.c_str(), "" ,nbins, minU, +maxU, nbins, minV, +maxV);
    _histoMap2D["hhitmap_sensor"]->SetXTitle("fit u [mm]");
    _histoMap2D["hhitmap_sensor"]->SetYTitle("fit v [mm]");
    _histoMap2D["hhitmap_sensor"]->SetZTitle("tracks");
    _histoMap2D["hhitmap_sensor"]->SetStats( false );
    
    histoName = "htrk_dudw_vs_u_sensor"+to_string( ipl );
    _profileMap[ "htrk_dudw_vs_u_sensor" ] = new TProfile(histoName.c_str(), "", nbins, minU, +maxU);
    _profileMap[ "htrk_dudw_vs_u_sensor" ]->SetYTitle("fit du/dw [rad]");
    _profileMap[ "htrk_dudw_vs_u_sensor" ]->SetXTitle("fit u [mm]");
     
    histoName = "htrk_dvdw_vs_u_sensor"+to_string( ipl );
    _profileMap[ "htrk_dvdw_vs_u_sensor" ] = new TProfile(histoName.c_str(), "", nbins, minU, +maxU);
    _profileMap[ "htrk_dvdw_vs_u_sensor" ]->SetYTitle("fit dv/dw [rad]");
    _profileMap[ "htrk_dvdw_vs_u_sensor" ]->SetXTitle("fit u [mm]");
    
    histoName = "htrk_dudw_vs_v_sensor"+to_string( ipl );
    _profileMap[ "htrk_dudw_vs_v_sensor" ] = new TProfile(histoName.c_str(), "", nbins, minV, +maxV);
    _profileMap[ "htrk_dudw_vs_v_sensor" ]->SetYTitle("fit du/dw [rad]");
    _profileMap[ "htrk_dudw_vs_v_sensor" ]->SetXTitle("fit v [mm]");
    
    histoName = "htrk_dvdw_vs_v_sensor"+to_string( ipl );
    _profileMap[ "htrk_dvdw_vs_v_sensor" ] = new TProfile(histoName.c_str(), "", nbins, minV, +maxV);
    _profileMap[ "htrk_dvdw_vs_v_sensor" ]->SetYTitle("dv/dw [rad]");
    _profileMap[ "htrk_dvdw_vs_v_sensor" ]->SetXTitle("v [mm]");
    
    histoName = "htrk_mom_vs_u_sensor"+to_string( ipl );
    _profileMap[ "htrk_mom_vs_u_sensor" ] = new TProfile(histoName.c_str(), "",  nbins, minU, +maxU);
    _profileMap[ "htrk_mom_vs_u_sensor" ]->SetYTitle("momentum [GeV]");
    _profileMap[ "htrk_mom_vs_u_sensor" ]->SetXTitle("u [mm]");
    
    histoName = "htrk_mom_vs_v_sensor"+to_string( ipl );
    _profileMap[ "htrk_mom_vs_v_sensor" ] = new TProfile(histoName.c_str(), "", nbins, minV, +maxV);
    _profileMap[ "htrk_mom_vs_v_sensor" ]->SetYTitle("momentum [GeV]");
    _profileMap[ "htrk_mom_vs_v_sensor" ]->SetXTitle("v [mm]");
    
    
    
    // Change current directory to root
    _rootFile->cd("");
  }
}

} // Namespace
