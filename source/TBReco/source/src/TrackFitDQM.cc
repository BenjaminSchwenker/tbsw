// TrackFitDQM implementation file
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// Local includes 
#include "TrackFitDQM.h"

// TBTools includes
#include "TBTrack.h"
#include "TrackInputProvider.h"
#include "GenericTrackFitter.h"
#include "PixelCluster.h"

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
using namespace CLHEP;

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
  
  registerProcessorParameter ("AlignmentDBFileName",
                             "This is the name of the LCIO file with the alignment constants (add .slcio)",
                             _alignmentDBFileName, 
                             static_cast< string > ( "alignmentDB.slcio" ) ); 
  
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
  
  // Read detector constants from gear file
  _detector.ReadGearConfiguration();    
      
  // Read alignment data base file 
  _detector.ReadAlignmentDB( _alignmentDBFileName );    
  
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
  
  GenericTrackFitter TrackFitter(_detector);
  TrackFitter.SetNumIterations(2); 
   
  LCCollection* inputCollection;
  try {
      inputCollection = evt->getCollection(_inputTrackCollectionName);
  } catch (DataNotAvailableException& e) {
      throw SkipEventException(this);
  }
  
  // Main loop over all tracks
  int nTracks = inputCollection->getNumberOfElements(); 
  _histoMap["hntracks"]->Fill(nTracks); 
  
  for (int itrk = 0; itrk < nTracks; itrk++) {
    
    // Retrieve track from LCIO 
    Track * inputtrack = dynamic_cast<Track*> (inputCollection->getElementAt(itrk));
    
    // Convert LCIO -> TB track  
    TBTrack track = TrackIO.MakeTBTrack( inputtrack, _detector );  
    
    // ReFit track 
    bool trkerr = TrackFitter.Fit(track);
    if ( trkerr ) {
      continue;
    }  
    
    // Fill DQM histos 
    // ===============

    //
    // Top level histograms 
    
    _histoMap["hnhits"]->Fill(track.GetNumHits());  
    _histoMap["htrkchi2"]->Fill(track.GetChiSqu());
    _histoMap["htrkchi2ndof"]->Fill(track.GetChiSqu()/track.GetNDF()); 
    _histoMap["hchi2prob"]->Fill(TMath::Prob(track.GetChiSqu(), track.GetNDF()));
    _histoMap["hmom"]->Fill(track.GetMomentum());  
    _histoMap["hcharge"]->Fill(track.GetCharge());  
      
    //
    // Sensor level histograms 
  
    // Get number of sensors
    int nSens = _detector.GetNSensors(); 
    std::string histoName; 

    for (int ipl= 0; ipl< nSens; ++ipl) {  
        
      // Get sensor data 
      //------------------------
      TBTrackElement& TE = track.GetTE(ipl);  
           
      // Get local track parameters 
      double trk_tu = TE.GetState().GetPars()[0][0];  // rad
      double trk_tv = TE.GetState().GetPars()[1][0];  // rad
      double trk_u = TE.GetState().GetPars()[2][0];   // mm
      double trk_v = TE.GetState().GetPars()[3][0];   // mm
      double trk_qp = TE.GetState().GetPars()[4][0];   // 1/GeV
         
      double trk_charge = track.GetCharge();
      double trk_mom = std::abs(trk_charge/trk_qp); 

      // Fill track parameter errors
      
      histoName = "hsigma2_tu_sensor"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TE.GetState().GetCov()[0][0] ) ; 
      
      histoName = "hsigma2_tv_sensor"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TE.GetState().GetCov()[1][1] ) ;  
      
      histoName = "hsigma2_u_sensor"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TE.GetState().GetCov()[2][2] ); 
      
      histoName = "hsigma2_v_sensor"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TE.GetState().GetCov()[3][3] );
      
      histoName = "hsigma2_qp_sensor"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TE.GetState().GetCov()[4][4] );
      
      //
      // Fill beam profile histograms 
      
      histoName = "htrk_dudw_vs_u_sensor"+to_string( ipl );
      _profileMap[ histoName ]->Fill( trk_u, trk_tu );
      
      histoName = "htrk_dvdw_vs_u_sensor"+to_string( ipl );
      _profileMap[ histoName ]->Fill( trk_u, trk_tv );
      
      histoName = "htrk_dudw_vs_v_sensor"+to_string( ipl );
      _profileMap[ histoName ]->Fill( trk_v, trk_tu );
       
      histoName = "htrk_dvdw_vs_v_sensor"+to_string( ipl );
      _profileMap[ histoName ]->Fill( trk_v, trk_tv );   
      
      histoName = "htrk_mom_vs_u_sensor"+to_string( ipl );
      _profileMap[ histoName ]->Fill( trk_u, trk_mom );
      
      histoName = "htrk_mom_vs_v_sensor"+to_string( ipl );
      _profileMap[ histoName ]->Fill( trk_v, trk_mom );   

      histoName = "htrk_u_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill(trk_u); 
      
      histoName = "htrk_v_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill(trk_v);
      
      histoName = "hhitmap_sensor"+to_string( ipl );
      _histoMap2D[ histoName  ]->Fill(trk_u,trk_v);   
      
      histoName = "htrk_tu_sensor"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(trk_tu);
       
      histoName = "htrk_tv_sensor"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(trk_tv);
      
      histoName = "htrk_mom_sensor"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(trk_mom);
       
      // Skip sensor w/o measurment
      if ( !TE.HasHit() ) continue;  
       
      // Get pixel residuals 
      double du = TE.GetHit().GetCoord()[0][0] - TE.GetState().GetPars()[2][0]; // mm 
      double dv = TE.GetHit().GetCoord()[1][0] - TE.GetState().GetPars()[3][0]; // mm
               
      double pull_u = du / TMath::Sqrt( TE.GetState().GetCov()[2][2] + TE.GetHit().GetCov()[0][0] ) ; 
      double pull_v = dv / TMath::Sqrt( TE.GetState().GetCov()[3][3] + TE.GetHit().GetCov()[1][1] ) ;  
      
      PixelCluster Cluster = TE.GetHit().GetCluster(); 
      
      string id = Cluster.getDigitalClusterID();
      
      // Register new cluster if needed
      if (_clusterSpectrumMap[ipl].find(id) == _clusterSpectrumMap[ipl].end() ) {
        //  Create a new folder for histos for this cluster id  
        TDirectory *newfolder = _clusterDirMap[ipl]->mkdir(id.c_str());    
        newfolder->cd(); 

        double maxU = 5.0*_detector.GetDet(ipl).GetPitchU(); 
        double maxV = 5.0*_detector.GetDet(ipl).GetPitchV(); 
        
        histoName = "hresU_sensor"+to_string( ipl )+"_"+id;
        _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 500, -maxU, +maxU);
        _histoMap[ histoName ]->SetXTitle("u residual [mm]"); 
        _histoMap[ histoName ]->SetYTitle("tracks"); 
        
        histoName = "hresV_sensor"+to_string( ipl )+"_"+id;
        _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 500, -maxV, +maxV);
        _histoMap[ histoName ]->SetXTitle("v residual [mm]"); 
        _histoMap[ histoName ]->SetYTitle("tracks"); 
        
        histoName = "hpull_resU_sensor"+to_string( ipl )+"_"+id;
        _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, -4, +4);
        _histoMap[ histoName ]->SetXTitle(" pull u residual"); 
        _histoMap[ histoName ]->SetYTitle(" tracks");    
        
        histoName = "hpull_resV_sensor"+to_string( ipl )+"_"+id;
        _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, -4, +4);
        _histoMap[ histoName ]->SetXTitle(" pull v residual"); 
        _histoMap[ histoName ]->SetYTitle(" tracks");          
        
        // Create a counter for occurences of this cluster id
        _clusterSpectrumMap[ipl][id] = 0;
         
        _rootFile->cd("");
      }
         
      // Increment the clusterID counter 
      _clusterSpectrumMap[ipl][id]++;
      
      histoName = "hresU_sensor"+to_string( ipl )+"_"+id;
      _histoMap[ histoName ]->Fill( du ); 

      histoName = "hresV_sensor"+to_string( ipl )+"_"+id;
      _histoMap[ histoName ]->Fill( dv ); 

      histoName = "hpull_resU_sensor"+to_string( ipl )+"_"+id;
      _histoMap[ histoName ]->Fill( pull_u );

      histoName = "hpull_resV_sensor"+to_string( ipl )+"_"+id;
      _histoMap[ histoName ]->Fill( pull_v );

      histoName = "hresU_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill( du ); 
       
      histoName = "hresV_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill( dv ); 
      
      histoName = "hpull_resU_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill( pull_u );
      
      histoName = "hpull_resV_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill( pull_v );
       
      // Fill alignment profile 
      
      histoName = "hduvsu_sensor"+to_string( ipl );
      _profileMap[ histoName ]->Fill( trk_u , du);
      
      histoName = "hdvvsv_sensor"+to_string( ipl );
      _profileMap[ histoName ]->Fill( trk_v , dv );
      
      histoName = "hduvsv_sensor"+to_string( ipl );
      _profileMap[ histoName ]->Fill( trk_v , du );
      
      histoName = "hdvvsu_sensor"+to_string( ipl );
      _profileMap[ histoName ]->Fill( trk_u , dv );
            
      histoName = "hduvsthetau_sensor"+to_string( ipl );
      _profileMap[ histoName ]->Fill( trk_tu , du );   
      
      histoName = "hdvvsthetav_sensor"+to_string( ipl );
      _profileMap[ histoName ]->Fill( trk_tv , dv );
                  
    } // end sensor loop   
    
                       
  } // end track loop 
}

//
// Method called after each event to check the data processed
//
void TrackFitDQM::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void TrackFitDQM::end()
{

  // Loop over all sensors
  for (int ipl=0 ; ipl < _detector.GetNSensors(); ipl++) {
    
    std::string histoName;

    // Fill summary histos on telescope resolution 
    
    histoName = "hsigma2_u_sensor"+to_string( ipl );
    double sigma_u = TMath::Sqrt(_histoMap[ histoName ]->GetMean());
    double sigma_error_u = 0; 
    if (sigma_u > 0) sigma_error_u = 0.5*_histoMap[ histoName ]->GetMeanError()/sigma_u;    

    histoName = "hsigma2_v_sensor"+to_string( ipl );
    double sigma_v = TMath::Sqrt(_histoMap[ histoName ]->GetMean());
    double sigma_error_v = 0; 
    if (sigma_v > 0) sigma_error_v = 0.5*_histoMap[ histoName ]->GetMeanError()/sigma_v;

    histoName = "hsigma2_tu_sensor"+to_string( ipl );
    double sigma_tu = TMath::Sqrt(_histoMap[ histoName ]->GetMean());
    double sigma_error_tu = 0; 
    if (sigma_tu > 0) sigma_error_tu = 0.5*_histoMap[ histoName ]->GetMeanError()/sigma_tu;
    
    histoName = "hsigma2_tv_sensor"+to_string( ipl );
    double sigma_tv = TMath::Sqrt(_histoMap[ histoName ]->GetMean());
    double sigma_error_tv = 0; 
    if (sigma_tv > 0) sigma_error_tv = 0.5*_histoMap[ histoName ]->GetMeanError()/sigma_tv;
    
    histoName = "hpull_resU_sensor"+to_string( ipl );
    double pull_rms_u = _histoMap[ histoName ]->GetRMS();
    double pull_rms_error_u = _histoMap[ histoName ]->GetRMSError(); 

    histoName = "hpull_resV_sensor"+to_string( ipl );
    double pull_rms_v = _histoMap[ histoName ]->GetRMS();
    double pull_rms_error_v = _histoMap[ histoName ]->GetRMSError();    

    histoName = "hresU_sensor"+to_string( ipl );
    double res_rms_u = _histoMap[ histoName ]->GetRMS();
    double res_rms_error_u = _histoMap[ histoName ]->GetRMSError();
    
    histoName = "hresV_sensor"+to_string( ipl );
    double res_rms_v = _histoMap[ histoName ]->GetRMS();
    double res_rms_error_v = _histoMap[ histoName ]->GetRMSError();

    _histoMap["hfit_sigma_u"]->SetBinContent(ipl+1, sigma_u ); 
    _histoMap["hfit_sigma_u"]->SetBinError(ipl+1, sigma_error_u ); 
    
    _histoMap["hfit_sigma_v"]->SetBinContent(ipl+1, sigma_v ); 
    _histoMap["hfit_sigma_v"]->SetBinError(ipl+1, sigma_error_v ); 

    _histoMap["hfit_sigma_tu"]->SetBinContent(ipl+1, sigma_tu ); 
    _histoMap["hfit_sigma_tu"]->SetBinError(ipl+1, sigma_error_tu ); 

    _histoMap["hfit_sigma_tv"]->SetBinContent(ipl+1, sigma_tv );  
    _histoMap["hfit_sigma_tv"]->SetBinError(ipl+1, sigma_error_tv ); 
    
    _histoMap["hfit_pull_rms_u"]->SetBinContent(ipl+1, pull_rms_u );
    _histoMap["hfit_pull_rms_u"]->SetBinError(ipl+1, pull_rms_error_u );  
    
    _histoMap["hfit_pull_rms_v"]->SetBinContent(ipl+1, pull_rms_v );   
    _histoMap["hfit_pull_rms_v"]->SetBinError(ipl+1, pull_rms_error_v );  
    
    _histoMap["hfit_res_rms_u"]->SetBinContent(ipl+1, res_rms_u ); 
    _histoMap["hfit_res_rms_u"]->SetBinError(ipl+1, res_rms_error_u );  
    
    _histoMap["hfit_res_rms_v"]->SetBinContent(ipl+1, res_rms_v );   
    _histoMap["hfit_res_rms_v"]->SetBinError(ipl+1, res_rms_error_v );  

    // Fill summary histos on cluster shapes
    
    // Change into cluster shape folder
    _clusterDirMap[ipl]->cd();    
    
    // Book histograms for clusterDB
    int NCLUSTERS = _clusterSpectrumMap[ipl].size();   
    
    // No clusterID means nothing to do
    if (NCLUSTERS == 0) continue;
    
    histoName = "hClusterID_fractions_sensor"+to_string( ipl );
    _histoMap[histoName] = new TH1D(histoName.c_str(),"",NCLUSTERS,0,NCLUSTERS);
    _histoMap[histoName]->SetYTitle("clusterID fraction"); 
    _histoMap[histoName]->SetStats( false );
    
    histoName = "hClusterID_rmsU_sensor"+to_string( ipl );
    _histoMap[histoName] = new TH1D(histoName.c_str(),"",NCLUSTERS,0,NCLUSTERS);
    _histoMap[histoName]->SetYTitle("RMS u residual [mm]"); 
    _histoMap[histoName]->SetStats( false );
    
    histoName = "hClusterID_rmsV_sensor"+to_string( ipl );
    _histoMap[histoName] = new TH1D(histoName.c_str(),"",NCLUSTERS,0,NCLUSTERS);
    _histoMap[histoName]->SetYTitle("RMS v residual [mm]"); 
    _histoMap[histoName]->SetStats( false );
    
    
    
    int i = 0; 
    
    // Go through all cluster shapes
    for (auto iter =_clusterSpectrumMap[ipl].begin(); iter!=_clusterSpectrumMap[ipl].end(); iter++ ) {
      int count = iter->second;  
      string id = iter->first;
      i++;  
       
      TH1D * tmpHistU = _histoMap["hresU_sensor"+to_string( ipl )+"_"+id];
      TH1D * tmpHistV = _histoMap["hresV_sensor"+to_string( ipl )+"_"+id];
      
      histoName = "hClusterID_fractions_sensor"+to_string( ipl );
      _histoMap[histoName]->SetBinContent( i, count );
      _histoMap[histoName]->SetBinError( i, TMath::Sqrt(count) );
      _histoMap[histoName]->GetXaxis()->SetBinLabel( i, id.c_str() );

      histoName = "hClusterID_rmsU_sensor"+to_string( ipl );
      _histoMap[histoName]->SetBinContent( i, tmpHistU->GetRMS() );
      _histoMap[histoName]->SetBinError( i, tmpHistU->GetRMSError() );
      _histoMap[histoName]->GetXaxis()->SetBinLabel( i, id.c_str() );
      
      histoName = "hClusterID_rmsV_sensor"+to_string( ipl );
      _histoMap[histoName]->SetBinContent( i, tmpHistV->GetRMS() );
      _histoMap[histoName]->SetBinError( i, tmpHistV->GetRMSError() );
      _histoMap[histoName]->GetXaxis()->SetBinLabel( i, id.c_str() );  
      
    }  
    
    _rootFile->cd("");
    
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
   
  _histoMap["hntracks"] = new TH1D("hntracks", "", 14, 0, 14); 
  _histoMap["hntracks"]->SetYTitle("events"); 
  _histoMap["hntracks"]->SetXTitle("number of tracks per event"); 

  _histoMap["hnhits"] = new TH1D("hnhits", "", 14, 0, 14); 
  _histoMap["hnhits"]->SetYTitle("tracks"); 
  _histoMap["hnhits"]->SetXTitle("number of hits per track"); 
  
  _histoMap["hmom"] = new TH1D("hmom", "", 100, 0, 0); 
  _histoMap["hmom"]->SetYTitle("tracks"); 
  _histoMap["hmom"]->SetXTitle("track momentum [GeV]");
  
  _histoMap["hcharge"] = new TH1D("hcharge", "", 10, -2, 2); 
  _histoMap["hcharge"]->SetYTitle("tracks"); 
  _histoMap["hcharge"]->SetXTitle("track charge [e]");
  
  _histoMap["htrkchi2"] = new TH1D("htrkchi2", "", 400, 0, 100 ); 
  _histoMap["htrkchi2"]->SetYTitle("tracks"); 
  _histoMap["htrkchi2"]->SetXTitle("track #chi^{2}");
   
  _histoMap["htrkchi2ndof"] = new TH1D("htrkchi2ndof", "", 200, 0, 10); 
  _histoMap["htrkchi2ndof"]->SetYTitle("tracks"); 
  _histoMap["htrkchi2ndof"]->SetXTitle("track #chi^{2}/ndof"); 

  _histoMap["hchi2prob"] = new TH1D("hchi2prob", "", 100, 0, 1);
  _histoMap["hchi2prob"]->SetXTitle("track p-value"); 
  _histoMap["hchi2prob"]->SetYTitle("tracks"); 
  _histoMap["hchi2prob"]->SetMinimum(0.);  
  
  // Get number of sensors
  int nSens = _detector.GetNSensors();
  
  _histoMap["hfit_sigma_u"] = new TH1D("hfit_sigma_u","",nSens,0,nSens);    
  _histoMap["hfit_sigma_u"]->SetStats( false );
  _histoMap["hfit_sigma_u"]->SetXTitle("plane number"); 
  _histoMap["hfit_sigma_u"]->SetYTitle("track fit sigma u [mm]"); 
  
  _histoMap["hfit_sigma_v"] = new TH1D("hfit_sigma_v","",nSens,0,nSens);    
  _histoMap["hfit_sigma_v"]->SetStats( false );
  _histoMap["hfit_sigma_v"]->SetXTitle("plane number"); 
  _histoMap["hfit_sigma_v"]->SetYTitle("track fit sigma v [mm]"); 

  _histoMap["hfit_sigma_tu"] = new TH1D("hfit_sigma_tu","",nSens,0,nSens);    
  _histoMap["hfit_sigma_tu"]->SetStats( false );
  _histoMap["hfit_sigma_tu"]->SetXTitle("plane number"); 
  _histoMap["hfit_sigma_tu"]->SetYTitle("track fit sigma dudw [rad]"); 
  
  _histoMap["hfit_sigma_tv"] = new TH1D("hfit_sigma_tv","",nSens,0,nSens);    
  _histoMap["hfit_sigma_tv"]->SetStats( false );
  _histoMap["hfit_sigma_tv"]->SetXTitle("plane number"); 
  _histoMap["hfit_sigma_tv"]->SetYTitle("track fit sigma dvdw [rad]"); 
   
  _histoMap["hfit_pull_rms_u"] = new TH1D("hfit_pull_rms_u","",nSens,0,nSens);    
  _histoMap["hfit_pull_rms_u"]->SetStats( false );
  _histoMap["hfit_pull_rms_u"]->SetXTitle("plane number"); 
  _histoMap["hfit_pull_rms_u"]->SetYTitle("RMS pull u residual"); 
  
  _histoMap["hfit_pull_rms_v"] = new TH1D("hfit_pull_rms_v","",nSens,0,nSens);    
  _histoMap["hfit_pull_rms_v"]->SetStats( false );
  _histoMap["hfit_pull_rms_v"]->SetXTitle("plane number"); 
  _histoMap["hfit_pull_rms_v"]->SetYTitle("RMS pull v residual");  

  _histoMap["hfit_res_rms_u"] = new TH1D("hfit_res_rms_u","",nSens,0,nSens);    
  _histoMap["hfit_res_rms_u"]->SetStats( false );
  _histoMap["hfit_res_rms_u"]->SetXTitle("plane number"); 
  _histoMap["hfit_res_rms_u"]->SetYTitle("RMS u residual"); 
  
  _histoMap["hfit_res_rms_v"] = new TH1D("hfit_res_rms_v","",nSens,0,nSens);    
  _histoMap["hfit_res_rms_v"]->SetStats( false );
  _histoMap["hfit_res_rms_v"]->SetXTitle("plane number"); 
  _histoMap["hfit_res_rms_v"]->SetYTitle("RMS v residual");  
   
  //
  // Sensor level histograms 
      
  // Loop over all sensors
  for (int ipl=0 ; ipl < nSens; ipl++) {
  
    // Create subdirs for sensors
    std::string dirName = "Sensor"+to_string( ipl );
    TDirectory *sensDir = _rootFile->mkdir(dirName.c_str());
    sensDir->cd();
    
    double max; 
    int nbins; 
    double safetyFactor = 1.1;
    std::string histoName;
      
    // Get handle to sensor data
    Det & Sensor = _detector.GetDet(ipl); 
    
    // Temporary histograms used to compute the mean variance
    // of track parameters. Histograms will not be added to the 
    // root output file.  
    
    histoName = "hsigma2_u_sensor"+to_string( ipl );
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 1, 0, 1); 
    _histoMap[ histoName  ]->StatOverflows(); 	
    _histoMap[ histoName  ]->SetDirectory(0);
    
    histoName = "hsigma2_v_sensor"+to_string( ipl );  
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 1, 0, 1);
    _histoMap[ histoName  ]->StatOverflows(); 	
    _histoMap[ histoName  ]->SetDirectory(0);
    
    histoName = "hsigma2_tu_sensor"+to_string( ipl );
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 1, 0, 1);
    _histoMap[ histoName  ]->StatOverflows(); 	
    _histoMap[ histoName  ]->SetDirectory(0);    

    histoName = "hsigma2_tv_sensor"+to_string( ipl );
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 1, 0, 1);  
    _histoMap[ histoName  ]->StatOverflows(); 	    
    _histoMap[ histoName  ]->SetDirectory(0);
    
    histoName = "hsigma2_qp_sensor"+to_string( ipl );
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 1, 0, 1);
    _histoMap[ histoName  ]->StatOverflows(); 	
    _histoMap[ histoName  ]->SetDirectory(0);     
    
    // Plot residuals U/V 
    
    histoName = "hresU_sensor"+to_string( ipl );
    max = 5*safetyFactor*Sensor.GetPitchU(); 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 500, -max, +max);
    _histoMap[ histoName ]->SetXTitle("u residual [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 
    
    histoName = "hresV_sensor"+to_string( ipl );
    max = 5*safetyFactor*Sensor.GetPitchV(); 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 500, -max, +max);
    _histoMap[ histoName ]->SetXTitle("v residual [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 
    
    // Plot residuals pulls  

    histoName = "hpull_resU_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, -4, +4);
    _histoMap[ histoName ]->SetXTitle(" pull u residual"); 
    _histoMap[ histoName ]->SetYTitle(" tracks");    
    
    histoName = "hpull_resV_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, -4, +4);
    _histoMap[ histoName ]->SetXTitle(" pull v residual"); 
    _histoMap[ histoName ]->SetYTitle(" tracks");
           
    // Plot residual profiles 
    
    histoName = "hduvsu_sensor"+to_string( ipl );
    max = safetyFactor*Sensor.GetSensitiveSizeU()/2; 
    nbins = 100; 
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", nbins, -max, +max);
    _profileMap[ histoName ]->SetXTitle("u [mm]"); 
    _profileMap[ histoName ]->SetYTitle("mean residual u [mm]");
      
    histoName = "hdvvsv_sensor"+to_string( ipl );
    max = safetyFactor*Sensor.GetSensitiveSizeV()/2;  
    nbins = 100;  
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", nbins, -max, +max); 
    _profileMap[ histoName ]->SetXTitle("v [mm]"); 
    _profileMap[ histoName ]->SetYTitle("mean residual v [mm]");
    
    histoName = "hduvsv_sensor"+to_string( ipl );
    max = safetyFactor*Sensor.GetSensitiveSizeV()/2; 
    nbins = 100;  
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", nbins, -max, +max); 
    _profileMap[ histoName ]->SetXTitle("v [mm]"); 
    _profileMap[ histoName ]->SetYTitle("mean residual u [mm]");    
       
    histoName = "hdvvsu_sensor"+to_string( ipl );
    max = safetyFactor*Sensor.GetSensitiveSizeU()/2;  
    nbins = 100; 
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", nbins, -max, +max); 
    _profileMap[ histoName ]->SetXTitle("u [mm]"); 
    _profileMap[ histoName ]->SetYTitle("mean residual v [mm]");
          
    histoName = "hduvsthetau_sensor"+to_string( ipl );
    max = 0;  // max track slope [rad]
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", 100, -max, +max); 
    _profileMap[ histoName ]->SetXTitle("du/dw [rad]"); 
    _profileMap[ histoName ]->SetYTitle("mean residual u [mm]");    
    
    histoName = "hdvvsthetav_sensor"+to_string( ipl );
    max = 0;  // max track slope [rad]
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", 100, -max, +max); 
    _profileMap[ histoName ]->SetXTitle("dv/dw [rad]"); 
    _profileMap[ histoName ]->SetYTitle("mean residual v [mm]");    

    // Plot beam profiles
    TDirectory *beamDir = sensDir->mkdir("BeamProfile");
    beamDir->cd(); 
    
    histoName = "htrk_mom_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, 0, 0); 
    _histoMap[ histoName ]->SetXTitle("track momentum [GeV]");
    _histoMap[ histoName ]->SetYTitle("tracks"); 
    
    histoName = "htrk_tu_sensor"+to_string( ipl ); 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, 0, 0); 
    _histoMap[ histoName ]->SetXTitle("fit du/dw [rad]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 
    
    histoName = "htrk_tv_sensor"+to_string( ipl ); 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, 0, 0); 
    _histoMap[ histoName ]->SetXTitle("fit dv/dw [rad]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 
      
    double maxU = safetyFactor*Sensor.GetSensitiveSizeU()/2;  
    double maxV = safetyFactor*Sensor.GetSensitiveSizeV()/2;  
    nbins = 100;           
    
    histoName = "htrk_u_sensor"+to_string( ipl ); 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", nbins, -maxU, +maxU); 
    _histoMap[ histoName ]->SetXTitle("fit u [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 

    histoName = "htrk_v_sensor"+to_string( ipl ); 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", nbins, -maxV, +maxV);
    _histoMap[ histoName ]->SetXTitle("fit v [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 

    histoName = "hhitmap_sensor"+to_string( ipl );
    _histoMap2D[ histoName] = new TH2D(histoName.c_str(), "" ,nbins, -maxU, +maxU, nbins, -maxV, +maxV);
    _histoMap2D[histoName]->SetXTitle("fit u [mm]"); 
    _histoMap2D[histoName]->SetYTitle("fit v [mm]");   
    _histoMap2D[histoName]->SetZTitle("tracks"); 
    _histoMap2D[histoName]->SetStats( false );
    
    histoName = "htrk_dudw_vs_u_sensor"+to_string( ipl );
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", nbins, -maxU, +maxU); 
    _profileMap[ histoName ]->SetYTitle("fit du/dw [rad]"); 
    _profileMap[ histoName ]->SetXTitle("fit u [mm]");    
     
    histoName = "htrk_dvdw_vs_u_sensor"+to_string( ipl );
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", nbins, -maxU, +maxU); 
    _profileMap[ histoName ]->SetYTitle("fit dv/dw [rad]"); 
    _profileMap[ histoName ]->SetXTitle("fit u [mm]");   
    
    histoName = "htrk_dudw_vs_v_sensor"+to_string( ipl );
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", nbins, -maxV, +maxV); 
    _profileMap[ histoName ]->SetYTitle("fit du/dw [rad]"); 
    _profileMap[ histoName ]->SetXTitle("fit v [mm]");    
    
    histoName = "htrk_dvdw_vs_v_sensor"+to_string( ipl );
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", nbins, -maxV, +maxV); 
    _profileMap[ histoName ]->SetYTitle("dv/dw [rad]"); 
    _profileMap[ histoName ]->SetXTitle("v [mm]");  
    
    histoName = "htrk_mom_vs_u_sensor"+to_string( ipl );
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "",  nbins, -maxU, +maxU); 
    _profileMap[ histoName ]->SetYTitle("momentum [GeV]"); 
    _profileMap[ histoName ]->SetXTitle("u [mm]");  
    
    histoName = "htrk_mom_vs_v_sensor"+to_string( ipl );
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", nbins, -maxV, +maxV); 
    _profileMap[ histoName ]->SetYTitle("momentum [GeV]"); 
    _profileMap[ histoName ]->SetXTitle("v [mm]");      
    
    // Create empty directory to store clusterID specific histos
    _clusterDirMap[ipl] = sensDir->mkdir("ClusterShapes");
    
    // Change current directory to root
    _rootFile->cd("");
  }
}

} // Namespace
