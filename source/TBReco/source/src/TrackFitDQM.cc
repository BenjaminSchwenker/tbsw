// TrackFitDQM implementation file
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// Local includes 
#include "TrackFitDQM.h"

// DEPFETTrackTools includes
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
  _description = "TrackFitDQM: DQM plots for monitor quality of track fitting";
   

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
      double trk_k = TE.GetState().GetPars()[4][0];   // 1/GeV
         
      double trk_charge = track.GetCharge();
      double trk_mom = std::abs(trk_charge/trk_k); 
       
      // Fill track parameter errors
      
      histoName = "hsigma_tu_sensor"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TMath::Sqrt( TE.GetState().GetCov()[0][0]) ) ; 
      
      histoName = "hsigma_tv_sensor"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TMath::Sqrt(TE.GetState().GetCov()[1][1]) ) ;  
      
      histoName = "hsigma_u_sensor"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TMath::Sqrt( TE.GetState().GetCov()[2][2]) ); 
      
      histoName = "hsigma_v_sensor"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( TMath::Sqrt( TE.GetState().GetCov()[3][3]) );

      histoName = "hsigma_mom_sensor"+to_string( ipl );
      _histoMap[ histoName  ]->Fill( trk_mom*trk_mom*TMath::Sqrt( TE.GetState().GetCov()[4][4] ) );

      // Fill track parameter histos 
      
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
              
      double hit_u = TE.GetHit().GetCoord()[0][0]; 
      double hit_v = TE.GetHit().GetCoord()[1][0]; 
       
      double hit_sigma_u = TMath::Sqrt(TE.GetHit().GetCov()[0][0]); 
      double hit_sigma_v = TMath::Sqrt(TE.GetHit().GetCov()[1][1]); 
      
      double pull_u = du / TMath::Sqrt( TE.GetState().GetCov()[2][2] + TE.GetHit().GetCov()[0][0] ) ; 
      double pull_v = dv / TMath::Sqrt( TE.GetState().GetCov()[3][3] + TE.GetHit().GetCov()[1][1] ) ;  
      
      cout << "benni get cluster " << endl; 

      PixelCluster Cluster = TE.GetHit().GetCluster();

      cout << "not seen benni get cluster " << endl; 

      histoName = "hcls_charge_sensor"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(Cluster.getCharge()); 
      
      histoName = "hseed_charge_sensor"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(Cluster.getSeedCharge()); 
      
      histoName = "hcls_type_sensor"+to_string( ipl );
      _histoMap[ histoName  ]->Fill(Cluster.getClusterType()); 
  
      histoName = "hsize_sensor"+to_string( ipl );   
      _histoMap[ histoName  ]->Fill(Cluster.getSize()); 
       
      histoName = "hsizeU_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill(Cluster.getUSize()); 
      
      histoName = "hsizeV_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill(Cluster.getVSize()); 
       
      if ( Cluster.getUSize() == 1 )  {
        histoName = "hpull_resU1_sensor"+to_string( ipl );
        _histoMap[ histoName ]->Fill( pull_u );  
      }  
        
      if ( Cluster.getUSize() == 2 )  {
        histoName = "hpull_resU2_sensor"+to_string( ipl );
        _histoMap[ histoName ]->Fill( pull_u );  
      }  
      
      if ( Cluster.getVSize() == 1 )  {
        histoName = "hpull_resV1_sensor"+to_string( ipl );
        _histoMap[ histoName ]->Fill( pull_v );  
      }  
      
      if ( Cluster.getVSize() == 2 )  {
        histoName = "hpull_resV2_sensor"+to_string( ipl );
        _histoMap[ histoName ]->Fill( pull_v );  
      }  
      
      histoName = "hsigma_clu_u_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill(hit_sigma_u); 
      
      histoName = "hsigma_clu_v_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill(hit_sigma_v); 

      histoName = "hhit_u_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill(hit_u); 
      
      histoName = "hhit_v_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill(hit_v);       
        
      histoName = "hresU_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill( du ); 
       
      histoName = "hresV_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill( dv ); 
      
      histoName = "hpull_resU_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill( pull_u );
      
      histoName = "hpull_resV_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill( pull_v );
     
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
    
    //
    // Beam profile histograms 

    // FIXME extrapolate track paramters to Z=0 frame
    HepMatrix p0 = track.GetTE(0).GetState().GetPars();
    double dxdz = p0[0][0];
    double dydz = p0[1][0];
    double x = p0[2][0];
    double y = p0[3][0];
    double kappa =  p0[4][0];
    
    double mom = std::abs(track.GetCharge()/kappa); 

    histoName = "hbeam_intensity"; 
    _histoMap2D[ histoName  ]->Fill(x,y); 

    histoName = "hbeam_dxdz";
    _histoMap[ histoName  ]->Fill(dxdz);
     
    histoName = "hbeam_dydz"; 
    _histoMap[ histoName  ]->Fill(dydz);
    
    histoName = "hbeam_dxdz_vs_x";
    _profileMap[ histoName ]->Fill( x, dxdz );
    
    histoName = "hbeam_dydz_vs_x";
    _profileMap[ histoName ]->Fill( x, dydz );
    
    histoName = "hbeam_dxdz_vs_y";
    _profileMap[ histoName ]->Fill( y, dxdz );
    
    histoName = "hbeam_dydz_vs_y";
    _profileMap[ histoName ]->Fill( y, dydz );   

    histoName = "hbeam_mom_vs_x";
    _profileMap[ histoName ]->Fill( x, mom );
    
    histoName = "hbeam_mom_vs_y";
    _profileMap[ histoName ]->Fill( y, mom );   
                       
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
  
  _histoMap["hmom"] = new TH1D("hmom", "", 10000, 0, 10); 
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
  
  //
  // Sensor level histograms 
  
  // Get number of sensors
  int nSens = _detector.GetNSensors();
  
  // Create subdirs for sensors
  std::string dirName; 
  std::string histoName;
  for (int ipl=0 ; ipl < nSens; ipl++) {
    std::string dirName = "Sensor"+to_string( ipl );
    _rootFile->mkdir(dirName.c_str());     
  }      

  // Loop over all sensors
  for (int ipl=0 ; ipl < nSens; ipl++) {
    
    dirName = "/Sensor"+to_string(ipl)+"/";
    _rootFile->cd(dirName.c_str());
    
    
    double max; 
    int nbins; 
    double safetyFactor = 1.1;
      
    // Get handle to sensor data
    Det & Sensor = _detector.GetDet(ipl); 
    
    // Plot tracking errors

    histoName = "hsigma_u_sensor"+to_string( ipl );
    max = 10*Sensor.GetResolutionU();
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 800, 0, max); 
    _histoMap[ histoName  ]->SetXTitle("fit sigma u [mm]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks");
    
    histoName = "hsigma_v_sensor"+to_string( ipl ); 
    max = 10*Sensor.GetResolutionV();
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 800, 0, max);
    _histoMap[ histoName  ]->SetXTitle("fit sigma v [mm]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
    
    histoName = "hsigma_tu_sensor"+to_string( ipl );
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 4000, 0, 0.01);
    _histoMap[ histoName  ]->SetXTitle("fit sigma du/dw [rad]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks");  
    
    histoName = "hsigma_tv_sensor"+to_string( ipl );
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 4000, 0, 0.01);  
    _histoMap[ histoName  ]->SetXTitle("fit sigma dv/dw [rad]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks");  
           
    histoName = "hsigma_clu_u_sensor"+to_string( ipl ); 
    max = 10*Sensor.GetResolutionU();
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 800, 0, max); 
    _histoMap[ histoName  ]->SetXTitle("cluster sigma u [mm]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks");
    
    histoName = "hsigma_clu_v_sensor"+to_string( ipl );
    max = 10*Sensor.GetResolutionV();
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 800, 0, max);
    _histoMap[ histoName  ]->SetXTitle("cluster sigma v [mm]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 
    
    histoName = "hsigma_mom_sensor"+to_string( ipl );
    _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 5000, 0, 10);
    _histoMap[ histoName  ]->SetXTitle("sigma momentum [GeV]"); 
    _histoMap[ histoName  ]->SetYTitle("tracks"); 

    // Local track parameters
    histoName = "htrk_mom_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 10000, 0, 10); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 
    _histoMap[ histoName ]->SetXTitle("track momentum [GeV]");

    histoName = "htrk_tu_sensor"+to_string( ipl ); 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100000, -1, 1); 
    _histoMap[ histoName ]->SetXTitle("fit du/dw [rad]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 

    histoName = "htrk_tv_sensor"+to_string( ipl ); 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100000, -1, 1); 
    _histoMap[ histoName ]->SetXTitle("fit dv/dw [rad]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 
      
    double  uBox = safetyFactor * 0.5 * Sensor.GetModuleBoxSizeU();   
    int uBins = Sensor.GetNColumns();            
    
    histoName = "htrk_u_sensor"+to_string( ipl ); 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", uBins, -uBox, +uBox); 
    _histoMap[ histoName ]->SetXTitle("fit u [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 

    histoName = "hhit_u_sensor"+to_string( ipl ); 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", uBins, -uBox, +uBox); 
    _histoMap[ histoName ]->SetXTitle("cluster u [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 
    
    double  vBox = safetyFactor * 0.5 * Sensor.GetModuleBoxSizeV();
    int vBins = Sensor.GetNRows(); 
    
    histoName = "htrk_v_sensor"+to_string( ipl ); 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", vBins, -vBox, +vBox);
    _histoMap[ histoName ]->SetXTitle("fit v [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 

    histoName = "hhit_v_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", vBins, -vBox, +vBox);
    _histoMap[ histoName ]->SetXTitle("cluster v [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks");     

    histoName = "hhitmap_sensor"+to_string( ipl );
    _histoMap2D[ histoName] = new TH2D(histoName.c_str(), "" ,uBins, -uBox, +uBox, vBins, -vBox, +vBox);
    _histoMap2D[histoName]->SetXTitle("fit u [mm]"); 
    _histoMap2D[histoName]->SetYTitle("fit v [mm]");    
    _histoMap2D[histoName]->SetStats( false );  
    
    // Plot residuals U/V 
    
    histoName = "hresU_sensor"+to_string( ipl );
    max = 100*safetyFactor*Sensor.GetPitchU(); 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 10000, -max, +max);
    _histoMap[ histoName ]->SetXTitle("u residual [mm]"); 
    _histoMap[ histoName ]->SetYTitle("tracks"); 
    
    histoName = "hresV_sensor"+to_string( ipl );
    max = 100*safetyFactor*Sensor.GetPitchV(); 
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 10000, -max, +max);
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

    // Plot cluster shape   

    histoName = "hcls_charge_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, 0, 100);
    _histoMap[ histoName ]->SetXTitle(" cluster charge [ADU]"); 
    _histoMap[ histoName ]->SetYTitle(" tracks");  

    histoName = "hseed_charge_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, 0, 100);
    _histoMap[ histoName ]->SetXTitle(" seed charge [ADU]"); 
    _histoMap[ histoName ]->SetYTitle(" tracks");      

    histoName = "hcls_type_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 10, 0, 10);
    _histoMap[ histoName ]->SetXTitle(" cluster type"); 
    _histoMap[ histoName ]->SetYTitle(" tracks");   

    histoName = "hsize_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 10, 0, 10);
    _histoMap[ histoName ]->SetXTitle(" cluster size [pixels]"); 
    _histoMap[ histoName ]->SetYTitle(" tracks");    
    
    histoName = "hsizeU_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 6, 0, 6);
    _histoMap[ histoName ]->SetXTitle(" cluster size u [cells]"); 
    _histoMap[ histoName ]->SetYTitle(" tracks");    
    
    histoName = "hsizeV_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 6, 0, 6);
    _histoMap[ histoName ]->SetXTitle(" cluster size v [cells]"); 
    _histoMap[ histoName ]->SetYTitle(" tracks");  

    histoName = "hpull_resU1_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, -4, +4);
    _histoMap[ histoName ]->SetXTitle(" pull u residual (one cell clusters)"); 
    _histoMap[ histoName ]->SetYTitle(" tracks");    
    
    histoName = "hpull_resU2_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, -4, +4);
    _histoMap[ histoName ]->SetXTitle(" pull u residual (two cell clusters)"); 
    _histoMap[ histoName ]->SetYTitle(" tracks");    

    histoName = "hpull_resV1_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, -4, +4);
    _histoMap[ histoName ]->SetXTitle(" pull v residual (one cell clusters)"); 
    _histoMap[ histoName ]->SetYTitle(" tracks");    
    
    histoName = "hpull_resV2_sensor"+to_string( ipl );
    _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 100, -4, +4);
    _histoMap[ histoName ]->SetXTitle(" pull v residual (two cell clusters)"); 
    _histoMap[ histoName ]->SetYTitle(" tracks");    
         
    // Plot residual profiles 
    
    histoName = "hduvsu_sensor"+to_string( ipl );
    max = safetyFactor*Sensor.GetSensitiveSizeU()/2; 
    nbins = Sensor.GetNColumns()/4;
    if ( nbins > 50 ) nbins = 50; 
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", nbins, -max, +max);
    _profileMap[ histoName ]->SetXTitle("u [mm]"); 
    _profileMap[ histoName ]->SetYTitle("mean residual u [mm]");
      
    histoName = "hdvvsv_sensor"+to_string( ipl );
    max = safetyFactor*Sensor.GetSensitiveSizeV()/2;  
    nbins = Sensor.GetNRows()/4;
    if ( nbins > 50 ) nbins = 50;  
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", nbins, -max, +max); 
    _profileMap[ histoName ]->SetXTitle("v [mm]"); 
    _profileMap[ histoName ]->SetYTitle("mean residual v [mm]");
    
    histoName = "hduvsv_sensor"+to_string( ipl );
    max = safetyFactor*Sensor.GetSensitiveSizeV()/2; 
    nbins = Sensor.GetNRows()/4;
    if ( nbins > 50 ) nbins = 50;   
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", nbins, -max, +max); 
    _profileMap[ histoName ]->SetXTitle("v [mm]"); 
    _profileMap[ histoName ]->SetYTitle("mean residual u [mm]");    
       
    histoName = "hdvvsu_sensor"+to_string( ipl );
    max = safetyFactor*Sensor.GetSensitiveSizeU()/2;  
    nbins = Sensor.GetNColumns()/4;
    if ( nbins > 50 ) nbins = 50;  
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", nbins, -max, +max); 
    _profileMap[ histoName ]->SetXTitle("u [mm]"); 
    _profileMap[ histoName ]->SetYTitle("mean residual v [mm]");
          
    histoName = "hduvsthetau_sensor"+to_string( ipl );
    max = 2;  // max track slope [rad]
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", 50000, -max, +max); 
    _profileMap[ histoName ]->SetXTitle("du/dw [rad]"); 
    _profileMap[ histoName ]->SetYTitle("mean residual u [mm]");    
    
    histoName = "hdvvsthetav_sensor"+to_string( ipl );
    max = 2;  // max track slope [rad]
    _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", 50000, -max, +max); 
    _profileMap[ histoName ]->SetXTitle("dv/dw [rad]"); 
    _profileMap[ histoName ]->SetYTitle("mean residual v [mm]");    

  }
  
  //
  // Beam profile histograms 
  _rootFile->mkdir("BeamProfile");           
  _rootFile->cd("BeamProfile");
  
  histoName = "hbeam_intensity";
  _histoMap2D[ histoName] = new TH2D(histoName.c_str(), "" ,576, -30, +30, 288, -15, +15);
  _histoMap2D[histoName]->SetXTitle("x [mm]"); 
  _histoMap2D[histoName]->SetYTitle("y [mm]"); 
  _histoMap2D[histoName]->SetZTitle("tracks");     
  _histoMap2D[histoName]->SetStats( false );  
  
  histoName = "hbeam_dxdz"; 
  _histoMap[ histoName ] = new TH1D(histoName.c_str(), "",  20000, -0.1, 0.1); 
  _histoMap[ histoName ]->SetXTitle("fit dx/dz [rad]"); 
  _histoMap[ histoName ]->SetYTitle("tracks");  

  histoName = "hbeam_dydz"; 
  _histoMap[ histoName ] = new TH1D(histoName.c_str(), "",  20000, -0.1, 0.1); 
  _histoMap[ histoName ]->SetXTitle("fit dy/dz [rad]"); 
  _histoMap[ histoName ]->SetYTitle("tracks");  
  
  histoName = "hbeam_dxdz_vs_x";
  _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", 576, -30, +30); 
  _profileMap[ histoName ]->SetYTitle("dx/dz [rad]"); 
  _profileMap[ histoName ]->SetXTitle("x [mm]");    

  histoName = "hbeam_dydz_vs_x";
  _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", 576, -30, +30); 
  _profileMap[ histoName ]->SetYTitle("dy/dz [rad]"); 
  _profileMap[ histoName ]->SetXTitle("x [mm]");   
    
  histoName = "hbeam_dxdz_vs_y";
  _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", 288, -15, +15); 
  _profileMap[ histoName ]->SetYTitle("dx/dz [rad]"); 
  _profileMap[ histoName ]->SetXTitle("y [mm]");    

  histoName = "hbeam_dydz_vs_y";
  _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", 288, -15, +15); 
  _profileMap[ histoName ]->SetYTitle("dy/dz [rad]"); 
  _profileMap[ histoName ]->SetXTitle("y [mm]");  

  histoName = "hbeam_mom_vs_x";
  _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", 576, -30, +30); 
  _profileMap[ histoName ]->SetYTitle("momentum [GeV]"); 
  _profileMap[ histoName ]->SetXTitle("x [mm]");  
    
  histoName = "hbeam_mom_vs_y";
  _profileMap[ histoName ] = new TProfile(histoName.c_str(), "", 288, -15, +15); 
  _profileMap[ histoName ]->SetYTitle("momentum [GeV]"); 
  _profileMap[ histoName ]->SetXTitle("y [mm]");          
  
  
}

} // Namespace
