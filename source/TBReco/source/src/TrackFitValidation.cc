// TrackFitValidation Processor  
// 
// See TrackFitValidation.h for full documentation of processor. 
// 
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


#include "TrackFitValidation.h"

// TBTools includes
#include "TBTrack.h"
#include "GenericTrackFitter.h"
#include "TrackInputProvider.h"
#include "PixelCluster.h"
#include "Det.h"
#include "Utilities.h" 

// Include basic C
#include <iostream>
#include <limits>
#include <iomanip>
#include <cmath>
#include <algorithm>

// Include LCIO classes
#include "lcio.h"
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>

// Include CLHEP classes
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Vector/ThreeVector.h>

// Used namespaces
using namespace std; 
using namespace CLHEP; 
using namespace lcio ;
using namespace marlin ;

namespace depfet {

  //
  // Instantiate this object
  //
  TrackFitValidation aTrackFitValidation ;
  
  //
  // Constructor
  //
  TrackFitValidation::TrackFitValidation() : Processor("TrackFitValidation")
  {
    
    // Processor description
    _description = "Produces Validation plots (pulls) for track fitting algorithms" ;
    
    //
    // Input collections 
    
    registerInputCollection( LCIO::SIMTRACKERHIT,
                             "SimTrackerHitCollection" ,
                             "Name of collection with simulated hits"  ,
                             _simHitColName ,
                             std::string("SimTrackerHits") ) ;
    
    registerInputCollection( LCIO::TRACK,
                           "TrackCollection" ,
                           "Name of collection with reco tracks "  ,
                           _trackColName ,
                           std::string("tracks") ) ;
    
    // Processor parameters:
    
    registerProcessorParameter ("AlignmentDBFileName",
                                "This is the name of the LCIO file with the alignment constants (add .slcio)",
                                _alignmentDBFileName, static_cast< string > ( "alignmentDB.slcio" ) );     
       
    registerProcessorParameter ("MaxResidualU",
                                "Maximum u residual for matching simHits to hits [mm]. Put -1 to deactivate cut.",
                                _maxResidualU,  static_cast < double > (0.2));
    
    registerProcessorParameter ("MaxResidualV",
                                "Maximum v residual for matching simHits to hits [mm]. Put -1 to deactivate cut.",
                                _maxResidualV,  static_cast < double > (0.2));
    
    registerProcessorParameter( "RootFileName",
                                "Output root file name",
                                _rootFileName, std::string("Validation.root"));  
    
    
  }
  
  //
  // Method called at the beginning of data processing
  //
  void TrackFitValidation::init() {
   
    // Initialize variables
    _nRun = 0 ;
    _nEvt = 0 ;
    _timeCPU = clock()/1000 ;
    
    // Print set parameters
    printProcessorParams();
    
    // Read detector constants from gear file
    _detector.ReadGearConfiguration();  
    
    // Read alignment data base file 
    _detector.ReadAlignmentDB( _alignmentDBFileName );    
    
    bookHistos(); 
  }
  
  //
  // Method called for each run
  //
  void TrackFitValidation::processRunHeader(LCRunHeader * run)
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
  void TrackFitValidation::processEvent(LCEvent * evt)
  {
    
    _nEvt ++ ;
     
    //
    // Get reco track collection
    //
     
    LCCollection* trackCol = NULL;
    int nTrack = 0;   
    try {
      trackCol = evt->getCollection( _trackColName ) ;
      nTrack = trackCol->getNumberOfElements();
    } catch (lcio::DataNotAvailableException& e) {
      streamlog_out(MESSAGE2) << "Not able to get collection "
                              << _trackColName
                              << " from event " << evt->getEventNumber()
                              << " in run " << evt->getRunNumber()  << endl << endl;   
    }  

    streamlog_out(MESSAGE2) << "Total of " << nTrack  << " track(s) in collection " << _trackColName << endl;
    
    
    // Read fitted reco tracks 
    // ----------------------------------
    std::vector<TBTrack> TrackStore; 

    // Configure Kalman track fitter
    GenericTrackFitter TrackFitter(_detector);
    TrackFitter.SetNumIterations(2); 
    
    TrackInputProvider TrackLCIOReader;  
    
    for(int itrk=0; itrk< nTrack ; itrk++)
    {  
      // Retrieve track from LCIO 
      Track * lciotrk = dynamic_cast<Track*> (trackCol->getElementAt(itrk));
      
      // Convert LCIO -> TB track  
      TBTrack trk = TrackLCIOReader.MakeTBTrack( lciotrk, _detector );  
      
      // Refit track in nominal alignment
      bool trkerr = TrackFitter.Fit(trk);
      if ( trkerr ) {
        streamlog_out ( MESSAGE3 ) << "Fit failed. Skipping track!" << endl;
        continue;
      } 
        
      TrackStore.push_back(trk);                                 
    } 
    
    //
    // Get simhit collection
    //
      
    LCCollection* simHitCol = NULL;
    int nSimHit = 0;   
    try {
      simHitCol = evt->getCollection( _simHitColName ) ;
      nSimHit = simHitCol->getNumberOfElements();
    } catch (lcio::DataNotAvailableException& e) {
      streamlog_out(MESSAGE2) << "Not able to get collection "
                              << _simHitColName
                              << " from event " << evt->getEventNumber()
                              << " in run " << evt->getRunNumber()  << endl << endl;   
    }  
      
    streamlog_out(MESSAGE2) << "Total of " << nSimHit << " simHits in collection " << _simHitColName << endl;  
    
    //
    // Fill track level histos 
    //
     
    for(int i=0;i<(int)TrackStore.size(); ++i) 
    {
      TBTrack& recoTrack = TrackStore[i];
      _histoMap["nhits"]->Fill(recoTrack.GetNumHits());    
      _histoMap["ndf"]->Fill(recoTrack.GetNDF());   
      _histoMap["chi2"]->Fill(recoTrack.GetChiSqu());
      _histoMap["chi2ndof"]->Fill(recoTrack.GetChiSqu()/recoTrack.GetNDF());
      _histoMap["chi2prob"]->Fill(TMath::Prob(recoTrack.GetChiSqu(), recoTrack.GetNDF()));
    }
    
    //
    // Main loop over all sensors (fill sensor level histograms) 
    // 
    
    for (int ipl=0 ; ipl < _detector.GetNSensors(); ipl++) {
                    
     
      
      // Read simHits for sensor 
      // -----------------------
      std::vector<SimTrackerHit*> SimHitStore;
      CellIDDecoder<SimTrackerHit> cellIDDec(simHitCol);   
       
      for(int i=0; i< nSimHit ; i++)
      {  
        // Retrieve simtrackerhit 
        SimTrackerHit * simHit = dynamic_cast<SimTrackerHit*> (simHitCol->getElementAt(i));
        
        // Set current - layer ID, ladder ID and sensor ID
        int sensorID = cellIDDec(simHit)["sensorID"];
         
        if( ipl == _detector.GetPlaneNumber(sensorID))
        {         
          streamlog_out(MESSAGE2) << " SimHit at plane " << ipl << " at: (" << simHit->getPosition()[0] << ", " << simHit->getPosition()[1] << ")" 
                                  << endl;
          
          SimHitStore.push_back(simHit);    
        }                                       
      } 
      
      // Record for each track a matched simhit
      vector<int> track2simhit(TrackStore.size(), -1);
      
      // Continue matching until all hits are matched 
      // or no hit is close enough!!
      
      { 
        double distmin=numeric_limits<double >::max();
        int bestsimhit=-1;   
        int besttrack=-1; 
      
        do{
          bestsimhit=-1;
          besttrack=-1; 
          distmin=numeric_limits<double >::max();
          
          // Find hit/simhit pair with minimum chi2 distance.  
          for(int i=0; i< (int)TrackStore.size() ; i++)
          {

            TBTrack& recoTrack = TrackStore[i];
            
            // If track is already matched, skip hit 
            if (track2simhit[i] >= 0) continue;

            // If track is not intersecting sensor, skip it
            if (!recoTrack.GetTE(ipl).IsCrossed()) continue; 
          
            for(int j=0;j<(int)SimHitStore.size(); j++)
            {
              SimTrackerHit * simHit = SimHitStore[j];
              double simHitPosU = simHit->getPosition()[0];
              double simHitPosV = simHit->getPosition()[1];
              
              
              double hitPosU = recoTrack.GetTE(ipl).GetState().GetPars()[2][0];
              double hitPosV = recoTrack.GetTE(ipl).GetState().GetPars()[3][0];
              
              // Skip all DUT hits with too large residuum 
              if ( std::abs(hitPosU-simHitPosU) >= _maxResidualU && _maxResidualU > 0 ) continue;  
              if ( std::abs(hitPosV-simHitPosV) >= _maxResidualV && _maxResidualV > 0 ) continue; 
              
              // Finally, we will use a simple 2D distance to select best matching hit
              double hitdist = 0; 
              if ( _maxResidualU > 0 )  hitdist += std::abs(hitPosU-simHitPosU); 
              if ( _maxResidualV > 0 )  hitdist += std::abs(hitPosV-simHitPosV); 
          
              if( hitdist<distmin )
              {
                distmin=hitdist;
                besttrack=i;
                bestsimhit=j;
              }
            }
          }
        
          streamlog_out(MESSAGE2) << "In matching loop: best track " << besttrack << " to simhit " << bestsimhit << endl; 
          streamlog_out(MESSAGE2) << "  distmin: " <<  distmin  << endl; 
          
          // Check if a match was found
          if( bestsimhit>-1 &&  besttrack>-1   )
          {   
            streamlog_out(MESSAGE2) << "  match found!!!"   << endl;
            track2simhit[besttrack] = bestsimhit;
          } 
            
        } // End of do loop of matching DUT hits to fitted positions
        while( bestsimhit>-1 &&  besttrack>-1);
      }
      
      streamlog_out(MESSAGE2) << "Start fill sensor histos" << endl; 
      
      for(int i=0;i<(int)TrackStore.size(); ++i)
      {
        if ( track2simhit[i] >= 0 ) {        
        
          std::string histoName;
          
          SimTrackerHit * simHit = SimHitStore[ track2simhit[i] ]; 
          Hep3Vector momentum(simHit->getMomentum()[0],simHit->getMomentum()[1],simHit->getMomentum()[2]);
          
          // Compute a local track state from simHit
          HepMatrix sim_x(5,1,0);
          sim_x[0][0] = momentum[0]/momentum[2];                              // du/dw   
          sim_x[1][0] = momentum[1]/momentum[2];                              // dv/dw      
          sim_x[2][0] = simHit->getPosition()[0];                             // u
          sim_x[3][0] = simHit->getPosition()[1];                             // v
          sim_x[4][0] = simHit->getMCParticle()->getCharge()/momentum.mag();  // q/p

          histoName = "htrk_dir_truth_det"+to_string( ipl ); 
          _histoMap2D[ histoName  ]->Fill(sim_x[0][0],sim_x[1][0]);  
          
          histoName = "htrk_tu_truth_det"+to_string( ipl );
          _histoMap[ histoName  ]->Fill(sim_x[0][0]);  
          
          histoName = "htrk_tv_truth_det"+to_string( ipl );
          _histoMap[ histoName  ]->Fill(sim_x[1][0]);  
          
          histoName = "htrk_u_truth_det"+to_string( ipl );
          _histoMap[ histoName  ]->Fill(sim_x[2][0]);  
          
          histoName = "htrk_v_truth_det"+to_string( ipl );
          _histoMap[ histoName  ]->Fill(sim_x[3][0]);  
          
          histoName = "htrk_mom_truth_det"+to_string( ipl );
          _histoMap[ histoName  ]->Fill(momentum.mag());
          
          
          TBTrack& recoTrack = TrackStore[i];
           
          HepMatrix rec_x = recoTrack.GetTE(ipl).GetState().GetPars();
          HepSymMatrix& rec_cov = recoTrack.GetTE(ipl).GetState().GetCov();  
          
          double rec_charge = recoTrack.GetCharge();       
          double rec_mom = std::abs(rec_charge/rec_x[4][0]);   

          // Fill reco track parameters
             
          histoName = "htrk_tu_det"+to_string( ipl );
          _histoMap[ histoName  ]->Fill(rec_x[0][0]);  
          
          histoName = "htrk_tv_det"+to_string( ipl );
          _histoMap[ histoName  ]->Fill(rec_x[1][0]); 
          
          histoName = "htrk_u_det"+to_string( ipl );
          _histoMap[ histoName  ]->Fill(rec_x[2][0]); 
          
          histoName = "htrk_v_det"+to_string( ipl );
          _histoMap[ histoName  ]->Fill(rec_x[3][0]); 
          
          histoName = "htrk_mom_det"+to_string( ipl );
          _histoMap[ histoName  ]->Fill( rec_mom ); 
          
          // Fill track parameter pulls
                     
          HepMatrix diff = rec_x-sim_x;
          
          // Case w/o magnetic field is spacial
          if ( ( std::abs(_detector.GetBx()) + std::abs(_detector.GetBy()) + std::abs(_detector.GetBz()) )  == 0 )  {
             
            int ierr; 
            HepMatrix jchisq = diff.sub(1,4,1,1).T()*rec_cov.sub(1,4).inverse(ierr)*diff.sub(1,4,1,1);
                        
            histoName = "hJ_det"+to_string( ipl );
            _histoMap[histoName]->Fill(jchisq[0][0]);  
            
            histoName = "hJp_det"+to_string( ipl );
            _histoMap[histoName]->Fill(TMath::Prob(jchisq[0][0], 4));
             
          } else {
            
            int ierr; 
            HepMatrix jchisq = diff.T()*rec_cov.inverse(ierr)*diff;
            
            histoName = "hJ_det"+to_string( ipl );
            _histoMap[histoName]->Fill(jchisq[0][0]);  
            
            histoName = "hJp_det"+to_string( ipl );
            _histoMap[histoName]->Fill(TMath::Prob(jchisq[0][0], 5));
            
          }
           
          
          histoName = "hp1_det"+to_string( ipl );
          _histoMap[ histoName ]->Fill(diff[0][0]/TMath::Sqrt(rec_cov[0][0]));
          
          histoName = "hp2_det"+to_string( ipl );
          _histoMap[ histoName ]->Fill(diff[1][0]/TMath::Sqrt(rec_cov[1][1]));
          
          histoName = "hp3_det"+to_string( ipl );
          _histoMap[ histoName ]->Fill(diff[2][0]/TMath::Sqrt(rec_cov[2][2]));
          
          histoName = "hp4_det"+to_string( ipl );
          _histoMap[ histoName ]->Fill(diff[3][0]/TMath::Sqrt(rec_cov[3][3]));  
          
          histoName = "hp5_det"+to_string( ipl );
          _histoMap[ histoName ]->Fill(diff[4][0]/TMath::Sqrt(rec_cov[4][4]));   
          
          // Fill track parameter errors
          
          histoName = "hsigma_tu_det"+to_string( ipl );
          _histoMap[ histoName  ]->Fill( TMath::Sqrt(rec_cov[0][0]) ) ; 
           
          histoName = "hsigma_tv_det"+to_string( ipl );
          _histoMap[ histoName  ]->Fill( TMath::Sqrt(rec_cov[1][1]) ) ;  
          
          histoName = "hsigma_u_det"+to_string( ipl );
          _histoMap[ histoName  ]->Fill( TMath::Sqrt( rec_cov[2][2]) ); 
          
          histoName = "hsigma_v_det"+to_string( ipl );
          _histoMap[ histoName  ]->Fill( TMath::Sqrt( rec_cov[3][3]) );
          
          histoName = "hsigma_mom_det"+to_string( ipl );
          _histoMap[ histoName  ]->Fill( rec_mom*rec_mom*TMath::Sqrt( rec_cov[4][4]) ); 
          
          // Hit variables  
          if ( recoTrack.GetTE(ipl).HasHit() ) {
          
            // Get local chi2   	     
            double  hitchi2 = recoTrack.GetTE(ipl).GetChiSqu();  
            
            // Get hit residual
            HepMatrix r = recoTrack.GetTE(ipl).GetHit().GetCoord();
            r -= recoTrack.GetTE(ipl).GetState().GetXCoord();
            
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
    return;
  }


  //
  // Method called after each event to check the data processed
  //
  void TrackFitValidation::check( LCEvent * evt ) {}

  //
  // Method called after all data processing
  //
  void TrackFitValidation::end()
  {
    
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

    // Close root  file
    _rootFile->Write();
    _rootFile->Close();
    delete _rootFile;
  }

  //
  // Method printing processor parameters
  //
  void TrackFitValidation::printProcessorParams() const 
  {
    
    streamlog_out(MESSAGE3)  << std::endl
                             << " "
                             << "TrackFitValidation Development Version, be carefull!!"
                             << " "
                             << std::endl  << std::endl;   
    
  }
  
  void TrackFitValidation::bookHistos()
  {   
     
    _rootFile = new TFile( _rootFileName.c_str(),"recreate");
    _rootFile->cd("");
     
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
    
    
    
    _histoMap["chi2ndof"] = new TH1D("hchi2ndof", "", 100, 0, 10); 
    _histoMap["chi2ndof"]->SetYTitle("tracks"); 
    _histoMap["chi2ndof"]->SetXTitle("#chi^{2}/ndof");  
    
    _histoMap["chi2prob"] = new TH1D("hchi2prob", "", 100, 0, 1);
    _histoMap["chi2prob"]->SetXTitle("p-value"); 
    _histoMap["chi2prob"]->SetYTitle("tracks");
    _histoMap["chi2prob"]->SetMinimum(0.);  
     
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
      _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, 0, 1.2*10); 
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
      _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, 0, 1.2*10); 
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
      _histoMap[ histoName  ] = new TH1D(histoName.c_str(), "", 4000, 0, 10);
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



