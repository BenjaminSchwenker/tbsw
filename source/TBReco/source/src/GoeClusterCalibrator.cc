// GoeClusterCalibrator Processor  
// 
// See GoeClusterCalibrator.h for full documentation of processor. 
// 
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


#include "GoeClusterCalibrator.h"

// TBTools includes
#include "TBTrack.h"
#include "TrackInputProvider.h"
#include "GenericTrackFitter.h"
#include "PixelCluster.h"

// Include basic C
#include <iostream>
#include <limits>
#include <iomanip>
#include <algorithm>
#include <sstream>


// Include LCIO classes
#include "lcio.h"
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <Exceptions.h>

// Include ROOT classes
#include <TMath.h>
#include <TFile.h>
#include <TVectorD.h>


// Used namespaces
using namespace std; 
using namespace CLHEP; 
using namespace lcio ;
using namespace marlin ;

namespace depfet {

  //
  // Instantiate this object
  //
  GoeClusterCalibrator aGoeClusterCalibrator ;
  
  //
  // Constructor
  //
  GoeClusterCalibrator::GoeClusterCalibrator() : Processor("GoeClusterCalibrator")
  {
    
    // Processor description
    _description = "GoeClusterCalibrator: Create clusterDB for clusters using reco tracks" ;
    
    //
    // Input collections 
    
    registerInputCollection(LCIO::TRACK,"InputTrackCollectionName",
                            "Track input collection",
                            _inputTrackCollectionName,std::string("tracks"));
    
    registerProcessorParameter( "ClusterDBFileName",
                                "Output clusterDB file name",
                                _clusterDBFileName, std::string("clusterDB.root"));  
     
    registerProcessorParameter ("MinClusters",
                                "Minimum number of cluster ID occurances for clusterDB",
                                _minClusters,  static_cast < int > (2000));
    
    std::vector<int> initSoftADCSteps;
    registerProcessorParameter ("SoftwareADC",
                                "List of steps for software ADC. An empty list gives a constant transfer curve (0bit limit)",
                                _swADCSteps, initSoftADCSteps);
    
    registerProcessorParameter ("MinVarianceU", "Minimum value of variance for u position measurement [mm^2]",
                                _minVarianceU, static_cast < float > (1.0E-6) );
    
    registerProcessorParameter ("MinVarianceV", "Minimum value of variance for v position measurement [mm^2]",
                                _minVarianceV, static_cast < float > (1.0E-6) );
     
    registerProcessorParameter ("AlignmentDBFileName",
                                "This is the name of the LCIO file with the alignment constants (add .slcio)",
                                _alignmentDBFileName, static_cast< string > ( "eudet-alignmentDB.slcio" ) ); 
    
    std::vector<int> initIgnoreIDVec;
    registerProcessorParameter ("IgnoreIDs",
                                "Ignore clusters from list of sensorIDs",
                                _ignoreIDVec, initIgnoreIDVec);
    
  }
  
  //
  // Method called at the beginning of data processing
  //
  void GoeClusterCalibrator::init() {
    
    // Initialize variables
    _nRun = 0 ;
    _nEvt = 0 ;
    _timeCPU = clock()/1000 ;
    
    if ( _minVarianceU <= 0 ) {  
      _minVarianceU = 1.0E-6;  
    } 
    
    if ( _minVarianceV <= 0 ) {
      _minVarianceV = 1.0E-6;  
    } 
    
    // Make sure adc steps are sorted
    std::sort(_swADCSteps.begin(), _swADCSteps.end());
    
    // Erase duplicated entries
    _swADCSteps.erase( std::unique( _swADCSteps.begin(), _swADCSteps.end() ), _swADCSteps.end() ); 
    
    for(auto step : _swADCSteps ) {
      streamlog_out( MESSAGE2 ) << " sw adc step "  << step << endl;
    }
    
    _trackVarUMap = new TH1F("","",1,0,1);
    _trackVarUMap->SetDirectory(0);
    _trackVarUMap->StatOverflows(); 	            
    
    _trackVarVMap = new TH1F("","",1,0,1); 
    _trackVarVMap->SetDirectory(0);
    _trackVarVMap->StatOverflows();         
    
    _trackCovUVMap = new TH1F("","",1,0,1);
    _trackCovUVMap->SetDirectory(0);
    _trackCovUVMap->StatOverflows(); 	
    
    // Print set parameters
    printProcessorParams();
    
    // Read detector constants from gear file
    _detector.ReadGearConfiguration();  
    
    // Read alignment data base file 
    _detector.ReadAlignmentDB( _alignmentDBFileName );     
  }
  
  //
  // Method called for each run
  //
  void GoeClusterCalibrator::processRunHeader(LCRunHeader * run)
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
  void GoeClusterCalibrator::processEvent(LCEvent * evt)
  {
    
    _nEvt ++ ;
    
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
    for (int itrk = 0; itrk < nTracks; itrk++) {   
      
      // Retrieve track from LCIO 
      Track * inputtrack = dynamic_cast<Track*> (inputCollection->getElementAt(itrk));
      
      // Convert LCIO -> TB track  
      TBTrack track = TrackIO.MakeTBTrack( inputtrack, _detector );  
      
      // Refit track 
      bool trkerr = TrackFitter.Fit(track);
      if ( trkerr ) {
        continue;
      } 
      
      //
      // Loop over all clusters in the track
      for (int ipl= 0; ipl< _detector.GetNSensors(); ++ipl) {   
        
        
        // Get sensor data 
        //------------------------
        TBTrackElement& TE = track.GetTE(ipl);  
        Det & Sensor = _detector.GetDet(ipl);  
        int sensorID = Sensor.GetDAQID();        
        
        bool ignoreID = false;
        for (auto id :  _ignoreIDVec)  {
          if  (id == sensorID) ignoreID = true; 
        }
         
        // Ignore track elements w/o measurment
        if ( TE.HasHit() && !ignoreID ) { 
          
          // Get local track parameters 
          double trk_tu = TE.GetState().GetPars()[0][0];  // rad
          double trk_tv = TE.GetState().GetPars()[1][0];  // rad
          double trk_u = TE.GetState().GetPars()[2][0];   // mm
          double trk_v = TE.GetState().GetPars()[3][0];   // mm
          double trk_qp = TE.GetState().GetPars()[4][0];  // 1/GeV
          double trk_charge = track.GetCharge();
          double trk_mom = std::abs(trk_charge/trk_qp); 
           
          double sigma2_u = TE.GetState().GetCov()[2][2]; 
          double sigma2_v = TE.GetState().GetCov()[3][3]; 
          double cov_uv   = TE.GetState().GetCov()[2][3];  
           
          // We use temporary histograms to compute an averaged 2x2 
          // covariance matrix for all reco tracks at given sensor.
             
          _trackVarUMap->Fill(sigma2_u);
          _trackVarVMap->Fill(sigma2_v); 
          _trackCovUVMap->Fill(cov_uv);  
          
          // Get cluster id  
          PixelCluster Cluster = TE.GetHit().GetCluster();  
          string id = Cluster.getLabel(_swADCSteps); 
          
          // Register new cluster if needed
          if (_sensorMap.find(id) == _sensorMap.end() ) {
            _sensorMap[id] = 0;
            
            _clusterUMap[id] = new TH1F("","",1,0,1);
            _clusterUMap[id]->SetDirectory(0);
            _clusterUMap[id]->StatOverflows(); 	
             
            _clusterVMap[id] = new TH1F("","",1,0,1);
            _clusterVMap[id]->SetDirectory(0); 
            _clusterVMap[id]->StatOverflows(); 	
            
            _clusterUVMap[id] = new TH2F("","",1,0,1,1,0,1);
            _clusterUVMap[id]->SetDirectory(0);
            _clusterUVMap[id]->StatOverflows(); 	
          }
          
          trk_u -= Sensor.GetPixelCenterCoordU( Cluster.getVStart(), Cluster.getUStart()); 
          trk_v -= Sensor.GetPixelCenterCoordV( Cluster.getVStart(), Cluster.getUStart()); 
          
          // Count how many times a label appears 
          _sensorMap[id]++;  
          _clusterUMap[id]->Fill( trk_u ); 
          _clusterVMap[id]->Fill( trk_v );     
          _clusterUVMap[id]->Fill( trk_u, trk_v ); 
          
          
        }
      }
    }  
    
    return;
  }
  
  
  //
  // Method called after each event to check the data processed
  //
  void GoeClusterCalibrator::check( LCEvent * evt ) {}
  
  //
  // Method called after all data processing
  //
  void GoeClusterCalibrator::end()
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
     
    // We remove all cluster IDs which have to few counts and cannot be 
    // calibrated :( 
    
    streamlog_out(MESSAGE3) << "Remove all cluster IDs with less than " << _minClusters  << " counts." << endl; 
    
    // Count all clusters
    double countAll = 0;
    std::map<std::string, int>  countAllMap;  
    
    // Count reject clusters
    int countReject = 0;
    std::map<std::string, int>  countRejectMap; 
    
    // Delete cluster ids with too small counter
    for(auto iter = _sensorMap.begin(); iter != _sensorMap.end(); ) {
      auto id = iter->first; 
      auto counter = iter->second;
      auto type = getClusterType(id); 
      
      // Add counter if type is new
      if (countAllMap.find(type) == countAllMap.end() ) {
        countAllMap[type] = 0;
        countRejectMap[type] = 0;
      }
        
      countAll += counter;
      countAllMap[type] += counter;
        
      if(counter < _minClusters ) {
        streamlog_out(MESSAGE3) << "  Deleting label:  " << id << " because too few counts (" << counter  << ")" << endl;
        iter = _sensorMap.erase(iter);
        countReject += counter;
        countRejectMap[type] += counter;
      } else {
        ++iter;
      }
    }
    
    streamlog_out(MESSAGE3) << "Number of rejected clusters is: " 
                            << countReject << " (" << 100.0*countReject/countAll  << "%)"  << endl;
    
    
    streamlog_out(MESSAGE3) << "Create the clusterDB ... " << endl; 
    
    TFile * _rootFile = new TFile( _clusterDBFileName.c_str(),"recreate");
    _rootFile->cd("");
    
    
    // Book histograms for clusterDB
    int NCLUSTERS = _sensorMap.size(); 
    string histoName;  
       
    _rootFile->cd("");
    
    histoName = "hDB_Coverage";
    _histoMap[histoName] = new TH1F(histoName.c_str(),"",1,0,1);
    _histoMap[histoName]->SetStats( false );
    _histoMap[histoName]->SetYTitle("coverage [%]");
    _histoMap[histoName]->SetBinContent( 1, 100.0 - 100.0*countReject/countAll );
    _histoMap[histoName]->GetXaxis()->SetBinLabel( 1, "cluster found in clusterDB" );

    histoName = "hDB_CoverageTypes";
    _histoMap[histoName] = new TH1F(histoName.c_str(),"",countAllMap.size(),0,countAllMap.size());
    _histoMap[histoName]->SetStats( false );
    _histoMap[histoName]->SetYTitle("coverage [%]");
    _histoMap[histoName]->SetXTitle("type");
    
    auto it1 = countAllMap.begin();
    auto it2 = countRejectMap.begin();
    for (int bin = 1; bin <= countAllMap.size(); ++bin)
    {
      string type = it1->first;
      int all = it1->second;  
      int reject = it2->second;  
      
      _histoMap[histoName]->SetBinContent( bin, 100.0 - 100.0*reject/all );
      _histoMap[histoName]->GetXaxis()->SetBinLabel( bin, type.c_str() ); 
      ++it1;
      ++it2;
    }
       
    histoName = "hDB_Weight";
    _histoMap[histoName] = new TH1F(histoName.c_str(),"",NCLUSTERS,0,NCLUSTERS);
    _histoMap[histoName]->SetStats( false );
    _histoMap[histoName]->SetYTitle("label weight");  
      
    histoName = "hDB_U";
    _histoMap[histoName] = new TH1F(histoName.c_str(),"",NCLUSTERS,0,NCLUSTERS);
    _histoMap[histoName]->SetStats( false );
    _histoMap[histoName]->SetYTitle("offset u [mm]");  
      
    histoName = "hDB_V";
    _histoMap[histoName] = new TH1F(histoName.c_str(),"",NCLUSTERS,0,NCLUSTERS);
    _histoMap[histoName]->SetStats( false );
    _histoMap[histoName]->SetYTitle("offset v [mm]");  
      
    histoName = "hDB_Sigma2_U";
    _histoMap[histoName] = new TH1F(histoName.c_str(),"",NCLUSTERS,0,NCLUSTERS);
    _histoMap[histoName]->SetStats( false );
    _histoMap[histoName]->SetYTitle("sigma2 offset u [mm^2]");        

    histoName = "hDB_Sigma2_V";
    _histoMap[histoName] = new TH1F(histoName.c_str(),"",NCLUSTERS,0,NCLUSTERS);
    _histoMap[histoName]->SetStats( false );
    _histoMap[histoName]->SetYTitle("sigma2 offset v [mm^2]"); 
    
    histoName = "hDB_Cov_UV";
    _histoMap[histoName] = new TH1F(histoName.c_str(),"",NCLUSTERS,0,NCLUSTERS);
    _histoMap[histoName]->SetStats( false );
    _histoMap[histoName]->SetYTitle("covariance u-v [mm^2]"); 
      
    // Compute the average 2x2 covariance matrix of reco tracks 
    double var_u = _trackVarUMap->GetMean();
    double var_u_error = _trackVarUMap->GetMeanError();      
    
    double var_v = _trackVarVMap->GetMean();  
    double var_v_error = _trackVarVMap->GetMeanError();       
        
    double cov_uv = _trackCovUVMap->GetMean();   
      
    // The estimated track variances may be negative too large. This can 
    // drive the cluster variance negative. To prevent this, we truncate
    // the average track variance.     
      
    double clu_minU = numeric_limits< double >::max();
    double clu_minV = numeric_limits< double >::max();
      
    for (auto iter =_sensorMap.begin(); iter!=_sensorMap.end(); iter++ ) {  
      string id = iter->first;
      
      double clu_rms2_u = std::pow(_clusterUMap[id]->GetRMS(),2);
      double clu_rms2_v = std::pow(_clusterVMap[id]->GetRMS(),2);
      
      if ( clu_rms2_u < clu_minU ) clu_minU = clu_rms2_u;
      if ( clu_rms2_v < clu_minV ) clu_minV = clu_rms2_v;
    }   
    
    streamlog_out(MESSAGE3) << std::endl
                            << "Apply minVariance cut at:   "
                            << std::setiosflags(std::ios::fixed | std::ios::internal )
                            << std::setprecision(10)
                            << _minVarianceU
                            << std::setprecision(3)
                            << std::endl;   
      
    if ( clu_minU - var_u < _minVarianceU ) {
      streamlog_out(MESSAGE3) << std::endl
                              << "original track variance is  "
                              << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(10)
                              << var_u
                              << std::setprecision(3)
                              << std::endl;
        
      var_u = clu_minU - _minVarianceU;
        
      streamlog_out(MESSAGE3) << std::endl
                              << "truncated track variance is  "
                              << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(10)
                              << var_u
                              << std::setprecision(3)
                              << std::endl;
    } 
    
    streamlog_out(MESSAGE3) << std::endl
                            << "Apply minVariance cut at:   "
                            << std::setiosflags(std::ios::fixed | std::ios::internal )
                            << std::setprecision(10)
                            << _minVarianceV
                            << std::setprecision(3)
                            << std::endl;   

    if ( clu_minV - var_v < _minVarianceV ) {
      streamlog_out(MESSAGE3) << std::endl
                              << "original track variance is  "
                              << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(10)
                              << var_v
                              << std::setprecision(3)
                              << std::endl;
      
      var_v = clu_minV - _minVarianceV;
        
      streamlog_out(MESSAGE3) << std::endl
                              << "truncated track variance is  "
                              << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(10)
                              << var_v
                              << std::setprecision(3)
                              << std::endl;
    } 
             
    // Go through all cluster shapes
    int i = 0; 
    for (auto iter =_sensorMap.begin(); iter!=_sensorMap.end(); iter++ ) {
      int counter = iter->second;  
      string id = iter->first;
      i++; 
        
      // Perform the calibration for u position 
      double clu_mean_u = _clusterUMap[id]->GetMean();
      double clu_mean_u_error = _clusterUMap[id]->GetMeanError();        
      double clu_rms_u = _clusterUMap[id]->GetRMS();
      double clu_rms_u_error = _clusterUMap[id]->GetRMSError();
        
      double clu_rms2_u = clu_rms_u*clu_rms_u;
      double clu_rms2_u_error = 2*clu_rms_u*clu_rms_u_error;
        
      // Subtract track fit variance
      clu_rms2_u -= var_u;
      clu_rms2_u_error += var_u_error;  
         
      // Perform the calibration for v position  
      double clu_mean_v = _clusterVMap[id]->GetMean();
      double clu_mean_v_error = _clusterVMap[id]->GetMeanError();   
      double clu_rms_v = _clusterVMap[id]->GetRMS();
      double clu_rms_v_error = _clusterVMap[id]->GetRMSError();   
        
      double clu_rms2_v = clu_rms_v*clu_rms_v;  
      double clu_rms2_v_error = 2*clu_rms_v*clu_rms_v_error;
        
      // Subtract track fit variance
      clu_rms2_v -= var_v;
      clu_rms2_v_error += var_v_error;  
         
      // Perform the calibration for uv covariance   
      double clu_cov_uv = _clusterUVMap[id]->GetCovariance();
      clu_cov_uv -= cov_uv;        
         
      // Store calibration result   
      histoName = "hDB_Weight";
      _histoMap[histoName]->SetBinContent( i, counter );
      _histoMap[histoName]->SetBinError( i, TMath::Sqrt(counter) );
      _histoMap[histoName]->GetXaxis()->SetBinLabel( i, id.c_str() );
        
      histoName = "hDB_U";
      _histoMap[histoName]->SetBinContent( i, clu_mean_u );
      _histoMap[histoName]->SetBinError( i, clu_mean_u_error );
      _histoMap[histoName]->GetXaxis()->SetBinLabel( i, id.c_str() );
        
      histoName = "hDB_V"; 
      _histoMap[histoName]->SetBinContent( i, clu_mean_v );
      _histoMap[histoName]->SetBinError( i, clu_mean_v_error );
      _histoMap[histoName]->GetXaxis()->SetBinLabel( i, id.c_str() );
        
      histoName = "hDB_Sigma2_U";
      _histoMap[histoName]->SetBinContent( i, clu_rms2_u );
      _histoMap[histoName]->SetBinError( i, clu_rms2_u_error );
      _histoMap[histoName]->GetXaxis()->SetBinLabel( i, id.c_str() );
         
      histoName = "hDB_Sigma2_V";
      _histoMap[histoName]->SetBinContent( i, clu_rms2_v );
      _histoMap[histoName]->SetBinError( i, clu_rms2_v_error );
      _histoMap[histoName]->GetXaxis()->SetBinLabel( i, id.c_str() );  
        
      histoName = "hDB_Cov_UV";
      _histoMap[histoName]->SetBinContent( i, clu_cov_uv );
      _histoMap[histoName]->SetBinError( i, 0 );
      _histoMap[histoName]->GetXaxis()->SetBinLabel( i, id.c_str() );
      
      streamlog_out(MESSAGE3) << "  Label:  " << id  << endl
                              << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(8)
                              << "  u: " << clu_mean_u  << ", sigma2: " << clu_rms2_u << endl
                              << "  v: " << clu_mean_v  << ", sigma2: " << clu_rms2_v << endl
                              << "  cov: " << clu_cov_uv
                              << std::setprecision(3)
                              << endl;
         
    }  
    
    // Finally, we must store the software ADC that was used to compute the 
    // cluster labels 
    TVectorD DB_swADCSteps(_swADCSteps.size());
    for(auto i = 0; i < _swADCSteps.size(); i++ ) {
      DB_swADCSteps[i] = _swADCSteps[i];
    }
    DB_swADCSteps.Write("DB_swADCSteps");
              
    streamlog_out(MESSAGE3) << "ClusterDB written to file " << _clusterDBFileName 
                            << endl; 
    
    // Close root  file
    _rootFile->Write();
    _rootFile->Close();
    delete _rootFile;
  }

  //
  // Method printing processor parameters
  //
  void GoeClusterCalibrator::printProcessorParams() const 
  {
    
    streamlog_out(MESSAGE3)  << std::endl
                             << " "
                             << "GoeClusterCalibrator Development Version, be carefull!!"
                             << " "
                             << std::endl  << std::endl;   
    
  }
  
  std::string GoeClusterCalibrator::getClusterType(std::string & id) 
  { 
    string type("");    
        
    istringstream label(id);
    string token;
        
    // Read number of digits in label  
    std::getline(label, token, 'D');      
    type += token;  
        
    // Read all digits in label 
    while (std::getline(label, token, 'D')) {
      istringstream codes(token);
      string number("");
      type += "D";
      std::getline(codes, number, '.');
      type += number + '.';
      std::getline(codes, number, '.');
      type += number; 
    }
    return type;
  }   
  

} // Namespace



