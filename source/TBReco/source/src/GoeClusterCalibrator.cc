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

// Include LCIO classes
#include "lcio.h"
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <Exceptions.h>

// Include ROOT classes
#include <TMath.h>
#include <TFile.h>


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
    
    registerProcessorParameter ("MinVarianceU", "Minimum value of variance for u position measurement [mm^2]",
                                _minVarianceU, static_cast < float > (1.0E-6) );
    
    registerProcessorParameter ("MinVarianceV", "Minimum value of variance for v position measurement [mm^2]",
                                _minVarianceV, static_cast < float > (1.0E-6) );
     
    registerProcessorParameter ("AlignmentDBFileName",
                                "This is the name of the LCIO file with the alignment constants (add .slcio)",
                                _alignmentDBFileName, static_cast< string > ( "eudet-alignmentDB.slcio" ) ); 
    
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
        
        TBTrackElement& TE = track.GetTE(ipl);  
        
        // Ignore track elements w/o measurment
        if ( TE.HasHit() ) { 
          
          // Get sensor data 
          //------------------------
          Det & Sensor = _detector.GetDet(ipl);  
          int sensorID = Sensor.GetDAQID(); 
          
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
          
          if (_trackVarUMap.find(sensorID) == _trackVarUMap.end() ) {

            _trackVarUMap[sensorID] = new TH1F("","",1,0,1);
            _trackVarUMap[sensorID]->SetDirectory(0);
            _trackVarUMap[sensorID]->StatOverflows(); 	            

            _trackVarVMap[sensorID] = new TH1F("","",1,0,1); 
            _trackVarVMap[sensorID]->SetDirectory(0);
            _trackVarVMap[sensorID]->StatOverflows();         

            _trackCovUVMap[sensorID] = new TH1F("","",1,0,1);
            _trackCovUVMap[sensorID]->SetDirectory(0);
            _trackCovUVMap[sensorID]->StatOverflows(); 	
            
          }
          
          _trackVarUMap[sensorID]->Fill(sigma2_u);
          _trackVarVMap[sensorID]->Fill(sigma2_v); 
          _trackCovUVMap[sensorID]->Fill(cov_uv);  
          
           
          // Get cluster id  
          PixelCluster Cluster = TE.GetHit().GetCluster(); 
          string id = Cluster.getClusterID(); 
          
          // Register new cluster if needed
          if (_sensorMap[sensorID].find(id) == _sensorMap[sensorID].end() ) {
            _sensorMap[sensorID][id] = 0;
            
            _clusterUMap[sensorID][id] = new TH1F("","",1,0,1);
            _clusterUMap[sensorID][id]->SetDirectory(0);
            _clusterUMap[sensorID][id]->StatOverflows(); 	
             
            _clusterVMap[sensorID][id] = new TH1F("","",1,0,1);
            _clusterVMap[sensorID][id]->SetDirectory(0); 
            _clusterVMap[sensorID][id]->StatOverflows(); 	
            
            _clusterUVMap[sensorID][id] = new TH2F("","",1,0,1,1,0,1);
            _clusterUVMap[sensorID][id]->SetDirectory(0);
            _clusterUVMap[sensorID][id]->StatOverflows(); 	
          }
          
          trk_u -= Sensor.GetPixelCenterCoordU( Cluster.getVStart(), Cluster.getUStart()); 
          trk_v -= Sensor.GetPixelCenterCoordV( Cluster.getVStart(), Cluster.getUStart()); 
          
          // Count how many times a clusterID appear 
          _sensorMap[sensorID][id]++;  
          _clusterUMap[sensorID][id]->Fill( trk_u ); 
          _clusterVMap[sensorID][id]->Fill( trk_v );     
          _clusterUVMap[sensorID][id]->Fill( trk_u, trk_v ); 
          
          
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
    
    for(auto it = _sensorMap.begin(); it != _sensorMap.end(); it++) {
      auto sensorID = it->first;
      auto&& clusterMap = it->second;
    
      // Count all clusters
      double countAll = 0;
    
      // Count reject clusters
      int countReject = 0;
  
      // Delete cluster ids with too small counter
      for(auto iter = clusterMap.begin(); iter != clusterMap.end(); ) {
        auto id = iter->first; 
        auto counter = iter->second;
        
        countAll += counter;
        
        if(counter < _minClusters ) {
          streamlog_out(MESSAGE2) << "  Deleting clusterId:  " << id << " because too few counts (" << counter << ") on sensorID " << sensorID 
                                  << endl;
          iter = clusterMap.erase(iter);
          countReject += counter;
        } else {
          ++iter;
        }
      }
      streamlog_out(MESSAGE3) << "Number of rejected clusters on sensorID " << sensorID << " is: " 
                              << countReject << " (" << 100.0*countReject/countAll  << "%)"  << endl;
    }
    
    streamlog_out(MESSAGE3) << "Create the clusterDB ... " << endl; 
    
    TFile * _rootFile = new TFile( _clusterDBFileName.c_str(),"recreate");
    _rootFile->cd("");
    
    // Loop over all registered sensors 
    for(auto it = _sensorMap.begin(); it != _sensorMap.end(); it++) {
      auto sensorID = it->first;
      auto&& clusterMap = it->second;
      
      // Book histograms for clusterDB
      int NCLUSTERS = clusterMap.size(); 
      string histoName;  
       
      _rootFile->cd("");
      
      histoName = "hDB_sensor" + to_string(sensorID) + "_ID";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",NCLUSTERS,0,NCLUSTERS);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("clusterID fraction");  
      
      histoName = "hDB_sensor" + to_string(sensorID) + "_U";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",NCLUSTERS,0,NCLUSTERS);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("offset u [mm]");  
      
      histoName = "hDB_sensor" + to_string(sensorID) + "_V";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",NCLUSTERS,0,NCLUSTERS);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("offset v [mm]");  
      
      histoName = "hDB_sensor" + to_string(sensorID) + "_Sigma2_U";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",NCLUSTERS,0,NCLUSTERS);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("sigma2 offset u [mm^2]");        

      histoName = "hDB_sensor" + to_string(sensorID) + "_Sigma2_V";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",NCLUSTERS,0,NCLUSTERS);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("sigma2 offset v [mm^2]"); 

      histoName = "hDB_sensor" + to_string(sensorID) + "_Cov_UV";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",NCLUSTERS,0,NCLUSTERS);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("covariance u-v [mm^2]"); 
      
      // Compute the average 2x2 covariance matrix of reco tracks 
      double var_u = _trackVarUMap[sensorID]->GetMean();
      double var_u_error = _trackVarUMap[sensorID]->GetMeanError();      

      double var_v = _trackVarVMap[sensorID]->GetMean();  
      double var_v_error = _trackVarVMap[sensorID]->GetMeanError();       
        
      double cov_uv = _trackCovUVMap[sensorID]->GetMean();   
      
      // The estimated track variances may be negative too large. This can 
      // drive the cluster variance negative. To prevent this, we truncate
      // the average track variance.     
      
      double clu_minU = numeric_limits< double >::max();
      double clu_minV = numeric_limits< double >::max();
      
      for (auto iter =clusterMap.begin(); iter!=clusterMap.end(); iter++ ) {  
        string id = iter->first;

        double clu_rms2_u = std::pow(_clusterUMap[sensorID][id]->GetRMS(),2);
        double clu_rms2_v = std::pow(_clusterVMap[sensorID][id]->GetRMS(),2);

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
      for (auto iter =clusterMap.begin(); iter!=clusterMap.end(); iter++ ) {
        int counter = iter->second;  
        string id = iter->first;
        i++; 
        
        // Perform the calibration for u position 
        double clu_mean_u = _clusterUMap[sensorID][id]->GetMean();
        double clu_mean_u_error = _clusterUMap[sensorID][id]->GetMeanError();        
        double clu_rms_u = _clusterUMap[sensorID][id]->GetRMS();
        double clu_rms_u_error = _clusterUMap[sensorID][id]->GetRMSError();
        
        double clu_rms2_u = clu_rms_u*clu_rms_u;
        double clu_rms2_u_error = 2*clu_rms_u*clu_rms_u_error;
        
        // Subtract track fit variance
        clu_rms2_u -= var_u;
        clu_rms2_u_error += var_u_error;  
         
        // Perform the calibration for v position  
        double clu_mean_v = _clusterVMap[sensorID][id]->GetMean();
        double clu_mean_v_error = _clusterVMap[sensorID][id]->GetMeanError();   
        double clu_rms_v = _clusterVMap[sensorID][id]->GetRMS();
        double clu_rms_v_error = _clusterVMap[sensorID][id]->GetRMSError();   
        
        double clu_rms2_v = clu_rms_v*clu_rms_v;  
        double clu_rms2_v_error = 2*clu_rms_v*clu_rms_v_error;
        
        // Subtract track fit variance
        clu_rms2_v -= var_v;
        clu_rms2_v_error += var_v_error;  
         
        // Perform the calibration for uv covariance   
        double clu_cov_uv = _clusterUVMap[sensorID][id]->GetCovariance();
        clu_cov_uv -= cov_uv;        
         
        // Store calibration result   
        histoName = "hDB_sensor" + to_string(sensorID) + "_ID";
        _histoMap[histoName]->SetBinContent( i, counter );
        _histoMap[histoName]->SetBinError( i, TMath::Sqrt(counter) );
        _histoMap[histoName]->GetXaxis()->SetBinLabel( i, id.c_str() );
        
        histoName = "hDB_sensor" + to_string(sensorID) + "_U";
        _histoMap[histoName]->SetBinContent( i, clu_mean_u );
        _histoMap[histoName]->SetBinError( i, clu_mean_u_error );
        _histoMap[histoName]->GetXaxis()->SetBinLabel( i, id.c_str() );
        
        histoName = "hDB_sensor" + to_string(sensorID) + "_V"; 
        _histoMap[histoName]->SetBinContent( i, clu_mean_v );
        _histoMap[histoName]->SetBinError( i, clu_mean_v_error );
        _histoMap[histoName]->GetXaxis()->SetBinLabel( i, id.c_str() );
        
        histoName = "hDB_sensor" + to_string(sensorID) + "_Sigma2_U";
        _histoMap[histoName]->SetBinContent( i, clu_rms2_u );
        _histoMap[histoName]->SetBinError( i, clu_rms2_u_error );
        _histoMap[histoName]->GetXaxis()->SetBinLabel( i, id.c_str() );
         
        histoName = "hDB_sensor" + to_string(sensorID) + "_Sigma2_V";
        _histoMap[histoName]->SetBinContent( i, clu_rms2_v );
        _histoMap[histoName]->SetBinError( i, clu_rms2_v_error );
        _histoMap[histoName]->GetXaxis()->SetBinLabel( i, id.c_str() );  
        
        histoName = "hDB_sensor" + to_string(sensorID) + "_Cov_UV";
        _histoMap[histoName]->SetBinContent( i, clu_cov_uv );
        _histoMap[histoName]->SetBinError( i, 0 );
        _histoMap[histoName]->GetXaxis()->SetBinLabel( i, id.c_str() );

        streamlog_out(MESSAGE3) << "  ClusterId:  " << id << " sensorID: " << sensorID << endl
                                << std::setiosflags(std::ios::fixed | std::ios::internal )
                                << std::setprecision(8)
                                << "  u: " << clu_mean_u  << ", sigma2: " << clu_rms2_u << endl
                                << "  v: " << clu_mean_v  << ", sigma2: " << clu_rms2_v << endl
                                << "  cov: " << clu_cov_uv
                                << std::setprecision(3)
                                << endl;
         
      }  
            
      // Normalaize the cluster ID spectrum to unit area for 
      // better comparison between data sets 
      
      histoName = "hDB_sensor" + to_string(sensorID) + "_ID";
      double normID = _histoMap[histoName]->Integral();
      if (normID > 0 ) _histoMap[histoName]->Scale(1.0/normID, "width");
      
    }

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
  
  

} // Namespace



