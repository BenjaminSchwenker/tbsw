// GoeClusterCalibratorForMC Processor  
// 
// See GoeClusterCalibratorForMC.h for full documentation of processor. 
// 
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


#include "GoeClusterCalibratorForMC.h"

// TBTools includes
#include "TBHit.h"
#include "PixelCluster.h"
#include "Det.h"
#include "Utilities.h" 


// Include basic C
#include <iostream>
#include <limits>
#include <iomanip>
#include <algorithm>

// Include LCIO classes
#include "lcio.h"
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>

// Include CLHEP classes
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Vector/ThreeVector.h>

// Include ROOT classes
#include <TFile.h>
#include <TMath.h>
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
  GoeClusterCalibratorForMC aGoeClusterCalibratorForMC ;
  
  //
  // Constructor
  //
  GoeClusterCalibratorForMC::GoeClusterCalibratorForMC() : Processor("GoeClusterCalibratorForMC")
  {
    
    // Processor description
    _description = "GoeClusterCalibratorForMC: Create clusterDB for clusters using simhits" ;
    
    //
    // Input collections 
    
    registerInputCollection( LCIO::SIMTRACKERHIT,
                             "SimTrackerHitCollection" ,
                             "Name of collection with simulated hits"  ,
                             _simHitColName ,
                             std::string("SimTrackerHits") ) ;
    
    registerInputCollection( LCIO::TRACKERHIT,
                             "HitCollection" ,
                             "Name of hit collection"  ,
                             _hitColName ,
                             std::string("hit") ) ;
     
    registerProcessorParameter ("MaxResidualU",
                                "Maximum u residual for matching simHits to hits [mm]. Put -1 to deactivate cut.",
                                _maxResidualU,  static_cast < double > (0.2));
    
    registerProcessorParameter ("MaxResidualV",
                                "Maximum v residual for matching simHits to hits [mm]. Put -1 to deactivate cut.",
                                _maxResidualV,  static_cast < double > (0.2));
    
    registerProcessorParameter( "ClusterDBFileName",
                                "Output clusterDB file name",
                                _clusterDBFileName, std::string("clusterDB.root"));  
     
    registerProcessorParameter ("MinClusters",
                                "Minimum number of cluster ID occurances for clusterDB",
                                _minClusters,  static_cast < int > (2000));
    
    std::vector<int> initIgnoreIDVec;
    registerProcessorParameter ("IgnoreIDs",
                                "Ignore clusters from list of sensorIDs",
                                _ignoreIDVec, initIgnoreIDVec);
    
    std::vector<int> initSoftADCSteps;
    registerProcessorParameter ("SoftwareADC",
                                "List of steps for software ADC. An empty list gives a constant transfer curve (0bit limit)",
                                _swADCSteps, initSoftADCSteps);
        
  }
  
  //
  // Method called at the beginning of data processing
  //
  void GoeClusterCalibratorForMC::init() {
    
    // Initialize variables
    _nRun = 0 ;
    _nEvt = 0 ;
    _timeCPU = clock()/1000 ;
    
    // Make sure adc steps are sorted
    std::sort(_swADCSteps.begin(), _swADCSteps.end());
    
    // Erase duplicated entries
    _swADCSteps.erase( std::unique( _swADCSteps.begin(), _swADCSteps.end() ), _swADCSteps.end() ); 
    
    for(auto step : _swADCSteps ) {
      streamlog_out( MESSAGE2 ) << " sw adc step "  << step << endl;
    }
    
    // Print set parameters
    printProcessorParams();
    
    // Read detector constants from gear file
    _detector.ReadGearConfiguration();  
  }
  
  //
  // Method called for each run
  //
  void GoeClusterCalibratorForMC::processRunHeader(LCRunHeader * run)
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
  void GoeClusterCalibratorForMC::processEvent(LCEvent * evt)
  {
    
    _nEvt ++ ;
     
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
    // Get hit collection 
    // 
    
    LCCollection* hitCol = 0;
    int nHit = 0; 
    try {
      hitCol = evt->getCollection( _hitColName ) ;
      nHit = hitCol->getNumberOfElements();
    } catch (lcio::DataNotAvailableException& e) {
      streamlog_out(MESSAGE2) << "Not able to get collection "
                              << _hitColName
                              << " from event " << evt->getEventNumber()
                              << " in run " << evt->getRunNumber() << endl << endl;
    } 
    
    streamlog_out(MESSAGE2) << "Total of " << nHit << " hit(s) in collection " << _hitColName << endl;
    
    // Read hits and simHits 
    // ---------------------
    
    std::map<int, std::vector<SimTrackerHit*> > SimHitMap;
    std::map<int, std::vector<TBHit> > RecoHitMap;    
    
    CellIDDecoder<SimTrackerHit> cellIDDec(simHitCol);
     
    for(int i=0; i< nSimHit ; i++)
    {  
      // Retrieve simtrackerhit 
      SimTrackerHit * simHit = dynamic_cast<SimTrackerHit*> (simHitCol->getElementAt(i));
      int sensorID = cellIDDec(simHit)["sensorID"];
      
      streamlog_out(MESSAGE2) << " SimHit with sensorID " << sensorID << " at: (" << simHit->getPosition()[0] << ", " << simHit->getPosition()[1] << ")" 
                              << endl;
          
      bool ignoreID = false;
      for (auto id :  _ignoreIDVec)  {
        if  (id == sensorID) ignoreID = true; 
      }
      
      if ( !ignoreID ) SimHitMap[sensorID].push_back(simHit);            
    } 
    
    
    for(int i=0; i< nHit ; i++)
    {
      // Built a TBHit
      TrackerHitImpl * lciohit = dynamic_cast<TrackerHitImpl*>( hitCol->getElementAt(i) ) ;
      TBHit recoHit ( lciohit  );        
      int sensorID = recoHit.GetDAQID();      
                
      streamlog_out(MESSAGE2) << " RecoHit with sensorID " << sensorID << " at: (" << recoHit.GetCoord()[0][0] << ", " << recoHit.GetCoord()[1][0] << ")" 
                              << endl;
        
      bool ignoreID = false;
      for (auto id :  _ignoreIDVec)  {
        if  (id == sensorID) ignoreID = true; 
      }

      if ( !ignoreID ) RecoHitMap[sensorID].push_back( recoHit );    
    } 
    
    // Go through all sensors 
    for (auto it=RecoHitMap.begin(); it!=RecoHitMap.end(); it++) {    
        
      auto sensorID = it->first; 
      
      // Go to next sensor, if there are no simHits 
      if (SimHitMap.find(sensorID) == SimHitMap.end() ) continue;
      
      // Vector with all recoHits on current sensor
      auto RecoHits = it->second;
      
      // Vector with all simHits on current sensor
      auto SimHits = SimHitMap[sensorID];
      
      // Record for each hit a matched simhit
      vector<int> hit2simhit(RecoHits.size(), -1);
      
      // Continue matching until all hits are matched 
      // or no hit is close enough!!
      
      { 
        double distmin=numeric_limits<double >::max();
        int bestsimhit=-1;   
        int besthit=-1; 
        
        do{
          bestsimhit=-1;
          besthit=-1; 
          distmin=numeric_limits<double >::max();
           
          // Find hit/simhit pair with minimum chi2 distance.  
          for(int i=0; i< (int)RecoHits.size() ; i++)
          {
            
            // If matched, skip hit 
            if (hit2simhit[i] >= 0) continue;
          
            for(int j=0;j<(int)SimHits.size(); j++)
            {
              SimTrackerHit * simHit = SimHits[j];
              double simHitPosU = simHit->getPosition()[0];
              double simHitPosV = simHit->getPosition()[1];
               
              TBHit& recoHit = RecoHits[i];
              double hitPosU = recoHit.GetCoord()[0][0];
              double hitPosV = recoHit.GetCoord()[1][0];
               
              // Skip all hits with too large residuum 
              if ( std::abs(hitPosU-simHitPosU) >= _maxResidualU && _maxResidualU > 0 ) continue;  
              if ( std::abs(hitPosV-simHitPosV) >= _maxResidualV && _maxResidualV > 0 ) continue; 
              
              // Finally, we will use a simple 2D distance to select best matching hit
              double hitdist = 0; 
              if ( _maxResidualU > 0 )  hitdist += std::abs(hitPosU-simHitPosU); 
              if ( _maxResidualV > 0 )  hitdist += std::abs(hitPosV-simHitPosV); 
               
              if( hitdist<distmin )
              {
                distmin=hitdist;
                besthit=i;
                bestsimhit=j;
              }
            }
          }
          
          streamlog_out(MESSAGE2) << "In matching loop: best hit " << besthit << " to simhit " << bestsimhit << endl; 
          streamlog_out(MESSAGE2) << "  distmin: " <<  distmin  << endl; 
          
          // Check if a match was found
          if( bestsimhit>-1 &&  besthit>-1   )
          {   
            streamlog_out(MESSAGE2) << "  match found!!!"   << endl;
            hit2simhit[besthit] = bestsimhit;
          } 
          
        } while( bestsimhit>-1 &&  besthit>-1);
      }
      
      for(int ihit=0;ihit<(int)RecoHits.size(); ++ihit)
      {
          
        if ( hit2simhit[ihit] >= 0 ) {        
          
          TBHit& hit = RecoHits[ihit];     
          int sensorID = hit.GetDAQID();     
          int ipl = _detector.GetPlaneNumber(sensorID);
          Det & Sensor = _detector.GetDet(ipl);   
           
          SimTrackerHit * simHit = SimHits[ hit2simhit[ihit] ]; 
          Hep3Vector momentum(simHit->getMomentum()[0],simHit->getMomentum()[1],simHit->getMomentum()[2]);
          
          // Get local track parameters 
          double trk_tu = momentum[0]/momentum[2];    // rad
          double trk_tv = momentum[1]/momentum[2];    // rad
          double trk_u = simHit->getPosition()[0];    // mm
          double trk_v = simHit->getPosition()[1];    // mm
          double trk_mom = momentum.mag();            // GeV
          
          // Get cluster label  
          PixelCluster Cluster = hit.GetCluster(); 
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
          
          // Count how many times a label appear 
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
  void GoeClusterCalibratorForMC::check( LCEvent * evt ) {}
  
  //
  // Method called after all data processing
  //
  void GoeClusterCalibratorForMC::end()
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
    
    // Count all clusters
    double countAll = 0;
    
    // Count reject clusters
    int countReject = 0;
    
    for(auto iter = _sensorMap.begin(); iter != _sensorMap.end(); ) {
      auto id = iter->first; 
      auto counter = iter->second;
      
      countAll += counter;
        
      if(counter < _minClusters ) {
        streamlog_out(MESSAGE3) << "  Deleting label:  " << id << " because too few counts (" << counter << ")" << endl;
        iter = _sensorMap.erase(iter);
        countReject += counter;
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
         
      // Perform the calibration for v position  
      double clu_mean_v = _clusterVMap[id]->GetMean();
      double clu_mean_v_error = _clusterVMap[id]->GetMeanError();   
      double clu_rms_v = _clusterVMap[id]->GetRMS();
      double clu_rms_v_error = _clusterVMap[id]->GetRMSError();   
        
      double clu_rms2_v = clu_rms_v*clu_rms_v;  
      double clu_rms2_v_error = 2*clu_rms_v*clu_rms_v_error;
         
      // Perform the calibration for uv covariance   
      double clu_cov_uv = _clusterUVMap[id]->GetCovariance();
                 
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
  void GoeClusterCalibratorForMC::printProcessorParams() const 
  {
    
    streamlog_out(MESSAGE3)  << std::endl
                             << " "
                             << "GoeClusterCalibratorForMC Development Version, be carefull!!"
                             << " "
                             << std::endl  << std::endl;   
    
  }
  
  

} // Namespace



