// GoeClusterCalibratorFromMC Processor  
// 
// See GoeClusterCalibratorFromMC.h for full documentation of processor. 
// 
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


#include "GoeClusterCalibratorFromMC.h"

// TBTools includes
#include "TBHit.h"
#include "PixelCluster.h"
#include "Det.h"
#include "Utilities.h" 
#include "PhysicalConstants.h"


// Include basic C
#include <iostream>
#include <limits>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <string>
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



// Used namespaces
using namespace std; 
using namespace CLHEP; 
using namespace lcio ;
using namespace marlin ;

namespace depfet {

  //
  // Instantiate this object
  //
  GoeClusterCalibratorFromMC aGoeClusterCalibratorFromMC ;
  
  //
  // Constructor
  //
  GoeClusterCalibratorFromMC::GoeClusterCalibratorFromMC() : Processor("GoeClusterCalibratorFromMC")
  {
    
    // Processor description
    _description = "GoeClusterCalibratorFromMC: Create clusterDB for clusters using simhits" ;
    
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
                                _minClusters,  static_cast < int > (100));
    
  }
  
  //
  // Method called at the beginning of data processing
  //
  void GoeClusterCalibratorFromMC::init() {
    
    // Initialize variables
    _nRun = 0 ;
    _nEvt = 0 ;
    _timeCPU = clock()/1000 ;
    
    // Print set parameters
    printProcessorParams();
    
    // Read detector constants from gear file
    _detector.ReadGearConfiguration();  
      
    bookHistos(); 
  }
  
  //
  // Method called for each run
  //
  void GoeClusterCalibratorFromMC::processRunHeader(LCRunHeader * run)
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
  void GoeClusterCalibratorFromMC::processEvent(LCEvent * evt)
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
          
      SimHitMap[sensorID].push_back(simHit);            
    } 
    
    
    for(int i=0; i< nHit ; i++)
    {
      // Built a TBHit
      TrackerHitImpl * lciohit = dynamic_cast<TrackerHitImpl*>( hitCol->getElementAt(i) ) ;
      TBHit recoHit ( lciohit  );        
      int sensorID = recoHit.GetDAQID();      
                
      streamlog_out(MESSAGE2) << " RecoHit with sensorID " << sensorID << " at: (" << recoHit.GetCoord()[0][0] << ", " << recoHit.GetCoord()[1][0] << ")" 
                              << endl;
        
      RecoHitMap[sensorID].push_back( recoHit );    
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
        // Provide default values for all branches
        _rootClusterID = "";                              
        _rootTrackMomentum = 0;               
        _rootTrackPosU = 0;           
        _rootTrackPosV = 0;                          
        _rootTrackdUdW = 0;           
        _rootTrackdVdW = 0;     
        
        if ( hit2simhit[ihit] >= 0 ) {        
          
          TBHit& hit = RecoHits[ihit];     
          int sensorID = hit.GetDAQID();     
          int ipl = _detector.GetPlaneNumber(sensorID);
           
          PixelCluster Cluster = hit.GetCluster();
          Det & Sensor = _detector.GetDet(ipl);   
          
          double originPosU = Sensor.GetPixelCenterCoordU( Cluster.getVStart(), Cluster.getUStart()); 
          double originPosV = Sensor.GetPixelCenterCoordV( Cluster.getVStart(), Cluster.getUStart()); 
          
          SimTrackerHit * simHit = SimHits[ hit2simhit[ihit] ]; 
          Hep3Vector momentum(simHit->getMomentum()[0],simHit->getMomentum()[1],simHit->getMomentum()[2]);
          
          _rootClusterID = (string) Cluster; 
          _rootTrackMomentum = momentum.mag();             
          _rootTrackPosU = simHit->getPosition()[0] - originPosU;          
          _rootTrackPosV = simHit->getPosition()[1] - originPosV;                           
          _rootTrackdUdW = momentum[0]/momentum[2];        
          _rootTrackdVdW = momentum[1]/momentum[2]; 
          
          // Fill tree with set variables 
          _rootFile->cd("");
          _rootClusterTree->Fill();
          
          // Register new cluster if needed
          if (_clusterMap.find(Cluster) == _clusterMap.end() ) {
            _clusterMap[Cluster] = 0;
          }
          
          // Count how many times a clusterID appear 
          _clusterMap[Cluster]++;  
        }
      }
    }
    return;
  }
  
  
  //
  // Method called after each event to check the data processed
  //
  void GoeClusterCalibratorFromMC::check( LCEvent * evt ) {}
  
  //
  // Method called after all data processing
  //
  void GoeClusterCalibratorFromMC::end()
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
    
    streamlog_out(MESSAGE3) << "Create the clusterDB ... " << endl; 
     
    // Erase all clusterIDs with too few entries  
    for(auto it = _clusterMap.begin(); it != _clusterMap.end(); ) {
      if(it->second < _minClusters ) {
        streamlog_out(MESSAGE3) << "  Deleting clusterId:  " << it->first << " because too few counts (" << it->second << ")." 
                                << endl;
        it = _clusterMap.erase(it);
      } else {
        ++it;
      }
    }
    
    // The histograms are used to compute unbinned moments of
    // cluster positions stored in "Cluster" tree. They will 
    // never be used to display data and are not added to 
    // the clusterDB.  
    for (auto iter =_clusterMap.begin(); iter!=_clusterMap.end(); iter++ ) {
      string id = iter->first;    
      string histoName;
    
      // TODO: the range of the histograms should not be hard coded. Better to 
      // use the size of the sensitive area insted. 
      histoName = id + "_U";
      _histoMapU[id] = new TH1F(histoName.c_str(), histoName.c_str(),1,-1,+1);
      _histoMapU[id]->SetDirectory(0);
          
      histoName = id + "_V";
      _histoMapV[id] = new TH1F(histoName.c_str(), histoName.c_str(),1,-1,+1);
      _histoMapV[id]->SetDirectory(0);
      
      histoName = id + "_UV";
      _histoMapUV[id] = new TH2F( histoName.c_str(), histoName.c_str(),1,-1,+1,1,-1,+1);
      _histoMapUV[id]->SetDirectory(0);          
    }
    
    // Now fill the histograms by looping over the matched simhits
    _rootClusterID = "";                              
    _rootTrackMomentum = 0;               
    _rootTrackPosU = 0;           
    _rootTrackPosV = 0;                          
    _rootTrackdUdW = 0;           
    _rootTrackdVdW = 0;     
    
    for (int ii = 0; ii < _rootClusterTree->GetEntries(); ii++) {
      _rootClusterTree->GetEntry(ii);
      
      std::string id = _rootClusterID;
      
      // Ignore simhits creating clusters with too few entries
      if (_clusterMap.find(id) == _clusterMap.end() ) {
        streamlog_out(MESSAGE2) << "Ignore data for clusterID " << id   
                                << endl;  
         
      } else {
      
        _histoMapU[id]->Fill( _rootTrackPosU  );
        _histoMapV[id]->Fill( _rootTrackPosV  );
        _histoMapUV[id]->Fill( _rootTrackPosU   , _rootTrackPosV);  
        
      } 
    }
     
    // Book histograms for clusterDB
    int NBINS = _clusterMap.size();   
    
    _rootFile->cd("");
    TH1F *hDB_ID         = new TH1F("hDB_ID","",NBINS,0,NBINS);
    TH1F *hDB_U          = new TH1F("hDB_U","",NBINS,0,NBINS);
    TH1F *hDB_V          = new TH1F("hDB_V","",NBINS,0,NBINS);   
    TH1F *hDB_Sigma_U    = new TH1F("hDB_Sigma_U","",NBINS,0,NBINS);  
    TH1F *hDB_Sigma_V    = new TH1F("hDB_Sigma_V","",NBINS,0,NBINS);
    TH1F *hDB_Cov_UV     = new TH1F("hDB_Cov_UV","",NBINS,0,NBINS);  
    
    hDB_ID->SetStats( false );
    hDB_ID->SetYTitle("clusterID fraction");  
    hDB_U->SetStats( false );
    hDB_U->SetYTitle("offset u [mm]");   
    hDB_V->SetStats( false );
    hDB_V->SetYTitle("offset v [mm]"); 
    hDB_Sigma_U->SetStats( false );
    hDB_Sigma_U->SetYTitle("sigma offset u [mm]");  
    hDB_Sigma_V->SetStats( false );
    hDB_Sigma_V->SetYTitle("sigma offset v [mm]"); 
    hDB_Cov_UV->SetStats( false );
    hDB_Cov_UV->SetYTitle("covariance u-v"); 
    
    int i = 0; 
    
    // Go through all cluster shapes
    for (auto iter =_clusterMap.begin(); iter!=_clusterMap.end(); iter++ ) {
      int count = iter->second;  
      string id = iter->first;
      i++;  
       
      hDB_ID->SetBinContent( i, count );
      hDB_ID->SetBinError( i, TMath::Sqrt(count) );
      hDB_ID->GetXaxis()->SetBinLabel( i, id.c_str() );

      hDB_U->SetBinContent( i, _histoMapU[id]->GetMean() );
      hDB_U->SetBinError( i, _histoMapU[id]->GetMeanError() );
      hDB_U->GetXaxis()->SetBinLabel( i, id.c_str() );
      
      hDB_V->SetBinContent( i, _histoMapV[id]->GetMean() );
      hDB_V->SetBinError( i, _histoMapV[id]->GetMeanError() );
      hDB_V->GetXaxis()->SetBinLabel( i, id.c_str() );
       
      hDB_Sigma_U->SetBinContent( i, _histoMapU[id]->GetRMS() );
      hDB_Sigma_U->SetBinError( i, _histoMapU[id]->GetRMSError() );
      hDB_Sigma_U->GetXaxis()->SetBinLabel( i, id.c_str() );
      
      hDB_Sigma_V->SetBinContent( i, _histoMapV[id]->GetRMS() );
      hDB_Sigma_V->SetBinError( i, _histoMapV[id]->GetRMSError() );
      hDB_Sigma_V->GetXaxis()->SetBinLabel( i, id.c_str() );  

      hDB_Cov_UV->SetBinContent( i, _histoMapUV[id]->GetCovariance() );
      hDB_Cov_UV->SetBinError( i, 0 );
      hDB_Cov_UV->GetXaxis()->SetBinLabel( i, id.c_str() );
      
      streamlog_out(MESSAGE3) << "  ClusterId:  " << id << " count: " << count << endl
                              << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(8)
                              << "  u: " << _histoMapU[id]->GetMean() << ", sigma: " << _histoMapU[id]->GetRMS() << endl
                              << "  v: " << _histoMapV[id]->GetMean() << ", sigma: " << _histoMapV[id]->GetRMS() << endl
                              << "  corr: " << _histoMapUV[id]->GetCorrelationFactor()
                              << std::setprecision(3)
                              << endl;
    
    }  
    
    streamlog_out(MESSAGE3) << "ClusterDB with " << _clusterMap.size() << " cluster shapes written to file "
                            << _clusterDBFileName << endl; 
     
    double normID = hDB_ID->Integral();
    if (normID > 0 ) hDB_ID->Scale(1.0/normID, "width");
    
    
    // Close root  file
    _rootFile->Write();
    _rootFile->Close();
    delete _rootFile;
  }

  //
  // Method printing processor parameters
  //
  void GoeClusterCalibratorFromMC::printProcessorParams() const 
  {
    
    streamlog_out(MESSAGE3)  << std::endl
                             << " "
                             << "GoeClusterCalibratorFromMC Development Version, be carefull!!"
                             << " "
                             << std::endl  << std::endl;   
    
  }
  
  void GoeClusterCalibratorFromMC::bookHistos()
  {   
     
    _rootFile = new TFile( _clusterDBFileName.c_str(),"recreate");
    _rootFile->cd("");
      
    // 
    // Cluster Tree 
    _rootClusterTree = new TTree("Cluster","Cluster info");
    _rootClusterTree->Branch("clusterID"             ,&_rootClusterID              );
    _rootClusterTree->Branch("trackMomentum"         ,&_rootTrackMomentum          ,"trackMomentum/D"); 
    _rootClusterTree->Branch("trackPosU"             ,&_rootTrackPosU              ,"trackPosU/D");
    _rootClusterTree->Branch("trackPosV"             ,&_rootTrackPosV              ,"trackPosV/D"); 
    _rootClusterTree->Branch("trackDuDw"             ,&_rootTrackdUdW              ,"trackDuDw/D");
    _rootClusterTree->Branch("trackDvDw"             ,&_rootTrackdVdW              ,"trackDvDw/D");   
    
    // Auto save every 5 MB
    _rootClusterTree->SetAutoSave(5000000);     
  }

} // Namespace



