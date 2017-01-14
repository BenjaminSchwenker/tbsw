// MCTreeProducer Processor  
// 
// See MCTreeProducer.h for full documentation of processor. 
// 
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


#include "MCTreeProducer.h"

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
  MCTreeProducer aMCTreeProducer ;
  
  //
  // Constructor
  //
  MCTreeProducer::MCTreeProducer() : Processor("MCTreeProducer")
  {
    
    // Processor description
    _description = "MCTreeProducer: Matching SimtrackerHits to clusters and moving data into output TTree(s)" ;
    
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
     
    
    // Processor parameters: 
    registerProcessorParameter ("DUTPlane",
                                "Plane number of DUT along the beam line (-1 for all sensors in input collection)",
                                _idut,  static_cast < int > (3));
       
    registerProcessorParameter ("MaxResidualU",
                                "Maximum u residual for matching simHits to hits [mm]. Put -1 to deactivate cut.",
                                _maxResidualU,  static_cast < double > (0.2));
    
    registerProcessorParameter ("MaxResidualV",
                                "Maximum v residual for matching simHits to hits [mm]. Put -1 to deactivate cut.",
                                _maxResidualV,  static_cast < double > (0.2));
    
    registerProcessorParameter( "RootFileName",
                                "Output root file name",
                                _rootFileName, std::string("histos.root"));  

    registerProcessorParameter ("MinClusters",
                                "Minimum number of cluster ID occurances for clusterDB",
                                _minClusters,  static_cast < int > (100));

    registerProcessorParameter ("CreateDBFromHitMaker",
                                "Create clusterDB from hitmaker flag",
                                _createDBFromHitMaker,  static_cast < bool > (false));

        
    std::vector<float> initBiningPos;
    initBiningPos.push_back(100);
    initBiningPos.push_back(-0.02);
    initBiningPos.push_back(0.02);
    registerProcessorParameter ("BinningU", "Parameters for histogram binning: NBins, Min [mm], Max [mm]",
                                _binningU, initBiningPos );
     
    registerProcessorParameter ("BinningV", "Parameters for histogram binning: NBins, Min [mm], Max [mm]",
                                _binningV, initBiningPos );
    
    std::vector<float> initBiningDir;
    initBiningDir.push_back(100);
    initBiningDir.push_back(-1);
    initBiningDir.push_back(+1);
    registerProcessorParameter ("BinningTu", "Parameters for histogram binning: NBins, Min [rad], Max [rad]",
                                _binningTu, initBiningDir );
     
    registerProcessorParameter ("BinningTv", "Parameters for histogram binning: NBins, Min [rad], Max [rad]",
                                _binningTv, initBiningDir );

    std::vector<float> initBiningMom;
    initBiningMom.push_back(100);
    initBiningMom.push_back(0);
    initBiningMom.push_back(+10);
    registerProcessorParameter ("BinningMom", "Parameters for histogram binning: NBins, Min [GeV], Max [GeV]",
                                _binningMom, initBiningMom );
    
    
  }
  
  //
  // Method called at the beginning of data processing
  //
  void MCTreeProducer::init() {
   
    // Initialize variables
    _nRun = 0 ;
    _nEvt = 0 ;
    _timeCPU = clock()/1000 ;
    
    if ( (int)_binningU.size() != 3 ) {
      streamlog_out ( MESSAGE3 ) <<  "Invalid u axis binning parameters: Using default values" << endl;  
      _binningU.clear();
      _binningU.push_back(100);
      _binningU.push_back(-0.1);
      _binningU.push_back(+0.1);
    } 
    
    if ( (int)_binningV.size() != 3 ) {
      streamlog_out ( MESSAGE3 ) <<  "Invalid v axis binning parameters: Using default values" << endl;  
      _binningV.clear();
      _binningV.push_back(100);
      _binningV.push_back(-0.1);
      _binningV.push_back(+0.1);
    }  
    
    if ( (int)_binningTu.size() != 3 ) {
      streamlog_out ( MESSAGE3 ) <<  "Invalid du/dw axis binning parameters: Using default values" << endl;  
      _binningTu.clear();
      _binningTu.push_back(100);
      _binningTu.push_back(-1);
      _binningTu.push_back(+1);
    }  
    
    if ( (int)_binningTv.size() != 3 ) {
      streamlog_out ( MESSAGE3 ) <<  "Invalid dv/dw axis binning parameters: Using default values" << endl;  
      _binningTv.clear();
      _binningTv.push_back(100);
      _binningTv.push_back(-1);
      _binningTv.push_back(+1);
    }  

    if ( (int)_binningMom.size() != 3 ) {
      streamlog_out ( MESSAGE3 ) <<  "Invalid momentum axis binning parameters: Using default values" << endl;  
      _binningMom.clear();
      _binningMom.push_back(100);
      _binningMom.push_back(-1);
      _binningMom.push_back(+1);
    }  
    
    // Print set parameters
    printProcessorParams();
    
    // Read detector constants from gear file
    _detector.ReadGearConfiguration();  
      
    bookHistos(); 
  }
  
  //
  // Method called for each run
  //
  void MCTreeProducer::processRunHeader(LCRunHeader * run)
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
  void MCTreeProducer::processEvent(LCEvent * evt)
  {
    
    _nEvt ++ ;
     
    // Load DUT module    
    Det & dut = _detector.GetDet(_idut);   
        
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
    // ----------------------------------
    
    std::vector<SimTrackerHit*> SimHitStore;
    std::vector<TBHit> HitStore;    
 
    CellIDDecoder<SimTrackerHit> cellIDDec(simHitCol);
     
    for(int i=0; i< nSimHit ; i++)
    {  
      // Retrieve simtrackerhit 
      SimTrackerHit * simHit = dynamic_cast<SimTrackerHit*> (simHitCol->getElementAt(i));
      
      // Set current - layer ID, ladder ID and sensor ID
      int sensorID = cellIDDec(simHit)["sensorID"];
      int ipl = _detector.GetPlaneNumber(sensorID);
      
      // Store all hits on the DUT module  
      if( dut.GetPlaneNumber() == ipl )
      {         
        streamlog_out(MESSAGE2) << " SimHit at plane " << ipl << " at: (" << simHit->getPosition()[0] << ", " << simHit->getPosition()[1] << ")" 
                                << endl;
          
        SimHitStore.push_back(simHit);    
      }                            
    } 
    
    
    for(int i=0; i< nHit ; i++)
    {
      // Built a TBHit
      TrackerHitImpl * lciohit = dynamic_cast<TrackerHitImpl*>( hitCol->getElementAt(i) ) ;
      TBHit RecoHit ( lciohit  );        
         
      // We have to find plane number of the hit 
      int sensorID = RecoHit.GetDAQID();      
      int ipl = _detector.GetPlaneNumber(sensorID);  
      
      // Store all hits on the DUT module  
      if( dut.GetPlaneNumber() == ipl )
      {         
        streamlog_out(MESSAGE2) << " RecoHit at plane " << ipl << " at: (" << RecoHit.GetCoord()[0][0] << ", " << RecoHit.GetCoord()[1][0] << ")" 
                                << endl;
        
        HitStore.push_back( RecoHit );  
      }
    } 
    
    streamlog_out(MESSAGE2) << "Total of " << SimHitStore.size() << " simHit(s) on DUT" << endl; 
    streamlog_out(MESSAGE2) << "Total of " << HitStore.size() << " hit(s) on DUT" << endl; 
    
    
    
    // Record for each hit a matched simhit
    vector<int> hit2simhit(HitStore.size(), -1);
    
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
        for(int i=0; i< (int)HitStore.size() ; i++)
        {
            
          // If matched, skip hit 
          if (hit2simhit[i] >= 0) continue;
          
          for(int j=0;j<(int)SimHitStore.size(); j++)
          {
            SimTrackerHit * simHit = SimHitStore[j];
            double simHitPosU = simHit->getPosition()[0];
            double simHitPosV = simHit->getPosition()[1];
            
            TBHit& recoHit = HitStore[i];
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
  
      } // End of do loop of matching  hits to fitted positions
      while( bestsimhit>-1 &&  besthit>-1);
    }
    
    // Fill hit tree   
    streamlog_out(MESSAGE2) << "Start fill hit tree" << endl; 
    
    for(int ihit=0;ihit<(int)HitStore.size(); ++ihit)
    {
      
      // Provide default values for all branches
      _rootClusterID = "";                              
      _rootTrackMomentum = 0;               
      _rootTrackPosU = 0;           
      _rootTrackPosV = 0;                          
      _rootTrackdUdW = 0;           
      _rootTrackdVdW = 0;     
      
      if ( hit2simhit[ihit] >= 0 ) {        
        
        TBHit& hit = HitStore[ihit];           
        PixelCluster Cluster = hit.GetCluster();
        
        _rootClusterID = (string) Cluster; 
        
          
        double originPosU = dut.GetPixelCenterCoordU( Cluster.getVStart(), Cluster.getUStart()); 
        double originPosV = dut.GetPixelCenterCoordV( Cluster.getVStart(), Cluster.getUStart()); 
        
        SimTrackerHit * simHit = SimHitStore[ hit2simhit[ihit] ]; 
        Hep3Vector momentum(simHit->getMomentum()[0],simHit->getMomentum()[1],simHit->getMomentum()[2]);
        
        _rootTrackMomentum = momentum.mag();             
        _rootTrackPosU = simHit->getPosition()[0] - originPosU;          
        _rootTrackPosV = simHit->getPosition()[1] - originPosV;                           
        _rootTrackdUdW = momentum[0]/momentum[2];        
        _rootTrackdVdW = momentum[1]/momentum[2]; 
        
        if ( _createDBFromHitMaker ) {
          _rootTrackMomentum = 0;         
          _rootTrackPosU = hit.GetCoord()[0][0] - originPosU + gRandom->Gaus(0, TMath::Sqrt( hit.GetCov()[0][0]) );          
          _rootTrackPosV = hit.GetCoord()[1][0] - originPosV + gRandom->Gaus(0, TMath::Sqrt( hit.GetCov()[1][1]) );                           
          _rootTrackdUdW = 0;        
          _rootTrackdVdW = 0; 
        }
        
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
    
    return;
  }


  //
  // Method called after each event to check the data processed
  //
  void MCTreeProducer::check( LCEvent * evt ) {}

  //
  // Method called after all data processing
  //
  void MCTreeProducer::end()
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
    
    // Book histograms for valid clusterIDs
    // Used to visualize data and for computing moements
    for (auto iter =_clusterMap.begin(); iter!=_clusterMap.end(); iter++ ) {
      string id = iter->first;
      
      _rootFile->mkdir(id.c_str());  
      _rootFile->cd(id.c_str());
          
      string histoName;
          
      histoName = id + "_U";
      _histoMapU[id] = new TH1F(histoName.c_str(), histoName.c_str(),(int)_binningU[0],_binningU[1],_binningU[2]);
      _histoMapU[id]->SetXTitle("position u [mm]");  
      _histoMapU[id]->SetYTitle("number of clusters");  
          
      histoName = id + "_V";
      _histoMapV[id] = new TH1F(histoName.c_str(), histoName.c_str(),(int)_binningV[0],_binningV[1],_binningV[2]);
      _histoMapV[id]->SetXTitle("position v [mm]");  
      _histoMapV[id]->SetYTitle("number of clusters");  
      
      histoName = id + "_Tu";
      _histoMapTu[id] = new TH1F(histoName.c_str(), histoName.c_str(),(int)_binningTu[0],_binningTu[1],_binningTu[2]);
      _histoMapTu[id]->SetXTitle("du/dw  [rad]");  
      _histoMapTu[id]->SetYTitle("number of clusters");  
          
      histoName = id + "_Tv";
      _histoMapTv[id] = new TH1F(histoName.c_str(), histoName.c_str(),(int)_binningTv[0],_binningTv[1],_binningTv[2]);
      _histoMapTv[id]->SetXTitle("dv/dw  [rad]");  
      _histoMapTv[id]->SetYTitle("number of clusters");  
          
      histoName = id + "_Mom";
      _histoMapMom[id] = new TH1F(histoName.c_str(), histoName.c_str(),(int)_binningMom[0],_binningMom[1],_binningMom[2]);
      _histoMapMom[id]->SetXTitle("momentum  [GeV]");  
      _histoMapMom[id]->SetYTitle("number of clusters");  
           
      histoName = id + "_U_V";
      _histoMapU_V[id] = new TH2F( histoName.c_str(), histoName.c_str(),(int)_binningU[0],_binningU[1],_binningU[2],(int)_binningV[0],_binningV[1],_binningV[2]);
      _histoMapU_V[id]->SetXTitle("position u [mm]"); 
      _histoMapU_V[id]->SetYTitle("position v [mm]");  
      _histoMapU_V[id]->SetZTitle("number of clusters");
          
      histoName = id + "_U_Tu";
      _histoMapU_Tu[id] = new TH2F( histoName.c_str(), histoName.c_str(),(int)_binningU[0],_binningU[1],_binningU[2],(int)_binningTu[0],_binningTu[1],_binningTu[2]);
      _histoMapU_Tu[id]->SetXTitle("position u [mm]"); 
      _histoMapU_Tu[id]->SetYTitle("du/dw [rad]");  
      _histoMapU_Tu[id]->SetZTitle("number of clusters");
          
      histoName = id + "_U_Tv";
      _histoMapU_Tv[id] = new TH2F( histoName.c_str(), histoName.c_str(),(int)_binningU[0],_binningU[1],_binningU[2],(int)_binningTv[0],_binningTv[1],_binningTv[2]);
      _histoMapU_Tv[id]->SetXTitle("position u [mm]"); 
      _histoMapU_Tv[id]->SetYTitle("dv/dw [rad]");  
      _histoMapU_Tv[id]->SetZTitle("number of clusters");      
          
      histoName = id + "_U_Mom";
      _histoMapU_Mom[id] = new TH2F( histoName.c_str(), histoName.c_str(),(int)_binningU[0],_binningU[1],_binningU[2],(int)_binningMom[0],_binningMom[1],_binningMom[2]);
      _histoMapU_Mom[id]->SetXTitle("position u [mm]"); 
      _histoMapU_Mom[id]->SetYTitle("momentum [GeV]");  
      _histoMapU_Mom[id]->SetZTitle("number of clusters");   
          
      histoName = id + "_V_Tu";
      _histoMapV_Tu[id] = new TH2F( histoName.c_str(), histoName.c_str(),(int)_binningV[0],_binningV[1],_binningV[2],(int)_binningTu[0],_binningTu[1],_binningTu[2]);
      _histoMapV_Tu[id]->SetXTitle("position v [mm]"); 
      _histoMapV_Tu[id]->SetYTitle("du/dw [rad]");  
      _histoMapV_Tu[id]->SetZTitle("number of clusters");            
          
      histoName = id + "_V_Tv";
      _histoMapV_Tv[id] = new TH2F( histoName.c_str(), histoName.c_str(),(int)_binningV[0],_binningV[1],_binningV[2],(int)_binningTv[0],_binningTv[1],_binningTv[2]);
      _histoMapV_Tv[id]->SetXTitle("position v [mm]"); 
      _histoMapV_Tv[id]->SetYTitle("dv/dw [rad]");  
      _histoMapV_Tv[id]->SetZTitle("number of clusters");    
          
      histoName = id + "_V_Mom";
      _histoMapV_Mom[id] = new TH2F( histoName.c_str(), histoName.c_str(),(int)_binningV[0],_binningV[1],_binningV[2],(int)_binningMom[0],_binningMom[1],_binningMom[2]);
      _histoMapV_Mom[id]->SetXTitle("position v [mm]"); 
      _histoMapV_Mom[id]->SetYTitle("momentum [GeV]");  
      _histoMapV_Mom[id]->SetZTitle("number of clusters");    
          
      histoName = id + "_Tv_Tu";
      _histoMapTv_Tu[id] = new TH2F( histoName.c_str(), histoName.c_str(),(int)_binningTv[0],_binningTv[1],_binningTv[2],(int)_binningTu[0],_binningTu[1],_binningTu[2]);
      _histoMapTv_Tu[id]->SetXTitle("dv/dw [rad]"); 
      _histoMapTv_Tu[id]->SetYTitle("du/dw [rad]");  
      _histoMapTv_Tu[id]->SetZTitle("number of clusters");   
          
      histoName = id + "_Tv_Mom";
      _histoMapTv_Mom[id] = new TH2F( histoName.c_str(), histoName.c_str(),(int)_binningTv[0],_binningTv[1],_binningTv[2],(int)_binningMom[0],_binningMom[1],_binningMom[2]);
      _histoMapTv_Mom[id]->SetXTitle("dv/dw [rad]"); 
      _histoMapTv_Mom[id]->SetYTitle("momentum [GeV]");  
      _histoMapTv_Mom[id]->SetZTitle("number of clusters");
          
      histoName = id + "_Tu_Mom";
      _histoMapTu_Mom[id] = new TH2F( histoName.c_str(), histoName.c_str(),(int)_binningTu[0],_binningTu[1],_binningTu[2],(int)_binningMom[0],_binningMom[1],_binningMom[2]);
      _histoMapTu_Mom[id]->SetXTitle("du/dw [rad]"); 
      _histoMapTu_Mom[id]->SetYTitle("momentum [GeV]");  
      _histoMapTu_Mom[id]->SetZTitle("number of clusters");
    
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
        _histoMapTu[id]->Fill( _rootTrackdUdW  );
        _histoMapTv[id]->Fill( _rootTrackdVdW  );   
        _histoMapMom[id]->Fill( _rootTrackMomentum  );   
      
        _histoMapU_V[id]->Fill( _rootTrackPosU   , _rootTrackPosV);  
        _histoMapU_Tu[id]->Fill( _rootTrackPosU   , _rootTrackdUdW);  
        _histoMapU_Tv[id]->Fill( _rootTrackPosU   , _rootTrackdVdW);
        _histoMapU_Mom[id]->Fill( _rootTrackPosU   , _rootTrackMomentum);
        _histoMapV_Tu[id]->Fill( _rootTrackPosV   , _rootTrackdUdW);  
        _histoMapV_Tv[id]->Fill( _rootTrackPosV   , _rootTrackdVdW);
        _histoMapV_Mom[id]->Fill( _rootTrackPosV   , _rootTrackMomentum);
        _histoMapTv_Tu[id]->Fill( _rootTrackdVdW   , _rootTrackdUdW);
        _histoMapTv_Mom[id]->Fill( _rootTrackdVdW   , _rootTrackMomentum);   
        _histoMapTu_Mom[id]->Fill( _rootTrackdUdW   , _rootTrackMomentum); 
      } 
    }
    
     
    // Book histograms for clusterDB
    int NBINS = _clusterMap.size();   
    
    _rootFile->cd("");
    TH1F *hDB_ID         = new TH1F("hDB_ID","",NBINS,0,NBINS);
    TH1F *hDB_U          = new TH1F("hDB_U","",NBINS,0,NBINS);
    TH1F *hDB_V          = new TH1F("hDB_V","",NBINS,0,NBINS);  
    TH1F *hDB_Tu         = new TH1F("hDB_Tu","",NBINS,0,NBINS);
    TH1F *hDB_Tv         = new TH1F("hDB_Tv","",NBINS,0,NBINS);    
    TH1F *hDB_Mom        = new TH1F("hDB_Mom","",NBINS,0,NBINS);
    TH1F *hDB_Sigma_U    = new TH1F("hDB_Sigma_U","",NBINS,0,NBINS);  
    TH1F *hDB_Sigma_V    = new TH1F("hDB_Sigma_V","",NBINS,0,NBINS);
    TH1F *hDB_Sigma_Tu   = new TH1F("hDB_Sigma_Tu","",NBINS,0,NBINS); 
    TH1F *hDB_Sigma_Tv   = new TH1F("hDB_Sigma_Tv","",NBINS,0,NBINS);
    TH1F *hDB_Sigma_Mom  = new TH1F("hDB_Sigma_Mom","",NBINS,0,NBINS);
    TH1F *hDB_Cov_UV     = new TH1F("hDB_Cov_UV","",NBINS,0,NBINS);  
    TH1F *hDB_Corr_UV    = new TH1F("hDB_Corr_UV","",NBINS,0,NBINS);
    TH1F *hDB_Corr_UTu   = new TH1F("hDB_Corr_UTu","",NBINS,0,NBINS);    
    TH1F *hDB_Corr_UTv   = new TH1F("hDB_Corr_UTv","",NBINS,0,NBINS);
    TH1F *hDB_Corr_UMom  = new TH1F("hDB_Corr_UMom","",NBINS,0,NBINS);  
    TH1F *hDB_Corr_VTu   = new TH1F("hDB_Corr_VTu","",NBINS,0,NBINS);
    TH1F *hDB_Corr_VTv   = new TH1F("hDB_Corr_VTv","",NBINS,0,NBINS);  
    TH1F *hDB_Corr_VMom  = new TH1F("hDB_Corr_VMom","",NBINS,0,NBINS);
    TH1F *hDB_Corr_TvTu  = new TH1F("hDB_Corr_TvTu","",NBINS,0,NBINS);  
    TH1F *hDB_Corr_TvMom = new TH1F("hDB_Corr_TvMom","",NBINS,0,NBINS);
    TH1F *hDB_Corr_TuMom = new TH1F("hDB_Corr_TuMom","",NBINS,0,NBINS);  
    
    hDB_ID->SetStats( false );
    hDB_ID->SetYTitle("clusterID fraction");  
    hDB_U->SetStats( false );
    hDB_U->SetYTitle("offset u [mm]");   
    hDB_V->SetStats( false );
    hDB_V->SetYTitle("offset v [mm]"); 
    hDB_Tu->SetStats( false );
    hDB_Tu->SetYTitle("du/dw");     
    hDB_Tv->SetStats( false );
    hDB_Tv->SetYTitle("dv/dw");
    hDB_Mom->SetStats( false );
    hDB_Mom->SetYTitle("momentum [GeV]");   
    hDB_Sigma_U->SetStats( false );
    hDB_Sigma_U->SetYTitle("sigma offset u [mm]");  
    hDB_Sigma_V->SetStats( false );
    hDB_Sigma_V->SetYTitle("sigma offset v [mm]"); 
    hDB_Sigma_Tu->SetStats( false );
    hDB_Sigma_Tu->SetYTitle("sigma du/dw");   
    hDB_Sigma_Tv->SetStats( false );
    hDB_Sigma_Tv->SetYTitle("sigma dv/dw");   
    hDB_Sigma_Mom->SetStats( false );
    hDB_Sigma_Mom->SetYTitle("sigma momentum [GeV]");  
    hDB_Cov_UV->SetStats( false );
    hDB_Cov_UV->SetYTitle("covariance u-v"); 
    hDB_Corr_UV->SetStats( false );
    hDB_Corr_UV->SetYTitle("correlation factor u-v");   
    hDB_Corr_UTv->SetStats( false );
    hDB_Corr_UTv->SetYTitle("correlation factor u-dv/dw");      
    hDB_Corr_UTu->SetStats( false );
    hDB_Corr_UTu->SetYTitle("correlation factor u-du/dw");      
    hDB_Corr_UMom->SetStats( false );
    hDB_Corr_UMom->SetYTitle("correlation factor u-mom");      
    hDB_Corr_VTv->SetStats( false );
    hDB_Corr_VTv->SetYTitle("correlation factor v-dv/dw");      
    hDB_Corr_VTu->SetStats( false );
    hDB_Corr_VTu->SetYTitle("correlation factor v-du/dw");      
    hDB_Corr_VMom->SetStats( false );
    hDB_Corr_VMom->SetYTitle("correlation factor v-mom");    
    hDB_Corr_TvTu->SetStats( false );
    hDB_Corr_TvTu->SetYTitle("correlation factor dv/dw-du/dw");      
    hDB_Corr_TvMom->SetStats( false );
    hDB_Corr_TvMom->SetYTitle("correlation factor dv/dw-mom");   
    hDB_Corr_TuMom->SetStats( false );
    hDB_Corr_TuMom->SetYTitle("correlation factor du/dw-mom");   
    
    
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
       
      hDB_Tu->SetBinContent( i, _histoMapTu[id]->GetMean() );
      hDB_Tu->SetBinError( i, _histoMapTu[id]->GetMeanError() );
      hDB_Tu->GetXaxis()->SetBinLabel( i, id.c_str() );   
   
      hDB_Tv->SetBinContent( i, _histoMapTv[id]->GetMean() );
      hDB_Tv->SetBinError( i, _histoMapTv[id]->GetMeanError() );
      hDB_Tv->GetXaxis()->SetBinLabel( i, id.c_str() );
      
      hDB_Mom->SetBinContent( i, _histoMapMom[id]->GetMean() );
      hDB_Mom->SetBinError( i, _histoMapMom[id]->GetMeanError() );
      hDB_Mom->GetXaxis()->SetBinLabel( i, id.c_str() );
      
      hDB_Sigma_U->SetBinContent( i, _histoMapU[id]->GetRMS() );
      hDB_Sigma_U->SetBinError( i, _histoMapU[id]->GetRMSError() );
      hDB_Sigma_U->GetXaxis()->SetBinLabel( i, id.c_str() );
      
      hDB_Sigma_V->SetBinContent( i, _histoMapV[id]->GetRMS() );
      hDB_Sigma_V->SetBinError( i, _histoMapV[id]->GetRMSError() );
      hDB_Sigma_V->GetXaxis()->SetBinLabel( i, id.c_str() );  

      hDB_Sigma_Tu->SetBinContent( i, _histoMapTu[id]->GetRMS() );
      hDB_Sigma_Tu->SetBinError( i, _histoMapTu[id]->GetRMSError() );
      hDB_Sigma_Tu->GetXaxis()->SetBinLabel( i, id.c_str() );

      hDB_Sigma_Tv->SetBinContent( i, _histoMapTv[id]->GetRMS() );
      hDB_Sigma_Tv->SetBinError( i, _histoMapTv[id]->GetRMSError() );
      hDB_Sigma_Tv->GetXaxis()->SetBinLabel( i, id.c_str() );

      hDB_Sigma_Mom->SetBinContent( i, _histoMapMom[id]->GetRMS() );
      hDB_Sigma_Mom->SetBinError( i, _histoMapMom[id]->GetRMSError() );
      hDB_Sigma_Mom->GetXaxis()->SetBinLabel( i, id.c_str() );

      hDB_Cov_UV->SetBinContent( i, _histoMapU_V[id]->GetCovariance() );
      hDB_Cov_UV->SetBinError( i, 0 );
      hDB_Cov_UV->GetXaxis()->SetBinLabel( i, id.c_str() );

      hDB_Corr_UV->SetBinContent( i, _histoMapU_V[id]->GetCorrelationFactor() );
      hDB_Corr_UV->SetBinError( i, 0 );
      hDB_Corr_UV->GetXaxis()->SetBinLabel( i, id.c_str() );

      hDB_Corr_UTu->SetBinContent( i, _histoMapU_Tu[id]->GetCorrelationFactor() );
      hDB_Corr_UTu->SetBinError( i, 0 );
      hDB_Corr_UTu->GetXaxis()->SetBinLabel( i, id.c_str() );
      
      hDB_Corr_UTv->SetBinContent( i, _histoMapU_Tv[id]->GetCorrelationFactor() );
      hDB_Corr_UTv->SetBinError( i, 0 );
      hDB_Corr_UTv->GetXaxis()->SetBinLabel( i, id.c_str() );

      hDB_Corr_UMom->SetBinContent( i, _histoMapU_Mom[id]->GetCorrelationFactor() );
      hDB_Corr_UMom->SetBinError( i, 0 );
      hDB_Corr_UMom->GetXaxis()->SetBinLabel( i, id.c_str() );

      hDB_Corr_VTu->SetBinContent( i, _histoMapV_Tu[id]->GetCorrelationFactor() );
      hDB_Corr_VTu->SetBinError( i, 0 );
      hDB_Corr_VTu->GetXaxis()->SetBinLabel( i, id.c_str() );
      
      hDB_Corr_VTv->SetBinContent( i, _histoMapV_Tv[id]->GetCorrelationFactor() );
      hDB_Corr_VTv->SetBinError( i, 0 );
      hDB_Corr_VTv->GetXaxis()->SetBinLabel( i, id.c_str() );

      hDB_Corr_VMom->SetBinContent( i, _histoMapV_Mom[id]->GetCorrelationFactor() );
      hDB_Corr_VMom->SetBinError( i, 0 );
      hDB_Corr_VMom->GetXaxis()->SetBinLabel( i, id.c_str() );

      hDB_Corr_TvTu->SetBinContent( i, _histoMapTv_Tu[id]->GetCorrelationFactor() );
      hDB_Corr_TvTu->SetBinError( i, 0 );
      hDB_Corr_TvTu->GetXaxis()->SetBinLabel( i, id.c_str() );

      hDB_Corr_TvMom->SetBinContent( i, _histoMapTv_Mom[id]->GetCorrelationFactor() );
      hDB_Corr_TvMom->SetBinError( i, 0 );
      hDB_Corr_TvMom->GetXaxis()->SetBinLabel( i, id.c_str() );
      
      hDB_Corr_TuMom->SetBinContent( i, _histoMapTu_Mom[id]->GetCorrelationFactor() );
      hDB_Corr_TuMom->SetBinError( i, 0 );
      hDB_Corr_TuMom->GetXaxis()->SetBinLabel( i, id.c_str() );
            
      streamlog_out(MESSAGE3) << "  ClusterId:  " << id << " count: " << count << endl
                              << std::setiosflags(std::ios::fixed | std::ios::internal )
                              << std::setprecision(8)
                              << "  u: " << _histoMapU[id]->GetMean() << ", sigma: " << _histoMapU[id]->GetRMS() << endl
                              << "  v: " << _histoMapV[id]->GetMean() << ", sigma: " << _histoMapV[id]->GetRMS() << endl
                              << "  corr: " << _histoMapU_V[id]->GetCorrelationFactor() 
                              << std::setprecision(3)
                              << endl;
    
    }  
    
    streamlog_out(MESSAGE3) << "ClusterDB with " << _clusterMap.size() << " cluster shapes written to file "
                            << _rootFileName << endl; 
     
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
  void MCTreeProducer::printProcessorParams() const 
  {
    
    streamlog_out(MESSAGE3)  << std::endl
                             << " "
                             << "MCTreeProducer Development Version, be carefull!!"
                             << " "
                             << std::endl  << std::endl;   
    
  }
  
  void MCTreeProducer::bookHistos()
  {   
     
    _rootFile = new TFile( _rootFileName.c_str(),"recreate");
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



