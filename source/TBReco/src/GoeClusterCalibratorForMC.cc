// GoeClusterCalibratorForMC Processor  
// 
// See GoeClusterCalibratorForMC.h for full documentation of processor. 
// 
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


#include "GoeClusterCalibratorForMC.h"

// TBTools includes
#include "TBDetector.h"
#include "TBHit.h"
#include "PixelCluster.h"
#include "PolyClusterDescriptor.h"
#include "Utilities.h"

// Include basic C
#include <iostream>
#include <limits>
#include <iomanip>
#include <algorithm>
#include <sstream>

// Include LCIO classes
#include "lcio.h"
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>


// Include ROOT classes
#include <TMath.h>
#include <TFile.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>


// Used namespaces
using namespace std; 
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
    
    registerProcessorParameter ("MaxEtaBins",
                                "Maximum number of eta bins for clusterDB",
                                _maxEtaBins,  static_cast < int > (1));
    
    registerProcessorParameter ("vCellPeriod",
                                "Periodicity for vCells used for clusterDB",
                                _vCellPeriod,  static_cast < int > (1));
     
    registerProcessorParameter ("uCellPeriod",
                                "Periodicity for uCells used for clusterDB",
                                _uCellPeriod,  static_cast < int > (1));
    
    
    std::vector<int> initSelectIDVec;
    registerProcessorParameter ("SelectPlanes",
                                "Select clusters from list of planes",
                                _selectIDVec, initSelectIDVec);
    
  }
  
  //
  // Method called at the beginning of data processing
  //
  void GoeClusterCalibratorForMC::init() {
    
    // Initialize variables
    _nRun = 0 ;
    _nEvt = 0 ;
    _timeCPU = clock()/1000 ;
    
    if ( _maxEtaBins<=0) _maxEtaBins=1;
    
    // Create a useful name for a file that should be deleted after Marlin is finished
    // There should be one collector output for each clusterDB
    // Add prefix "tmp" to indicate that the file is temporary 
    
    std::set<char> delims{'/'};
    _collectorOutputFileName = "tmpCollectorOutputFor_" + splitpath(_clusterDBFileName, delims).back();
    
    _rootCollectorOutputFile = new TFile( _collectorOutputFileName.c_str(),"recreate");
    _rootCollectorOutputFile->cd("");
    
    _trackVarUHisto = new TH1F("trackVarUHisto","trackVarUHisto",1,0,1);
    _trackVarUHisto->StatOverflows(); 	            
    
    _trackVarVHisto = new TH1F("trackVarVHisto","trackVarVHisto",1,0,1); 
    _trackVarVHisto->StatOverflows();         
    
    _trackCovUVHisto = new TH1F("trackCovUVHisto","trackCovUVHisto",1,0,1);
    _trackCovUVHisto->StatOverflows(); 	

    _trackDuDwHisto = new TH1F("trackDuDwHisto","trackDuDwHisto",1,0,1);
    _trackDuDwHisto->StatOverflows(); 
     
    _trackDvDwHisto = new TH1F("trackDvDwHisto","trackDvDwHisto",1,0,1);
    _trackDuDwHisto->StatOverflows();      
    
    m_rootTree = new TTree("tree","Cluster info");
    m_rootTree->Branch<string>("TypeName", &m_typeName);
    m_rootTree->Branch<float>("ClusterEtaPP", &m_clusterEtaPP);
    m_rootTree->Branch<float>("ClusterEtaNP", &m_clusterEtaNP);
    m_rootTree->Branch<float>("ClusterEtaPN", &m_clusterEtaPN);
    m_rootTree->Branch<float>("ClusterEtaNN", &m_clusterEtaNN);
    m_rootTree->Branch<float>("OffsetU", &m_positionOffsetU);
    m_rootTree->Branch<float>("OffsetV", &m_positionOffsetV);
    
    // Print set parameters
    printProcessorParams();
    
   
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
  
    LCCollection* simHitCol = nullptr;
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
      int ipl = TBDetector::GetInstance().GetPlaneNumber(sensorID);
      
      streamlog_out(MESSAGE2) << " SimHit with sensorID " << sensorID << " at: (" << simHit->getPosition()[0] << ", " << simHit->getPosition()[1] << ")" 
                              << endl;
          
      
      bool selectID = false;
      for (auto id :  _selectIDVec)  {
        if  (id == ipl) selectID = true; 
      }
      
      if ( selectID ) SimHitMap[sensorID].push_back(simHit);            
    } 
    
    
    for(int i=0; i< nHit ; i++)
    {
      // Built a TBHit
      TrackerHitImpl * lciohit = dynamic_cast<TrackerHitImpl*>( hitCol->getElementAt(i) ) ;
      TBHit recoHit ( lciohit  );        
      int sensorID = recoHit.GetSensorID();      
      int ipl = TBDetector::GetInstance().GetPlaneNumber(sensorID);
                
      streamlog_out(MESSAGE2) << " RecoHit with sensorID " << sensorID << " at: (" << recoHit.GetCoord()[0] << ", " << recoHit.GetCoord()[1] << ")" 
                              << endl;
        
      bool selectID = false;
      for (auto id :  _selectIDVec)  {
        if  (id == ipl) selectID = true; 
      }

      if ( selectID ) RecoHitMap[sensorID].push_back( recoHit );    
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
              double hitPosU = recoHit.GetCoord()[0];
              double hitPosV = recoHit.GetCoord()[1];
               
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
          int sensorID = hit.GetSensorID();     
          int ipl = TBDetector::GetInstance().GetPlaneNumber(sensorID);
          const Det & Sensor = TBDetector::Get(ipl);   

          // This is the plane number of one plane to which 
          // the clusterDB would be applied
          _setOfPlaneNumbers.insert(ipl);
           
          SimTrackerHit * simHit = SimHits[ hit2simhit[ihit] ]; 
          Vector3d momentum;
          momentum << simHit->getMomentum()[0], simHit->getMomentum()[1],simHit->getMomentum()[2];
          
          // Get local track parameters 
          double trk_tu = momentum[0]/momentum[2];    // rad
          double trk_tv = momentum[1]/momentum[2];    // rad
          double trk_u = simHit->getPosition()[0];    // mm
          double trk_v = simHit->getPosition()[1];    // mm
          //double trk_mom = momentum.norm();            // GeV
           
          _trackDuDwHisto->Fill(trk_tu);  
          _trackDvDwHisto->Fill(trk_tv);  
           
          // FIXME: This cluster is only needed to be handed to the ctor of PolyClusterDescriptor. 
          // Check if its creation can be avoided. 
          PixelCluster Cluster = hit.GetCluster();  
          PolyClusterDescriptor Descriptor(Cluster, Sensor);
          
          // Fill collector output
          
          // A string to identify the cluster type, it quantifies the configuration of firing pixels 
          // but does not using the measured pixel signals. Details depend on the implementation of 
          // the cluster descriptor. 
          m_typeName = Descriptor.getType(_vCellPeriod, _uCellPeriod);
          
          // The eta value is a scalar computed from the pixel charges. It value may depend on the sign of
          // the incidence angle of the beam into the sensor. But details depend on the implementation of 
          // the cluster descriptor. 
          m_clusterEtaPP = Descriptor.computeEta(+1, +1);
          m_clusterEtaPN = Descriptor.computeEta(+1, -1);
          m_clusterEtaNP = Descriptor.computeEta(-1, +1);
          m_clusterEtaNN = Descriptor.computeEta(-1, -1); 
           
          m_positionOffsetU = trk_u - Descriptor.getOriginU();  
          m_positionOffsetV = trk_v - Descriptor.getOriginV();  
          m_rootTree->Fill(); 
        }
      }
    }  
    
    return;
  }
  
  
  //
  // Method called after each event to check the data processed
  //
  void GoeClusterCalibratorForMC::check( LCEvent *  ) {}
  
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
     
    // Make sure that all data was collected for a set of planes 
    // have the same protopixels. 
    bool consistencyTestPassed = true;
    std::vector<int> planeNumbersVec(_setOfPlaneNumbers.begin(), _setOfPlaneNumbers.end()); 
    for(auto i : planeNumbersVec ) {
      for(auto j : planeNumbersVec ) {    
        if ( TBDetector::Get(i).GetProtopixels() ==  TBDetector::Get(j).GetProtopixels() ) {
          streamlog_out(MESSAGE3) << "Plane " << i << " and plane " << j << " have identical map of protopixels." << std::endl;   
        } else {
          streamlog_out(MESSAGE3) << "Plane " << i << " and plane " << j << " have different map of protopixels!" << std::endl;   
          consistencyTestPassed = false;
        } 
      }  
    }  
     
    // Compute the average incidence angles into the sensor
    // We assume a strongly collimated beam and the rms of the 
    // should be small (< few mrad) 
    double thetaU = _trackDuDwHisto->GetMean();
    double thetaV = _trackDvDwHisto->GetMean();
    
    // Enumerate all types by unique name and count their
    // occurence in training data.
    vector< pair<string, float> > typeList;
    
    const auto nEntries = m_rootTree->GetEntries();
    for (int i = 0; i < nEntries; ++i) {
      m_rootTree->GetEntry(i);
      
      auto it = std::find_if(typeList.begin(), typeList.end(),
                             [&](const pair<string, float>& element) { return element.first == m_typeName;});
      
      //Shape name exists in vector
      if (it != typeList.end()) {
        //increment key in map
        it->second++;
      }
      //Shape name does not exist
      else {
        //Not found, insert in vector
        typeList.push_back(pair<string, int>(m_typeName, 1));
      }
    }
    
    // Loop over typeList to select types with enough data for
    // next calibration step
    
    // Vector with eta histograms for selected shapes
    vector< pair<string, TH1D> > etaHistos;
 
    // Coverage of position offsets on training data
    double coverage = 0.0;
    
    for (auto iter : typeList) {
      auto name = iter.first;
      auto counter = iter.second;
      if (counter >=  _minClusters) {
        coverage += counter / nEntries;
        string etaname = "eta_" + name;     
        TH1D etaHisto(etaname.c_str(), etaname.c_str(), 301, 0, 1);
        etaHisto.SetDirectory(0);
        etaHistos.push_back(pair<string, TH1D>(name, etaHisto));
      } else {
        streamlog_out(MESSAGE3) << "  Unable to calibrate cluster type:  " << name << " because too few counts (" << counter  << ")" << endl;
      } 
    }     
    
    // Loop over the tree is to fill the eta histograms for
    // selected shapes.
    for (int i = 0; i < nEntries; ++i) {
      m_rootTree->GetEntry(i);
      auto it = std::find_if(etaHistos.begin(), etaHistos.end(),
                    [&](const pair<string, TH1D>& element) { return element.first == m_typeName;});
      //Item exists in map
      if (it != etaHistos.end()) {
        // increment key in map
        auto clusterEta = m_clusterEtaPP;
        if (thetaU > 0 && thetaV < 0) {clusterEta = m_clusterEtaPN;}
        else if (thetaU < 0 && thetaV > 0) {clusterEta = m_clusterEtaNP;}
        else if (thetaU < 0 && thetaV < 0) {clusterEta = m_clusterEtaNN;}   
        it->second.Fill(clusterEta);
      }
    }
    
    // Vector for eta bin edges stored by type name 
    vector< pair< string, vector<double> > > etaBinEdgesVec;
    
    // Vector for offset histograms storing pairs of typename and a vector of 2d histos for different eta bins
    vector< pair< string, vector<TH2D> > > offsetHistosVec;
    
    for (auto iter : etaHistos) {
      auto name = iter.first;
      auto& histo = iter.second;
      int nClusters = histo.GetEntries();
      
      streamlog_out(MESSAGE3) << "Eta histogram " << name << " has " << nClusters << " entries" << std::endl; 
      
      // Try to split clusters into n bins with _minClusters clusters
      int nEtaBins  = std::max(int(nClusters / _minClusters), 1);
      nEtaBins =  std::min(nEtaBins, _maxEtaBins);  
      
      // We have to check for singular cases where eta distribution is concentrated is concentrated in 
      // less bins than required number of eta bins 
      int nFilledBins = 0;  
      for (auto ibin = 1; ibin <= histo.GetXaxis()->GetNbins(); ibin++) {  
        if (histo.GetBinContent(ibin) > 0) nFilledBins++;  
      }
      if (nFilledBins <= nEtaBins) {
        nEtaBins = 1; 
        streamlog_out(MESSAGE3) << "Eta histogram " << name << " is a delta spike" << std::endl;   
      }
      
      vector<double> etaBinEdges;
      vector< TH2D > offsetHistos;
      
      for (int i = 0; i < nEtaBins; i++) {
        // Position where to compute the quantiles in [0,1]
        double xq = double(i) / nEtaBins;
        // Double to contain the quantile
        double yq = 0;
        histo.GetQuantiles(1, &yq, &xq);
        streamlog_out(MESSAGE3) << " Quantile at xq =" << xq << " is yq=" << yq << std::endl;
        etaBinEdges.push_back(yq);
        
        string offsetname = "E" + std::to_string(i) + name;         
        TH2D offsetHisto(offsetname.c_str(), offsetname.c_str(), 1, 0, 1, 1, 0, 1);
        offsetHisto.StatOverflows();
        offsetHistos.push_back(offsetHisto);
      }
      etaBinEdgesVec.push_back(pair< string, vector<double> >(name, etaBinEdges));
      offsetHistosVec.push_back(pair< string, vector<TH2D> >(name, offsetHistos));
    }
    
    // Loop over the tree is to fill offset histograms
    for (int i = 0; i < nEntries; ++i) {
      m_rootTree->GetEntry(i); 
      
      auto it = std::find_if(offsetHistosVec.begin(), offsetHistosVec.end(),
                        [&](const pair<string, vector<TH2D>>& element) { return element.first == m_typeName;});

      auto it2 = std::find_if(etaBinEdgesVec.begin(), etaBinEdgesVec.end(),
                        [&](const pair<string, vector<double>>& element) { return element.first == m_typeName;});

      //Item exists in maps
      if (it != offsetHistosVec.end()  && it2 != etaBinEdgesVec.end() ) {
        auto clusterEta = m_clusterEtaPP;
        if (thetaU > 0 && thetaV < 0) {clusterEta = m_clusterEtaPN;}
        else if (thetaU < 0 && thetaV > 0) {clusterEta = m_clusterEtaNP;}
        else if (thetaU < 0 && thetaV < 0) {clusterEta = m_clusterEtaNN;}   
        // FIXME add a switch to select which descriptor to use
        auto etaBin = PolyClusterDescriptor::computeEtaBin(clusterEta, it2->second);
        it->second.at(etaBin).Fill(m_positionOffsetU, m_positionOffsetV);
      }
    }
    
    // Count total number of types
    int nTypes = 0;
 
    // Count total number of shapes 
    int nShapes = 0; 
      
    // Compute the moments of the offset histograms 
    for (auto iter : offsetHistosVec) {
      nShapes += iter.second.size();
      nTypes += 1;
    }
        
    // Close collector file
    _rootCollectorOutputFile->Write();
    _rootCollectorOutputFile->Close();
    delete _rootCollectorOutputFile;   
    
    if (nShapes > 0 and consistencyTestPassed) {
      
      streamlog_out(MESSAGE3) << "Create the clusterDB ... " << endl; 
      
      TFile * _rootClusterDBFile = new TFile( _clusterDBFileName.c_str(),"recreate");
      _rootClusterDBFile->cd("");
      
      // Book histograms for clusterDB
      
      string histoName;  
      
      histoName = "hDB_Coverage";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",1,0,1);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("coverage [%]");
      _histoMap[histoName]->SetBinContent( 1, 100 * coverage );
      _histoMap[histoName]->GetXaxis()->SetBinLabel( 1, "cluster found in clusterDB" );
      
      histoName = "hDB_Types";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",nTypes,0,nTypes);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("type");  
      
      histoName = "hDB_Weight";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",nShapes,0,nShapes);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("label weight");  
      
      histoName = "hDB_U";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",nShapes,0,nShapes);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("offset u [mm]");  
      
      histoName = "hDB_V";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",nShapes,0,nShapes);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("offset v [mm]");  
      
      histoName = "hDB_Sigma2_U";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",nShapes,0,nShapes);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("sigma2 offset u [mm^2]");        
      
      histoName = "hDB_Sigma2_V";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",nShapes,0,nShapes);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("sigma2 offset v [mm^2]"); 
      
      histoName = "hDB_Cov_UV";
      _histoMap[histoName] = new TH1F(histoName.c_str(),"",nShapes,0,nShapes);
      _histoMap[histoName]->SetStats( false );
      _histoMap[histoName]->SetYTitle("covariance u-v [mm^2]"); 
      
      
      // Compute the moments of the offset histograms 
      int offsetBin = 0;
      int typeBin = 0; 
      for (auto iter : offsetHistosVec) {
        typeBin++;
        auto name = iter.first;
        auto& histovec = iter.second;
        int typeCounter = 0; 
        
        // Loop over eta bins
        int etaBin = -1; 
        for (auto& histo : histovec) {
          // Compute offset moments
          offsetBin++;
          etaBin++; 
          auto shapeName = "E" + std::to_string(etaBin) + name; 
          int counter = histo.GetEntries();
          typeCounter += counter;
          double offsetU = histo.GetMean(1);
          double offsetV = histo.GetMean(2);
          double covUV = histo.GetCovariance();
          double covU = pow(histo.GetRMS(1), 2);
          double covV = pow(histo.GetRMS(2), 2);
          
          streamlog_out(MESSAGE3) << "Name " << shapeName  << " entries=" << counter << ", posU=" << offsetU << ", posV=" << offsetV 
                                  << ", sigmaU=" << sqrt(covU) << ", sigmaV=" << sqrt(covV) << ", corrUV=" << covUV/sqrt(covU)/sqrt(covV) << std::endl;
          
          TMatrixDSym HitCov(2);
          HitCov(0, 0) = covU;
          HitCov(1, 0) = covUV;
          HitCov(0, 1) = covUV;
          HitCov(1, 1) = covV;
          
          TMatrixDSymEigen HitCovE(HitCov);
          TVectorD eigenval = HitCovE.GetEigenValues();
          if (eigenval(0) <= 0 || eigenval(1) <= 0) {
            streamlog_out(MESSAGE3) << "Estimated covariance matrix not positive definite." << std::endl;
          }
          
          // Store calibration result   
          histoName = "hDB_Weight";
          _histoMap[histoName]->SetBinContent( offsetBin, counter );
          _histoMap[histoName]->GetXaxis()->SetBinLabel( offsetBin, shapeName.c_str() );
          
          histoName = "hDB_U";
          _histoMap[histoName]->SetBinContent( offsetBin, offsetU );
          _histoMap[histoName]->SetBinError( offsetBin, histo.GetMeanError(1) );
          _histoMap[histoName]->GetXaxis()->SetBinLabel( offsetBin, shapeName.c_str() );
          
          histoName = "hDB_V"; 
          _histoMap[histoName]->SetBinContent( offsetBin, offsetV );
          _histoMap[histoName]->SetBinError( offsetBin, histo.GetMeanError(2) );
          _histoMap[histoName]->GetXaxis()->SetBinLabel( offsetBin, shapeName.c_str() );
               
          histoName = "hDB_Sigma2_U";
          _histoMap[histoName]->SetBinContent( offsetBin, covU );
          _histoMap[histoName]->SetBinError( offsetBin, 2 * histo.GetRMS(1) * histo.GetRMSError(1) );
          _histoMap[histoName]->GetXaxis()->SetBinLabel( offsetBin, shapeName.c_str() );
          
          histoName = "hDB_Sigma2_V";
          _histoMap[histoName]->SetBinContent( offsetBin, covV );
          _histoMap[histoName]->SetBinError( offsetBin, 2 * histo.GetRMS(2) * histo.GetRMSError(2) );
          _histoMap[histoName]->GetXaxis()->SetBinLabel( offsetBin, shapeName.c_str() );  
        
          histoName = "hDB_Cov_UV";
          _histoMap[histoName]->SetBinContent( offsetBin, covUV );
          _histoMap[histoName]->SetBinError( offsetBin, 0 );
          _histoMap[histoName]->GetXaxis()->SetBinLabel( offsetBin, shapeName.c_str() );  
        }

        // Add cluster type 
        histoName = "hDB_Types";
        _histoMap[histoName]->SetBinContent( typeBin, typeCounter );
        _histoMap[histoName]->GetXaxis()->SetBinLabel( typeBin, name.c_str() ); 
      }
      
      for (auto iter : etaBinEdgesVec) {
        auto name = iter.first;
        auto& etaBinEdges = iter.second;
        
        // Add eta bin edges for type
        TVectorD DB_etaBinEdges( etaBinEdges.size() );
        for ( size_t iBin=0; iBin < etaBinEdges.size(); iBin++) {
          DB_etaBinEdges[iBin] = etaBinEdges[iBin];
        }
        DB_etaBinEdges.Write(  string("DB_etaBinEdges_"+name).c_str() );
      }
               
      TVectorD DB_periods( 2 );
      DB_periods[0] = _vCellPeriod;
      DB_periods[1] = _uCellPeriod;
      DB_periods.Write("DB_periods");
      
      TVectorD DB_angles( 2 );
      DB_angles[0] = thetaU;
      DB_angles[1] = thetaV;
      DB_angles.Write("DB_angles");
      
      TVectorD DB_telcov( 3 );
      DB_telcov[0] = 0.001;
      DB_telcov[1] = 0.001;
      DB_telcov[2] = 0;
      DB_telcov.Write("DB_telcov");
      
      if (not planeNumbersVec.empty()) {
        // Loop over all protopixel for one sensor to which the cluseterDB is going to 
        // be used. 
        for (auto protopixel : TBDetector::Get(planeNumbersVec[0]).GetProtopixels())  {
          auto pixeltype = protopixel.first;
          auto points = protopixel.second;
          
          // Add eta bin edges for type
          TVectorD DB_protopixel( 2*points.size() );
          int iBin = 0;
          for (auto point: points){
            DB_protopixel[iBin] = std::get<0>(point);
            DB_protopixel[iBin+1] = std::get<1>(point);
            iBin+=2;
          }
          DB_protopixel.Write( string("DB_protopixel_"+std::to_string(pixeltype)).c_str() ); 
        }
      } else {
        streamlog_out(ERROR3) << "No information about protopixels is available. This is strange! "  << endl;   
      }
      
      streamlog_out(MESSAGE3) << "Created clusterDB with coverage " << 100 * coverage << " percent on training data sample." << std::endl; 
      
      streamlog_out(MESSAGE3) << "ClusterDB written to file " << _clusterDBFileName << std::endl; 
      
      // Close clusterDB  file
      _rootClusterDBFile->Write();
      _rootClusterDBFile->Close();
      delete _rootClusterDBFile;
      
    }
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



