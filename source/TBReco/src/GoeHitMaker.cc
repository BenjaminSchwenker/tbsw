// GoeHitMaker - Marlin Processor
// 
// Compute hits from pixel clusters using clusterDB
//
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 
#include "GoeHitMaker.h"
#include "TBHit.h"
#include "TBDetector.h"
#include "PolyClusterDescriptor.h"

// Include basic C
#include <limits>
#include <cmath>
#include <iomanip>

// ROOT includes
#include <TMath.h>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerPulseImpl.h>

// Include ROOT classes
#include <TFile.h>
#include <TVectorD.h>
using namespace std::string_literals;

// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;

namespace depfet {

  //
  // Instantiate this object
  //
  GoeHitMaker aGoeHitMaker ;
  
  //
  // Constructor
  //
  GoeHitMaker::GoeHitMaker() : Processor("GoeHitMaker"),_inputDecodeHelper("")
  {
    
    // Processor description
    _description = "GoeHitMaker:  Compute hits from pixel clusters using a clusterDB";
    
    //   
    // First of all we need to register the input/output collections
    
    registerInputCollection( LCIO::TRACKERPULSE,
                             "ClusterCollection" ,
                             "Name of cluster collection"  ,
                             _clusterCollectionName , std::string("cluster") ) ;
    
    registerOutputCollection(LCIO::TRACKERHIT, "HitCollectionName",
                             "Name of hit collection",
                             _hitCollectionName, 
                             string("hit"));
    
    registerProcessorParameter("ClusterDBFileName",
                               "This is the name of the ROOT file with the cluster constants (add .root)",
                               _clusterDBFileName, static_cast< string > ( "ClusterDB.root" ) ); 
    
    registerProcessorParameter("UseCenterOfGravityFallback",
                               "Set true to use center of gravity hit when no correction available from clusterDB",
                               _useCoGFallback, static_cast< bool > ( true ) ); 
     
    std::vector<float> initSigmaUCorrrections;
    initSigmaUCorrrections.push_back(1.0);
    registerProcessorParameter ("SigmaUCorrections", 
                                "List of correction factors for sigma U for sizeU=1,2,... . Sigma U will be computed as factor*uLength/sqrt(12). Defaults to factor=1.",
                                _sigmaUCorrections, initSigmaUCorrrections); 
    
    std::vector<float> initSigmaVCorrections;
    initSigmaVCorrections.push_back(1.0);
    registerProcessorParameter ("SigmaVCorrections", 
                                "List of correction factors for sigma V for sizeV=1,2,... . Sigma V will be computed as factor*vLength/sqrt(12). Defaults to factor=1.",
                                _sigmaVCorrections, initSigmaVCorrections); 
     
  }
  
  //
  // Method called at the beginning of data processing
  //
  void GoeHitMaker::init() {
    
    // Initialize variables
    _nRun = 0 ;
    _nEvt = 0 ;
    _timeCPU = clock()/1000;
                 
    // Print set parameters
    printProcessorParams();
    
    
    
    // Open clusterDB file 
    TFile * clusterDBFile = new TFile(_clusterDBFileName.c_str(), "READ");
    
    if (clusterDBFile->IsOpen()) {  
      // Read calibration data  
      string histoName;
      
      histoName = "hDB_Weight";
      if ( (TH1F *) clusterDBFile->Get(histoName.c_str()) != nullptr) {
        m_DB_Weight = (TH1F *) clusterDBFile->Get(histoName.c_str());  
        m_DB_Weight->SetDirectory(nullptr);
      } 
       
      histoName = "hDB_U";
      if ((TH1F *) clusterDBFile->Get(histoName.c_str()) != nullptr) {
        m_DB_U = (TH1F *) clusterDBFile->Get(histoName.c_str());
        m_DB_U->SetDirectory(nullptr);
      } 
      
      histoName = "hDB_V"; 
      if ((TH1F *) clusterDBFile->Get(histoName.c_str()) != nullptr) {
        m_DB_V = (TH1F *) clusterDBFile->Get(histoName.c_str());
        m_DB_V->SetDirectory(nullptr);
      } 
      
      histoName = "hDB_Sigma2_U";
      if ((TH1F *) clusterDBFile->Get(histoName.c_str()) != nullptr) {
        m_DB_Sigma2_U = (TH1F *) clusterDBFile->Get(histoName.c_str());
        m_DB_Sigma2_U->SetDirectory(nullptr);
      } 
      
      histoName = "hDB_Sigma2_V";
      if ((TH1F *) clusterDBFile->Get(histoName.c_str()) != nullptr) {
        m_DB_Sigma2_V = (TH1F *) clusterDBFile->Get(histoName.c_str());
        m_DB_Sigma2_V->SetDirectory(nullptr);
      } 
      
      histoName = "hDB_Cov_UV";
      if ((TH1F *) clusterDBFile->Get(histoName.c_str()) != nullptr) {
        m_DB_Cov_UV = (TH1F *) clusterDBFile->Get(histoName.c_str());
        m_DB_Cov_UV->SetDirectory(nullptr);
      } 
    
      histoName = "DB_periods";
      if ((TVectorD*) clusterDBFile->Get(histoName.c_str()) != nullptr) {
        TVectorD *DB_periods = (TVectorD*) clusterDBFile->Get( histoName.c_str() );
        _vCellPeriod  = (*DB_periods)[0] ; 
        _uCellPeriod  = (*DB_periods)[1] ; 
      } 
      
      histoName = "DB_angles";
      if ((TVectorD*) clusterDBFile->Get(histoName.c_str()) != nullptr) {
        TVectorD *DB_angles = (TVectorD*) clusterDBFile->Get( histoName.c_str() );
        _thetaU  = (*DB_angles)[0]; 
        _thetaV  = (*DB_angles)[1]; 
      } 

      histoName = "hDB_Types";
      if ((TH1F *) clusterDBFile->Get(histoName.c_str()) != nullptr) {
        m_DB_Types = (TH1F *) clusterDBFile->Get(histoName.c_str());
        m_DB_Types->SetDirectory(nullptr);
      } 
    
      // We need the eta bin edges for all cluster types  
      for (auto i = 1; i <= m_DB_Types->GetXaxis()->GetNbins(); i++) {    
        string typeName =  m_DB_Types->GetXaxis()->GetBinLabel(i);
        histoName = "DB_etaBinEdges_" + typeName;
        if ((TVectorD*) clusterDBFile->Get(histoName.c_str()) != nullptr) {
          TVectorD *DB_edges = (TVectorD*) clusterDBFile->Get( histoName.c_str() );
          for(auto i = 0; i < DB_edges->GetNrows(); i++ ) {
            m_etaBinEdgesMap[typeName].push_back( (*DB_edges)[i] ); 
          }  
        }    
      }
    }
    
    // Close root  file
    clusterDBFile->Close();
    // delete pointer
    delete clusterDBFile;
    
    for(int ipl=0;ipl<TBDetector::GetInstance().GetNSensors();ipl++)  { 
      int sensorID = TBDetector::Get(ipl).GetSensorID();
      _countAllMap[sensorID] = 0;   
      _countCalMap[sensorID] = 0;  
    }
  }
  
  //
  // Method called for each run
  //
  void GoeHitMaker::processRunHeader(LCRunHeader * run)
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
  void GoeHitMaker::processEvent(LCEvent * evt)
  {

    _nEvt ++ ;
    
    // Print event number
    if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                << (evt->getEventNumber())
                                                                << std::endl << std::endl;
    
    

    //
    // Open collections
    try {

      LCCollectionVec* clusterCollection = dynamic_cast < LCCollectionVec * >  ( evt->getCollection(_clusterCollectionName) );
      
      // Create hit collection  
      LCCollectionVec * hitCollection = new LCCollectionVec(LCIO::TRACKERHIT) ;

      // Helper class for decoding cluster ID's
      CellIDDecoder<TrackerPulseImpl> clusterDecoder( clusterCollection,&_inputDecodeHelper);

      // Loop on all clusters 
      for (unsigned int iClu = 0; iClu < clusterCollection->size(); iClu++) 
      {

        // Read cluster header


        TrackerPulseImpl* cluster = dynamic_cast<TrackerPulseImpl* > ( clusterCollection->getElementAt(iClu) )  ;       

        static auto idx_clust_sensorID=clusterDecoder(cluster).index("sensorID"s); //find the address ONCE.
        int sensorID = clusterDecoder(cluster)[idx_clust_sensorID];
        int ipl = TBDetector::GetInstance().GetPlaneNumber(sensorID);
        const Det& Det = TBDetector::Get(ipl);

        
        // Increment the cluster counter
        _countAllMap[sensorID]++;
       
        PixelCluster aCluster(cluster->getTrackerData());  

        // Compute the cluster ID string
        // FIXME add a switch to choose among descriptors
        PolyClusterDescriptor Descriptor(aCluster, Det);
        
        // A string to identify the cluster type, it quantifies the configuration of firing pixels 
        // but does not using the measured pixel signals. Details depend on the implementation of 
        // the cluster descriptor. 
        string typeName = Descriptor.getType(_vCellPeriod, _uCellPeriod);
        
        // The eta value is a scalar computed from the pixel charges. It value may depend on the sign 
        // of the incidence angles of the beam into the sensor. But details depend on the implementation of 
        // the cluster descriptor. 
        double eta = Descriptor.computeEta(_thetaU, _thetaV);
        
        // The eta values are quantized. The optimal bin edges have been computed during the cluster 
        // calibration and are read back from the clusterDB.  
        int etaBin = PolyClusterDescriptor::computeEtaBin(eta, m_etaBinEdgesMap[typeName]);

        // A string to identify the cluster shape, including the information from analog pixel charges. 
        string shapeName = Descriptor.getEtaBinString(etaBin)+typeName;
        
        streamlog_out(MESSAGE2) << "Processing cluster on sensorID " << sensorID << " with shape " << shapeName << endl; 
        
        // Compute position measurement and its 2x2 covariance matrix   
        double u{0.0}, v{0.0}, sig2_u{0.0}, sig2_v{0.0}, cov_uv{0.0};
        int quality = 0; 
        
        bool found = searchDB(sensorID, shapeName, u, v, sig2_u, sig2_v, cov_uv); 
        if (found) {
          // Count matched clusters
          _countCalMap[sensorID]++; 
          // Shift position into local sensor coordinates
          u += Descriptor.getOriginU();  
          v += Descriptor.getOriginV();  
          
          streamlog_out(MESSAGE2) << "  Shape " << shapeName << " found: " << endl
                                << std::setprecision(8)
                                << "  u: " << u << ", sigmaU: " << TMath::Sqrt(sig2_u) << endl
                                << "  v: " << v << ", sigmaV: " << TMath::Sqrt(sig2_v) << endl
                                << "  corr(u,v): " << cov_uv/TMath::Sqrt(sig2_u)/TMath::Sqrt(sig2_v)
                                << std::setprecision(3)
                                << endl; 
        
          // Make TBHit 
          TBHit hit(sensorID, u, v, sig2_u, sig2_v, cov_uv, quality);
        
          // Make LCIO TrackerHit
          TrackerHitImpl * trackerhit = hit.MakeLCIOHit();  
            
          // Add link to full cluster data 
          LCObjectVec clusterVec;
          clusterVec.push_back( cluster->getTrackerData() );
          trackerhit->rawHits() = clusterVec;
          
          //  Add hit to the hit collection
          hitCollection->push_back( trackerhit );
             
        } else if (_useCoGFallback)  {
           
          aCluster.getCenterOfGravity(Det, u, v, sig2_u, sig2_v, cov_uv); 
          
          // Override sigma u from user input
          if ( aCluster.getUSize()-1 < int(_sigmaUCorrections.size()) ) {
            sig2_u *= pow(_sigmaUCorrections[aCluster.getUSize()-1],2);  
          }
           
          // Override sigma v from user input    
          if ( aCluster.getVSize()-1 < int(_sigmaVCorrections.size()) ) {
            sig2_v *= pow(_sigmaVCorrections[aCluster.getVSize()-1],2); 
          }        
          
          // Mark as fallback
          quality = 1;  
          
          // Make TBHit 
          TBHit hit(sensorID, u, v, sig2_u, sig2_v, cov_uv, quality);
          
          // Make LCIO TrackerHit
          TrackerHitImpl * trackerhit = hit.MakeLCIOHit();  
            
          // Add link to full cluster data 
          LCObjectVec clusterVec;
          clusterVec.push_back( cluster->getTrackerData() );
          trackerhit->rawHits() = clusterVec;
          
          //  Add hit to the hit collection
          hitCollection->push_back( trackerhit );
        }
               
      } // End cluster loop 
       
      // Store hitCollection in LCIO file
      evt->addCollection( hitCollection, _hitCollectionName );     
       
    } catch (DataNotAvailableException& e) {
      streamlog_out(MESSAGE2) << "Missing cluster collection in event: " << evt->getEventNumber() << std::endl;
    }
  }

  //
  // Method called after each event to check the data processed
  //
  void GoeHitMaker::check( LCEvent * )
  {
  }

  //
  // Method called after all data processing
  //
  void GoeHitMaker::end()
  {
    
    for(int ipl=0;ipl<TBDetector::GetInstance().GetNSensors();ipl++)  { 
      int sensorID = TBDetector::Get(ipl).GetSensorID();
      
      float coverage_efficiency = 100.0*((float)_countCalMap[sensorID]/_countAllMap[sensorID]);
      

      // Print message
      streamlog_out(MESSAGE3) << std::endl
                              << " "
                              << "ClusterDB coverage efficiency on sensorID " << sensorID << ": "
                              << std::setprecision(4)
                              <<  coverage_efficiency
                              << " %"
                              << std::endl
                              << std::endl;
    }
     
    // CPU time end
    _timeCPU = clock()/1000 - _timeCPU;
     
    // Print message
    streamlog_out(MESSAGE) << std::endl
                            << " "
                            << "Time per event: "
                            << std::setiosflags(std::ios::fixed | std::ios::internal )
                            << std::setprecision(3)
                            << _timeCPU/_nEvt
                            << " ms"
                            << std::setprecision(3)
                            << std::endl
                            << " "
                            << "Processor succesfully finished!"
                            << std::endl;
  }
  
  bool GoeHitMaker::searchDB(int /*sensorID*/, string id, double& u, double& v, double& sig2_u, double& sig2_v, double& cov_uv)
  {
   
    if ( m_DB_Weight == nullptr ) {
      return false;  
    }
     
    int bin = m_DB_Weight->GetXaxis()->FindFixBin(id.c_str());
    if (bin == -1) {
      return false;
    }
    
    // Get calibrated measurement 
    u      =  m_DB_U->GetBinContent(bin);
    v      =  m_DB_V->GetBinContent(bin);
    sig2_u =  m_DB_Sigma2_U->GetBinContent(bin);
    sig2_v =  m_DB_Sigma2_V->GetBinContent(bin);
    cov_uv =  m_DB_Cov_UV->GetBinContent(bin);
    
    if (sig2_u <=0 || sig2_v <= 0) {
      return false;
    } 
    
    return true;
  }
  
  //
  // Method printing processor parameters
  //
  void GoeHitMaker::printProcessorParams() const
  {
  
    streamlog_out(MESSAGE3)  << std::endl
                             << " "
                             << "GoeHitMaker Development Version, be carefull!!"
                             << " "
                             << std::endl  << std::endl;   
  }

} // Namespace

