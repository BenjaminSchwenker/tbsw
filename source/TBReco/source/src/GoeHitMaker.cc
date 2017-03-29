// GoeHitMaker - Marlin Processor
// 
// Compute hits from pixel clusters using clusterDB
//
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 
#include "GoeHitMaker.h"
#include "TBHit.h"

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
  GoeHitMaker::GoeHitMaker() : Processor("GoeHitMaker")
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
    
    registerProcessorParameter ("MaxSizeU", 
                                "Clusters having more u cells marked as bad [and can be filterd in track finder]",
                                _maxSizeU,  
                                static_cast < int > (3));
    
    registerProcessorParameter ("MaxSizeV", 
                                "Clusters having more v cells marked as bad [and can be filterd in track finder]",
                                _maxSizeV,  
                                static_cast < int > (3));
     
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
    
    // Read detector constants from gear file
    _detector.ReadGearConfiguration();    
    
    // Open clusterDB file 
    TFile * clusterDBFile = new TFile(_clusterDBFileName.c_str(), "READ");
    
    // Read calibration data  
    string histoName;

    histoName = "hDB_ID";
    if ( (TH1F *) clusterDBFile->Get(histoName.c_str()) != nullptr) {
      m_DB_Weight = (TH1F *) clusterDBFile->Get(histoName.c_str());  
      m_DB_Weight->SetDirectory(0);
    } 
       
    histoName = "hDB_U";
    if ((TH1F *) clusterDBFile->Get(histoName.c_str()) != nullptr) {
      m_DB_U = (TH1F *) clusterDBFile->Get(histoName.c_str());
      m_DB_U->SetDirectory(0);
    } 
      
    histoName = "hDB_V"; 
    if ((TH1F *) clusterDBFile->Get(histoName.c_str()) != nullptr) {
      m_DB_V = (TH1F *) clusterDBFile->Get(histoName.c_str());
      m_DB_V->SetDirectory(0);
    } 
      
    histoName = "hDB_Sigma2_U";
    if ((TH1F *) clusterDBFile->Get(histoName.c_str()) != nullptr) {
      m_DB_Sigma2_U = (TH1F *) clusterDBFile->Get(histoName.c_str());
      m_DB_Sigma2_U->SetDirectory(0);
    } 
      
    histoName = "hDB_Sigma2_V";
    if ((TH1F *) clusterDBFile->Get(histoName.c_str()) != nullptr) {
      m_DB_Sigma2_V = (TH1F *) clusterDBFile->Get(histoName.c_str());
      m_DB_Sigma2_V->SetDirectory(0);
    } 
      
    histoName = "hDB_Cov_UV";
    if ((TH1F *) clusterDBFile->Get(histoName.c_str()) != nullptr) {
      m_DB_Cov_UV = (TH1F *) clusterDBFile->Get(histoName.c_str());
      m_DB_Cov_UV->SetDirectory(0);
    } 
    
    histoName = "DB_Scale";
    if ((TVectorD*) clusterDBFile->Get(histoName.c_str()) != nullptr) {
      TVectorD *DB_Scale = (TVectorD*) clusterDBFile->Get( histoName.c_str() );
      m_scale = (*DB_Scale)[0];  
    } else {
      m_scale = 1;
    } 
    
    // Close root  file
    clusterDBFile->Close();
    delete clusterDBFile;

    for(int ipl=0;ipl<_detector.GetNSensors();ipl++)  { 
      int sensorID = _detector.GetDet(ipl).GetDAQID();
      _countAllMap[sensorID] = 0;   
      _countBadMap[sensorID] = 0;   
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
    LCCollectionVec* clusterCollection;
    try {
      clusterCollection = dynamic_cast < LCCollectionVec * >  ( evt->getCollection(_clusterCollectionName) );
    } catch (DataNotAvailableException& e) {
      throw SkipEventException(this);
    }
        
    // Create hit collection  
    LCCollectionVec * hitCollection = new LCCollectionVec(LCIO::TRACKERHIT) ;
         
    // Loop on all clusters 
    for (unsigned int iClu = 0; iClu < clusterCollection->size(); iClu++) 
    {
      // Helper class for decoding cluster ID's 
      CellIDDecoder<TrackerPulseImpl> clusterDecoder( clusterCollection ); 
      
      // Read cluster header
      TrackerPulseImpl* cluster = dynamic_cast<TrackerPulseImpl* > ( clusterCollection->getElementAt(iClu) )  ;       
      int sensorID = clusterDecoder(cluster)["sensorID"]; 
      int ipl = _detector.GetPlaneNumber(sensorID);
      Det& Det = _detector.GetDet(ipl);
      
      // Increment the cluster counter
      _countAllMap[sensorID]++;
      
      // Compute the cluster ID string
      PixelCluster aCluster(cluster->getTrackerData());    
      string id = aCluster.getLabel(m_scale); 
      
      streamlog_out(MESSAGE2) << "Processing cluster on sensorID " << sensorID << " with label " << id << endl; 
      
      // Compute position measurement and its 2x2 covariance matrix   
      double u{0.0}, v{0.0}, sig2_u{0.0}, sig2_v{0.0}, cov_uv{0.0};
      int quality = 0; 
       
      bool found = searchDB(sensorID, id, u, v, sig2_u, sig2_v, cov_uv); 
      if (found) {
        // Count matched clusters
        _countCalMap[sensorID]++; 
        // Shift position into local sensor coordinates
        u += Det.GetPixelCenterCoordU( aCluster.getVStart(), aCluster.getUStart()); 
        v += Det.GetPixelCenterCoordV( aCluster.getVStart(), aCluster.getUStart()); 
      } else { 
        // Fallback for computing position measurement
        double ustart = Det.GetPixelCenterCoordU( aCluster.getVStart(), aCluster.getUStart());
        double ustop  = Det.GetPixelCenterCoordU( aCluster.getVStart(), aCluster.getUStart() + aCluster.getUSize()-1);
        u       = 0.5*( ustart + ustop );
        double vstart = Det.GetPixelCenterCoordV( aCluster.getVStart(), aCluster.getUStart());
        double vstop  = Det.GetPixelCenterCoordV( aCluster.getVStart() + aCluster.getVSize()-1 , aCluster.getUStart());
        v       = 0.5*( vstart + vstop ); 
        sig2_u  = std::pow(aCluster.getUSize()*Det.GetPitchU()/TMath::Sqrt(12),2);  
        sig2_v  = std::pow(aCluster.getVSize()*Det.GetPitchV()/TMath::Sqrt(12),2);   
        cov_uv  = 0; 
      }
      
      // Too large clusters are marked as bad and can be ignored 
      // in the track finder stage
      if ( aCluster.getUSize() > _maxSizeU  ||  aCluster.getVSize() > _maxSizeV ) {
        quality = 1; 
        _countBadMap[sensorID]++;   
      }
      
      // Make TBHit 
      TBHit hit(sensorID, u, v, sig2_u, sig2_v, cov_uv, quality);
      
      // Make LCIO TrackerHit
      TrackerHitImpl * trackerhit = hit.MakeLCIOHit();  
            
      // Add link to full cluster data 
      LCObjectVec clusterVec;
      clusterVec.push_back( cluster->getTrackerData() );
      trackerhit->rawHits() = clusterVec;
      
      // Add hit to the hit collection
      hitCollection->push_back( trackerhit );
      
      streamlog_out(MESSAGE2) << "  Label " << id << " found: " << endl
                                << std::setiosflags(std::ios::fixed | std::ios::internal )
                                << std::setprecision(8)
                                << "  offset u: " << u << ", sigma: " << TMath::Sqrt(sig2_u) << endl
                                << "  offset v: " << v << ", sigma: " << TMath::Sqrt(sig2_v) << endl
                                << "  cov(u,v): " << cov_uv
                                << std::setprecision(3)
                                << endl;   
          
    } // End cluster loop 
      
    // Store hitCollection in LCIO file
    evt->addCollection( hitCollection, _hitCollectionName );      
  }

  //
  // Method called after each event to check the data processed
  //
  void GoeHitMaker::check( LCEvent * evt )
  {
  }

  //
  // Method called after all data processing
  //
  void GoeHitMaker::end()
  {
    
    for(int ipl=0;ipl<_detector.GetNSensors();ipl++)  { 
      int sensorID = _detector.GetDet(ipl).GetDAQID();
      
      // Print message
      streamlog_out(MESSAGE3) << std::endl
                              << " "
                              << "ClusterDB coverage on sensorID " << sensorID << ": "
                              << std::setprecision(4)
                              << 100.0*((float)_countCalMap[sensorID]/_countAllMap[sensorID])
                              << " %"
                              << std::endl
                              << "Fraction of bad clusters on sensorID " << sensorID << ": "
                              << std::setprecision(4)
                              << 100.0*((float)_countBadMap[sensorID]/_countAllMap[sensorID])
                              << " %"
                              << std::endl    
                              << std::endl;
    }
     
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
                            << std::setprecision(3)
                            << std::endl
                            << " "
                            << "Processor succesfully finished!"
                            << std::endl;
  }
  
  bool GoeHitMaker::searchDB(int sensorID, string id, double& u, double& v, double& sig2_u, double& sig2_v, double& cov_uv)
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

