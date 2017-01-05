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
                            _clusterCollectionName ,
                            std::string("cluster") ) ;
    
    registerOutputCollection (LCIO::TRACKERHIT, "HitCollectionName",
                             "Name of hit collection",
                             _hitCollectionName, 
                             string("hit"));
    
    registerProcessorParameter ("ClusterDBFileName",
                                "This is the name of the LCIO file with the alignment constants (add .slcio)",
                                _clusterDBFileName, static_cast< string > ( "ClusterDB.root" ) ); 
   
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
    _clusterDBFile = new TFile(_clusterDBFileName.c_str(), "UPDATE");
    if (_clusterDBFile == 0 || _clusterDBFile->IsOpen() != kTRUE) {
      streamlog_out ( ERROR4) << "Could not open clusterDB file!!" << endl;
      exit(-1);
    }
    
    // Get histograms storing position offsets 
    _DB_U  = (TH1F *) _clusterDBFile->Get("hDB_U");
    if (_DB_U == 0) {
      streamlog_out ( ERROR4) << "Could not find cluster offset histograms in file " << _clusterDBFileName << endl;
      exit(-1);  
    }
    _DB_V  = (TH1F *) _clusterDBFile->Get("hDB_V");
    if (_DB_U == 0) {
      streamlog_out ( ERROR4) << "Could not find cluster offset histograms in file " << _clusterDBFileName << endl;
      exit(-1);  
    }
   
    // Get histograms storing components of covariance matrix  
    _DB_Sigma_U  = (TH1F *) _clusterDBFile->Get("hDB_Sigma_U");
    if (_DB_Sigma_U == 0) {
      streamlog_out ( ERROR4) << "Could not find cluster covariance histograms in file " << _clusterDBFileName << endl;
      exit(-1);  
    }
    _DB_Sigma_V  = (TH1F *) _clusterDBFile->Get("hDB_Sigma_V");
    if (_DB_Sigma_V == 0) {
      streamlog_out ( ERROR4) << "Could not find cluster covariance histograms in file " << _clusterDBFileName << endl;
      exit(-1);  
    }
    _DB_Cov_UV  = (TH1F *) _clusterDBFile->Get("hDB_Cov_UV");
    if (_DB_Cov_UV == 0) {
      streamlog_out ( ERROR4) << "Could not find cluster covariance histograms in file " << _clusterDBFileName << endl;
      exit(-1);  
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
      
      // Compute the cluster ID string
      PixelCluster aCluster(cluster->getTrackerData());   
      string id = (string) aCluster; 
      
      streamlog_out(MESSAGE2) << "Processing cluster on sensorID " << sensorID << " with clusterID " << id << endl; 
      
      // Compute hit position and 2x2 covariance matrix 
      double originPosU = Det.GetPixelCenterCoordU( aCluster.getVStart(), aCluster.getUStart()); 
      double originPosV = Det.GetPixelCenterCoordV( aCluster.getVStart(), aCluster.getUStart()); 
      double u{originPosU}, v{originPosV}, sig2_u{0.0}, sig2_v{0.0}, cov_uv{0.0}; 
         
      // Lookup hash for clusterID
      int bin = _DB_U->GetXaxis()->FindFixBin(id.c_str());
      
      if (bin == -1) {
        streamlog_out(MESSAGE2) << "  ClusterId " << id << " not found in clusterDB. Ignoring cluster." << endl;
        continue;
      }
      
      streamlog_out(MESSAGE2) << "  ClusterId " << id << " found at bin " << bin << endl
                                << std::setiosflags(std::ios::fixed | std::ios::internal )
                                << std::setprecision(8)
                                << "  offset u: " << _DB_U->GetBinContent(bin) << ", sigma: " << _DB_Sigma_U->GetBinContent(bin) << endl
                                << "  offset v: " << _DB_V->GetBinContent(bin) << ", sigma: " << _DB_Sigma_V->GetBinContent(bin) << endl
                                << "  cov(u,v): " << _DB_Cov_UV->GetBinContent(bin)
                                << std::setprecision(3)
                                << endl;   
                            
      // Lookup offsets from clusterDB 
      u += _DB_U->GetBinContent(bin); 
      v += _DB_V->GetBinContent(bin); 
          
      // Lookup sigmas and covariance from clusterDB
      sig2_u = pow(_DB_Sigma_U->GetBinContent(bin),2);     
      sig2_v = pow(_DB_Sigma_V->GetBinContent(bin),2);    
      cov_uv = _DB_Cov_UV->GetBinContent(bin);
      
      TBHit hit(sensorID, u, v, sig2_u, sig2_v, cov_uv, 0);
      
      // Make LCIO TrackerHit
      TrackerHitImpl * trackerhit = hit.MakeLCIOHit();  
            
      // Add link to full cluster data 
      LCObjectVec clusterVec;
      clusterVec.push_back( cluster->getTrackerData() );
      trackerhit->rawHits() = clusterVec;
      
      // Add hit to the hit collection
      hitCollection->push_back( trackerhit );
          
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
  void GoeHitMaker::printProcessorParams() const
  {
  
    streamlog_out(MESSAGE3)  << std::endl
                             << " "
                             << "GoeHitMaker Development Version, be carefull!!"
                             << " "
                             << std::endl  << std::endl;   

  }

} // Namespace

