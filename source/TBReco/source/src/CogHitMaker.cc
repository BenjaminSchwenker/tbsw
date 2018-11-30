// CogHitMaker - Marlin Processor
// 
// Compute center of gravity hits from clusters 
//
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 
#include "CogHitMaker.h"
#include "TBHit.h"

// Include basic C
#include <limits>
#include <cmath>
#include <iomanip>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerPulseImpl.h>

// ROOT includes
#include <TMath.h>



// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;

namespace depfet {

//
// Instantiate this object
//
CogHitMaker aCogHitMaker ;

//
// Constructor
//
CogHitMaker::CogHitMaker() : Processor("CogHitMaker")
{

// Processor description
   _description = "CogHitMaker:  Compute center of gravity hits from pixel clusters ";
   
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
void CogHitMaker::init() {
    
   // Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
   _timeCPU = clock()/1000;
                 
   // Print set parameters
   printProcessorParams();
   
   // Read detector constants from gear file
   _detector.ReadGearConfiguration();    
    
}

//
// Method called for each run
//
void CogHitMaker::processRunHeader(LCRunHeader * run)
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
void CogHitMaker::processEvent(LCEvent * evt)
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
        
        streamlog_out(MESSAGE2) << "Processing cluster on sensorID " << sensorID  << endl; 
        
        PixelCluster myCluster(cluster->getTrackerData());  
        
        // Calculate center of gravity hit coord in local frame in mm
        double u{0.0}, v{0.0}, sig2_u{0.0}, sig2_v{0.0}, cov_uv{0.0};
        myCluster.getCenterOfGravity(Det, u, v, sig2_u, sig2_v, cov_uv); 
        
        // Override sigmas from user input
        if ( myCluster.getUSize()-1 < _sigmaUCorrections.size() ) {
          sig2_u *= pow(_sigmaUCorrections[myCluster.getUSize()-1],2);  
        }
        
        if ( myCluster.getVSize()-1 < _sigmaVCorrections.size() ) {
          sig2_v *= pow(_sigmaVCorrections[myCluster.getVSize()-1],2); 
        }        
        
        // != 0 means the cluster is marked bad 
        int clsType = 0; 
                       
        streamlog_out(MESSAGE2) << "Stored cluster on sensorID " << sensorID << " at u=" << u << " v=" << v << endl; 
        
        TBHit hit(sensorID, u, v, sig2_u, sig2_v, cov_uv, clsType);
        
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
  
  } catch(DataNotAvailableException &e){
    streamlog_out(MESSAGE2) << "Missing cluster collection in event: " << evt->getEventNumber() << std::endl;
  }   
   
}


//
// Method called after each event to check the data processed
//
void CogHitMaker::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void CogHitMaker::end()
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
void CogHitMaker::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "CogHitMaker Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}



} // Namespace

