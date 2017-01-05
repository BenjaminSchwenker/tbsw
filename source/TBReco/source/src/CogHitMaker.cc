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

   registerProcessorParameter ("SigmaU1", 
                               "SigmaU for cluster with one uCell contributing [mm]",
                               _SigmaU1,  
                               static_cast < double > (1)); 

   registerProcessorParameter ("SigmaU2", 
                               "SigmaU for cluster with two uCells contributing [mm]",
                               _SigmaU2,  
                               static_cast < double > (1));  

   registerProcessorParameter ("SigmaU3", 
                               "SigmaU for cluster with three or more uCells contributing [mm]",
                               _SigmaU3,  
                               static_cast < double > (1));  

   registerProcessorParameter ("SigmaV1", 
                               "SigmaV for cluster with one vCell contributing [mm]",
                               _SigmaV1,  
                               static_cast < double > (1)); 

   registerProcessorParameter ("SigmaV2", 
                               "SigmaV for cluster with two vCells contributing [mm]",
                               _SigmaV2,  
                               static_cast < double > (1));  
    
   registerProcessorParameter ("SigmaV3", 
                               "SigmaV for cluster with three or more vCells contributing [mm]",
                               _SigmaV3,  
                               static_cast < double > (1));

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
    
    streamlog_out(MESSAGE2) << "Processing cluster on sensorID " << sensorID  << endl; 
    
    PixelCluster myCluster(cluster->getTrackerData());  
       
    // Calculate hit coord in local frame in mm
    float u = 0; 
    float v = 0;  
    float total = 0; 
      
    // Loop over digits and compute hit coordinates
    TrackerData * clusterDigits =  cluster->getTrackerData();
    FloatVec sparsePixels = clusterDigits->getChargeValues();
    int nPixels = sparsePixels.size()/3; 
       
    for ( int index=0; index<nPixels;  index++) { 
         
      int col = static_cast<int> (sparsePixels[index * 3]);
      int row = static_cast<int> (sparsePixels[index * 3 + 1]);
      float chargeValue =  sparsePixels[index * 3 + 2];
        
      total += chargeValue;
      u += Det.GetPixelCenterCoordU(row, col)*chargeValue;
      v += Det.GetPixelCenterCoordV(row, col)*chargeValue;       
    }
      
    if ( total > 0)  {
      u /= total;
      v /= total; 
    } 
      
    double cov_u{0.0} , cov_v {0.0}, cov_uv{0.0};    

    if ( myCluster.getUSize() == 1) cov_u = pow(_SigmaU1,2);
    else if ( myCluster.getUSize() == 2) cov_u = pow(_SigmaU2,2); 
    else if ( myCluster.getUSize() == 3) cov_u = pow(_SigmaU3,2); 
    else cov_u = pow(myCluster.getUSize()*Det.GetPitchU()/TMath::Sqrt(12),2);  
     
    if ( myCluster.getVSize() == 1) cov_v = pow(_SigmaV1,2);
    else if ( myCluster.getVSize() == 2) cov_v = pow(_SigmaV2,2);
    else if ( myCluster.getVSize() == 3) cov_v = pow(_SigmaV3,2);
    else cov_v = pow(myCluster.getVSize()*Det.GetPitchV()/TMath::Sqrt(12),2);
     
    // = 1 means the cluster is marked bad 
    unsigned short clsType = 0; 

    if ( myCluster.getUSize() > _maxSizeU ) clsType = 1; 
    if ( myCluster.getVSize() > _maxSizeV ) clsType = 1;        

    TBHit hit(sensorID, u, v, cov_u, cov_v, cov_uv, clsType);
      
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

