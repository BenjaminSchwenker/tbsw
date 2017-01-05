// StripHitMaker - Marlin Processor
// 
// Compute center of gravity hits from clusters 
//
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 
#include "StripHitMaker.h"
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



// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;

namespace depfet {

//
// Instantiate this object
//
StripHitMaker aStripHitMaker ;

//
// Constructor
//
StripHitMaker::StripHitMaker() : Processor("StripHitMaker")
{

// Processor description
   _description = "StripHitMaker:  Compute center of gravity hits from strip clusters ";
   
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
   
   registerProcessorParameter("ClusterQualitySelection",
                              "To use only kGoodQuality write 0 here",
                              _clusterQualitySelect, 
                              static_cast<int> ( 0 ));
   
}

//
// Method called at the beginning of data processing
//
void StripHitMaker::init() {
    
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
void StripHitMaker::processRunHeader(LCRunHeader * run)
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
void StripHitMaker::processEvent(LCEvent * evt)
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
        
  // Helper class for decoding cluster ID's 
  CellIDDecoder<TrackerPulseImpl> clusterDecoder( clusterCollection ); 
         
  // Create intermediate hit collections  
  std::vector<TBHit> HitStoreU;   
  std::vector<TBHit> HitStoreV;   
  
  // Loop on all clusters 
  for (unsigned int iClu = 0; iClu < clusterCollection->size(); iClu++) 
  {   
    // Read cluster header
    TrackerPulseImpl* cluster = dynamic_cast<TrackerPulseImpl* > ( clusterCollection->getElementAt(iClu) )  ;       
    int sensorID = clusterDecoder(cluster)["sensorID"]; 
    int ipl = _detector.GetPlaneNumber(sensorID);
    Det& Sensor = _detector.GetDet(ipl);
    
    streamlog_out(MESSAGE2) << "Processing cluster on sensorID " << sensorID  << endl; 
      
    // Calculate hit coord in local frame in mm
    float u = 0;
    float v = 0;    
    float total = 0; 
    int isU = 0; 
    
    // Loop over digits and compute hit coordinates
    TrackerData * clusterDigits =  cluster->getTrackerData();
    FloatVec rawDigits = clusterDigits->getChargeValues();
    int nDigits = rawDigits.size()/3; 
       
    for ( int i=0; i<nDigits;  i++) { 
         
      isU = static_cast<int> (rawDigits[i * 3]);
      int cell = static_cast<int> (rawDigits[i * 3 + 1]);
      float signal = rawDigits[i * 3 + 2]; 
      
      total += signal;
      
      if (isU) 
        u += Sensor.GetPixelCenterCoordU(0, cell)*signal;  
      else 
        v += Sensor.GetPixelCenterCoordV(cell, 0)*signal;
              
    }
      
    if ( total > 0)  {
      u /= total;
      v /= total; 
    } 
      
    double cov_v = pow(Sensor.GetResolutionV(),2); 
    double cov_u = pow(Sensor.GetResolutionU(),2);
    double cov_uv = 0;         

    TBHit hit(sensorID, u, v, cov_u, cov_v, cov_uv, 0);
    hit.SetUniqueID(iClu);      
    
    if (isU)
      HitStoreU.push_back(hit);  
    else  
      HitStoreV.push_back(hit);  
         
  } // End cluster loop 
  
  // Now, try to merge clusters along u and v 

  LCCollectionVec * hitCollection = new LCCollectionVec(LCIO::TRACKERHIT) ;      
  
  for (int iU = 0; iU < (int) HitStoreU.size(); iU++ ) {
    for (int iV = 0; iV < (int) HitStoreV.size(); iV++ ) {
       
      TBHit hitU = HitStoreU[iU]; 
      TBHit hitV = HitStoreV[iV];
      
      if (hitU.GetDAQID() == hitV.GetDAQID() ) {
        
        // Make LCIO TrackerHit 
        TBHit hit(hitU.GetDAQID(), hitU.GetCoord()[0][0], hitV.GetCoord()[1][0], hitU.GetCov()[0][0], hitV.GetCov()[1][1], 0, 0);
        TrackerHitImpl * trackerhit = hit.MakeLCIOHit();  
            
        // Add link to full cluster data 
        TrackerPulseImpl* clusterU = dynamic_cast<TrackerPulseImpl* > ( clusterCollection->getElementAt(hitU.GetUniqueID()) );
        TrackerPulseImpl* clusterV = dynamic_cast<TrackerPulseImpl* > ( clusterCollection->getElementAt(hitV.GetUniqueID()) )  ;       
   
        LCObjectVec clusterVec;
        clusterVec.push_back( clusterU->getTrackerData() );
        clusterVec.push_back( clusterV->getTrackerData() );
        trackerhit->rawHits() = clusterVec;

        hitCollection->push_back( trackerhit );
        
      }

    }
  }          


  // Store hitCollection in LCIO file
  evt->addCollection( hitCollection, _hitCollectionName );     
     
}


//
// Method called after each event to check the data processed
//
void StripHitMaker::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void StripHitMaker::end()
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
void StripHitMaker::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "StripHitMaker Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}



} // Namespace

