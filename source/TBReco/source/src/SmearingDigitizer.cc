// SmearingDigitizer
//                       
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "SmearingDigitizer.h"
#include "TBHit.h"

// Include TBTools  
#include "PhysicalConstants.h"


// Include basic C
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <string>
#include <map>

// Include LCIO classes
#include <lcio.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

#include <TRandom.h>
#include <TRandom3.h>
#include <TMath.h>



// Used namespaces
using namespace CLHEP;
using namespace lcio ;
using namespace marlin ;

namespace depfet {

  //
  // Instantiate this object
  //
  SmearingDigitizer aSmearingDigitizer ;
  
  //
  // Constructor
  //
  SmearingDigitizer::SmearingDigitizer() : Processor("SmearingDigitizer")
  {
    
    // Processor description
    _description = "Smearing digitizer creates TrackerHits by applying noise to SimTrackerHits" ;
    
    //
    // Input collections  
    registerInputCollection (LCIO::SIMTRACKERHIT, "SimTrackerHitCollectionName",
                             "Collection name for SimTrackerHits",
                             m_SimTrackerHitCollectionName, std::string ("SimTrackerHits") );
      
    registerOutputCollection(LCIO::TRACKERHIT, "HitCollectionName",
                             "Collection name for hits",
                              m_hitCollectionName, std::string("hits"));
    
    std::vector<int> initFilterIDs;
    registerProcessorParameter ("FilterIDs",
                                "Apply digitization only to SimTrackerHits for sensors having DAQ IDs in this list",
                                m_filterIDs, initFilterIDs);
    
    registerProcessorParameter( "TanLorentz",
                                "Tangent of Lorentz angle",
                                m_tanLorentzAngle,
                                double(0.0));
    
    registerProcessorParameter( "IntegrationWindow",
                                "Use integration window?",
                                m_integrationWindow,
                                bool(true));
     
    registerProcessorParameter( "StartIntegration",
                                "Only Simulated hits after the StartIntegration time in ns will be digitized",
                                m_startIntegration,
                                double(0.0));
   
    registerProcessorParameter( "StopIntegration",
                                "Only Simulated hits before the StopIntegration time in ns will be digitized",
                                m_stopIntegration,
                                double(20000.0));
    
    registerProcessorParameter( "ClusterSigmaU",
                                "Sigma for cluster covariance matrix [mm]",
                                m_sigmaU,
                                double(0.0032));
    
    registerProcessorParameter( "ClusterSigmaV",
                                "Sigma for cluster covariance matrix [mm]",
                                m_sigmaV,
                                double(0.0032));


   
  }

  //
  // Method called at the beginning of data processing
  //
  void SmearingDigitizer::init() {
    
    // Initialize variables
    m_nRun = 0 ;
    m_nEvt = 0 ;
    
    // Set variables in appropriate physical units
    m_startIntegration     *= ns;
    m_stopIntegration      *= ns;
    
    // Read detector constants from gear file
    m_detector.ReadGearConfiguration();      
    
    // Print set parameters
    printProcessorParams();
    
    // CPU time start
    m_timeCPU = clock()*ms/1000;
  }

  //
  // Method called for each run
  //
  void SmearingDigitizer::processRunHeader(LCRunHeader * run)
  {
  
    // Print run number
    streamlog_out(MESSAGE3) << "Processing run: "
                            << (run->getRunNumber())
                            << std::endl << std::endl;
    
    m_nRun++ ;
  }
  
  //
  // Method called for each event
  //
  void SmearingDigitizer::processEvent(LCEvent * evt)
  {
   
    //
    // Open collections
    try {
      
      // Maybe this event got recorded because of a fake trigger 
      bool isFakeTrigger = false;
      try {
        evt->getCollection( "FakeTrigger" );
        isFakeTrigger = true;
      } catch (lcio::DataNotAvailableException& e) {}  
       
      // Open SimTrackerHit collection
      LCCollection * simHitCol = evt->getCollection( m_SimTrackerHitCollectionName );
      
      // Set collection decoder
      CellIDDecoder<SimTrackerHit> cellIDDec(simHitCol);
      
      // Number of SimTrackerHits
      int nSimHits = simHitCol->getNumberOfElements();
      
      // Create hit collection  
      LCCollectionVec * hitCollection = new LCCollectionVec(LCIO::TRACKERHIT) ;      
      
      //
      // Loop over SimTracker hits;
      streamlog_out(MESSAGE2) << " Producing all SimTrackerHits ..." << std::endl;
      
      for (int i=0; i<nSimHits; ++i) {
        
        SimTrackerHit * simTrkHit = dynamic_cast<SimTrackerHit*>(simHitCol->getElementAt(i));
        
        // Set current - layer ID, ladder ID and sensor ID
        int sensorID = cellIDDec(simTrkHit)["sensorID"];
        
         
        streamlog_out(MESSAGE1) << " Found SimTrackerHit with sensorID: " << sensorID  
                                << std::setprecision(10)
                                << " at time[s]:" << simTrkHit->getTime()
                                << std::setprecision(0) << std::endl;

        streamlog_out(MESSAGE1) << std::setiosflags(std::ios::fixed | std::ios::internal )
                                << std::setprecision(10)
                                << "Integration start/stop[s]: " << m_startIntegration << "/" << m_stopIntegration
                                << std::setprecision(0) << std::endl;
                              
         
        // Cut on simHit creation time --> simulate integration time of a sensor (if option switched on))
	    if ((simTrkHit != 0) && (m_integrationWindow)) {
	      if (simTrkHit->getTime() < m_startIntegration || simTrkHit->getTime() > m_stopIntegration) {
            streamlog_out(MESSAGE1) << " Skipped simHit out of integration time" << std::endl; 
	        continue;
	      }		
	    }
        
        // Cut on simHit creation time --> events with fake trigger never have simhits at t=0 
	    if ((simTrkHit != 0) && (isFakeTrigger)) {
	      if (simTrkHit->getTime() == 0 ) {
            streamlog_out(MESSAGE1) << " Skipped simHit at t=0 because of a fake trigger" << std::endl; 
	        continue;
	      }		
	    }
        
        if ( std::find(m_filterIDs.begin(), m_filterIDs.end(), sensorID) == m_filterIDs.end() ) {
          streamlog_out(MESSAGE2) << " Ignore SimTrackerHit with sensorID: "  << sensorID << std::endl;
          continue;
        }
        
        unsigned short clsType = 0; 
        double u = simTrkHit->getPosition()[0]*mm + gRandom->Gaus(0, m_sigmaU); 
        double v = simTrkHit->getPosition()[1]*mm + gRandom->Gaus(0, m_sigmaV); 
        double cov_u = std::pow(m_sigmaU,2);
        double cov_v = std::pow(m_sigmaV,2); 
        double cov_uv = 0; 
        
        // Create TrackerHit    
        TBHit hit(sensorID, u, v, cov_u, cov_v, cov_uv, clsType);
        
        // Make LCIO TrackerHit
        TrackerHitImpl * trackerhit = hit.MakeLCIOHit();  
              
        // Add hit to the hit collection
        hitCollection->push_back( trackerhit );
     
      } // Loop over simTrkHits
      
      // Store hitCollection in LCIO file
      evt->addCollection( hitCollection, m_hitCollectionName ); 
      
    }
    catch(DataNotAvailableException &e){}
    
    m_nEvt ++ ;
  }

  //
  // Method called after each event to check the data processed
  //
  void SmearingDigitizer::check( LCEvent * evt ) { }
  
  //
  // Method called after all data processing
  //
  void SmearingDigitizer::end()
  {
    // CPU time end
    m_timeCPU = clock()*ms/1000 - m_timeCPU;
    
    // Print message
    streamlog_out(MESSAGE3) << std::endl
                            << " "
                            << "Time per event: "
                            << std::setiosflags(std::ios::fixed | std::ios::internal )
                            << std::setprecision(3)
                            << m_timeCPU/m_nEvt/ms
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
  void SmearingDigitizer::printProcessorParams() const
  {
     
    streamlog_out(MESSAGE3)  << std::endl
                             << " "
                             << "SmearingDigitizer development Version, be carefull!!"
                             << " "
                             << std::endl  << std::endl;   
    
  }

} // Namespace

