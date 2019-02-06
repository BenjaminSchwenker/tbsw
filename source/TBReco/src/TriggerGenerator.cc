// TriggerGenerator
//                       
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// user includes
#include "TriggerGenerator.h"

// C++ includes
#include <iostream>
#include <iomanip>


// Include LCIO classes
#include <lcio.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/SimTrackerHitImpl.h> 
#include <UTIL/CellIDDecoder.h>


// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;

namespace depfet {

  //
  // Instantiate this object
  //
  TriggerGenerator aTriggerGenerator ;
  
  bool TrgScinti::isPointInSensor( double u, double v) 
  {
    // Boundary set
    if ( u > m_minU && u < m_maxU && v > m_minV && v < m_maxV ) return true; 
    else return false; 
  }


  //
  // Constructor
  //
  TriggerGenerator::TriggerGenerator() : Processor("TriggerGenerator")
  {
    
    // Processor description
    _description = "Simulates the trigger signal from scinti coincidence for test beam experiments";
    
    //
    // Output collections  
    registerOutputCollection(LCIO::SIMTRACKERHIT,"SimTrackerHitCollectionName",
                             "Collection name for SimTrackerHits",
                             m_SimTrackerHitCollectionName, string ("SimTrackerHits"));
    
    registerProcessorParameter ("FakeTriggerPeriod",
                             "Every nth simulated event will be accepted regardless of trigger coincidence",
                              m_fakeTriggerPeriod, static_cast< int > ( 0 ) );  
    
    std::vector<float> initScinti;
    registerProcessorParameter ("ScinitNo1", "Scinti parameters: DAQID, Umin[mm], Vmin[mm], Umax[mm], Vmax[mm] (leave empty to deactivate)",
                              _scintiNo1, initScinti );
    
    registerProcessorParameter ("ScinitNo2", "Scinti parameters: DAQID, Umin[mm], Vmin[mm], Umax[mm], Vmax[mm] (leave empty to deactivate)",
                              _scintiNo2, initScinti );
    
    registerProcessorParameter ("ScinitNo3", "Scinti parameters: DAQID, Umin[mm], Vmin[mm], Umax[mm], Vmax[mm] (leave empty to deactivate)",
                              _scintiNo3, initScinti );
    
    registerProcessorParameter ("ScinitNo4", "Scinti parameters: DAQID, Umin[mm], Vmin[mm], Umax[mm], Vmax[mm] (leave empty to deactivate)",
                              _scintiNo4, initScinti );
                                 
  }
  
  void TriggerGenerator::init () {
  
    // Initialize variables
    m_nRun = 0 ;
    m_nEvt = 0 ;
  
    // Print set parameters
    printProcessorParams();
  
    // CPU time start
    m_timeCPU = clock()/1000;

    // Read detector constants from gear file
    m_detector.ReadGearConfiguration();  
    
    if ( _scintiNo1.size() >= 5 ) {
      int sensorID = (int) _scintiNo1[0];
      int ipl = m_detector.GetPlaneNumber(sensorID);
      float minU = _scintiNo1[1];
      float maxU = _scintiNo1[3];
      float minV = _scintiNo1[2];
      float maxV = _scintiNo1[4];
      
      if ( minU < maxU && minV < maxV && ipl >= 0 && ipl < m_detector.GetNSensors() ) {
        TrgScinti  sct(sensorID, minU, maxU, minV, maxV);
        m_scintiVec.push_back(sct);
      }
    }  
    
    
    if ( _scintiNo2.size() >= 5 ) {
      int sensorID = (int) _scintiNo2[0];
      int ipl = m_detector.GetPlaneNumber(sensorID);
      float minU = _scintiNo2[1];
      float maxU = _scintiNo2[3];
      float minV = _scintiNo2[2];
      float maxV = _scintiNo2[4];
      
      if ( minU < maxU && minV < maxV && ipl >= 0 && ipl < m_detector.GetNSensors() ) {
        TrgScinti  sct(sensorID, minU, maxU, minV, maxV);
        m_scintiVec.push_back(sct);
      }
    }  

    if ( _scintiNo3.size() >= 5 ) {
      int sensorID = (int) _scintiNo3[0];
      int ipl = m_detector.GetPlaneNumber(sensorID);
      float minU = _scintiNo3[1];
      float maxU = _scintiNo3[3];
      float minV = _scintiNo3[2];
      float maxV = _scintiNo3[4];
      
      if ( minU < maxU && minV < maxV && ipl >= 0 && ipl < m_detector.GetNSensors() ) {
        TrgScinti  sct(sensorID, minU, maxU, minV, maxV);
        m_scintiVec.push_back(sct);
      }
    }  
    
    
    if ( _scintiNo4.size() >= 5 ) {
      int sensorID = (int) _scintiNo4[0];
      int ipl = m_detector.GetPlaneNumber(sensorID);
      float minU = _scintiNo4[1];
      float maxU = _scintiNo4[3];
      float minV = _scintiNo4[2];
      float maxV = _scintiNo4[4];
      
      if ( minU < maxU && minV < maxV && ipl >= 0 && ipl < m_detector.GetNSensors() ) {
        TrgScinti  sct(sensorID, minU, maxU, minV, maxV);
        m_scintiVec.push_back(sct);
      }
    }

    streamlog_out(MESSAGE3) << "Number of active scintis is "
                            << m_scintiVec.size()
                            << std::endl << std::endl;
    
  }
  
  //
  // Method called for each run
  //
  void TriggerGenerator::processRunHeader(LCRunHeader * run)
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
  void TriggerGenerator::processEvent(LCEvent * evt)
  {
     
    //
    // Open collections
    try {
      
      // Open SimTrackerHit collection
      LCCollection * simHitCol = evt->getCollection( m_SimTrackerHitCollectionName );
      
      // Set collection decoder
      CellIDDecoder<SimTrackerHit> cellIDDec(simHitCol);
      
      // Number of SimTrackerHits
      int nSimHits = simHitCol->getNumberOfElements(); 
      
      // State for active scintis 
      std::vector< bool > sctStateVec( m_scintiVec.size(), false);  
      
      for (int i=0; i<nSimHits; ++i) {
         
        SimTrackerHit * simTrkHit = dynamic_cast<SimTrackerHit*>(simHitCol->getElementAt(i));
        
        // Set current - layer ID, ladder ID and sensor ID
        int sensorID = cellIDDec(simTrkHit)["sensorID"];
        float t = simTrkHit->getTime();
        float u = simTrkHit->getPosition()[0]; 
        float v = simTrkHit->getPosition()[1];
        
        // Only t=0 simhits can raise triggers
        if (t > 0) continue; 
        
        // Check if simhit can trigger a scinti 
        for ( size_t i = 0; i < m_scintiVec.size(); i++ )  {
          bool signal = ( m_scintiVec[i].GetDAQID() == sensorID ) && m_scintiVec[i].isPointInSensor(u,v);         
          sctStateVec[i] = sctStateVec[i] || signal; 
        }
      }
      
      // If true, we keep the simulated event for digitization
      bool trg = true; 
       
      // Compute the coincidence of scinti signals 
      for ( auto isHigh : sctStateVec )  trg = trg && isHigh ; 
       
      // Maybe it is a event with a fake trigger
      if ( m_fakeTriggerPeriod > 0 ) { 
        if ((evt->getEventNumber())%m_fakeTriggerPeriod == 0) {
          trg = true;
          // We signal a fake trigger by inserting an empty collection called "FakeTrigger" into the event
          LCCollectionVec * triggerVec = new LCCollectionVec(LCIO::SIMTRACKERHIT) ;
          evt->addCollection(triggerVec, "FakeTrigger"); 
        }

      }

      

        
      if (!trg) throw( marlin::SkipEventException(this) );
      
    } catch(DataNotAvailableException &e){}
    
    m_nEvt ++ ; 
    
  }

  //
  // Method called after each event to check the data processed
  //
  void TriggerGenerator::check( LCEvent * evt ) {}

  //
  // Method called after all data processing
  //
  void TriggerGenerator::end()
  {
   
    streamlog_out ( MESSAGE3 ) << endl;
    streamlog_out ( MESSAGE3 ) << "Successfully finished" << endl;
    
    // CPU time end
    m_timeCPU = clock()/1000 - m_timeCPU;
   
    // Print message
    streamlog_out(MESSAGE3) << std::endl
                             << " "
                             << "Time per event: "
                             << std::setiosflags(std::ios::fixed | std::ios::internal )
                             << std::setprecision(3)
                             << m_timeCPU/m_nEvt
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
  void TriggerGenerator::printProcessorParams() const
  {

    streamlog_out(MESSAGE3)  << std::endl
                              << " "
                              << "TriggerGenerator development Version, be carefull!!"
                              << " "
                              << std::endl  << std::endl;   


  }

} // Namespace


