/*=====================================================================*/
/*          EventInfoSetter                                            */
/*                                                                     */
/*          Author: Benjamin Schwenker                                 */
/*                (benjamin.schwenker@phys.uni-goettingen.de)          */
/*                                                                     */
/*=====================================================================*/

// user includes
#include "EventInfoSetter.h"

// marlin includes
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"

// lcio includes
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTime.h>

// system includes
#include <string.h>


using namespace std;
using namespace marlin;
using namespace depfet;

EventInfoSetter::EventInfoSetter ():DataSourceProcessor  ("EventInfoSetter") {
  
  _description =
    "Creates a lcio file with run header and empty events\n"
    "Make sure to not specify any LCIOInputFiles in the steering in order";
    
  
  registerProcessorParameter ("RunNumber", "Set run number for lcio file",
                              m_runNumber,static_cast < int >(0)); 
  
  registerProcessorParameter( "DetectorName", "Set name of the detector",
                              m_detectorName, std::string("EUDET"));
   
}


EventInfoSetter * EventInfoSetter::newProcessor () {  
  return new EventInfoSetter;
}


void EventInfoSetter::init () {
  printParameters ();
}


void EventInfoSetter::readDataSource (int Ntrig) {
   
  // Create LCIO run header here 
  IMPL::LCRunHeaderImpl* lcHeader = new IMPL::LCRunHeaderImpl;
  lcHeader->setDescription("Simulated events");
  lcHeader->setRunNumber (m_runNumber);
  lcHeader->setDetectorName(m_detectorName);
         
  // Add run header to LCIO file
  ProcessorMgr::instance ()->processRunHeader ( static_cast<lcio::LCRunHeader*> (lcHeader) );
  delete lcHeader; 
   
  int fakeEventNumber = -1; 
   
   
  // Start reading raw data stream  
  while( fakeEventNumber < Ntrig){
     
    ++fakeEventNumber; 
     
    // Print event number 
    if ( fakeEventNumber%100 == 0 ) 
      streamlog_out(MESSAGE3) << "Events processed: " << fakeEventNumber
                               << std::endl << std::endl;
     
    // Copy to LCIO event 
    LCEventImpl* evt = new LCEventImpl;
     
    // Copy genaral info 
    evt->setDetectorName(m_detectorName);
    evt->setRunNumber(m_runNumber); 
    evt->setEventNumber(fakeEventNumber); 
    LCTime * now = new LCTime;
    evt->setTimeStamp(now->timeStamp());
    delete now;
      
    // Add event to LCIO file
    ProcessorMgr::instance ()->processEvent (static_cast<LCEventImpl*> (evt));
    delete evt;
     
  } // End loop over events 
         
}


void EventInfoSetter::end () 
{ }



