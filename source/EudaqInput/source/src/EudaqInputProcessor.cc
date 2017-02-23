/*=====================================================================*/
/*          EudaqInputProcessor                                        */
/*                                                                     */
/*          Author: Benjamin Schwenker                                 */
/*                (benjamin.schwenker@phys.uni-goettingen.de)          */
/*                                                                     */
/*=====================================================================*/

// user includes
#include "EudaqInputProcessor.h"

#include "FileReader.hh"
#include "DetectorEvent.hh"
#include "RawDataEvent.hh"
#include "PluginManager.hh"


// marlin includes
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"

// lcio includes
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTime.h>

// system includes
#include <string>
#include <iomanip>
#include <iostream>


using namespace std;
using namespace marlin;
using namespace eudaqinput;

namespace {
  //static const unsigned TLUID = eudaqinput::Event::str2id("_TLU");
  static const unsigned IDMASK = 0x7fff;
}

EudaqInputProcessor::EudaqInputProcessor ():DataSourceProcessor  ("EudaqInputProcessor") {
  
  _description =
    "Reads a EUDAQ .raw file and converts it event by event into LCIO format\n"
    "Make sure to not specify any LCIOInputFiles in the steering to read raw files.";
  
  registerProcessorParameter ("FileName", "Input file",
                              m_filename, std::string ("input.raw"));
  
  registerProcessorParameter ("SyncTriggerId", "Syncronize subevents by trigger id",
                              m_synctriggerid, static_cast < bool > (false));

  registerProcessorParameter( "DetectorName", "Set name of the detector",
                              m_detectorName, std::string("EUTelescope"));
    
}


EudaqInputProcessor * EudaqInputProcessor::newProcessor () {  
  return new EudaqInputProcessor;
}


void EudaqInputProcessor::init () {
  printParameters ();
  
  // Initialize some counter
  m_ndata = 0;
  m_nbore = 0; 
  m_neore = 0;     
}


void EudaqInputProcessor::readDataSource (int Ntrig) {
  
  // ------------------------------------------------------- // 
  // Create a file reader for raw data  
  eudaqinput::FileReader reader(m_filename, "", m_synctriggerid);
  
  // ------------------------------------------------------- // 
  // Read rawdata file   
  do {
    const eudaqinput::Event & ev = reader.GetEvent();
        
    if (ev.IsBORE()) {
      m_nbore++;
      if (m_nbore == 1) {
        // Process BORE event 
        const eudaqinput::DetectorEvent * dev = dynamic_cast<const eudaqinput::DetectorEvent *>(&ev);
         
        // Initialize pugin manager
        eudaqinput::PluginManager::Initialize(*dev);
        
        // Write LCIO run header here 
        IMPL::LCRunHeaderImpl* lcHeader = new IMPL::LCRunHeaderImpl;
        lcHeader->setDescription(" Reading from file " + reader.Filename());
        lcHeader->setRunNumber((*dev).GetRunNumber());
        lcHeader->setDetectorName(m_detectorName);
        
        // Add run header to LCIO file
        ProcessorMgr::instance ()->processRunHeader ( static_cast<lcio::LCRunHeader*> (lcHeader) );
        delete lcHeader; 
        
      } else { 
         streamlog_out ( MESSAGE3 ) << "Multiple BOREs " << to_string(m_nbore) << std::endl ;
      }    
    } else if (ev.IsEORE()) {
      m_neore++;
      if (m_neore > 1) {
        streamlog_out ( MESSAGE3 ) << "Warning: Multiple EOREs: " << to_string(m_neore) << std::endl;
      }    
    } else {
      m_ndata++;
                      
      const eudaqinput::DetectorEvent * dev = dynamic_cast<const eudaqinput::DetectorEvent *>(&ev);
                           
      for (size_t i = 0; i < (*dev).NumEvents(); ++i) {
        const eudaqinput::Event* subev = (*dev).GetEvent(i);           
               
        streamlog_out ( MESSAGE2 ) << "  TrigID  " << (eudaqinput::PluginManager::GetTriggerID(*subev) & IDMASK) 
                                   << "  (" << subev->GetSubType() << ")"         
                                   << endl;  
      } 
      
      // This actually converts the detector raw data 
      // FIXME store raw data and move unpacking of raw data into seperate processors ('unpackers')             
      LCEvent * lcEvent = eudaqinput::PluginManager::ConvertToLCIO(*dev);
            
      ProcessorMgr::instance ()->processEvent (static_cast<LCEventImpl*> (lcEvent));
      delete lcEvent;        
    }       
  } while (reader.NextEvent());
      
}


void EudaqInputProcessor::end () 
{ 
  streamlog_out ( MESSAGE3 ) << "Number of data events: " << to_string(m_ndata) << std::endl;
  if (!m_nbore) {
    streamlog_out ( MESSAGE3 ) << "Warning: No BORE found" << std::endl;
  } else if (m_nbore > 1) {
    streamlog_out ( MESSAGE3 ) << "Warning: Multiple BOREs found: " << m_nbore << std::endl;
  }
            
  if (!m_neore) {
    streamlog_out ( MESSAGE3 ) << "Warning: No EORE found, possibly truncated file." << std::endl;
  } else if (m_neore > 1) {
    streamlog_out ( MESSAGE3 ) << "Warning: Multiple EOREs found: " << m_nbore << std::endl;
  }
}



