/*=====================================================================*/
/*          EudaqInputProcessor                                        */
/*                                                                     */
/*          Author: Benjamin Schwenker                                 */
/*                (benjamin.schwenker@phys.uni-goettingen.de)          */
/*                                                                     */
/*  New version stores raw data into lcio::LCEvent and delegates       */
/*  unpacking (converting) of raw data to other downstream Marlin      */
/*  processors called unpackers.                                        */
/*=====================================================================*/

// user includes
#include "EudaqInputProcessor.h"

#include "FileReader.hh"
#include "DetectorEvent.hh"
#include "RawDataEvent.hh"

// marlin includes
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"

// lcio includes
#include "lcio.h"
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <UTIL/CellIDEncoder.h>

// system includes
#include <string>
#include <iomanip>
#include <iostream>


using namespace std;
using namespace marlin;
using namespace eudaqinput;

namespace {
  typedef std::vector<unsigned char> datavect;
  typedef std::vector<unsigned char>::const_iterator datait;
}

EudaqInputProcessor::EudaqInputProcessor ():DataSourceProcessor  ("EudaqInputProcessor") {
  
  _description =
    "Reads a EUDAQ .raw file event by event and stores raw data into lcio::LCEvents.\n"
    "Make sure to not specify any LCIOInputFiles.";
  
  registerProcessorParameter ("FileName", "Input file",
                              m_filename, std::string ("input.raw"));
  
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
  eudaqinput::FileReader reader(m_filename, "", false);
  
  // ------------------------------------------------------- // 
  // Read rawdata file   
  do {
    const eudaqinput::Event & ev = reader.GetEvent();
        
    if (ev.IsBORE()) {
      m_nbore++;
      
      if (m_nbore == 1) {
        // Process BORE event 

        const eudaqinput::DetectorEvent & dev = *dynamic_cast<const eudaqinput::DetectorEvent *>(&ev);
            
        // Write LCIO run header here 
        IMPL::LCRunHeaderImpl* lcHeader = new IMPL::LCRunHeaderImpl;
        lcHeader->setDescription(" Reading from file " + reader.Filename());
        lcHeader->setRunNumber(dev.GetRunNumber());
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
                      
      // Process data event 
      const eudaqinput::DetectorEvent & dev = *dynamic_cast<const eudaqinput::DetectorEvent *>(&ev);
       
      lcio::LCEventImpl * lcEvent = new lcio::LCEventImpl;
      lcEvent->setEventNumber(dev.GetEventNumber());
      lcEvent->setRunNumber(dev.GetRunNumber());
      lcEvent->setTimeStamp(dev.GetTimestamp());
       
      for (size_t i = 0; i < dev.NumEvents(); ++i) {
          
        const eudaqinput::Event& subev = *dev.GetEvent(i);    
        
        // process RawDataEvent
        if (const RawDataEvent * rawev = dynamic_cast<const RawDataEvent *>(&subev)) { 
          
          // this string will be the name of the raw data collection
          std::string type = rawev->GetSubType();
          
          // prepare the collections for the raw data 
          LCCollectionVec * rawDataCollection = new LCCollectionVec( lcio::LCIO::TRACKERRAWDATA );
                    
          // set the proper cell encoder
          CellIDEncoder<TrackerRawDataImpl> lcEncoder( "blockID:6", rawDataCollection  ); 
          
          // loop over all data blocks in sub event
          for (size_t j = 0; j < rawev->NumBlocks(); ++j) {
          
            // get a handle to data block
            const datavect & data = rawev->GetBlock(j);
            // prepare lcio::TrackerRawData
            lcio::TrackerRawDataImpl* lcBlock =  new lcio::TrackerRawDataImpl;
            lcEncoder["blockID"] = rawev->GetID(j); 
            lcEncoder.setCellID( lcBlock );
            
            // loop over the data block
            for (size_t it = 0; it < data.size(); ++it) {
              lcBlock->adcValues().push_back( data[it] );
            }
            
            // add the lcio::TrackerRawData to the collection
            rawDataCollection->push_back( lcBlock );  
          } 
          
          // add collection to the event
          lcEvent->addCollection( rawDataCollection, type.c_str() );
        }        
      }
            
      ProcessorMgr::instance ()->processEvent (static_cast<LCEventImpl*> (lcEvent));
      delete lcEvent;        
    }       
  } while (reader.NextEvent()  &&  m_ndata < Ntrig );
      
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



