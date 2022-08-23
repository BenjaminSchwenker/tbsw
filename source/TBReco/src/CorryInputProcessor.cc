/*=====================================================================*/
/*          CorryInputProcessor                                        */
/*                                                                     */
/*          Author: Benjamin Schwenker                                 */
/*                (benjamin.schwenker@phys.uni-goettingen.de)          */
/*                                                                     */
/*  Processers converts raw digits in hit tables from one or more      */ 
/*  ascii input files into series of lcio::LCEvents. In case multiple  */
/*  ascii input files are are supplied, raw digits are aligned using   */
/*  the event number.                                                  */
/*=====================================================================*/

// user includes
#include "CorryInputProcessor.h"

// marlin includes
#include "marlin/ProcessorMgr.h"

// lcio includes
#include "lcio.h"
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>


// system includes
#include <iomanip>
#include <iostream>
#include <assert.h>


using namespace std;
using namespace marlin;
using namespace depfet;
using namespace std::string_literals;



CorryInputProcessor::CorryInputProcessor ():DataSourceProcessor("CorryInputProcessor") {
  
  _description =
    "Reads hit table from one or more ascii files. Events are formed and hits are written to lcio::TrackerData collections.\n"
    "Make sure to not specify any LCIOInputFiles.";
  
  std::vector< string > inputFileNamesVecExample;
  inputFileNamesVecExample.push_back("input.txt");
  registerProcessorParameter ("FileNames", "Whitespace seperated list of input file names to read",
                              m_fileNameVec, inputFileNamesVecExample);

  std::vector< string > inputSensorNamesVecExample;
  inputSensorNamesVecExample.push_back("dummy");
  registerProcessorParameter ("SensorNames", "Whitespace seperated list of sensor names",
                              m_sensorNamesVec, inputSensorNamesVecExample);

  std::vector< int > inputSensorIDsVecExample;
  inputSensorIDsVecExample.push_back(-1);
  registerProcessorParameter ("SensorIDs", "Whitespace seperated list of sensorIDs",
                              m_sensorIDsVec, inputSensorIDsVecExample);
  
  registerProcessorParameter( "DetectorName", "Set name of the detector. Needs to be the same name as in gear file.",
                              m_detectorName, std::string("EUTelescope"));

  registerProcessorParameter( "RawHitCollectionName", "Set name of the raw hit collection.",
                              m_rawHitCollectionName, std::string("zsdata"));

  registerProcessorParameter ("RunNumber", "Set run number.",
                              m_runNumber,  static_cast < int > (0));
}


CorryInputProcessor * CorryInputProcessor::newProcessor () {  
  return new CorryInputProcessor;
}


void CorryInputProcessor::init () {
  m_chunkSize = 1000;
  m_ntriggers = 0;

  assert(m_sensorNamesVec.size() == m_sensorIDsVec.size());
  for (size_t i = 0; i < m_sensorNamesVec.size(); ++i) {
    streamlog_out(MESSAGE3) <<  m_sensorNamesVec[i] << ":"<<  m_sensorIDsVec[i] << std::endl;
    m_sensor_lookup[m_sensorNamesVec[i]] = m_sensorIDsVec[i];
  }

  printParameters ();
}


void CorryInputProcessor::readDataSource (int Ntrig) {
  
  if (m_fileNameVec.size() == 0) return;
  string filename = m_fileNameVec[0];

  streamlog_out ( MESSAGE3 ) << "Start reading Ntrig=" << Ntrig << " from raw file " << filename << std::endl;
   
  // Write LCIO run header here 
  IMPL::LCRunHeaderImpl* lcHeader = new IMPL::LCRunHeaderImpl;
  lcHeader->setDescription(" Reading from file " + filename);
  lcHeader->setRunNumber( m_runNumber );  
  lcHeader->setDetectorName(m_detectorName);
          
  // Add run header to LCIO file
  ProcessorMgr::instance ()->processRunHeader ( static_cast<lcio::LCRunHeader*> (lcHeader) );
  delete lcHeader; 
  
  streamlog_out ( MESSAGE2 ) << "Start reading next chunk from raw file " << filename << std::endl;
     
  // Open file stream
  fstream fin(filename);
  assert (!fin.fail( ));   
  
  // Stores the curent read position to other files 
  vector<streampos> other_pos;
  for (size_t iFile=1; iFile<m_fileNameVec.size(); iFile++) {
    string other_filename = m_fileNameVec[iFile];
    // Open file stream
    fstream other_fin(other_filename);  
    assert (!other_fin.fail( ));
    // Remember current position 
    other_pos.push_back(other_fin.tellg());     
  }
  
  bool readNextChunk=true;
  while(readNextChunk and m_ntriggers < Ntrig) {  
    
    // Read chunk of data from the first file 
    vector< vector<int> > hits;
    vector< vector<int> > events;
    
    streamlog_out ( MESSAGE2 ) << "Start reading next chunk from raw file " << filename << std::endl;
    readNextChunk = getAllEventFromFile(fin, hits, events, m_chunkSize, -1);  
    
    // Read data from all other files 
    vector< vector< vector<int> > > all_other_hits;
    vector< vector< vector<int> > > all_other_events;
    
    // now add data from other files  
    for (size_t iFile=1; iFile< m_fileNameVec.size(); iFile++) {
      string other_filename = m_fileNameVec[iFile];
      
      streamlog_out ( MESSAGE2 ) << "Start reading next chunk from raw file " << other_filename << std::endl;
      
      vector< vector<int> > other_hits;
      vector< vector<int> > other_events;
      
      // Open file stream
      fstream other_fin(other_filename);  
      assert (!other_fin.fail( ));         
      
      // Jump to position where reading stopped after last chunk
      other_fin.seekg ( other_pos[iFile-1] );
      
      // Read next chunk from file   
      getAllEventFromFile(other_fin, other_hits, other_events, -1, events[events.size()-1][0]); 
      
      // Remember current reading position for next chunk 
      other_pos[iFile-1] = other_fin.tellg();
      
      all_other_hits.push_back(other_hits);
      all_other_events.push_back(other_events);   
    }       
    
    // Process all events from curent chunk 
    for (size_t iEvt = 0; iEvt < events.size(); iEvt++) {
      // Count number of LCIO events   
      m_ntriggers++; 
      
      // Process data event 
      lcio::LCEventImpl * lcEvent = new lcio::LCEventImpl;
      lcEvent->setEventNumber(events[iEvt][0]);
      lcEvent->setRunNumber(m_runNumber); 
      lcEvent->setTimeStamp(0);  
      
      // Prepare the collections for the raw data 
      LCCollectionVec * rawDataCollection = new LCCollectionVec( lcio::LCIO::TRACKERDATA );           
       
      // Prepare a new lcio::TrackerData for the ZS data
      lcio::TrackerDataImpl* zsFrame =  new lcio::TrackerDataImpl;
      
      auto& chargevec=zsFrame->chargeValues();
      chargevec.reserve(chargevec.size()+size_t(60)  );
      
      for (int iHit=events[iEvt][1]; iHit<=events[iEvt][2]; iHit++ ) {
        chargevec.push_back( hits[iHit][1] );
        chargevec.push_back( hits[iHit][2] );
        chargevec.push_back( hits[iHit][3] );
        chargevec.push_back( hits[iHit][4] );
        chargevec.push_back( hits[iHit][5] );
      }
      
      // now add data from other file 
      for (size_t iFile=1; iFile<m_fileNameVec.size(); iFile++) {   
        for (size_t oEvt=0; oEvt<all_other_events[iFile-1].size(); oEvt++) {  
          if (all_other_events[iFile-1][oEvt][0] == events[iEvt][0]) { 
            for (int iHit = all_other_events[iFile-1][oEvt][1]; iHit<=all_other_events[iFile-1][oEvt][2]; iHit++ ) {  
              chargevec.push_back( all_other_hits[iFile-1][iHit][1] );
              chargevec.push_back( all_other_hits[iFile-1][iHit][2] );
              chargevec.push_back( all_other_hits[iFile-1][iHit][3] );
              chargevec.push_back( all_other_hits[iFile-1][iHit][4] );
              chargevec.push_back( all_other_hits[iFile-1][iHit][5] );
            } 
          }
        }
      }
      
      // Now add the TrackerData to the collection
      rawDataCollection->push_back( zsFrame );
      
      // Add output collection to event
      lcEvent->addCollection( rawDataCollection, m_rawHitCollectionName ); 
      
      ProcessorMgr::instance ()->processEvent (static_cast<LCEventImpl*> (lcEvent));
      delete lcEvent;  
    }
  }
}

bool CorryInputProcessor::getAllEventFromFile(fstream& fin, vector< vector<int> >& hits, std::vector< std::vector<int> >& events, int max_events, int max_event_number) 
{
  
  string line, word; 
  int event_number=-1;  
  int sensorID=-1; 
  int column=-1; 
  int row=-1; 
  int charge=-1;   
  int time=-1;
  int firstHit = 0;
  int current_event = -1; 
  bool validHit = false;

  // stores the position
  streampos oldpos = fin.tellg();
  
  while (getline(fin, line)) {
         
    if ( line.substr(0, 3) == "===") {
      // Start of event: Set event number
      stringstream sd(line); 
      getline(sd, word, ' ');   
      getline(sd, word, ' ');   
      event_number=stoi(word);
      validHit = false;
    } else if ( line.substr(0, 3) == "---") {
      // Start of sensor data: Set sensorID 
      stringstream sd(line); 
      getline(sd, word, ' ');   
      getline(sd, word, ' ');   
      auto it = m_sensor_lookup.find( word );
      if (it == m_sensor_lookup.end()) {
        sensorID = -1;
      }
      sensorID = it->second;  
      validHit = false;
    } else if ( line.substr(0, 5) == "Pixel" ) {
      // Pixel data: Read pixel hit
      stringstream sd(line);
      getline(sd, word, ' ');
      getline(sd, word, ',');
      column=stoi(word);
      getline(sd, word, ',');
      row=stoi(word);
      getline(sd, word, ','); 
      charge=stoi(word);     
      time=0;
      if (event_number >= 0 && sensorID >= 0) {
        validHit = true;
      }
    } else {
      // Something else
      validHit = false;
    }

    if (!validHit) {
      continue;
    }
 
    if (current_event == -1) {
      //cout << "START CHUNK WITH EVT=" << event_number << endl;
      current_event = event_number;
    }
    
    if (current_event !=  event_number ) {
      //cout << "FINISH EVT=" << current_event << " firstHit=" << firstHit << " lastHit=" << hits.size()-1  << endl;
      events.push_back( vector<int>{ current_event, firstHit, (int)hits.size()-1 } );
      
      if (event_number > max_event_number and max_event_number >=0) {
        //cout << "CURRENT EVT=" << event_number << " BIGGER THAN MAX EVT=" << max_event_number << ". LAST EVT FOR CURRENT CHUNK." << endl;
        // Get back to the position before reading the current hit 
        fin.seekg (oldpos);  
        return true;
      }
      
      if ( events.size() >= max_events) {
        // cout << "LAST EVT FOR CURRENT CHUNK." << endl;
        // Get back to the position before reading the current hit 
        fin.seekg (oldpos);  
        return true;
      }
      
      //cout << "START NEW EVENT!!" << endl;
      firstHit=hits.size();
    }     
    
    // Cache last event number
    current_event = event_number;
      
    // Buffer event data 
    hits.push_back(vector<int>{event_number, sensorID, column, row, charge, time} );
    
    // Stores the position in file stream
    oldpos = fin.tellg();  
  }
  
  // Reached end of file, nothing left to read 
  return false;
}


void CorryInputProcessor::end () 
{ 
  
  // Print message
  streamlog_out(MESSAGE3) << "Number of data events: " << to_string(m_ntriggers) << std::endl;
  streamlog_out(MESSAGE3) << std::endl
                          << "Processor succesfully finished!"
                          << std::endl;
}



