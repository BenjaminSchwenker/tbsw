/*=====================================================================*/
/*          AsciiInputProcessor                                        */
/*                                                                     */
/*          Author: Benjamin Schwenker                                 */
/*                (benjamin.schwenker@phys.uni-goettingen.de)          */
/*                                                                     */
/*  Processers converts raw digits in Hits tables from one or more     */ 
/*  ascii input files into series of lcio::LCEvents. In case multiple   */
/*  ascii input files are are supplied, raw digits are aligned using    */
/*  the event number.                                                  */
/*=====================================================================*/

// user includes
#include "AsciiInputProcessor.h"

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



AsciiInputProcessor::AsciiInputProcessor ():DataSourceProcessor("AsciiInputProcessor") {
  
  _description =
    "Reads Hits table from one or more ascii files. Events are formed and hits are written to TrackerData collections.\n"
    "Make sure to not specify any LCIOInputFiles.";
  
  std::vector< string > inputFileNameVecExample;
  inputFileNameVecExample.push_back("input.txt");
  registerProcessorParameter ("FileNames", "Whitespace seperated list of input file names to read",
                              m_fileNameVec, inputFileNameVecExample);
  
  registerProcessorParameter( "DetectorName", "Set name of the detector. Needs to be the same name as in gear file.",
                              m_detectorName, std::string("EUTelescope"));

  registerProcessorParameter( "RawHitCollectionName", "Set name of the raw hit collection.",
                              m_rawHitCollectionName, std::string("zsdata"));

  registerProcessorParameter ("RunNumber", "Set run number.",
                              m_runNumber,  static_cast < int > (0));
    
}


AsciiInputProcessor * AsciiInputProcessor::newProcessor () {  
  return new AsciiInputProcessor;
}


void AsciiInputProcessor::init () {
  m_chunkSize = 1000;
  printParameters ();
}


void AsciiInputProcessor::readDataSource (int Ntrig) {
  
  if (m_fileNameVec.size() == 0) return;
  string filename = m_fileNameVec[0];
   
  // Write LCIO run header here 
  IMPL::LCRunHeaderImpl* lcHeader = new IMPL::LCRunHeaderImpl;
  lcHeader->setDescription(" Reading from file " + filename);
  lcHeader->setRunNumber( m_runNumber );  
  lcHeader->setDetectorName(m_detectorName);
          
  // Add run header to LCIO file
  ProcessorMgr::instance ()->processRunHeader ( static_cast<lcio::LCRunHeader*> (lcHeader) );
  delete lcHeader; 
  
  streamlog_out ( MESSAGE3 ) << "Start reading next chunk from raw file " << filename << std::endl;
     
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
  while(readNextChunk) {  
    
    // Read chunk of data from the first file 
    vector< vector<int> > hits;
    vector< vector<int> > events;
    readNextChunk = getAllEventFromFile(fin, hits, events, m_chunkSize, -1);  
    
    // Read data from all other files 
    vector< vector< vector<int> > > all_other_hits;
    vector< vector< vector<int> > > all_other_events;
    
    // now add data from other file  
    for (size_t iFile=1; iFile< m_fileNameVec.size(); iFile++) {
      string other_filename = m_fileNameVec[iFile];
      
      streamlog_out ( MESSAGE3 ) << "Start reading next chunk from raw file " << other_filename << std::endl;
      
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
        
      //cout << "evt=" << events[iEvt][0] << endl;
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
        
        //cout << "   ipl=" << hits[iHit][1] << endl;
        //cout << "   col=" << hits[iHit][2] << endl;
        //cout << "   row=" << hits[iHit][3] << endl;
        //cout << "   charge=" << hits[iHit][4] << endl;
        
        chargevec.push_back( hits[iHit][1] );
        chargevec.push_back( hits[iHit][2] );
        chargevec.push_back( hits[iHit][3] );
        chargevec.push_back( hits[iHit][4] );
      }
      
      // now add data from other file 
      for (size_t iFile=1; iFile<m_fileNameVec.size(); iFile++) {   
        for (size_t oEvt=0; oEvt<all_other_events[iFile-1].size(); oEvt++) {  
          if (all_other_events[iFile-1][oEvt][0] == events[iEvt][0]) { 
            for (int iHit = all_other_events[iFile-1][oEvt][1]; iHit<=all_other_events[iFile-1][oEvt][2]; iHit++ ) {  
              //cout << "   ipl=" << all_other_hits[iFile-1][iHit][1] << endl;
              //cout << "   col=" << all_other_hits[iFile-1][iHit][2] << endl;
              //cout << "   row=" << all_other_hits[iFile-1][iHit][3] << endl;
              //cout << "   charge=" << all_other_hits[iFile-1][iHit][4] << endl;
              
              chargevec.push_back( all_other_hits[iFile-1][iHit][1] );
              chargevec.push_back( all_other_hits[iFile-1][iHit][2] );
              chargevec.push_back( all_other_hits[iFile-1][iHit][3] );
              chargevec.push_back( all_other_hits[iFile-1][iHit][4] );
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

bool AsciiInputProcessor::getAllEventFromFile(fstream& fin, vector< vector<int> >& hits, std::vector< std::vector<int> >& events, int max_events, int max_event_number) 
{
  
  string line, word; 
  int event_number=-1;  
  int plane=-1; 
  int column=-1; 
  int row=-1; 
  int charge=-1;   
  int firstHit = 0;
  int current_event = -1; 
  
  // stores the position
  streampos oldpos = fin.tellg();
  
  while (getline(fin, line)) {
             
    // used for breaking words 
    stringstream sd(line); 
     
    getline(sd, word, ',');    
    event_number=stoi(word);
    getline(sd, word, ',');
    plane=stoi(word);
    getline(sd, word, ',');
    column=stoi(word);
    getline(sd, word, ',');
    row=stoi(word);
    getline(sd, word, ','); 
    charge=stoi(word);       
    
    if (current_event == -1) {
      //cout << "START CHUNK WITH EVT=" << event_number << endl;
      current_event = event_number;
    }
    
    if (current_event !=  event_number ) {
      //cout << "FINISH EVT=" << current_event << " firstHit=" << firstHit << " lastHit=" << hits.size()-1  << endl;
      events.push_back( vector<int>{ current_event, firstHit, hits.size()-1 } );
      
      if (event_number > max_event_number and max_event_number >=0) {
        //cout << "CURRENT EVT=" << event_number << " BIGGER THAN MAX EVT=" << max_event_number << ". LAST EVT FOR CURRENT CHUNK." << endl;
        // Get back to the position before reading the current hit 
        fin.seekg (oldpos);  
        return true;
      }
      
      if (events.size() >= max_events) {
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
    hits.push_back(vector<int>{event_number, plane, column, row, charge} );
    
    // Stores the position in file stream
    oldpos = fin.tellg();  
  }
  
  // Reached end of file, nothing left to read 
  return false;
}


void AsciiInputProcessor::end () 
{ 
  // Print message
  streamlog_out(MESSAGE3) << std::endl
                          << "Processor succesfully finished!"
                          << std::endl;
}



