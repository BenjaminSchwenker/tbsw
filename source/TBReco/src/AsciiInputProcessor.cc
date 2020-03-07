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
#include <fstream>

using namespace std;
using namespace marlin;
using namespace depfet;
using namespace std::string_literals;



AsciiInputProcessor::AsciiInputProcessor ():DataSourceProcessor("AsciiInputProcessor"), _outputEncoderHelper( "sensorID:6,sparsePixelType:5") {
  
  _description =
    "Reads Hits table from one or more ascii files. Events are formed and hits are written to TrackerData collections.\n"
    "Make sure to not specify any LCIOInputFiles.";
  
  std::vector< string > inputFileNameVecExample;
  inputFileNameVecExample.push_back("input.txt");
  registerProcessorParameter ("FileNames", "Whitespace seperated list of input file names to read",
                              m_fileNameVec, inputFileNameVecExample);
  
  registerProcessorParameter( "DetectorName", "Set name of the detector. Needs to be the same name as in gear file.",
                              m_detectorName, std::string("EUTelescope"));
    
}


AsciiInputProcessor * AsciiInputProcessor::newProcessor () {  
  return new AsciiInputProcessor;
}


void AsciiInputProcessor::init () {
  printParameters ();
}


void AsciiInputProcessor::readDataSource (int Ntrig) {
  
  for (auto filename : m_fileNameVec) {

    streamlog_out ( MESSAGE3 ) << "Start reading raw file " << filename << std::endl;
     
    // Open file stream
    fstream fin(filename);
    
    // Read the Data from the file 
    // as String Vector 
    vector<int> row; 
    vector<string> col_names;
    string line, word, temp; 
    
    // Read first line with column names     
    getline(fin, line); 
  
    // used for breaking words 
    stringstream s(line); 
         
    // read every column name of first row and 
    // store it in a string variable, 'word' 
    while (getline(s, word, ',')) { 
  
      // add all the column data 
      // of a row to a vector 
      col_names.push_back(word); 
    }  
    
    for (auto name : col_names) cout << name << endl;
    size_t num_cols = col_names.size();
      
    // a table for all hits 
    vector< vector<int> > hits;
     
    while (getline(fin, line)) {
      
      // clear last row   
      row.clear(); 
           
      // used for breaking words 
      stringstream sd(line); 
         
      // read every column data of a row and 
      // store it in a string variable, 'word' 
      while (getline(sd, word, ',')) { 
  
        // add all the column data 
        // of a row to a vector 
        row.push_back(stoi(word)); 
      } 
         
      hits.push_back(row);     
    }
    
    vector<int> evr_eventNumber(hits.size());
    vector<int> evr_firstHit(hits.size());
    vector<int> evr_lastHit(hits.size());
    
    // Assume that first hit starts first event
    //cout << "START NEW EVENT!!" << endl;
    int currEvt = hits[0][0]; // TODO assuming event number is first
    int firstHit = 0;
    int nevt = 0;
    
    
    for (int i = 0; i<(int)hits.size(); i++) { 
       
      if (currEvt !=  hits[i][0] ) {
        
        evr_eventNumber[nevt] = currEvt;
        evr_firstHit[nevt] = firstHit;
        evr_lastHit[nevt] = i-1; 
        //cout << "FINISH EVENT: evt=" << evr_eventNumber[nevt] << " firstHit=" << evr_firstHit[nevt] << " lastHit=" << evr_lastHit[nevt] << endl;
        //cout << "START NEW EVENT!!" << endl;
        firstHit=i;
        nevt++;
        
      }      
      
      // Cache last event number
      currEvt = hits[i][0]; 
    }
      
    // Write LCIO run header here 
    IMPL::LCRunHeaderImpl* lcHeader = new IMPL::LCRunHeaderImpl;
    lcHeader->setDescription(" Reading from file " + filename);
    lcHeader->setRunNumber( 0 );  // TODO: how to know the run number???
    lcHeader->setDetectorName(m_detectorName);
          
    // Add run header to LCIO file
    ProcessorMgr::instance ()->processRunHeader ( static_cast<lcio::LCRunHeader*> (lcHeader) );
    delete lcHeader; 
    
    for (long i = 0; i< nevt; i++) {
      
      //cout << "benni 2 evt=" << evr_eventNumber[i] << endl;
      // Process data event 
      lcio::LCEventImpl * lcEvent = new lcio::LCEventImpl;
      lcEvent->setEventNumber(evr_eventNumber[i]);
      lcEvent->setRunNumber(0);  // TODO Again run number issue
      lcEvent->setTimeStamp(0);  // TODO how to know the timestamp?? not used anyway
      
      // prepare the collections for the raw data 
      LCCollectionVec * rawDataCollection = new LCCollectionVec( lcio::LCIO::TRACKERDATA );           
      // set the proper cell encoder
      CellIDEncoder<TrackerDataImpl> lcEncoder( "sensorID:6,sparsePixelType:5", rawDataCollection , &_outputEncoderHelper );
       
      // Prepare a new lcio::TrackerData for the ZS data
      lcio::TrackerDataImpl* zsFrame =  new lcio::TrackerDataImpl;
      lcEncoder["sensorID"s] = 0;  // TODO how to get sensorID???
      lcEncoder["sparsePixelType"s] = 0;
      lcEncoder.setCellID( zsFrame );
      
      auto& chargevec=zsFrame->chargeValues();
      chargevec.reserve(chargevec.size()+size_t(1.2* (evr_lastHit[i] - evr_firstHit[i] + 1 )  ));
      
      for (long j=evr_firstHit[i]; j<=evr_lastHit[i]; j++ ) {
 
        //cout << "   ipl=" << hits[j][1] << endl;
        //cout << "   col=" << hits[j][2] << endl;
        //cout << "   row=" << hits[j][3] << endl;
        //cout << "   charge=" << hits[j][4] << endl;
        
        chargevec.push_back( hits[j][1] );
        chargevec.push_back( hits[j][2] );
        chargevec.push_back( hits[j][3] );
        chargevec.push_back( hits[j][4] );
      }
      
      // Now add the TrackerData to the collection
      rawDataCollection->push_back( zsFrame );
      
      // Add output collection to event
      lcEvent->addCollection( rawDataCollection, "zsdata" ); // TODO how to get the collection name
      
      ProcessorMgr::instance ()->processEvent (static_cast<LCEventImpl*> (lcEvent));
      delete lcEvent;
      
    }
    
    
  } // end loop over raw files    
}


void AsciiInputProcessor::end () 
{ 
  /*
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
  */
}



