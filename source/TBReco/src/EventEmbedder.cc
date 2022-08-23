// EventEmbedder Processor  
// 
// See EventEmbedder.h for full documentation of processor. 
// 
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "EventEmbedder.h"

// Include basic C
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>

// Include LCIO classes
#include <lcio.h>
#include "IO/LCReader.h"


// Used namespaces
using namespace std; 
using namespace lcio ;
using namespace marlin ;

namespace depfet {

//
// Instantiate this object
//
EventEmbedder aEventEmbedder ;

//
// Constructor
//
EventEmbedder::EventEmbedder() : Processor("EventEmbedder"), _inputDecodeHelper(""), _outputEncoderHelper( "sensorID:6,sparsePixelType:5")
{

// Processor description
   _description = "EventEmbedder: Embed event data from other file into current LCEvent" ;

//
// Processor parameters

  registerProcessorParameter( "OtherFileName",
                               "Other file name, from where to add new collections",
                               _otherFileName,
                               std::string("otherfile.slcio"));
   
  registerProcessorParameter( "TiggerFileName",
                               "Input file with triggered event numbers",
                               _triggerFileName,
                               std::string("dummy.txt"));

  registerProcessorParameter( "PixelDataCollection",
                               "Name of pixel data collection to copy",
                               _DigitCollectionName,
                               std::string("zsdata"));  
  
}

//
// Method called at the beginning of data processing
//
void EventEmbedder::init() {

  // Initialize variables
  _nRun = 0 ;
  _nEvt = 0 ;
   
  // Print set parameters
  printProcessorParams();
   
   
  // CPU time start
  _timeCPU = clock()/1000;
   
  // Read trigger file 
  _triggers.clear();  
       
  string line;
  ifstream triggerFile( _triggerFileName.c_str() );
  if ( triggerFile.is_open() ) { 
       
    streamlog_out(MESSAGE3) << "Reading trigger file ... "<< endl; 
    while (! triggerFile.eof() ) {
      getline(triggerFile ,line);
      streamlog_out(MESSAGE3) << "reading line: "<< line; 

      int i = -1;   
      if( sscanf( line.c_str() , "%d", &i)  == 1 && i>-1 ){
        streamlog_out(MESSAGE2) << "  converted eventID "<< i << endl;  
        _triggers.push_back(i);
      } else {
        streamlog_out(MESSAGE2) << "  bad eventID, skip " << endl;  
      }
    }
    triggerFile.close();
       
    // sort triggers to improve search speed
    sort (_triggers.begin(), _triggers.end());
       
  } else {
    streamlog_out(ERROR) << "Unable to open trigger file. " << endl;  
    _triggers.push_back(0);
  }    
}

//
// Method called for each run
//
void EventEmbedder::processRunHeader(LCRunHeader * run)
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
void EventEmbedder::processEvent(LCEvent * evt)
{
   
  // Print event number
  if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                << (evt->getEventNumber())
                                                                << std::endl << std::endl;

  // create lcio reader 
  LCReader* lcReader = LCFactory::getInstance()->createLCReader() ;

  // open file to read data
  lcReader->open(_otherFileName); 
  LCEvent * otherEvent = lcReader->readNextEvent();

  // copy digit collection into current event
  LCCollectionVec * copiedDigitCollection = new LCCollectionVec(LCIO::TRACKERDATA);
  copy_digits(otherEvent, copiedDigitCollection);    
  evt->addCollection( copiedDigitCollection, _DigitCollectionName );

  // close reader
  lcReader->close();
  delete lcReader;

  
  _nEvt ++ ;
}

//
// Method called after each event to check the data processed
//
void EventEmbedder::check( LCEvent * )
{
}

//
// Method called after all data processing
//
void EventEmbedder::end()
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
void EventEmbedder::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "EventEmbedder Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}


void EventEmbedder::copy_digits( LCEvent * evt , LCCollectionVec * digitCollection) 
{  
  // Open zero suppressed pixel data  
  LCCollectionVec * Pix_collection = dynamic_cast < LCCollectionVec * > (evt->getCollection(_DigitCollectionName)); 
  // Helper class for decoding pixel data 
  CellIDDecoder<TrackerDataImpl> PixelID( Pix_collection,&_inputDecodeHelper );
  // Helper class for encoding copied pixel data
  CellIDEncoder<TrackerDataImpl> digitEncoder( "sensorID:6,sparsePixelType:5", digitCollection, &_outputEncoderHelper);

  // Loop over pixel detectors 
  for (unsigned int iDet = 0; iDet < Pix_collection->size(); iDet++) { 
    // Get zs pixels from next pixel detector   
    TrackerDataImpl * pixModule = dynamic_cast<TrackerDataImpl* > ( Pix_collection->getElementAt(iDet) );

    // Sensor ID for pixel detector
    int sensorID = PixelID( pixModule ) ["sensorID"s];

    // List of firing pixels. Each pixel has a iU, iV, charge and time 
    FloatVec pixVector = pixModule->getChargeValues();

    // Copy into digitCollection ...
    TrackerDataImpl* zsFrame = new TrackerDataImpl; 

    float telEvents = 99605;
    float dutEvents = 1863;
    float ratio = dutEvents/telEvents;
    int pivot_index = int(ratio*_nEvt);

    cout << "curent event " << _nEvt << " pivot index " << pivot_index << endl;

    int nDigits = pixVector.size()/4; 
    for (int iDigit = 0; iDigit < nDigits; iDigit++) { 

      if ( (iDigit > pivot_index - 160) && (iDigit < pivot_index +160)  ) {
        zsFrame->chargeValues().push_back( pixVector[iDigit * 4] );
        zsFrame->chargeValues().push_back( pixVector[iDigit * 4 +1] );
        zsFrame->chargeValues().push_back( pixVector[iDigit * 4 +2] );
        zsFrame->chargeValues().push_back( pixVector[iDigit * 4 +3] );   
      }
    }

    // Add copied digits to digit collection
    digitEncoder["sensorID"s] = sensorID;
    digitEncoder["sparsePixelType"s] = 0;
    digitEncoder.setCellID(zsFrame);
    digitCollection->push_back(zsFrame);
  }

}

} // Namespace
