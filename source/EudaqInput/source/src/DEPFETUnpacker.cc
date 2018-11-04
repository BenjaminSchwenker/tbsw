// DEPFETHybrid5Unpacker implementation file
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "DEPFETUnpacker.h"
#include "Utils.hh"

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <stdint.h>


// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;





namespace eudaqinput {

  //
  // Instantiate this object
  //
  static DEPFETUnpacker aDEPFETUnpacker ;

  //
  // Constructor
  //
  DEPFETUnpacker::DEPFETUnpacker() : Processor("DEPFETUnpacker")
  {
    
    // Processor description
    _description = "DEPFETUnpacker: Unpacker for DEPFET detectors (small and large Belle2 PXD9 sensor) raw data" ;
    
    //   
    // First of all, we need to register the input/output collections
    
    registerInputCollection (LCIO::TRACKERRAWDATA, "InputCollectionName",
                             "Name of collection containing raw data",
                             _inputCollectionName, string("DEPFET"));
    
    registerOutputCollection (LCIO::TRACKERDATA, "OutputCollectionName",
                              "Name of the output digit collection",
                              _outputCollectionName, string("zsdata_dep"));
    
    
    registerProcessorParameter ("OverrideSensorID", "Override sensorID. Put -1 to keep sensorID from rawdata",
                                _resetSensorID,  static_cast < int > (-1));
    
    registerProcessorParameter ("moduleID", "ModuleID which should be read. Put -1 to keep all",
                                _modID,  static_cast < int > (-1));
    registerProcessorParameter ("dhpID", "DHP ID which should be read. Put -1 to keep all",
                                _dhpID,  static_cast < int > (-1));
    registerProcessorParameter ("dheID", "DHE ID which should be read. Put -1 to keep all",
                                _dheID,  static_cast < int > (-1));
    registerProcessorParameter ("Mapping", "Sensor type for mapping. Can be PXD9_[OF,IF,OB,IB], HYBRID5, automatic or none. Put 'automatic' to determine mapping from dheID. Case insensitive.",
                                _mappingString, static_cast< string > ( "automatic" ));
    registerProcessorParameter ("SwapAxis", "Swap rows and columns",
                                _swapAxes,  static_cast < bool > (false));

   
    
  }

  //
  // Method called at the beginning of data processing
  //
  void DEPFETUnpacker::init() {
    streamlog_out(MESSAGE3) << "Building interpreter for mapping: "
                                << _mappingString<< std::endl;
    interpreter=DepfetInterpreter(_modID,_dhpID,_dheID);
    //interpreter.setDebug(true);
    interpreter.setMapping(mappingFromString(_mappingString)); 
    interpreter.setSkipRaw(true);
    interpreter.swapRowCol(_swapAxes);
     
    // Initialize variables;
    _nRun = 0 ;
    _nEvt = 0 ;
                 
    // Print set parameters
    printProcessorParams();
     
    // CPU time start
    _timeCPU = clock()/1000;
  }

  //
  // Method called for each run
  //
  void DEPFETUnpacker::processRunHeader(LCRunHeader * run) {
    
    // Print run number
    streamlog_out(MESSAGE3) << "Processing run: "
                            << (run->getRunNumber())
                            << std::endl << std::endl;
    
    _nRun++ ;  
  }
  
  //
  // Method called for each event
  //
  void DEPFETUnpacker::processEvent(LCEvent * evt)
  {
    
    // Print event number
    if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                  << (evt->getEventNumber())
                                                                  << std::endl << std::endl;
     
    // More detailed event numbering for testing
    streamlog_out(MESSAGE2) << std::endl << "Starting with Event Number " << evt->getEventNumber()  << std::endl;  
    
    //
    // Open collections
    try {
      
      // Open raw data collection   
      LCCollectionVec * source = dynamic_cast < LCCollectionVec * > (evt->getCollection(_inputCollectionName)); 
       
      // Helper class for decoding raw data 
      CellIDDecoder<TrackerRawDataImpl> inputDecoder( source );  
      
      //-----------------------------------------------
      // Decode event data to a DEPFETEvent format   
      size_t numplanes =  source->size();  
      std::vector<depfet_event> all_events;
      for (size_t iPlane = 0; iPlane < numplanes; ++iPlane) { 
      
        TrackerRawDataImpl* rawDataBlock = dynamic_cast<TrackerRawDataImpl* > ( source->getElementAt(iPlane) );
        
        // Get blockID 
        int blockID = inputDecoder(rawDataBlock)["blockID"]; 

        // FIXME: This is only needed, because i do not know how to put 
        // vector<unsigned char> directly into the lcio::LCEvent. 
        datavect  char_data;
        ShortVec& short_data = rawDataBlock->adcValues();
        std::copy(std::begin(short_data),std::end(short_data),std::back_inserter(char_data));

        streamlog_out(MESSAGE2) << "ConvertDEPFETEvent() of block id " << blockID << " of length " << char_data.size() << std::endl;
        auto start=all_events.size();
        interpreter.Interprete(all_events,&(char_data[0]),char_data.size());
        for(auto i=start;i<all_events.size();i++){
			all_events[i].modID=blockID;		
		}
      }
      
      // Create DEPFET_EVENT_INFO if it does not exist
      try {
        LCCollectionVec* eventinfo = dynamic_cast < LCCollectionVec * > (evt->getCollection( "DEPFET_EVENT_INFO" )) ;
      } catch (lcio::DataNotAvailableException& e) {
        LCCollectionVec* depfet_info = new LCCollectionVec( LCIO::LCGENERICOBJECT )  ;
        
        // Decode DEPFET event info to LCIO format      
        for(const auto & event:all_events){
            LCGenericObjectImpl* metaobj = new LCGenericObjectImpl(2,0,0);
            metaobj->setIntVal(0,event.isGood);
            metaobj->setIntVal(1,event.triggerOffset);
        }
        evt->addCollection( depfet_info , "DEPFET_EVENT_INFO" ) ;
      }   
      
      streamlog_out(MESSAGE2) << "Create LCIO::TrackerData collection for DEPFET." << std::endl;
      
      // Prepare collection for unpacked digits 
      LCCollectionVec * result = new LCCollectionVec(LCIO::TRACKERDATA);

      // Helper class for encoding digit data 
      CellIDEncoder<TrackerDataImpl> outputEncoder( "sensorID:6,sparsePixelType:5", result  );
      
      //-----------------------------------------------
      // Decode event data to a LCIO format
      
      for(const auto & event:all_events){
          // Prepare a TrackerData to store zs pixels
          TrackerDataImpl* zsFrame = new TrackerDataImpl;
          // Set description for zsFrame
          if (_resetSensorID >= 0) {
            outputEncoder["sensorID"] = _resetSensorID;
          } else {
            outputEncoder["sensorID"] = event.modID;
          }
          
          //streamlog_out(MESSAGE2) << "Write data for moduleID " << outputEncoder["sensorID"] << " to lcio collection " <<  _outputCollectionName << "." << std::endl;
            
          outputEncoder["sparsePixelType"] = 0;
          outputEncoder.setCellID( zsFrame );
          auto & charge_values=zsFrame->chargeValues();
          charge_values.reserve(charge_values.size()+3*event.zs_data.size());
          for ( const auto & hit : event.zs_data ) {
              // Store pixel data int lcio format
              charge_values.push_back( hit.col );
              charge_values.push_back( hit.row );
              charge_values.push_back( hit.val );
              //streamlog_out(MESSAGE2) << "Hit at col= " << hit.col << ", row=" <<  hit.row << ", val=" << hit.val << std::endl;
          }
          result->push_back(zsFrame);
      }

      // Add output collection to event
      evt->addCollection( result, _outputCollectionName );
             	    
    } catch(DataNotAvailableException &e){}  
    _nEvt ++ ;
  }

  // Method called after each event to check the data processed
  //
  void DEPFETUnpacker::check( LCEvent * evt )
  {
  }
  
  //
  // Method called after all data processing
  //
  void DEPFETUnpacker::end()
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
  void DEPFETUnpacker::printProcessorParams() const
  { 
    streamlog_out(MESSAGE3)  << std::endl
                             << " "
                             << "DEPFETUnpacker Development Version, be carefull!!"
                             << " "
                             << std::endl  <<"Parameters:"<< std::endl;

    
     streamlog_out(MESSAGE3) << std::endl;
  }
} // Namespace



