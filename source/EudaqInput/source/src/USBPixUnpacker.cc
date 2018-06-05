// USBPixUnpacker implementation file
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "USBPixUnpacker.h"
#include <iomanip>
#include "Utils.hh"


#define GET(d, i) getlittleendian<unsigned>(&(d)[(i)*4])

// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;


namespace eudaqinput {

  //
  // Instantiate this object
  //
  USBPixUnpacker aUSBPixUnpacker ;

  //
  // Constructor
  //
  USBPixUnpacker::USBPixUnpacker() : Processor("USBPixUnpacker")
  {
    
    // Processor description
    _description = "USBPixUnpacker: Unpacker for USBPix sensors in EUDET/AIDA telescope raw data" ;
    
    //   
    // First of all, we need to register the input/output collections
    
    registerInputCollection (LCIO::TRACKERRAWDATA, "InputCollectionName",
                             "Name of collection containing raw data",
                             _inputCollectionName, string("USBPIX_GEN3"));
    
    registerOutputCollection (LCIO::TRACKERDATA, "OutputCollectionName",
                              "Name of the output digit collection",
                              _outputCollectionName, string("zsdata_usbpix"));
  
  }
  
  //
  // Method called at the beginning of data processing
  //
  void USBPixUnpacker::init() {
     
    // Initialize variables
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
  void USBPixUnpacker::processRunHeader(LCRunHeader * run) {
    
    // Print run number
    streamlog_out(MESSAGE3) << "Processing run: "
                            << (run->getRunNumber())
                            << std::endl << std::endl;
    
    _nRun++ ;  
  }
  
  //
  // Method called for each event
  //
  void USBPixUnpacker::processEvent(LCEvent * evt)
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
      LCCollectionVec * inputCollection = dynamic_cast < LCCollectionVec * > (evt->getCollection(_inputCollectionName)); 
      
      // Prepare collection for unpacked digits 
      LCCollectionVec * outputCollection = new LCCollectionVec(LCIO::TRACKERDATA);
      
      // Unpack raw data collection
      UnpackRawCollection(outputCollection, inputCollection);
      
      // Add output collection to event
      evt->addCollection( outputCollection, _outputCollectionName );
             	    
    } catch(DataNotAvailableException &e){}  
    _nEvt ++ ;
  }
  
  
  bool USBPixUnpacker::UnpackRawCollection(LCCollectionVec * result, LCCollectionVec * source){
     
    // Helper class for decoding raw data 
    CellIDDecoder<TrackerRawDataImpl> inputDecoder( source );  
    
    // Helper class for encoding digit data 
    CellIDEncoder<TrackerDataImpl> outputEncoder( "sensorID:6,sparsePixelType:5", result  );
    
    // Check number of frames 
    if ( source->size() != 1 ) {
      streamlog_out(MESSAGE3) << "Ignoring bad event " << _nEvt << " having size " << source->size() << std::endl;
      return false;
    }
    
    // Read first data block 
    datavect  data0;
    ShortVec& short_data0 = dynamic_cast<TrackerRawDataImpl* > ( source->getElementAt(0) )->adcValues();
    
    for (size_t j = 0; j < short_data0.size(); j++) 
    {   
      data0.push_back( static_cast<unsigned char> (short_data0[j]) ); 
    }
    
	//getChannels will determine all the channels from a board, making the assumption that every channel (i.e. FrontEnd)
	//wrote date into the data block. This holds true if the FE is responding. Then for every trigger there will be
	//data headers (DHs) in the data stream
	auto channels = getChannels(data0);
    
    for (auto channel : channels) {
      std::cout << "found channel  " << channel << std::endl;
    }
     
    std::map<int, std::unique_ptr<lcio::TrackerDataImpl>> frameMap;
	std::map<int, std::unique_ptr<eutelescope::EUTelTrackerDataInterfacerImpl<eutelescope::EUTelGenericSparsePixel>>> frameInterfaceMap;
	 
    /*
    for(auto channel: boardChannels.at(boardID)){
		auto frame = std::unique_ptr<lcio::TrackerDataImpl>(new lcio::TrackerDataImpl);
		auto sensorID = 20 + channel;
		cellIDEncoder["sensorID"] = sensorID;
		cellIDEncoder.setCellID(frame.get());
		auto frameInterface = 	std::unique_ptr<eutelescope::EUTelTrackerDataInterfacerImpl<eutelescope::EUTelGenericSparsePixel>>( 
						new eutelescope::EUTelTrackerDataInterfacerImpl<eutelescope::EUTelGenericSparsePixel>(frame.get())
					);
		frameInterfaceMap[channel] = std::move(frameInterface);
		frameMap[channel] = std::move(frame);
	}
    */

	
    auto pixelVec = decodeFEI4Data(data0);    
    for(auto& hitPixel: pixelVec) {
		frameInterfaceMap[hitPixel.channel]->emplace_back(hitPixel.x, hitPixel.y, hitPixel.tot+1+hitDiscConf, hitPixel.lv1);
	}

	for(auto& framePair: frameMap){
		dataCollection->push_back( framePair.second.release() );
	}    
    

    // Prepare a new lcio::TrackerData for the ZS data
    lcio::TrackerDataImpl* zsFrame =  new lcio::TrackerDataImpl;
    outputEncoder["sensorID"] = id;
    outputEncoder["sparsePixelType"] = 0;
    outputEncoder.setCellID( zsFrame );
      
    // Fill ZS data into new lcio::TrackerData object
    /*  
    for (unsigned j = 0; j < num+1; ++j) {
      zsFrame->chargeValues().push_back( column+j );
      zsFrame->chargeValues().push_back( row );
      zsFrame->chargeValues().push_back( 1 );       
    }
    */
      
    // Now add the TrackerData to the collection
    result->push_back( zsFrame );  
     
    return true;
  }
  
  std::vector<int> USBPixUnpacker::getChannels(std::vector<unsigned char> const & data) 
  {
    std::vector<int> channels;
    for(size_t index = 0;  index < data.size(); index+=4) {
      uint32_t i =( static_cast<uint32_t>(data[index+3]) << 24 ) | 
                  ( static_cast<uint32_t>(data[index+2]) << 16 ) | 
				  ( static_cast<uint32_t>(data[index+1]) << 8 ) | 
				  ( static_cast<uint32_t>(data[index+0]) );
        
	  uint8_t channel = selectBits(i, 24, 8);
        
      if(!(channel >> 7)) { //Trigger
	    if(std::find(channels.begin(), channels.end(), static_cast<int>(channel)) == channels.end()) {
	      channels.emplace_back(static_cast<int>(channel));
		}
	  } 
	}
	return channels;
  }
  
  std::vector<APIXPix> USBPixUnpacker::decodeFEI4Data(std::vector<unsigned char> const & data) {
	std::vector<APIXPix> result;
	std::array<size_t,8> no_data_headers;
    
	for(size_t index = 0;  index < data.size(); index+=4) {
        uint32_t i =( static_cast<uint32_t>(data[index+3]) << 24 ) | 
                	( static_cast<uint32_t>(data[index+2]) << 16 ) | 
					( static_cast<uint32_t>(data[index+1]) << 8 ) | 
					( static_cast<uint32_t>(data[index+0]) );

		uint8_t channel = selectBits(i, 24, 8);

		if(channel >> 7) { //Trigger
		
			no_data_headers.fill(0);
			uint32_t trigger_number = selectBits(i, 0, 31); //testme
		} else {
			uint8_t type = selectBits(i, 16, 8);
			switch(type) {
				case data_header: {
					//int bcid = selectBits(i, 0, 10);
					//int lv1id = selectBits(i, 10, 5);
					//int flag = selectBits(i, 15, 1);
					no_data_headers.at(channel)++;
					break;
				}
				case address_record: {
					//int type = selectBits(i, 15, 1);
					//int address = selectBits(i, 0, 15);
					break;
				}
				case value_record: {
					//int value = selectBits(i, 0, 16);
					break;
				}
				case service_record: {
					//int code = selectBits(i, 10, 6);
					//int count = selectBits(i, 0, 10);
					break;
				}
				default: { //data record

					unsigned tot2 = selectBits(i, 0, 4);
					unsigned tot1 = selectBits(i, 4, 4);
					unsigned row = selectBits(i, 8, 9) - 1;
					unsigned column = selectBits(i, 17, 7) - 1;

					if(column < 80 && ((tot2 == 0xF && row < 336) || (tot2 < 0xF && row < 335) )) {
						//If tot2 != 0b1111 (0xF) then the tot2 is the tot code for pixel (col, row+1)
                        auto lv1 = static_cast<unsigned>(no_data_headers.at(channel)-1);
                        if(lv1 > 15) lv1 = 15;
						if( tot2 != 0xF) result.emplace_back(column, row+1, tot2, lv1, channel);
						result.emplace_back(column, row, tot1, lv1, channel);
					} else { // invalid data record
						//invalid_dr ++;
					}
					break;
				}
			}
		}
	}
	return result;
  }
 
  //
  // Method called after each event to check the data processed
  //
  void USBPixUnpacker::check( LCEvent * evt )
  {
  }
  
  //
  // Method called after all data processing
  //
  void USBPixUnpacker::end()
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
  void USBPixUnpacker::printProcessorParams() const
  { 
    streamlog_out(MESSAGE3)  << std::endl
                             << " "
                             << "USBPixUnpacker Development Version, be carefull!!"
                             << " "
                             << std::endl  << std::endl;   
  
  }
} // Namespace



