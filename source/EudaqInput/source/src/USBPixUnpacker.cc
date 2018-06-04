// USBPixUnpacker implementation file
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "USBPixUnpacker.h"
#include <iomanip>
#include "Utils.hh"


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
    if ( source->size() != 2 ) {
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
    
    // Read second data block 
    datavect data1;
    ShortVec& short_data1 = dynamic_cast<TrackerRawDataImpl* > ( source->getElementAt(1) )->adcValues();
    
    for (size_t j = 0; j < short_data1.size(); j++) 
    {   
      data1.push_back( static_cast<unsigned char> (short_data1[j]) ); 
    }
    
    while (it0 < data0.end() && it1 < data1.end()) {
    
      // Prepare a new lcio::TrackerData for the ZS data
      lcio::TrackerDataImpl* zsFrame =  new lcio::TrackerDataImpl;
      outputEncoder["sensorID"] = id;
      outputEncoder["sparsePixelType"] = 0;
      outputEncoder.setCellID( zsFrame );
      
      // Fill ZS data into new lcio::TrackerData object
      DecodeFrame(zsFrame, len0, it0+8, 0, pivotpixel);
      DecodeFrame(zsFrame, len1, it1+8, 1, pivotpixel);
      
      // Now add the TrackerData to the collection
      result->push_back( zsFrame );  
      
    } 
    
     
    
   
     
    return true;
  }
  
  void USBPixUnpacker::DecodeFrame(TrackerDataImpl* zsFrame, size_t len, datait it, int frame, unsigned pivotpixel) 
  {
    std::vector<unsigned short> vec;
    for (size_t i = 0; i < len; ++i) {
      unsigned v = GET(it, i);
      vec.push_back(v & 0xffff);
      vec.push_back(v >> 16);
    }
      
    unsigned npixels = 0;
    for (size_t i = 0; i < vec.size(); ++i) {
      //  std::cout << "  " << i << " : " << hexdec(vec[i]) << std::endl;
      if (i == vec.size() - 1) break;
      unsigned numstates = vec[i] & 0xf;
      unsigned row = vec[i]>>4 & 0x7ff;
      if (numstates+1 > vec.size()-i) {
        // Ignoring bad line
        //std::cout << "Ignoring bad line " << row << " (too many states)" << std::endl;
        break;
      }
        
      bool pivot = (row >= (pivotpixel / 16));
      for (unsigned s = 0; s < numstates; ++s) {
        unsigned v = vec.at(++i);
        unsigned column = v>>2 & 0x7ff;
        unsigned num = v & 3;
          
        for (unsigned j = 0; j < num+1; ++j) {
            zsFrame->chargeValues().push_back( column+j );
            zsFrame->chargeValues().push_back( row );
            zsFrame->chargeValues().push_back( 1 );       
        }
        npixels += num + 1;
      }
       
    }
    //streamlog_out(MESSAGE2) << "Total pixels " << frame << " = " << npixels << std::endl;
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



