// NIUnpacker implementation file
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "NIUnpacker.h"
#include <iomanip>
#include "Utils.hh"

#define GET(d, i) getlittleendian<unsigned>(&(d)[(i)*4])

// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;

namespace {
  int PIVOTPIXELOFFSET = 64;
}

namespace eudaqinput {

  //
  // Instantiate this object
  //
  NIUnpacker aNIUnpacker ;

  //
  // Constructor
  //
  NIUnpacker::NIUnpacker() : Processor("NIUnpacker")
  {
    
    // Processor description
    _description = "NIUnpacker: Unpacker for M26 sensors in EUDET/AIDA telescope raw data" ;
    
    //   
    // First of all, we need to register the input/output collections
    
    registerInputCollection (LCIO::TRACKERRAWDATA, "InputCollectionName",
                             "Name of collection containing raw data",
                             _inputCollectionName, string("NI"));
    
    registerOutputCollection (LCIO::TRACKERDATA, "OutputCollectionName",
                              "Name of the output digit collection",
                              _outputCollectionName, string("zsdata_m26"));
  
  }
  
  //
  // Method called at the beginning of data processing
  //
  void NIUnpacker::init() {
     
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
  void NIUnpacker::processRunHeader(LCRunHeader * run) {
    
    // Print run number
    streamlog_out(MESSAGE3) << "Processing run: "
                            << (run->getRunNumber())
                            << std::endl << std::endl;
    
    _nRun++ ;  
  }
  
  //
  // Method called for each event
  //
  void NIUnpacker::processEvent(LCEvent * evt)
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
  
  
  bool NIUnpacker::UnpackRawCollection(LCCollectionVec * result, LCCollectionVec * source){
     
    // Helper class for decoding raw data 
    CellIDDecoder<TrackerRawDataImpl> inputDecoder( source );  
    
    // Helper class for encoding digit data 
    CellIDEncoder<TrackerDataImpl> outputEncoder( "sensorID:6,sparsePixelType:5", result  );
    
    // Check number of frames 
    if ( source->size() != 2 ) {
      streamlog_out(MESSAGE3) << "Ignoring bad event " << _nEvt << std::endl;
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
    
    // to understand if we have problem with de-syncronisation, let
    // me prepare a vector of size_t to contain the
    // pivot pixel position
    std::vector<size_t > pivotPixelPosVec;
    
    unsigned header0 = GET(data0, 0);
    unsigned header1 = GET(data1, 0);
    unsigned pivot = GET(data0, 1) & 0xffff;
    datait it0 = data0.begin() + 8;
    datait it1 = data1.begin() + 8;
    unsigned board = 0;
    while (it0 < data0.end() && it1 < data1.end()) {
      unsigned id = board;
      if (it0 + 2 >= data0.end()) {
        std::cout << "Trailing rubbish in first frame" << std::endl;
        break;
      }
      if (it1 + 2 >= data1.end()) {
        std::cout << "Trailing rubbish in second frame" << std::endl;
        break;
      }
      unsigned len0 = GET(it0, 1);
      if ((len0 & 0xffff) != (len0 >> 16)) {
        std::cout << "Mismatched lengths decoding first frame (" +
                   to_string(len0 & 0xffff) + ", " + to_string(len0 >> 16) + ")" << std::endl;
        len0 = std::max(len0 & 0xffff, len0 >> 16);
      }
      len0 &= 0xffff;
      unsigned len1 = GET(it1, 1);
      if ((len1 & 0xffff) != (len1 >> 16)) {
        std::cout << "Mismatched lengths decoding second frame (" +
                   to_string(len1 & 0xffff) + ", " + to_string(len1 >> 16) + ")" << std::endl;
        len1 = std::max(len1 & 0xffff, len1 >> 16);
      }
      len1 &= 0xffff;
      if (it0 + len0*4 + 12 > data0.end()) {
        std::cout << "Bad length in first frame" << std::endl;
        break;
      }
      if (it1 + len1*4 + 12 > data1.end()) {
        std::cout << "Bad length in second frame" << std::endl;
        break;
      }

      unsigned pivotpixel = (9216 + pivot + PIVOTPIXELOFFSET) % 9216;
      pivotPixelPosVec.push_back( pivotpixel );
       
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
      
      it0 += len0*4 + 16;
      it1 += len1*4 + 16;
      if (it0 <= data0.end()) header0 = GET(it0, -1);
      if (it1 <= data1.end()) header1 = GET(it1, -1);
      ++board;
    }
    
     
    // check if all the boards where running in synchronous mode or
    // not. Remember that the last pivot pixel entry is the one of the
    // master board.
    bool outOfSyncFlag = false;
    std::vector<size_t >::iterator masterBoardPivotAddress = pivotPixelPosVec.end() - 1;
    std::vector<size_t >::iterator slaveBoardPivotAddress  = pivotPixelPosVec.begin();
    while ( slaveBoardPivotAddress < masterBoardPivotAddress ) {
      if ( *slaveBoardPivotAddress - *masterBoardPivotAddress >= 2 ) {
        outOfSyncFlag = true;

        // we don't need to continue looping over all boards if one of
        // them is already out of sync
        break;
      }
      ++slaveBoardPivotAddress;
    }
    if ( outOfSyncFlag ) {

      if ( _nEvt  < 20 ) {
        // in this case we have the responsibility to tell the user that
        // the event was out of sync
        std::cout << "Event number " << _nEvt << " seems to be out of sync" << std::endl;
        std::vector<size_t >::iterator masterBoardPivotAddress = pivotPixelPosVec.end() - 1;
        std::vector<size_t >::iterator slaveBoardPivotAddress  = pivotPixelPosVec.begin();
        while ( slaveBoardPivotAddress < masterBoardPivotAddress ) {
          // print out all the slave boards first
          std::cout << " --> Board (S) " <<  std::setw(3) << setiosflags( std::ios::right )
                    << slaveBoardPivotAddress - pivotPixelPosVec.begin() << resetiosflags( std::ios::right )
                    << " = " << std::setw(15) << setiosflags( std::ios::right )
                    << (*slaveBoardPivotAddress) << resetiosflags( std::ios::right )
                    << " (" << std::setw(15) << setiosflags( std::ios::right )
                    << (signed) (*masterBoardPivotAddress) - (signed) (*slaveBoardPivotAddress) << resetiosflags( std::ios::right)
                    << ")" << std::endl;
          ++slaveBoardPivotAddress;
        }
        // print out also the master. It is impossible that the master
        // is out of sync with respect to itself, but for completeness...
        std::cout  << " --> Board (M) "  <<  std::setw(3) << setiosflags( std::ios::right )
                   << slaveBoardPivotAddress - pivotPixelPosVec.begin() << resetiosflags( std::ios::right )
                   << " = " << std::setw(15) << setiosflags( std::ios::right )
                   << (*slaveBoardPivotAddress) << resetiosflags( std::ios::right )
                   << " (" << std::setw(15)  << setiosflags( std::ios::right )
                   << (signed) (*masterBoardPivotAddress) - (signed) (*slaveBoardPivotAddress) << resetiosflags( std::ios::right)
                   << ")" << std::endl;

      } else if ( _nEvt  == 20 ) {
        // if the number of consecutive warnings is equal to the maximum
        // allowed, don't bother the user anymore with this message,
        // because it's very likely the run was taken unsynchronized on
        // purpose
        std::cout << "The maximum number of unsychronized events has been reached." << std::endl
                  << "Assuming the run was taken in asynchronous mode" << std::endl;
      }
    }
     
    return true;
  }
  
  void NIUnpacker::DecodeFrame(TrackerDataImpl* zsFrame, size_t len, datait it, int frame, unsigned pivotpixel) 
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
  void NIUnpacker::check( LCEvent * evt )
  {
  }
  
  //
  // Method called after all data processing
  //
  void NIUnpacker::end()
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
  void NIUnpacker::printProcessorParams() const
  { 
    streamlog_out(MESSAGE3)  << std::endl
                             << " "
                             << "NIUnpacker Development Version, be carefull!!"
                             << " "
                             << std::endl  << std::endl;   
  
  }
} // Namespace



