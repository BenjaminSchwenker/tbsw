// PyBARUnpacker implementation file
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "PyBARUnpacker.h"
#include <iomanip>
#include "Utils.hh"
#include "PyBAR.h"


// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;


namespace eudaqinput {

  //
  // Instantiate this object
  //
  PyBARUnpacker aPyBARUnpacker ;

  //
  // Constructor
  //
  PyBARUnpacker::PyBARUnpacker() : Processor("PyBARUnpacker")
  {
    
    // Processor description
    _description = "PyBARUnpacker: Unpacker for USBPix sensors in EUDET/AIDA telescope raw data" ;
    
    //   
    // First of all, we need to register the input/output collections
    
    registerInputCollection (LCIO::TRACKERRAWDATA, "InputCollectionName",
                             "Name of collection containing raw data",
                             _inputCollectionName, string("USBPIX_GEN3"));
    
    registerOutputCollection (LCIO::TRACKERDATA, "OutputCollectionName",
                              "Name of the output digit collection",
                              _outputCollectionName, string("zsdata_usbpix"));

    registerProcessorParameter( "SensorID","Set SensorID for data",
                               _sensorID, static_cast<int > (21));
     
  }
  
  //
  // Method called at the beginning of data processing
  //
  void PyBARUnpacker::init() {
     
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
  void PyBARUnpacker::processRunHeader(LCRunHeader * run) {
    
    // Print run number
    streamlog_out(MESSAGE3) << "Processing run: "
                            << (run->getRunNumber())
                            << std::endl << std::endl;
    
    _nRun++ ;  
  }
  
  //
  // Method called for each event
  //
  void PyBARUnpacker::processEvent(LCEvent * evt)
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
  
  
  bool PyBARUnpacker::UnpackRawCollection(LCCollectionVec * result, LCCollectionVec * source){
     
    // Helper class for encoding digit data 
    CellIDEncoder<TrackerDataImpl> outputEncoder( "sensorID:6,sparsePixelType:5", result  );
    
    // Check number of frames 
    if ( source->size() != 1 ) {
      streamlog_out(MESSAGE3) << "Ignoring bad event " << _nEvt << " having size " << source->size() << std::endl;
      return false;
    }
    
    // Read first data block 
    ShortVec& short_data0 = dynamic_cast<TrackerRawDataImpl* > ( source->getElementAt(0) )->adcValues();
    
    datavect  data0; 
    for (size_t j = 0; j < short_data0.size(); j++) 
    {   
      data0.push_back( static_cast<unsigned char> (short_data0[j]) ); 
    }
    
    // Prepare a new lcio::TrackerData for the ZS data
    lcio::TrackerDataImpl* frame =  new lcio::TrackerDataImpl;
    outputEncoder["sensorID"] = _sensorID;
    outputEncoder["sparsePixelType"] = 0;
    outputEncoder.setCellID( frame );
    
    // Decode data
    auto pixelVec = decodeFEI4Data(data0);   
    
    // Fill ZS data into new lcio::TrackerData object 
    for(auto& hitPixel: pixelVec) {
      frame->chargeValues().push_back( hitPixel.x );
      frame->chargeValues().push_back( hitPixel.y );
      frame->chargeValues().push_back( hitPixel.tot);    
	}
      
    // Now add the TrackerData to the collection
    result->push_back( frame );
       
    return true;
  }
  
  std::vector<PyBARUnpacker::APIXPix> PyBARUnpacker::decodeFEI4Data(std::vector<unsigned char> const & data) {
	std::vector<PyBARUnpacker::APIXPix> result;
      
    // check for consistency
    //bool valid=isEventValid(data);
    unsigned int ToT = 0;
    unsigned int Col = 0;
    unsigned int Row = 0;
    // FE-I4: DH with lv1 before Data Record
    unsigned int lvl1 = 0;
    
    for (unsigned int i=4; i < data.size(); i += 4) {
      unsigned int Word = (((unsigned int)data[i]) << 24) | (((unsigned int)data[i+1]) << 16) | (((unsigned int)data[i+2]) << 8) | (unsigned int)data[i+3];
      if (PYBAR_DATA_HEADER_MACRO(Word)) {
        lvl1++;
      } else {
        // First Hit
        if (getHitData(Word, false, Col, Row, ToT)) {
          result.emplace_back(Col, Row, ToT, lvl1-1);
        }
        // Second Hit
        if (getHitData(Word, true, Col, Row, ToT)) {
          result.emplace_back(Col, Row, ToT, lvl1-1);
        }
      }
    }
	return result;
  }
  
  bool PyBARUnpacker::getHitData(unsigned int &Word, bool second_hit, unsigned int &Col, unsigned int &Row, unsigned int &ToT) 
  {
    static const unsigned int CHIP_MIN_COL = 1;
    static const unsigned int CHIP_MAX_COL = 80;
    static const unsigned int CHIP_MIN_ROW = 1;
    static const unsigned int CHIP_MAX_ROW = 336;
    static const unsigned int CHIP_MAX_ROW_NORM = CHIP_MAX_ROW - CHIP_MIN_ROW;
    static const unsigned int CHIP_MAX_COL_NORM = CHIP_MAX_COL - CHIP_MIN_COL;
    static const int chip_id_offset = 20;
    
    if ( !PYBAR_DATA_RECORD_MACRO(Word) ) return false;	// No Data Record
     
    unsigned int t_Col=0;
    unsigned int t_Row=0;
    unsigned int t_ToT=15;

    if (!second_hit) {
      t_ToT = PYBAR_DATA_RECORD_TOT1_MACRO(Word);
      t_Col = PYBAR_DATA_RECORD_COLUMN1_MACRO(Word);
      t_Row = PYBAR_DATA_RECORD_ROW1_MACRO(Word);
    } else {
      t_ToT = PYBAR_DATA_RECORD_TOT2_MACRO(Word);
      t_Col = PYBAR_DATA_RECORD_COLUMN2_MACRO(Word);
      t_Row = PYBAR_DATA_RECORD_ROW2_MACRO(Word);
    }
    
    // translate FE-I4 ToT code into tot
    if (tot_mode==1) {
      if (t_ToT==15) return false;
      if (t_ToT==14) ToT = 1;
      else ToT = t_ToT + 2;
    } else if (tot_mode==2) {
      // No tot = 2 ?
      if (t_ToT==15) return false;
      if (t_ToT==14) ToT = 1;
      else ToT = t_ToT + 3;
    } else {
      // 0
      if (t_ToT==14 || t_ToT==15) return false;
      ToT = t_ToT + 1;
    }
    
    if (t_Row > CHIP_MAX_ROW || t_Row < CHIP_MIN_ROW) {
      std::cout << "Invalid row: " << t_Row << std::endl;
      return false;
    }
    if (t_Col > CHIP_MAX_COL || t_Col < CHIP_MIN_COL) {
      std::cout << "Invalid col: " << t_Col << std::endl;
      return false;
    }
     
    // Normalize Pixelpositions
    t_Col -= CHIP_MIN_COL;
    t_Row -= CHIP_MIN_ROW;
    Col = t_Col;
    Row = t_Row;
     
    return true;
  }
  
  //
  // Method called after each event to check the data processed
  //
  void PyBARUnpacker::check( LCEvent * evt )
  {
  }
  
  //
  // Method called after all data processing
  //
  void PyBARUnpacker::end()
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
  void PyBARUnpacker::printProcessorParams() const
  { 
    streamlog_out(MESSAGE3)  << std::endl
                             << " "
                             << "PyBARUnpacker Development Version, be carefull!!"
                             << " "
                             << std::endl  << std::endl;   
  
  }
} // Namespace



