// DEPFETHybrid6Unpacker implementation file
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "DEPFETHybrid6Unpacker.h"
#include "Utils.hh"
#include "DEPFETADCValues.hh"
#include "DEPFETEvent.hh"

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

namespace {


struct EvtHeader {
    unsigned int    EventSize: 20;
    unsigned short   flag0: 1;
    unsigned short   flag1: 1;
    unsigned short  EventType: 2;
    unsigned short  ModuleNo: 4;
    unsigned short  DeviceType: 4;
    unsigned int    Triggernumber;
};

struct InfoWord  {
    unsigned int framecnt: 10; // number of Bits
    unsigned int startgate: 10; //jf new for 128x128
    unsigned int zerosupp: 1;
    unsigned int startgate_ver: 1;
    unsigned int temperature: 10;
};


struct DHH_Header_t{
    unsigned  flag:1;
    unsigned  DataType:3;
    unsigned  FrameNr:4;
    unsigned  DHH_ID:6;
    unsigned  Chip_ID:2;
    unsigned short TriggerNr;
};

struct Onsen_Header_t {
    unsigned int magic;
    unsigned int size;
    unsigned int dummy1;
    unsigned int dummy2;
    unsigned int DHH_Header;
};
  
}

namespace eudaqinput {

  //
  // Instantiate this object
  //
  DEPFETHybrid6Unpacker aDEPFETHybrid6Unpacker ;

  //
  // Constructor
  //
  DEPFETHybrid6Unpacker::DEPFETHybrid6Unpacker() : Processor("DEPFETHybrid6Unpacker")
  {
    
    // Processor description
    _description = "DEPFETHybrid6Unpacker: Unpacker for DEPFET Hybrid6 detectors (large Belle2 PXD6 sensor) raw data" ;
    
    //   
    // First of all, we need to register the input/output collections
    
    registerInputCollection (LCIO::TRACKERRAWDATA, "InputCollectionName",
                             "Name of collection containing raw data",
                             _inputCollectionName, string("DEPFET"));
    
    registerOutputCollection (LCIO::TRACKERDATA, "OutputCollectionName",
                              "Name of the output digit collection",
                              _outputCollectionName, string("zsdata_dep"));
  
  }
  
  //
  // Method called at the beginning of data processing
  //
  void DEPFETHybrid6Unpacker::init() {
     
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
  void DEPFETHybrid6Unpacker::processRunHeader(LCRunHeader * run) {
    
    // Print run number
    streamlog_out(MESSAGE3) << "Processing run: "
                            << (run->getRunNumber())
                            << std::endl << std::endl;
    
    _nRun++ ;  
  }
  
  //
  // Method called for each event
  //
  void DEPFETHybrid6Unpacker::processEvent(LCEvent * evt)
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
      depfet::DEPFETEvent m_event; 
      
      for (size_t iPlane = 0; iPlane < numplanes; ++iPlane) { 
      
        TrackerRawDataImpl* rawDataBlock = dynamic_cast<TrackerRawDataImpl* > ( source->getElementAt(iPlane) );
        
        // Get blockID 
        int blockID = inputDecoder(rawDataBlock)["blockID"]; 
        
        // FIXME: This is only needed, because i do not know how to put 
        // vector<unsigned char> directly into the lcio::LCEvent. 
        datavect  char_data;
        ShortVec& short_data = rawDataBlock->adcValues();
        for (size_t j = 0; j < short_data.size(); j++) 
        {   
          char_data.push_back( static_cast<unsigned char> (short_data[j]) ); 
        }
        
        // Get read fully encoded subevent
        depfet::DEPFETADCValues m_subevent = ConvertDEPFETEvent(char_data, blockID);    
        m_event.push_back(m_subevent);
      }
      
      LCCollectionVec* depfet_info = new LCCollectionVec( LCIO::LCGENERICOBJECT )  ;
      
      //-----------------------------------------------
      // Decode DEPFET event info to LCIO format      
      for (size_t iframe=0;iframe<m_event.size();iframe++) {
        
        // Get read fully decoded data  
        const depfet::DEPFETADCValues& data = m_event[iframe];
        
        LCGenericObjectImpl* metaobj = new LCGenericObjectImpl(2,0,0);
        metaobj->setIntVal(0,data.getGoodEvent());
        metaobj->setIntVal(1,data.getStartGate());
        
        depfet_info->addElement( metaobj) ;   
      }
      evt->addCollection( depfet_info , "DEPFET_EVENT_INFO" ) ;
      
      // Prepare collection for unpacked digits 
      LCCollectionVec * result = new LCCollectionVec(LCIO::TRACKERDATA);
      
      // Helper class for encoding digit data 
      CellIDEncoder<TrackerDataImpl> outputEncoder( "sensorID:6,sparsePixelType:5", result  );
      
      //-----------------------------------------------
      // Decode event data to a LCIO format      
      for (size_t iframe=0;iframe<m_event.size();iframe++) {
      
        // Get read fully decoded data  
        const depfet::DEPFETADCValues& data = m_event[iframe];
        
        // Prepare a TrackerData to store zs pixels
        TrackerDataImpl* zsFrame = new TrackerDataImpl;
      
        // Set description for zsFrame
        outputEncoder["sensorID"] = data.getModuleNr();
        outputEncoder["sparsePixelType"] = 0;
        outputEncoder.setCellID( zsFrame );
        
        int nPixel = (int) data.at(0).size(); 
           
        for ( int iPixel = 0; iPixel < nPixel; iPixel++ ) {
          int val = data.at(2).at(iPixel);
          int col = data.at(1).at(iPixel);   
          int row = data.at(0).at(iPixel);  
          //int cm = data.at(3).at(iPixel);
          
          if (row%8<4) {  
            // remap data 
            row = 4*(row/8) + row%4;   
            zsFrame->chargeValues().push_back( col );
            zsFrame->chargeValues().push_back( row );
            zsFrame->chargeValues().push_back( val );   
          }
        
        }
      
        // Add event to LCIO collection 
        result->push_back(zsFrame);
      }  
      
      // Add output collection to event
      evt->addCollection( result, _outputCollectionName );
      
      // Prepare collection for unpacked digits  (2nd sample) 
      LCCollectionVec * result2 = new LCCollectionVec(LCIO::TRACKERDATA);
      
      // Helper class for encoding digit data 
      CellIDEncoder<TrackerDataImpl> outputEncoder2( "sensorID:6,sparsePixelType:5", result2  );
      
      //-----------------------------------------------
      // Decode event data to a LCIO format      
      for (size_t iframe=0;iframe<m_event.size();iframe++) {
      
        // Get read fully decoded data  
        const depfet::DEPFETADCValues& data = m_event[iframe];
        
        // Prepare a TrackerData to store zs pixels
        TrackerDataImpl* zsFrame = new TrackerDataImpl;
      
        // Set description for zsFrame
        outputEncoder2["sensorID"] = data.getModuleNr();
        outputEncoder2["sparsePixelType"] = 0;
        outputEncoder2.setCellID( zsFrame );
        
        int nPixel = (int) data.at(0).size(); 
           
        for ( int iPixel = 0; iPixel < nPixel; iPixel++ ) {
          int val = data.at(2).at(iPixel);
          int col = data.at(1).at(iPixel);   
          int row = data.at(0).at(iPixel);  
          //int cm = data.at(3).at(iPixel);
           
          if (row%8>3) {  
            // remap data 
            row = 4*(row/8) + row%4;   
            zsFrame->chargeValues().push_back( col );
            zsFrame->chargeValues().push_back( row );
            zsFrame->chargeValues().push_back( val );   
          }
        }
        
        // Add event to LCIO collection 
        result2->push_back(zsFrame);
      }  
      
      // Add output collection to event
      evt->addCollection( result2, _outputCollectionName + to_string(2) );
      
    } catch(DataNotAvailableException &e){}  
    _nEvt ++ ;
  }
  
  
  depfet::DEPFETADCValues DEPFETHybrid6Unpacker::ConvertDEPFETEvent(const std::vector<unsigned char> & data, unsigned id)
  {
    //---------------------------------------------------------------------------------------------      
    depfet::DEPFETADCValues m_event; 
    
    bool debug = false; // false; 
    long int rc = 0; 
    
    streamlog_out(MESSAGE2) << "ConvertDEPFETEvent() of block id " << id << " of length " << data.size() << std::endl;  
    
    //---------------------------------------------------------------------------------------------
    struct EvtHeader GROUP_Header = *(struct EvtHeader *) &data[rc];
    rc += sizeof(struct EvtHeader); 
     
    if(debug) printf(" read Header: DevTyp=%2d   EvtType=%d  Trig=%#10X  Evt_size=%6d  ModID=%2d \n"
                        ,GROUP_Header.DeviceType
                        ,GROUP_Header.EventType
                        ,GROUP_Header.Triggernumber
                        ,GROUP_Header.EventSize
                        ,GROUP_Header.ModuleNo
                        );
    
    int ModID=GROUP_Header.ModuleNo;
    m_event.setModuleNr(ModID);
    m_event.setTriggerNr(GROUP_Header.Triggernumber);
    
    //---------------------------------------------------------------------------------------------
    if (GROUP_Header.EventType==0 && GROUP_Header.DeviceType==14) {  //-- INFO event
      if(debug)printf("     Info Event\n     Run Number=%d \n",GROUP_Header.Triggernumber);  
    }
    
    //---------------------------------------------------------------------------------------------
    if (GROUP_Header.EventType==2 && GROUP_Header.DeviceType==13) {  //-- TLU event 
      if(debug)printf(" TLU Event found!! This will probably crash!!"); 
    }
    
    //---------------------------------------------------------------------------------------------
    //              check DATA event
    //---------------------------------------------------------------------------------------------
    if (GROUP_Header.EventType==2
      && (GROUP_Header.DeviceType==2  || GROUP_Header.DeviceType==3 || GROUP_Header.DeviceType==4 || GROUP_Header.DeviceType==11)
      &&  ModID>0) {
             if(debug)printf(" Data event from DevType=%d: this would be Hybrid 4 or 4.1 Geant, etc...\n", GROUP_Header.DeviceType); 
    }
    
    
        
    if (GROUP_Header.EventType==2
        && (GROUP_Header.DeviceType==7)   // This is the DHH/Onsen system; Should be DeviceType==7
        &&  ModID>0) {
         
      // The InfoWord is a 32 bit frame header
      // behind the Bonn header
      //struct InfoWord frameInfoData = *(struct InfoWord *) &data[rc];
      rc += sizeof(struct InfoWord);
         
      
      uint16_t tempZSFrameData[131080];
         
      // Npixel counts how many byte in payload
      // EventSize is size of payload + header in 32bit integers
      int Npixels=(GROUP_Header.EventSize-3)*4;  
      
      // Start looping over DHH/Onsen frames in the event
         
      unsigned int nextHeaderIndex=rc;
      unsigned int lastByte=Npixels+rc; 
      bool finished = false; 
	  bool first_pixel=true;
            
      while( (!finished) && (nextHeaderIndex<lastByte) ){ 
          
	    //char dhhType = (data[nextHeaderIndex + 17] & 0x70) >> 4;
        
        //printf("dhhType = %d\n", dhhType);
            
        // Extract ONSEN and DHH Header (==Onsen_Header_t )
        // DHH header is last integer in Onsen_Header_t 
        // firstHeader.size is size of payload (without Onsen header but with DHH header) in byte
        struct Onsen_Header_t firstHeader = *(struct Onsen_Header_t * ) &data[nextHeaderIndex];
        rc += sizeof(struct Onsen_Header_t);
           
        if(firstHeader.magic!=0xCAFEBABE){
          printf("\nERROR, no magic detected! got 0x%08x, expected 0x%08x \n\n",firstHeader.magic,0xCAFEBABE);
          return m_event;
        }
               
        unsigned int currentPayloadBeginIndex=nextHeaderIndex+sizeof(Onsen_Header_t); // in byte
        unsigned int payloadLength=firstHeader.size-sizeof(unsigned int);   // in byte
            
        nextHeaderIndex+= firstHeader.size + sizeof(Onsen_Header_t) - sizeof(unsigned int);
        if((firstHeader.size + sizeof(Onsen_Header_t) - sizeof(unsigned int))%4!=0){
          nextHeaderIndex+=2;
        }
             
        //Correct DHH Header Endianess:
        unsigned int DHHHeaderNotCorrected= firstHeader.DHH_Header;
           
        DHH_Header_t  DHH_Header;
        DHH_Header.flag=(DHHHeaderNotCorrected &0x8000)>>15;
        DHH_Header.DataType=(DHHHeaderNotCorrected &0x7800)>>11;
        DHH_Header.FrameNr=(DHHHeaderNotCorrected &0x0f00)>>8;
        DHH_Header.DHH_ID=(DHHHeaderNotCorrected &0x00fc)>>2;
        DHH_Header.Chip_ID=(DHHHeaderNotCorrected &0x0003);
        DHH_Header.TriggerNr=(DHHHeaderNotCorrected)>>16;
           
        if(debug) printf("ErrorFlag %d, DataType %d, FrameNr %d, DHH_ID %d, Chip_ID %d, TriggerNr 0x%04x\n",
                           DHH_Header.flag,
                           DHH_Header.DataType,
                           DHH_Header.FrameNr,
                           DHH_Header.DHH_ID,
                           DHH_Header.Chip_ID,
                           DHH_Header.TriggerNr);
           
        if(false) {
          for (unsigned int index =0;index<payloadLength;index++){
            printf(" 0x%02x ", data[currentPayloadBeginIndex+index]);
          }
          printf("\n");
        }
            
        if(DHH_Header.DataType==3){
          if(debug)printf("StartOfFrame arrived\n");
                    
          // This is data
          rc = currentPayloadBeginIndex; 
          for (unsigned int index =0;index<payloadLength/2;index++){
            tempZSFrameData[index] = *(uint16_t *) &data[rc];
            rc += sizeof(uint16_t);  
            if(debug)printf(" 0x%02x ", tempZSFrameData[index]);
          }
                     
          int delay_frsyc2trg = ((tempZSFrameData[3] & 0x3ff));
          int lastdhpframe = (tempZSFrameData[3] & 0xfc00) >> 10;
                    
          if(debug) std::cout << " delay " << delay_frsyc2trg << std::endl;
          if(debug) std::cout << " last dhp frame id " << lastdhpframe << std::endl; 
                     
          continue;
        }
        if(DHH_Header.DataType==4){
          if(debug)printf("EndOfFrame arrived\n");
          finished=true;
          continue;
        }
        
        if(DHH_Header.DataType==2){
          if(debug)printf("GhostFrame arrived\n");
          // if this arrived, the DHE responded!
          m_event.setGoodEvent(true);
          continue;
        }
           
        if(DHH_Header.DataType==5){  
             
          // This is ZS data from DHH/Onsen
          // if this arrived, the DHE responded!
          m_event.setGoodEvent(true);
          rc = currentPayloadBeginIndex; 
          for (unsigned int index =0;index<payloadLength/2;index++){
            tempZSFrameData[index] = *(uint16_t *) &data[rc];
            rc += sizeof(uint16_t);  
          }
                     
          //if(debug)printf("ZSData arrived\n");
          uint16_t current_word,current_row_base,current_row_bit,current_val,current_col,current_CM,current_row;
          //uint16_t frame_type = ((tempZSFrameData[0]) >> 13);
          unsigned char mappedCol;
          uint16_t mappedRow;
                                      
          // Here: word means 16bit integer
          // Here: any header can be skipped. 
          // 
          // In case last word is 0000, it means last 
          // three words are footer (header+padding)
          // and should be skipped.
          // Otherwise, last two words are header. 
          // Just skip it. 
          // First two words are always header (skip
          // it). If a000 appears, new DHP reading 
          // frame starts. Skip this word and the next
          // word.
                        
          // Check for padding 
          int lastword; 
          if (tempZSFrameData[payloadLength/2-1] == 0x0000 )
            lastword = (payloadLength-6)/2; // skip last three words  
          else  
            lastword = (payloadLength-4)/2; // skip last two words 
                        
          //if(debug) std::cout << "DHP frame id " << tempZSFrameData[1] << std::endl; 
                  
          //lastWord = tempZSFrameData[payloadLength/2-3] ;
	       
          //printf("Drittes wort vom ende %04x\n",tempZSFrameData[payloadLength/2-3]); 
             
          unsigned short mydhpid = tempZSFrameData[0];
                          
          if ( (mydhpid & 0xFFFC) != 0xA000  ) {
            if(debug) std::cout << "EEEERRRRROOOOOOOOOOR" << std::endl;  
          }         
 
          mydhpid = mydhpid & 0x0003; 

          for (int ipix = 2; ipix < lastword; ipix++) {   
            //if(debug) printf("ZSData arrived %04x\n",tempZSFrameData[ipix]);
            current_word = tempZSFrameData[ipix];
                             
            // Detect start of new DHP reading frame 
            if ( current_word == 0xa000 ) { 
              // Skip this word and the next
              ipix++;
              continue;  
            }
                            
            if (((current_word) >> 15) == 0) { //row header
              current_row_base = 2 * ((current_word) >> 6);
              current_CM = (current_word) & 0x3f;
              //if(debug) printf("New row Header 0x%4x:  row base %2d,with common mode of:%d\n",current_word,current_row_base, current_CM);
            } else { //row content
              current_val = (current_word) & 0xFF;
              current_col = ((current_word) >> 8) & 0x3F;
              current_row_bit = (((current_word) >> 8) & 0x40) / 64;
              current_row = current_row_base + current_row_bit;
              //if(debug) printf("New hit token 0x%4x: col: %3d val: %2d]\n",current_word, current_col, current_val);
                 
              //mapping
              mappedCol= (current_col^0x3C) + 64*mydhpid; 
              mappedRow=current_row;   

              if(first_pixel){
                  first_pixel=false;
                  m_event.setStartGate(mappedRow/4);
              }        
              if (   (mappedRow < 960) && (mappedCol < 192)){     
                                    m_event.at(0).push_back(mappedRow);
                                    m_event.at(1).push_back(mappedCol);
                                    m_event.at(2).push_back(current_val);
                                    m_event.at(3).push_back(current_CM);
              } else printf("Warning: col/val: [%d|%d], in row  %d, out of range. Max Row=%d, max col=%d\n",mappedCol, current_val, mappedRow,480,192);
            }
                
          } //-- end dhp reading frame (hits) 
                               
        } //-- end zs data frame 
           
      } //-- end loop over DHH Frames             
           
    } //-- end DHH/Onsen DATA event 
            
    return m_event;
  }
  
  //
  // Method called after each event to check the data processed
  //
  void DEPFETHybrid6Unpacker::check( LCEvent * evt )
  {
  }
  
  //
  // Method called after all data processing
  //
  void DEPFETHybrid6Unpacker::end()
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
  void DEPFETHybrid6Unpacker::printProcessorParams() const
  { 
    streamlog_out(MESSAGE3)  << std::endl
                             << " "
                             << "DEPFETHybrid6Unpacker Development Version, be carefull!!"
                             << " "
                             << std::endl  << std::endl;   
  
  }
} // Namespace



