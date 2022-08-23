// HitsFilterProcessor implementation file
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "HitsFilterProcessor.h"
#include <iomanip>

// Include TBTools 
#include "TBDetector.h"


// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;
using namespace std::string_literals;
namespace depfet {

//
// Instantiate this object
//
HitsFilterProcessor aHitsFilterProcessor ;

//
// Constructor
//
HitsFilterProcessor::HitsFilterProcessor() : Processor("HitsFilterProcessor") , _outputEncoderHelper( "sensorID:6,sparsePixelType:5")
{

   // Processor description
   _description = "HitsFilterProcessor: Reformat telescope cluster data to raw data" ;

   //   
   // First of all, we need to register the input/output collections
   registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
                            "Name of collection containing original zero suppressed data",
                            _inputCollectionName, string("zsdata"));
   
   
   registerOutputCollection (LCIO::TRACKERDATA, "OutputCollectionName",
                            "Name of the output filtered collection",
                            _outputCollectionName, string("zsdata_filtered"));

   std::vector<int> initFilterIDs;
   registerProcessorParameter ("FilterIDs",
                              " Only use hits from these these sensorIDs",
                              _filterIDs, initFilterIDs);
   
}

//
// Method called at the beginning of data processing
//
void HitsFilterProcessor::init() {
   
   // Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
   _hitElements = 5;
   
   // Print set parameters
   printProcessorParams();
   
   // CPU time start
   _timeCPU = clock()/1000;
}

//
// Method called for each run
//
void HitsFilterProcessor::processRunHeader(LCRunHeader * run)
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
void HitsFilterProcessor::processEvent(LCEvent * evt)
{
   
   streamlog_out(MESSAGE0) << "Start HitsFilterProcessor " << endl; 

   // Print event number
   if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                   << (evt->getEventNumber())
                                                                   << std::endl << std::endl;
   
   // More detailed event numbering for testing
   streamlog_out(MESSAGE2) << std::endl << "Starting with Event Number " << evt->getEventNumber()  << std::endl;  
   
   //
   // Open collections
   try {
     
     // Open zero suppressed pixel data  
     LCCollectionVec * inputCollection = dynamic_cast < LCCollectionVec * > (evt->getCollection(_inputCollectionName)); 
      
     // The original data collection contains one or more frames with zs data
     LCCollectionVec * outputCollection = new LCCollectionVec(LCIO::TRACKERDATA);
     
     // Prepare a TrackerData to store reformatted raw data 
     std::map<int, TrackerDataImpl*> outputDigitsMap;
     for (auto id :  _filterIDs)  outputDigitsMap[id] = new TrackerDataImpl ; 
           
     // Frame loop
     for (size_t iFrame = 0; iFrame < inputCollection->size(); iFrame++) { 
       
       TrackerDataImpl * zsFrame = dynamic_cast<TrackerDataImpl* > ( inputCollection->getElementAt(iFrame) );
        
       
             
       // Loop over digits
       FloatVec rawData = zsFrame->getChargeValues();
       int nDigits = rawData.size()/_hitElements; 
         
       for (int iDigit = 0; iDigit < nDigits; iDigit++) 
       {   
        
        
        int sensorID = static_cast<int> (rawData[iDigit * _hitElements]);
        int col = static_cast<int> (rawData[iDigit * _hitElements + 1]);
        int row = static_cast<int> (rawData[iDigit * _hitElements + 2]);
        float charge =  rawData[iDigit * _hitElements + 3];    
        float time = rawData[iDigit * _hitElements + 4];
        
         
        // Print detailed pixel summary, for testing/debugging only !!! 
        streamlog_out(MESSAGE1) << "Digit Nr. " << iDigit << " on sensor " << sensorID  
                                << std::endl;  
        streamlog_out(MESSAGE1) << "   column:" << col << ", row:" << row
                                << ", charge:" << charge
                                << ", time:" << time
                                << std::endl;
        
        // Ignore digits from this sensor. 
        if ( outputDigitsMap.find(sensorID) == outputDigitsMap.end() ) continue;
         
        // Store raw digits in tbsw format 
        outputDigitsMap[sensorID]->chargeValues().push_back( col );
        outputDigitsMap[sensorID]->chargeValues().push_back( row );
        outputDigitsMap[sensorID]->chargeValues().push_back( charge );
        outputDigitsMap[sensorID]->chargeValues().push_back( time );   
       }  

     } // End frame loop 
       
     // Set the proper cell encoder
     CellIDEncoder<TrackerDataImpl> outputEncoder( "sensorID:6,sparsePixelType:5", outputCollection , &_outputEncoderHelper );
 
     // Add output digits to output collection
     for ( auto& cached :  outputDigitsMap ) 
     { 
       outputEncoder["sensorID"s] = cached.first;
       outputEncoder["sparsePixelType"s] = 0;
       outputEncoder.setCellID( cached.second );
       outputCollection->push_back( cached.second );
     }
     
     // Add output collection to event
     evt->addCollection( outputCollection, _outputCollectionName );
             	    
   } catch(DataNotAvailableException &e){}  
   _nEvt ++ ;
}


//
// Method called after each event to check the data processed
//
void HitsFilterProcessor::check( LCEvent * )
{
}

//
// Method called after all data processing
//
void HitsFilterProcessor::end()
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
void HitsFilterProcessor::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "HitsFilterProcessor Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}


} // Namespace



