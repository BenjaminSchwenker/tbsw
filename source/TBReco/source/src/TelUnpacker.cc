// TelUnpacker implementation file
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "TelUnpacker.h"

#include <iomanip>

// Include DEPFETTrackTools 
#include "DEPFET.h" 


// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;

namespace depfet {

//
// Instantiate this object
//
TelUnpacker aTelUnpacker ;

//
// Constructor
//
TelUnpacker::TelUnpacker() : Processor("TelUnpacker")
{

// Processor description
   _description = "TelUnpacker: Reformat telescope cluster data to raw data" ;

//   
// First of all, we need to register the input/output collections
   
   registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
                            "Name of collection containing telescope data",
                            _inputCollectionName, string("original_zsdata"));
   
   
   registerOutputCollection (LCIO::TRACKERDATA, "ClusterCollectionName",
                            "Name of the output cluster collection",
                            _outputCollectionName, string("zsdata_m26"));

   registerProcessorParameter( "Modulus","How many floats are used to encode a digit",
                               m_modulus, static_cast<int > (4));

   std::vector<int> initFilterIDs;
   registerProcessorParameter ("FilterIDs",
                              "Unpack digits from detectors having DAQ IDs in this list",
                              _filterIDs, initFilterIDs);

}

//
// Method called at the beginning of data processing
//
void TelUnpacker::init() {
   
   // Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
   
   // Read detector constants from gear file
   _detector.ReadGearConfiguration();    
               
   // Print set parameters
   printProcessorParams();
   
   // CPU time start
   _timeCPU = clock()/1000;
}

//
// Method called for each run
//
void TelUnpacker::processRunHeader(LCRunHeader * run)
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
void TelUnpacker::processEvent(LCEvent * evt)
{
   
   streamlog_out(MESSAGE0) << "start telunpacker " << endl; 

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
     // Helper class for decoding pixel data 
     CellIDDecoder<TrackerDataImpl> inputDecoder( inputCollection );  
     
     // The original data collection contains the sparse pixel data for 
     // each accpeted cluster.
     LCCollectionVec * outputCollection = new LCCollectionVec(LCIO::TRACKERDATA);
     
     // Prepare a TrackerData to store reformatted raw data 
     std::map<int, TrackerDataImpl*> outputDigitsMap;
     for (auto id :  _filterIDs)  outputDigitsMap[id] = new TrackerDataImpl ; 
           
     // Cluster loop
     for (size_t iClu = 0; iClu < inputCollection->size(); iClu++) { 
     
       TrackerDataImpl * cluster = dynamic_cast<TrackerDataImpl* > ( inputCollection->getElementAt(iClu) );
        
       // DAQ ID for pixel detector
       int sensorID = inputDecoder( cluster ) ["sensorID"];

       // Ignore digits from this sensor. 
       if ( outputDigitsMap.find(sensorID) == outputDigitsMap.end() ) continue;
             
       // Loop over digits
       FloatVec rawData = cluster->getChargeValues();
       int nDigits = rawData.size()/m_modulus; 
         
       for (int iDigit = 0; iDigit < nDigits; iDigit++) 
       {   
        
        int col = static_cast<int> (rawData[iDigit * m_modulus]);
        int row = static_cast<int> (rawData[iDigit * m_modulus + 1]);
        float charge =  rawData[iDigit * m_modulus + 2];     
         
        // Print detailed pixel summary, for testing/debugging only !!! 
        streamlog_out(MESSAGE1) << "Digit Nr. " << iDigit << " on sensor " << sensorID  
                                << std::endl;  
        streamlog_out(MESSAGE1) << "   column:" << col << ", row:" << row
                                << ", charge:" << charge
                                << std::endl;
        
        // Store raw digits in tbsw format 
        outputDigitsMap[sensorID]->chargeValues().push_back( col );
        outputDigitsMap[sensorID]->chargeValues().push_back( row );
        outputDigitsMap[sensorID]->chargeValues().push_back( charge );   
         
       }  

     } // End cluster loop 
     
     // CellID encoding string  
     CellIDEncoder<TrackerDataImpl> outputEncoder( "sensorID:6,sparsePixelType:5" , outputCollection ); 
     
     // Add output digits to output collection
     for ( auto& cached :  outputDigitsMap ) 
     { 
       outputEncoder["sensorID"] = cached.first;
       outputEncoder["sparsePixelType"] = 0;
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
void TelUnpacker::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void TelUnpacker::end()
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
void TelUnpacker::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "TelUnpacker Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}


} // Namespace



