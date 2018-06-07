// DiamondRawHitSorter implementation file
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "DiamondRawHitSorter.h"

#include <iomanip>
 
#include "DEPFET.h" 


// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;

namespace depfet {

//
// Instantiate this object
//
DiamondRawHitSorter aDiamondRawHitSorter ;

//
// Constructor
//
DiamondRawHitSorter::DiamondRawHitSorter() : Processor("DiamondRawHitSorter")
{

// Processor description
   _description = "DiamondRawHitSorter: Select a subset of hits on Diamond hits with a given pixel type" ;

//   
// First of all, we need to register the input/output collections
   
   registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
                            "Name of collection containing raw FEI4 hits",
                            _inputCollectionName, string("zsdata_raw_dia"));
   
   
   registerOutputCollection (LCIO::TRACKERDATA, "OutputCollectionName",
                            "Name of the output cluster collection",
                            _outputCollectionName, string("zsdata_dia"));

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
void DiamondRawHitSorter::init() {
   
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
void DiamondRawHitSorter::processRunHeader(LCRunHeader * run)
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
void DiamondRawHitSorter::processEvent(LCEvent * evt)
{
   
   streamlog_out(MESSAGE0) << "start DiamondRawHitSorter " << endl; 

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
void DiamondRawHitSorter::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void DiamondRawHitSorter::end()
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
void DiamondRawHitSorter::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "DiamondRawHitSorter Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}

//
// Function for type of pixel: 0 not on diamond, 1 FE only, 2 rectangular, 3 hexagonal
// 
int DiamondRawHitSorter::getPixelType(int x, int y){
  
  if (x >= layoutStartX && x < metalStartX && y >= layoutStartY && y <= layoutEndY){ // left free
    return 1;
  }
  else if (x >= metalStartX && x <= metalEndX && y >= metalEndY && y <= layoutEndY){ // botFree
    return 1;
  }
   
  else if (x >= metalStartX && x <= metalEndX && y >= metalStartY && y < hexStartY){ // rectangular
    return 2;
  }
	
  else if (x >= metalStartX && x <= metalEndX && y >= hexStartY && y < metalEndY){
    return 3;
  }
  return 0; // not on diamond
}

//
// Get the new grid coordinates without knowing the type
//
vector<int> DiamondRawHitSorter::getPixelCoordinates(int x, int y){
  int type = getPixelType(x, y);
  return getPixelCoordinates(x, y, type);
}

//
// Get the coordinates knowing the type
//
vector<int> DiamondRawHitSorter::getPixelCoordinates(int x, int y, int type){
  vector<int> coord;
  if (type == 0 || type == 1){
    coord.push_back(x);
    coord.push_back(y);
  }
  else if (type == 2){
    coord = getRectCoord(x, y);
  }
  else if (type == 3){
    coord = getHexCoord(x, y);
  }
  return coord;
}


vector<int> DiamondRawHitSorter::getRectCoord(int x, int y){
  vector<int> coord;
  int c = 0, r = 0;
  	
  // columns
  if (x %2 == 0){ // right column
    if ((x/2) %2 == 0){ // right column of even double column
      if (y %2 == 0){
        if (x == layoutStartY){
          c = (x-metalStartX)/2 *rectPeriod; // +0 in first row of pixels
        }
        else {
          c = (x-metalStartX)/2 *rectPeriod +2; // (x-...)=number of double column knowing it is the right column of even dcol, +2 for correct pixel in this row
        }
      }
      else {
        c = (x-metalStartX)/2 *rectPeriod +1; // +1 for correct pixel in this row
      }
    }
    else { // right column of odd double column
      if (y %2 == 0){
        if (x == layoutStartY){
          c = (x-metalStartX)/2 *rectPeriod; // +0 in first row of pixels
        }
        else {
          c = (x-metalStartX)/2 *rectPeriod +1; // last +1 for correct pixel
        }
      }
      else {
        c = (x-metalStartX)/2 *rectPeriod; // +0 for correct pixel
      }
    }
  }
  else { // left column
    if ((x+1)/2 % 2 == 0){ //left column of even double column
      if (y %2 == 0){
        if (x == layoutStartY){
          c = (x-metalStartX+1)/2 *rectPeriod -1; // +0 in first row of pixels
        }
        else {
          c = (x-metalStartX+1)/2 *rectPeriod -1 -2; // -1 because the grid period is oriented on the right column, -2 for correct pixel
        }
      }
      else {
        c = (x-metalStartX+1)/2 *rectPeriod -1 -1; // -1 for correct pixel position
      }
    }
    else { //left column of odd double column
      if (y %2 == 0){
        if (x == layoutStartY){
          c = (x-metalStartX+1)/2 *rectPeriod -1; // +0 in first row of pixels
        }
        else {
          c = (x-metalStartX+1)/2 *rectPeriod -1 -1;
        }
      }
      else {
        c = (x-metalStartX+1)/2 *rectPeriod -1;
      }
    }
  }
  coord.push_back(c);
	
  // rows
  if (y%2 == 0) {
    r = y/2;
  }
  else {
    r = (y+1)/2;
  }
  coord.push_back(r);
  
  return coord;
}


svector<int> DiamondRawHitSorter::getHexCoord(int x, int y){
  int c = 0, r  = 0;
  vector<int> coord;
  
  // columns
  // at the moment no special handling for last rows which are not connected to the expected pixels
  if (x%2 == 0){ // right column
    if ((y-hexStartY)%4 == 0){
      c = (x-metalStartX)/2 *hexPeriod +1; // +1
    }
    else if ((y-hexStartY)%4 == 1){
      c = (x-metalStartX)/2 *hexPeriod; // +0
    }
    else if ((y-hexStartY)%4 == 2){
      c = (x-metalStartX)/2 *hexPeriod +1; // +1
    }
    else {
      if ((x/2) %2 == 0){ // right column of even double column
        c = (x-metalStartX)/2 *hexPeriod; // +0
      }
      else { // right column of odd double column
        c = (x-metalStartX)/2 *hexPeriod +2; // +2
      }
    }
  }
  else { // left column
    if ((y-hexStartY)%4 == 0){
      c = (x-metalStartX+1)/2 *hexPeriod -1 -1; // -1
    }
    else if ((y-hexStartY)%4 == 1){
      c = (x-metalStartX+1)/2 *hexPeriod -1; // -0
    }
    else if ((y-hexStartY)%4 == 2){
      c = (x-metalStartX+1)/2 *hexPeriod -1 -1; // +1
    }
    else {
      if ((x+1)/2 %2 == 0){ // left column of even double column
        c = (x-metalStartX+1)/2 *hexPeriod -1; // -0
      }
      else { // left column of odd double column
        c = (x-metalStartX+1)/2 *hexPeriod -2; // -2
      }
    }
  }
  coord.push_back(c);
  //rows
  if ((y-hexStartY) %2 == 0){
    r = (y-hexStartY)/2;
  }
  else {
    r = (y-hexStartY-1)/2;
  }
  coord.push_back(r);
  return coord;
}


} // Namespace



