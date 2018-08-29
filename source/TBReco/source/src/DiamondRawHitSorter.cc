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
                            "Name of the output data collection",
                            _outputCollectionName, string("zsdata_dia"));

   registerProcessorParameter( "PixelType","Select type of diamond pixel: 0 not on diamond, 1 FE only, 2 rectangular, 3 hexagonal",
                                m_type, static_cast<int > (2));

   registerProcessorParameter( "Modulus","How many floats are used to encode a digit",
                               m_modulus, static_cast<int > (3));

}

//
// Method called at the beginning of data processing
//
void DiamondRawHitSorter::init() {
   
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
     
     // The output collection with filtered and mapped pixels 
     LCCollectionVec * outputCollection = new LCCollectionVec(LCIO::TRACKERDATA);

      // CellID encoding string  
     CellIDEncoder<TrackerDataImpl> outputEncoder( "sensorID:6,sparsePixelType:5" , outputCollection ); 
     
      
           
     // Loop over data blocks
     for (size_t iBlock = 0; iBlock < inputCollection->size(); iBlock++) { 
     
       // Prepare a TrackerData to store the mapped zs data frame 
       TrackerDataImpl* mapped_frame = new TrackerDataImpl ;
  
       // Get the unmapped data 
       TrackerDataImpl * frame = dynamic_cast<TrackerDataImpl* > ( inputCollection->getElementAt(iBlock) );
        
       int sensorID = inputDecoder( frame ) ["sensorID"];
       
       // Loop over digits and filter digits
       FloatVec rawData = frame->getChargeValues();
       int nDigits = rawData.size()/m_modulus; 
         
       for (int iDigit = 0; iDigit < nDigits; iDigit++) 
       {   
        
         int x = static_cast<int> (rawData[iDigit * m_modulus]);
         int y = static_cast<int> (rawData[iDigit * m_modulus + 1]);
         float charge =  rawData[iDigit * m_modulus + 2];     
         
         // Print detailed pixel summary, for testing/debugging only !!! 
         streamlog_out(MESSAGE1) << "Digit Nr. " << iDigit << " on sensor " << sensorID
                                 << std::endl;  
         streamlog_out(MESSAGE1) << "   column:" << x << ", row:" << y
                                 << ", charge:" << charge
                                 << std::endl;

         if ( getPixelType(x, y) == m_type) {
           vector<int> newcoord = getPixelCoordinates(x, y);
           streamlog_out(MESSAGE1) << "FE coord (" << x << ", " << y << "), new coord: (" << newcoord[0] << ", " << newcoord[1] << ")" << std::endl;
           // Store raw digits in tbsw format 
           mapped_frame->chargeValues().push_back( newcoord[0] );
           mapped_frame->chargeValues().push_back( newcoord[1] );
           mapped_frame->chargeValues().push_back( charge );   
         }
         
       }  
       
       outputEncoder["sensorID"] = sensorID;
       outputEncoder["sparsePixelType"] = 0;
       outputEncoder.setCellID( mapped_frame );
       outputCollection->push_back( mapped_frame );
       
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


vector<int> DiamondRawHitSorter::getHexCoord(int x, int y){
  int c = 0, r  = 0;
  std::vector<int> coord;
  	
  // column periodicity is 4 in one row of hex and 5 in the second if you count the small rectangle without dedicated connection in the odd d columns.
  // So do not count these small rectangles to get a 4 periodicity every where. This entails extra handling of some odd dcolumn rows and left columns
  
	// columns
	// at the moment no special handling for last rows which are not connected to the expected pixels
	if (x%2 == 0){ // right column
		if ((y-hexStartY)%4 == 0){
			c = (x-metalStartX)/2 *hexPeriod +1; // +1
		}
		else if ((y-hexStartY)%4 == 1){
			c = (x-metalStartX)/2 *hexPeriod; // +0
		}
		
		else {
			if ((x/2) %2 == 0){ // right column of even double column
				if ((y-hexStartY)%4 == 2){
					c = (x-metalStartX)/2 *hexPeriod +1; // +1
				}
				else if ((y-hexStartY)%4 == 3){
					c = (x-metalStartX)/2 *hexPeriod; // +0
				}
			}
			else { // right column of odd double column
				if ((y-hexStartY)%4 == 2){
					c = (x-metalStartX)/2 *hexPeriod; // +0 because small rectangle does not count
				}
				else if ((y-hexStartY)%4 == 3){
					c = (x-metalStartX)/2 *hexPeriod +1; // +1 because small rectangle
				}
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
		else {
			if ((x+1)/2 %2 == 0){ // left column of even double column
				if ((y-hexStartY)%4 == 2){
					c = (x-metalStartX+1)/2 *hexPeriod -1 -1; // -1 
				}
				else if ((y-hexStartY)%4 == 3){
					c = (x-metalStartX+1)/2 *hexPeriod -1; // -0
				}
			}
			else { // left column of odd double column
				if ((y-hexStartY)%4 == 2){
					c = (x-metalStartX+1)/2 *hexPeriod -1; // -0 because no small rect counting
				}
				else if ((y-hexStartY)%4 == 3){
					c = (x-metalStartX+1)/2 *hexPeriod -2; // -1
				}
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



