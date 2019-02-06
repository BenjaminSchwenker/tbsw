// RawDataSparsifier
// 
// See RawDataSparsifier.h for full documentation of processor. 
// 
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "RawDataSparsifier.h"

// Include DEPFETTrackTools 
#include "DEPFET.h" 
#include "MatrixDecoder.h"

// Include basic C
#include <cstdlib>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>


// Include LCIO classes
#include <lcio.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

 
// Used namespaces
using namespace std; 
using namespace lcio ;
using namespace marlin ;


namespace depfet {


//
// Instantiate this object
//
RawDataSparsifier aRawDataSparsifier ;

//
// Constructor
//
RawDataSparsifier::RawDataSparsifier() : Processor("RawDataSparsifier")
{

// Processor description
   _description = "RawDataSparsifier: Finding firing pixels in full analog pixel data" ;

//   
// First of all we need to register the input collection
   registerInputCollection (LCIO::TRACKERDATA, "PixelDataCollection",
                           "Input pixel data",
                           _DataCollectionName, string ("data"));

   registerInputCollection (LCIO::TRACKERDATA, "NoiseCollection",
                           "Input pixel noise data",
                           _NoiseCollectionName, string("noise"));

   registerInputCollection (LCIO::TRACKERRAWDATA, "StatusCollection",
                           "Input pixel status data",
                           _StatusCollectionName, string("status"));

   registerOutputCollection (LCIO::TRACKERDATA, "SparseDataCollection",
                            "Output 0-suppressed pixels",
                            _SparseDataCollectionName, string("zsdata"));
    
//
// Processor parameters
     
   vector<float > ThresholdVecExample;
   ThresholdVecExample.push_back(3.0);
   
   registerProcessorParameter("SignalThreshold","Signal threshold for firing pixels",
                             _ThresholdVec, ThresholdVecExample);
   
   registerProcessorParameter( "UseSNRCut","Iff true, cut on signal to noise ratio",
                               _useSNRCut, static_cast<bool > (false) );
   
}

//
// Method called at the beginning of data processing
//
void RawDataSparsifier::init() {
   
// Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
   _isFirstEvent = true;
                   
// Print set parameters
   printProcessorParams();
   
// CPU time start
   _timeCPU = clock()/1000;
      
}

//
// Method called for each run
//
void RawDataSparsifier::processRunHeader(LCRunHeader * run)
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
void RawDataSparsifier::processEvent(LCEvent * evt)
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
    
     
     LCCollectionVec * DataCollection = dynamic_cast < LCCollectionVec * > (evt->getCollection(_DataCollectionName));
     LCCollectionVec * NoiseCollection = dynamic_cast < LCCollectionVec * > (evt->getCollection(_NoiseCollectionName));
     LCCollectionVec * StatusCollection = dynamic_cast < LCCollectionVec * > (evt->getCollection(_StatusCollectionName));
     
     
     if (isFirstEvent()) {
        
       // this is the right place to cross check data to steering file
       int noOfDetector = DataCollection->getNumberOfElements();
  
       // let's check if the number of sigma cut components is the same of
       // the detector number.
       if ( (DataCollection->getNumberOfElements() != NoiseCollection->getNumberOfElements()) ) {
         streamlog_out(ERROR3) <<  "Input data and noise are incompatible" << endl
            << "full data collection has    " << DataCollection->getNumberOfElements()    << " detectors," << endl
            << "noise collection has " << NoiseCollection->getNumberOfElements() << " detectors." << endl;
         exit(-1); 
       }
       
       if ( (DataCollection->getNumberOfElements() != StatusCollection->getNumberOfElements()) ) {
         streamlog_out(ERROR3) << "Input data and status are incompatible" << endl
           << "full data collection has    " << DataCollection->getNumberOfElements()    << " detectors," << endl
           << "status collection has " << StatusCollection->getNumberOfElements() << " detectors." << endl;
         exit(-1); 
       } 
        
       if ( (unsigned) noOfDetector != _ThresholdVec.size() ) {
         streamlog_out( WARNING2 ) << "The number of values in the threshold vector does not match the number of detectors\n"
                                   << "Padding values!!" << endl;
         _ThresholdVec.resize(noOfDetector, _ThresholdVec.back());
       }
       
       for ( int iDetector = 0; iDetector < noOfDetector; iDetector++) {

         CellIDDecoder<TrackerDataImpl> idMatrixDecoder(DataCollection);
         CellIDDecoder<TrackerDataImpl> idNoiseDecoder(NoiseCollection);
         CellIDDecoder<TrackerRawDataImpl> idStatusDecoder(StatusCollection);   
              
          
         TrackerDataImpl    * matrix = dynamic_cast < TrackerDataImpl * >   (DataCollection->getElementAt(iDetector));
         TrackerDataImpl    * noise = dynamic_cast < TrackerDataImpl * >   (NoiseCollection->getElementAt(iDetector));
         TrackerRawDataImpl * status = dynamic_cast < TrackerRawDataImpl * >(StatusCollection->getElementAt(iDetector));
           
         int dataID = static_cast<int> ( idMatrixDecoder(matrix)["sensorID"] );
         int noiseID = static_cast<int> ( idNoiseDecoder(noise)["sensorID"] );
         int statusID = static_cast<int> ( idStatusDecoder(status)["sensorID"] );
         
         // We are assuming that the input matrix, noise and
         // status collections are aligned according to the sensorID. In
         // other words, we are assuming the element i-th in the all the
         // collections corresponds to the same sensorID.
         if (dataID != noiseID) {
           streamlog_out(ERROR3) << "noiseID inconsistent!" << std::endl << std::endl;      
           exit(-1); 
         }
         if (dataID != statusID) {
           streamlog_out(ERROR3) << "statusID inconsistent!" << std::endl << std::endl;      
           exit(-1); 
         }

         if (matrix->getChargeValues().size() != noise->getChargeValues().size()) {   
           streamlog_out(ERROR3) << "Input data and noise are incompatible" << endl
              << "Detector " << iDetector << " has " <<  matrix->getChargeValues().size() << " pixels in the full data " << endl
              << "while " << noise->getChargeValues().size() << " in the noise data " << endl;
           exit(-1); 
         }

         if (matrix->getChargeValues().size() != status->getADCValues().size()){   
           streamlog_out(ERROR3) << "Input data and status are incompatible" << endl
              << "Detector " << iDetector << " has " <<  matrix->getChargeValues().size() << " pixels in the full data " << endl
              << "while " << status->getADCValues().size() << " in the status data " << endl;
           exit(-1); 
         }  
          
         // Before continuing it is a better idea to verify that the
         // input data is really a full frame data set and not already
         // sparsified! To do that we try to get from the cellDecoder
         // some information stored only in full frame mode, like the
         // "xMin". If xMin it is not there, then an un-caught
         // exception will be thrown and the execution stopped here.
         idMatrixDecoder(matrix)["xMin"];
       }
       _isFirstEvent = false;
     }
        
     // Initialize sparse data collection
     LCCollectionVec * sparseDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);
     CellIDEncoder<TrackerDataImpl> sparseDataEncoder( "sensorID:6,sparsePixelType:5" , sparseDataCollection );  
     	    
     // Loop over all sensors 
     for (unsigned int iSensor = 0; iSensor < DataCollection->size(); iSensor++) { 
       
       CellIDDecoder<TrackerDataImpl> idMatrixDecoder(DataCollection);
       
       // we are assuming that the input matrix, noise and
       // status collections are aligned according to the sensorID. In
       // other words, we are assuming the element i-th in the all the
       // collections corresponds to the same sensorID.
       TrackerDataImpl * matrix = dynamic_cast<TrackerDataImpl* > (DataCollection->getElementAt(iSensor));
       TrackerDataImpl * noise = dynamic_cast < TrackerDataImpl * >   (NoiseCollection->getElementAt(iSensor));
       TrackerRawDataImpl * status = dynamic_cast < TrackerRawDataImpl * >(StatusCollection->getElementAt(iSensor));
       
       // 
       // Matrix encoding 
       MatrixDecoder matrixDecoder(idMatrixDecoder, matrix); 
       int sensorID = static_cast<int> ( idMatrixDecoder(matrix)["sensorID"] );
        
       // 
       // Access data vectors
       FloatVec matrixVec  = matrix->getChargeValues(); 
       FloatVec noiseVec   = noise->getChargeValues();
       ShortVec statusVec  = status->getADCValues();
       
       // Prepare a TrackerData to store zs pixels
       TrackerDataImpl* zspixels = new TrackerDataImpl;
       
       // Set description for zspixels 
       sparseDataEncoder["sensorID"] = sensorID;
       sparseDataEncoder["sparsePixelType"] = 0;
       sparseDataEncoder.setCellID( zspixels );
       
       for ( int iPixel = 0; iPixel < (int) matrixVec.size(); iPixel++ ) {
            
          if (  statusVec[iPixel]  == 0 ) {
            
            float data  = matrixVec[iPixel];
                 
            // Use a global ADU threshold 
            float threshold = _ThresholdVec[iSensor]; 
            
            // Use a global SNR threshold
            if ( _useSNRCut ) threshold *= noiseVec[iPixel]; 
             
            if ( data >= threshold  ) {
              
              int xPixel,yPixel; 
              matrixDecoder.getXYFromIndex(iPixel,  xPixel, yPixel);
              
              // Store pixel data int EUTelescope format 
              zspixels->chargeValues().push_back( xPixel );
              zspixels->chargeValues().push_back( yPixel );
              zspixels->chargeValues().push_back( data );    
                
              // Print detailed pixel summary 
              streamlog_out(MESSAGE1) << "Found pixel Nr. " << iPixel << " on sensor " << sensorID << endl
                                      << "   x:" << xPixel << ", y:" << yPixel << ", charge:" << data
                                      << ", quality:" << statusVec[iPixel] << endl;
               
              
               
            }  // Endif threshold 
          }  // Endif status
       } // Endfor ipixel
       
       // Add zs pixels to collection 
       sparseDataCollection->push_back( zspixels );
             
     } // end sensor loop

     evt->addCollection(sparseDataCollection, _SparseDataCollectionName );
     
   } catch(DataNotAvailableException &e){}  
   _nEvt ++ ;
}


//
// Method called after each event to check the data processed
//
void RawDataSparsifier::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void RawDataSparsifier::end()
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
void RawDataSparsifier::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "RawDataSparsifier Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}

} // Namespace
