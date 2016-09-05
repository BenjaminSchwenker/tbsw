// DEPFET HotPixelKiller Tool   
// 		
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "DEPFETHotPixelMasker.h"

// Include DEPFETTrackTools 
#include "DEPFET.h" 
#include "MatrixDecoder.h"
#include "TBHit.h"

// Include basic C
#include <cstdlib>
#include <iostream>
#include <memory>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/LCTime.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerHitImpl.h>


// Used namespaces
using namespace std; 
using namespace lcio ;
using namespace marlin ;

namespace depfet {

//
// Instantiate this object
//
DEPFETHotPixelMasker aDEPFETHotPixelMasker ;

//
// Constructor
//
DEPFETHotPixelMasker::DEPFETHotPixelMasker() : Processor("DEPFETHotPixelMasker")
{
   
   // Processor description
   _description = "DEPFETHotPixelMasker: Remove hits from hot pixels";
   
    
   //   
   // First of all, we need to register the input/output collections
   
   registerInputCollection (LCIO::TRACKERRAWDATA, "StatusCollectionName",
                           "Input status collection name",
                           _statusCollectionName, string("status")); 
    
   vector< string > inputHitCollectionNameVecExample;
   inputHitCollectionNameVecExample.push_back( "hit" );
   
   registerInputCollections (LCIO::TRACKERHIT, "InputHitCollectionNameVec",
                            "Hit collection names",
                            _inputHitCollectionNameVec, inputHitCollectionNameVecExample );

}

//
// Method called at the beginning of data processing
//
void DEPFETHotPixelMasker::init() {
   
   // Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
   _timeCPU = clock()/1000;
              
   // Read detector constants from gear file
   _detector.ReadGearConfiguration();    

  
   // Print set parameters
   printProcessorParams();
       
}

//
// Method called for each run
//
void DEPFETHotPixelMasker::processRunHeader(LCRunHeader * run)
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
void DEPFETHotPixelMasker::processEvent(LCEvent * evt)
{
   
  // Print event number
  if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                 << (evt->getEventNumber())
                                                                 << std::endl << std::endl;
   
  // More detailed event numbering for testing
  streamlog_out(MESSAGE2) << std::endl << "Starting with Event Number " << evt->getEventNumber()  << std::endl; 
    
  // Filter hits 
  // ===========================
    
  for ( size_t iCol = 0 ; iCol < _inputHitCollectionNameVec.size(); ++iCol ) {
     
    try {

       LCCollectionVec * statusCollection = dynamic_cast < LCCollectionVec * > (evt->getCollection(_statusCollectionName));
       
       LCCollectionVec * hitcol = dynamic_cast < LCCollectionVec * > (evt->getCollection( _inputHitCollectionNameVec.at( iCol ) ) );
         
       LCCollectionVec * new_hitcol = new LCCollectionVec(LCIO::TRACKERHIT) ;
              
       for ( int ihit = 0 ; ihit < (int) hitcol->size() ; ++ihit ) {
             
         // Built a TBHit     
         TrackerHitImpl * lciohit = dynamic_cast< TrackerHitImpl* > ( hitcol->getElementAt( ihit ) );
         TBHit aHit ( lciohit ); 
         
           
         // Read hit data 
         int sensorID = aHit.GetDAQID();
         int ipl = _detector.GetPlaneNumber(sensorID); 
         
         // skip hit from unregistered layers 
         if (ipl < 0) continue; 
           
         double um = aHit.GetCoord()[0][0];
         double vm = aHit.GetCoord()[1][0];   
          
         Det& adet = _detector.GetDet(ipl);
         int noOfXPixels = adet.GetNColumns(); 
         int noOfYPixels = adet.GetNRows();
         MatrixDecoder matrixDecoder( noOfXPixels, noOfYPixels); 
         
         int xPixel = adet.GetColumnFromCoord( um, vm );  
         int yPixel = adet.GetRowFromCoord( um, vm );  

         if ( xPixel >= noOfXPixels  || xPixel < 0 || yPixel >= noOfYPixels  || yPixel < 0) {
           streamlog_out(MESSAGE2) << " skipping bad hit!!!!! ---------------- " << std::endl; 
          continue;
         }
         
         int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);  

         // Get status collection  
         TrackerRawDataImpl * status = dynamic_cast<TrackerRawDataImpl*> (statusCollection->getElementAt( ipl ));  
          
         // Only keep hits from good channels 
         if ( status->getADCValues() [  iPixel ]   == 0 ) {

           // Make LCIO TrackerHit
           TrackerHitImpl * trackerhit = aHit.MakeLCIOHit();  
            
           // Add hit to the hit collection
           new_hitcol->push_back( trackerhit );
           
         }
                      
       } // End for hit loop   
       
       // Derive a new name for filtered hit collection
       std::string newHitCollectionName = _inputHitCollectionNameVec.at( iCol ) +"2";
        
       // Store hitCollection in LCIO file
       evt->addCollection( new_hitcol, newHitCollectionName );   
          
    } catch (lcio::DataNotAvailableException& e) {
       streamlog_out ( MESSAGE2 ) << "Not able to get collection "
                                  << _inputHitCollectionNameVec.at( iCol )
                                  << "\nfrom event " << evt->getEventNumber()
                                  << " in run " << evt->getRunNumber()  << endl;
       
    }  
  }   
  
   
  // increment the event number
  ++_nEvt;
     
   
}


//
// Method called after each event to check the data processed
//
void DEPFETHotPixelMasker::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void DEPFETHotPixelMasker::end()
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
void DEPFETHotPixelMasker::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "DEPFETHotPixelMasker Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   

}







} // Namespace

