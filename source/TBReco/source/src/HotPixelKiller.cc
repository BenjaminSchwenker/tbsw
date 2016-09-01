// HotPixelKiller processor   
// 		
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "HotPixelKiller.h"

// Include DEPFETTrackTools 
#include "DEPFET.h" 
#include "MatrixDecoder.h"

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


// Used namespaces
using namespace std; 
using namespace lcio ;
using namespace marlin ;

namespace depfet {

//
// Instantiate this object
//
HotPixelKiller aHotPixelKiller ;

//
// Constructor
//
HotPixelKiller::HotPixelKiller() : Processor("HotPixelKiller")
{
   
   // Processor description
   _description = "HotPixelKiller: Masking of hot pixels";
   
   
   //   
   // First of all, we need to register the input/output collections
   
   registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
                            "Name of zero suppressed data collection",
                            _zeroSuppressedDataCollectionName, string("zsdata"));

   registerProcessorParameter("StatusCollectionName",
                            "Name of status data collection",
                            _statusDataCollectionName, string("status"));
   
   registerProcessorParameter("OutputRootFileName",
                              "This is the name of the output root file",
                              _rootFileName, string("NoiseDB.root"));

   registerProcessorParameter("NoiseDBFileName",
                              "This is the name of the noise (hotpixel) data base",
                              _noiseDBFileName, string("NoiseDB.slcio"));

   registerProcessorParameter ("MaxOccupancy",
                              "Maximum occupancy for good pixels",
                              _maxOccupancy, static_cast < float >(0.01));

   registerProcessorParameter ("OfflineZSThreshold",
                              "Offline zero suppression threshold",
                              _offlineZSCut, static_cast < float >(0));


}

//
// Method called at the beginning of data processing
//
void HotPixelKiller::init() {
   
   // Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
   _timeCPU = clock()/1000;
              
   // Read detector constants from gear file
   _detector.ReadGearConfiguration();    
   
   // Print set parameters
   printProcessorParams();
   
   // ROOT Output 
   _rootFile = new TFile(_rootFileName.c_str(),"recreate");
   _rootFile->cd("");
     
   // Note: this tree contains snapshots for complete module, sampled every 5000 events
   _rootOccTree = new TTree("Occupancy","Occupancy info");
   _rootOccTree->Branch("det"             ,&_rootDetectorID     ,"det/I");
   _rootOccTree->Branch("px_x"            ,&_rootCol            ,"px_x/I");
   _rootOccTree->Branch("px_y"            ,&_rootRow            ,"px_y/I");
   _rootOccTree->Branch("status"          ,&_rootStatus         ,"status/I");
   _rootOccTree->Branch("hitFreq"         ,&_rootHitFrequency   ,"hitFreq/D");
   
   _rootEventTree = new TTree("Event","Event info");
   _rootEventTree->Branch("iEvt"            ,&_rootEventNumber    ,"iEvt/I");
   _rootEventTree->Branch("det"             ,&_rootDetectorID     ,"det/I");
   _rootEventTree->Branch("nhits"           ,&_rootNHits          ,"nhits/I");
   _rootEventTree->Branch("ngoodhits"       ,&_rootNGoodHits      ,"ngoodhits/I");
    
}

//
// Method called for each run
//
void HotPixelKiller::processRunHeader(LCRunHeader * run)
{
     
  // Print run number
  streamlog_out(MESSAGE3) << "Processing run: "
                          << (run->getRunNumber())
                          << std::endl << std::endl;
  _nRun++ ;
  
  // Write the current header to the output noise file
  LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
  
  try {
    lcWriter->open(_noiseDBFileName, LCIO::WRITE_NEW);
  } catch (IOException& e) {
    cerr << e.what() << endl;
    return;
  }
  
  LCRunHeaderImpl * lcHeader  = new LCRunHeaderImpl;
  lcHeader->setRunNumber( 0 );
   
  lcWriter->writeRunHeader(lcHeader);
  delete lcHeader;
  
  lcWriter->close();
  
}

//
// Method called for each event
//
void HotPixelKiller::processEvent(LCEvent * evt)
{
   
   // Print event number
   if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                 << (evt->getEventNumber())
                                                                 << std::endl << std::endl;
   
   // More detailed event numbering for testing
   streamlog_out(MESSAGE2) << std::endl << "Starting with Event Number " << evt->getEventNumber()  << std::endl; 
    
   
   if ( isFirstEvent() ) {
      
      // Initialize pedestal, commom mode and quality algorithms
      initializeAlgorithms(evt); 
      
      _isFirstEvent = false;
   }
    
   readEventData(evt);
   
   if ( _nEvt == 20000 ) {
     streamlog_out(MESSAGE4) << "Compute intermediate mask "
                             << (evt->getEventNumber())
                             << std::endl << std::endl;
      
     computeMask();
   }
   
}


//
// Method called after each event to check the data processed
//
void HotPixelKiller::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void HotPixelKiller::end()
{
  
  // Compute hit pixel mask   
  computeMask(); 
  
  // CPU time end
  _timeCPU = clock()/1000 - _timeCPU;
  
  streamlog_out ( MESSAGE3 ) << "Writing the noise data base file" << endl;
  
  LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();

  try {
     lcWriter->open(_noiseDBFileName,LCIO::WRITE_APPEND);
  } catch (IOException& e) {
     cerr << e.what() << endl;
     return;
  }

  LCEventImpl * event = new LCEventImpl();
  event->setRunNumber(_nRun);

  LCTime * now = new LCTime;
  event->setTimeStamp(now->timeStamp());
  delete now;

  LCCollectionVec * statusCollection   = new LCCollectionVec(LCIO::TRACKERRAWDATA);
  
  for ( int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
  
    streamlog_out( MESSAGE3 ) << "Writing status data for sensor " << iDetector << endl;
    
    // Read geometry info for sensor
    int ipl = _planeNumbers[iDetector];       
    Det& adet = _detector.GetDet(ipl);      
      
    int noOfXPixels = adet.GetNColumns(); 
    int noOfYPixels = adet.GetNRows(); 
    int sensorID = adet.GetDAQID(); 
    
    TrackerRawDataImpl * statusMatrix   = new TrackerRawDataImpl;
    
    CellIDEncoder<TrackerRawDataImpl> idStatusEncoder(DEPFET::MATRIXDEFAULTENCODING, statusCollection);
  
    idStatusEncoder["sensorID"]   = sensorID;
    idStatusEncoder["xMin"]       = 0;
    idStatusEncoder["xMax"]       = noOfXPixels-1;
    idStatusEncoder["yMin"]       = 0;
    idStatusEncoder["yMax"]       = noOfYPixels-1;
    idStatusEncoder.setCellID(statusMatrix);

    statusMatrix->setADCValues(_status[iDetector]);
    statusCollection->push_back(statusMatrix);
          
    // Loop over all pixels 
    MatrixDecoder matrixDecoder( noOfXPixels, noOfYPixels );  
    
    for (int yPixel = 0; yPixel < noOfYPixels; yPixel++) {
      for (int xPixel = 0; xPixel < noOfXPixels; xPixel++) {
             
        if (statusMatrix->getADCValues() [ matrixDecoder.getIndexFromXY(xPixel,yPixel) ]   != 0) {
          
            streamlog_out(MESSAGE3) << "Masking pixel on sensor " << iDetector 
                                  << "   x: " << xPixel << ", y: " << yPixel 
                                  << endl;
        
        }      
           
      }
    }  
  }
  
  event->addCollection(statusCollection, _statusDataCollectionName);

  lcWriter->writeEvent(event);
  delete event;
  lcWriter->close();
   
  streamlog_out(MESSAGE3) << std::endl << " " << "Writing ROOT file ..." << std::endl;   
     
  // Fill occupancy ntuple    
  fillOccupancyTuple();
  
  _rootFile->Write();
   
  streamlog_out(MESSAGE3) << std::endl << " " << "Closing ROOT file ..." << std::endl;   
   
  _rootFile->Close();
  
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
void HotPixelKiller::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "HotPixelKiller Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   

}


 

//
// Method doing sparse pixel conversion 
//
void HotPixelKiller::readEventData(LCEvent * evt) {

  try {
    
    // Open sparse pixel collection 
    LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection ( _zeroSuppressedDataCollectionName ));
          
    for ( size_t iDetector = 0; iDetector < collectionVec->size() ; iDetector++) {
           
      // Get the TrackerRawData object from the collection for this plane
      TrackerDataImpl *trackerData = dynamic_cast < TrackerDataImpl * >(collectionVec->getElementAt (iDetector)); 
      
      // Get DAQ ID of sensor 
      CellIDDecoder< TrackerDataImpl > DataDecoder( collectionVec );  
      int sensorID =  DataDecoder( trackerData ) ["sensorID"] ;
      
      // Read geometry info for sensor 
      int ipl = _detector.GetPlaneNumber(sensorID);      
      Det& adet = _detector.GetDet(ipl);
      
      int noOfXPixels = adet.GetNColumns(); 
      int noOfYPixels = adet.GetNRows();
      MatrixDecoder matrixDecoder( noOfXPixels, noOfYPixels); 
          
      FloatVec sparsePixels = trackerData->getChargeValues();
      int nPixels = sparsePixels.size()/3; 
      int nGoodPixels =  0;       


      streamlog_out(MESSAGE2) << "Module sensorID " << sensorID << " having pixels " <<  nPixels << std::endl; 
       
      for ( int index=0; index<nPixels;  index++) { 
            
        int xPixel = static_cast<int> (sparsePixels[index * 3]);
        int yPixel = static_cast<int> (sparsePixels[index * 3 + 1]);
        float chargeValue =  sparsePixels[index * 3 + 2]; 
        int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);  
        
        // If pixel signal below threshold, skip it
        if ( chargeValue < _offlineZSCut ) {
          streamlog_out(MESSAGE2) << "  Signal below ZS cut. Skipping it." << std::endl; 
          continue;
        }
        
        // Estimate firing rate for all pixels 
        _hitCounter[ iDetector ][ iPixel ]++; 
        
        // Count number of good hits per frame
        if (_status[iDetector][iPixel] == 0) ++nGoodPixels;
          
        // Print detailed pixel summary, for testing/debugging only !!! 
        streamlog_out(MESSAGE1) << "Pixel on sensor " << sensorID 
                                << std::endl;  
        streamlog_out(MESSAGE1) << "   x: " << xPixel << ", y: " << yPixel 
                                << ", charge: " << chargeValue << ", quality: " << _status[iDetector][iPixel]
                                << std::endl;
        
        // Check if pixel address is valid
          
        if ( xPixel < 0 || xPixel >= noOfXPixels || yPixel < 0 || yPixel >= noOfYPixels) 
        {   
          streamlog_out(MESSAGE4) << "Pixel on sensor " << sensorID 
                                    << "   x:" << xPixel << ", y:" << yPixel 
                                    << " in event " << evt->getEventNumber() 
                                    << " is out of range!! " << std::endl;    
          
        }
        
        // Check if pixel is a duplicate   
        
        for ( int id=0; id<index;  id++) {
            
          int xPixel2 = static_cast<int> (sparsePixels[id * 3]);
          int yPixel2 = static_cast<int> (sparsePixels[id * 3 + 1]);
          
          if ( xPixel2 == xPixel && yPixel2==yPixel ) {
            streamlog_out(MESSAGE2) << "Pixel on sensor " << sensorID 
                                    << "   x:" << xPixel << ", y:" << yPixel 
                                    << " in event " << evt->getEventNumber() 
                                    << " is duplicated!! " << std::endl; 
          }     
        }  
         
      } // End pixel loop
      
      // Time to review quality of event
      _rootEventNumber = evt->getEventNumber();
      _rootDetectorID = sensorID;
      _rootNHits = nPixels; 
      _rootNGoodHits = nGoodPixels;
              
      _rootFile->cd("");
      _rootEventTree->Fill();
        
    }  // End loop on detectors
    
    // increment the event number
    ++_nEvt;
    
  } catch (DataNotAvailableException& e) {
        streamlog_out ( WARNING2 ) << "No input collection " << _zeroSuppressedDataCollectionName << " is not available in the current event" << endl;
  }
  
}


void HotPixelKiller::initializeAlgorithms(LCEvent * evt) {
  
  // Clear all internal vectors
  _status.clear();
  _hitCounter.clear();
   
  try {
      
    LCCollectionVec * collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection ( _zeroSuppressedDataCollectionName ));
    _noOfDetector = (int)collectionVec->size();
     
    

    for ( int iDetector = 0 ; iDetector < _noOfDetector ; ++iDetector ) {
      
      // Get the TrackerRawData object from the collection for this detector
      TrackerData *trackerData = dynamic_cast < TrackerData * >(collectionVec->getElementAt (iDetector));
      
      // Get DAQ ID of sensor 
      CellIDDecoder< TrackerData > DataDecoder( collectionVec );
      int sensorID =  DataDecoder( trackerData ) ["sensorID"] ;

      // Read geometry info for sensor 
      int ipl = _detector.GetPlaneNumber(sensorID);      
      Det& adet = _detector.GetDet(ipl);
      
      _planeNumbers.push_back(ipl); 
      int nPixel = adet.GetNColumns() * adet.GetNRows() ;
       
      // Initialize all pixels as GOODPIXEL
      _status.push_back(ShortVec( nPixel, kGoodPixel ));
      
      // Initialize hit counter 
      _hitCounter.push_back(FloatVec( nPixel, 0.));
      
       
    } // End of detector loop
    
  } catch (DataNotAvailableException& e) {
        streamlog_out ( WARNING2 ) << "No input collection " << _zeroSuppressedDataCollectionName << " is not available in the current event" << endl;
  }
     
}

   
void HotPixelKiller::fillOccupancyTuple() {

  for ( int iDetector = 0; iDetector < _noOfDetector ; iDetector++) {
        
    // Read geometry info for sensor
    int ipl = _planeNumbers[iDetector];       
    Det& adet = _detector.GetDet(ipl);      
      
    int noOfXPixels = adet.GetNColumns(); 
    int noOfYPixels = adet.GetNRows(); 
    int sensorID = adet.GetDAQID(); 
      
    MatrixDecoder matrixDecoder( noOfXPixels, noOfYPixels); 
      
    // start looping on all pixels
    for (int yPixel = 0; yPixel < noOfYPixels; yPixel++) {
      for (int xPixel = 0; xPixel < noOfXPixels; xPixel++) {
          
        int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);
          
        // ROOT Output 
        _rootDetectorID = sensorID ;  
        _rootCol = xPixel;                
        _rootRow = yPixel; 
        _rootStatus = _status[iDetector][iPixel];     
        _rootHitFrequency = _hitCounter[iDetector][iPixel]/_nEvt;  
        _rootFile->cd("");
        _rootOccTree->Fill();
    
      } // end loop on xPixel
    } // end loop on yPixel
  }  // end loop on detectors

}


void HotPixelKiller::computeMask() {

    
  for (  int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
              
    // Read geometry info for sensor
    int ipl = _planeNumbers[iDetector];       
    Det& adet = _detector.GetDet(ipl);      
      
    int noOfXPixels = adet.GetNColumns(); 
    int noOfYPixels = adet.GetNRows(); 
    int sensorID = adet.GetDAQID(); 
        
    // Use standard matrix encoding 
    MatrixDecoder matrixDecoder( noOfXPixels, noOfYPixels );  
      
    // Loop over all pixels 
    for (int yPixel = 0; yPixel < noOfYPixels; yPixel++) {
      for (int xPixel = 0; xPixel < noOfXPixels; xPixel++) {
           
        int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel); 
          
        // Mask pixel with very high hit frequency -> hot pixel killer 
        double firingFreq =  _hitCounter[ iDetector ][ iPixel ] / _nEvt;
        if ( firingFreq  > _maxOccupancy ) {
           
          
          _status[iDetector][iPixel] = 1;
          
          streamlog_out(MESSAGE1) << "Mask pixel on sensor " << sensorID 
                                  << "   x: " << xPixel << ", y: " << yPixel 
                                  << "   (" << firingFreq <<  ")" << endl;
                                  

        }      
           
      }
    } // 2x pixel loop      
  } // iDetector loop   
  
}

} // Namespace

