// HotPixelKiller processor   
// 		
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "HotPixelKiller.h"

// Include DEPFETTrackTools 
#include "DEPFET.h" 
#include "MatrixDecoder.h"

// Include basic C
#include <iostream>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/LCTime.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <EVENT/LCParameters.h>
#include <IMPL/LCCollectionVec.h>


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
                              "This is the name of the root file with DQM histos",
                              _rootFileName, string("NoiseDB.root"));
   
   registerProcessorParameter("NoiseDBFileName",
                              "This is the name of the noiseDB file",
                              _noiseDBFileName, string("NoiseDB.slcio"));
   
   registerProcessorParameter ("MaxOccupancy",
                              "Maximum hit rate (hits/events) for normal pixels",
                              _maxOccupancy, static_cast < float >(0.01));
   
   registerProcessorParameter ("EventsForMask",
                              "Events used to compute hot pixel mask",
                              _eventsForMask, static_cast < int >(20000));
   
   registerProcessorParameter ("OfflineZSThreshold",
                              "Zero suppression threshold for digits [ADU]",
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

   // Book all needed histograms 
   bookHistos();
    
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
   
   
   if ( isFirstEvent() ) {
      
      // Initialize pedestal, commom mode and quality algorithms
      initializeAlgorithms(evt); 
      
      _isFirstEvent = false;
   }
    
   accumulateHits(evt);
   
   if ( _nEvt == _eventsForMask ) {
     streamlog_out(MESSAGE4) << "Compute hotpixel mask "
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
  
  // Compute final pixel mask   
  computeMask(); 
  
  // CPU time end
  _timeCPU = clock()/1000 - _timeCPU;
  
  streamlog_out ( MESSAGE3 ) << "Writing the noise DB lcio file" << endl;
  
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
    idStatusEncoder["uMax"]       = noOfXPixels-1;
    idStatusEncoder["vMax"]       = noOfYPixels-1;
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

  // Root DQM histograms 
 
  std::string histoName;   

  for ( int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
     
    // Read geometry info for sensor
    int ipl = _planeNumbers[iDetector];       
    Det& Sensor = _detector.GetDet(ipl);      
      
    int noOfXPixels = Sensor.GetNColumns(); 
    int noOfYPixels = Sensor.GetNRows(); 
    
    // Loop over all pixels 
    MatrixDecoder matrixDecoder( noOfXPixels, noOfYPixels );  
    
    // Count masked channels 
    double nmasked = 0;  
    
    for (int yPixel = 0; yPixel < noOfYPixels; yPixel++) {
      for (int xPixel = 0; xPixel < noOfXPixels; xPixel++) {
        
        int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);
          
        if ( _status[iDetector][iPixel] != 0 ) nmasked++;
            
        histoName = "h2mask_sensor"+to_string( ipl );
        _histoMap2D[histoName]->Fill(xPixel,yPixel, _status[iDetector][iPixel]  );
             
        histoName = "h2occ_sensor"+to_string( ipl );
        _histoMap2D[histoName]->Fill(xPixel,yPixel, _hitCounter[iDetector][iPixel]/_nEvt  ); 
    
        histoName = "hocc_sensor"+to_string( ipl );
        _histoMap[ histoName ]->Fill(_hitCounter[iDetector][iPixel]/_nEvt); 
           
      }
    }
    
    histoName = "hnhits_sensor"+to_string( ipl );
    double nhits = _histoMap[ histoName ]->GetMean();
      
    histoName = "hnhits_mask_sensor"+to_string( ipl );
    double ngoodhits = _histoMap[ histoName ]->GetMean();
    
    histoName = "hnhits";
    _histoMap[ histoName ]->SetBinContent(ipl+1, nhits);  
      
    histoName = "hnhits_mask";
    _histoMap[ histoName ]->SetBinContent(ipl+1, ngoodhits);
      
    histoName = "hnmasked";
    _histoMap[ histoName ]->SetBinContent(ipl+1, nmasked);
  
  }
   
  streamlog_out(MESSAGE3) << std::endl << " " << "Writing ROOT file ..." << std::endl;   
     
  _rootFile->Write(); 
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
void HotPixelKiller::accumulateHits(LCEvent * evt) {
  
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
          
      FloatVec sparseDigits = trackerData->getChargeValues();
      int nDigits = sparseDigits.size()/3; 
      int nGoodDigits =  0;       
      
      streamlog_out(MESSAGE2) << "Module sensorID " << sensorID << " having digits " <<  nDigits << std::endl; 
       
      for ( int index=0; index<nDigits;  index++) { 
            
        int xPixel = static_cast<int> (sparseDigits[index * 3]);
        int yPixel = static_cast<int> (sparseDigits[index * 3 + 1]);
        float chargeValue =  sparseDigits[index * 3 + 2]; 

        int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);  
        
        // If digit signal below threshold, skip it
        if ( chargeValue < _offlineZSCut ) {
          streamlog_out(MESSAGE2) << "  Signal below ZS cut. Skipping it." << std::endl; 
          continue;
        }
        
        // Check if digit cellIDs are valid
        if ( xPixel < 0 || xPixel >= noOfXPixels || yPixel < 0 || yPixel >= noOfYPixels) 
        {   
          streamlog_out(MESSAGE2) << "Digit on sensor " << sensorID 
                                    << "   x:" << xPixel << ", y:" << yPixel 
                                    << " in event " << evt->getEventNumber() 
                                    << " is out of range!! " << std::endl;    
          continue; 
        }
               
        // Check if digit is a duplicate   
        for ( int id=0; id<index;  id++) {
            
          int xPixel2 = static_cast<int> (sparseDigits[id * 3]);
          int yPixel2 = static_cast<int> (sparseDigits[id * 3 + 1]);
          
          if ( xPixel2 == xPixel && yPixel2==yPixel ) {
            streamlog_out(MESSAGE2) << "Digit on sensor " << sensorID 
                                    << "   x:" << xPixel << ", y:" << yPixel 
                                    << " in event " << evt->getEventNumber() 
                                    << " is duplicated!! " << std::endl; 
            continue; 
          }     
        }  
        
        // Estimate firing rate for pixel 
        _hitCounter[ iDetector ][ iPixel ]++; 
        
        // Count number of good hits per frame  
        if (_status[iDetector][iPixel] == 0 ) ++nGoodDigits;
          
        // Print detailed pixel summary, for testing/debugging only !!! 
        streamlog_out(MESSAGE1) << "Digit on sensor " << sensorID 
                                << std::endl;  
        streamlog_out(MESSAGE1) << "   x: " << xPixel << ", y: " << yPixel 
                                << ", charge: " << chargeValue << ", quality: " << _status[iDetector][iPixel]
                                << std::endl;
        
        
         
      } // End digit loop
      
      // Time to review quality of event
      _rootEventNumber = evt->getEventNumber();
      _rootDetectorID = sensorID;
      _rootNHits = nDigits; 
      _rootNGoodHits = nGoodDigits;
              
      _rootFile->cd("");
      _rootEventTree->Fill();
    
      std::string histoName; 
      
      histoName = "hnhits_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill(nDigits);
      
      if ( _nEvt > _eventsForMask ) {
        histoName = "hnhits_mask_sensor"+to_string( ipl );
        _histoMap[ histoName ]->Fill(nGoodDigits);
      }
        
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
      _status.push_back(ShortVec( nPixel, 0 ));
      
      // Initialize hit counter 
      _hitCounter.push_back(FloatVec( nPixel, 0.));
          
    } // End of detector loop
    
    
    // Create DQM histograms 
     
    for ( int iDetector = 0 ; iDetector < _noOfDetector ; ++iDetector ) {
      
      // Get the TrackerRawData object from the collection for this detector
      TrackerData *trackerData = dynamic_cast < TrackerData * >(collectionVec->getElementAt (iDetector));
       
      // Get DAQ ID of sensor 
      CellIDDecoder< TrackerData > DataDecoder( collectionVec );
      int sensorID =  DataDecoder( trackerData ) ["sensorID"] ;
       
      // Read geometry info for sensor 
      int ipl = _detector.GetPlaneNumber(sensorID);      
      Det& Sensor = _detector.GetDet(ipl);
      
      std::string dirName; 
      std::string histoName;
      
      dirName = "Sensor"+to_string( ipl );
      _rootFile->mkdir(dirName.c_str());    
      
      dirName = "/Sensor"+to_string(ipl)+"/";
      _rootFile->cd(dirName.c_str());
      
      int uBins = Sensor.GetNColumns();  
      int vBins = Sensor.GetNRows(); 
      
      histoName = "h2mask_sensor"+to_string( ipl );
      _histoMap2D[histoName] = new TH2D(histoName.c_str(), "" ,uBins, 0, uBins, vBins, 0, vBins);
      _histoMap2D[histoName]->SetXTitle("uCell [cellID]"); 
      _histoMap2D[histoName]->SetYTitle("vCell [cellID]"); 
      _histoMap2D[histoName]->SetZTitle("mask");     
      _histoMap2D[histoName]->SetStats( false );      
      
      histoName = "h2occ_sensor"+to_string( ipl );
      _histoMap2D[histoName] = new TH2D(histoName.c_str(), "" ,uBins, 0, uBins, vBins, 0, vBins);
      _histoMap2D[histoName]->SetXTitle("uCell [cellID]"); 
      _histoMap2D[histoName]->SetYTitle("vCell [cellID]"); 
      _histoMap2D[histoName]->SetZTitle("hit occupancy");       
      _histoMap2D[histoName]->SetStats( false );    
      
      histoName = "hocc_sensor"+to_string( ipl );
      _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 10000000, 0, 0.01);
      _histoMap[ histoName ]->SetXTitle("hit occupancy"); 
      _histoMap[ histoName ]->SetYTitle("pixels");   
    
      histoName = "hnhits_sensor"+to_string( ipl );
      _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, 0, 200);
      _histoMap[ histoName ]->SetXTitle("hits per event"); 
      _histoMap[ histoName ]->SetYTitle("events");   
      
      histoName = "hnhits_mask_sensor"+to_string( ipl );
      _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, 0, 200);
      _histoMap[ histoName ]->SetXTitle("hits per event (after mask)"); 
      _histoMap[ histoName ]->SetYTitle("events");
            
    }
   
  } catch (DataNotAvailableException& e) {
        streamlog_out ( WARNING2 ) << "No input collection " << _zeroSuppressedDataCollectionName << " is not available in the current event" << endl;
  }
     
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

void HotPixelKiller::bookHistos()
{   
  
  _rootFile = new TFile(_rootFileName.c_str(),"recreate");
  _rootFile->cd("");
   
  _rootEventTree = new TTree("Event","Event info");
  _rootEventTree->Branch("iEvt"            ,&_rootEventNumber    ,"iEvt/I");
  _rootEventTree->Branch("det"             ,&_rootDetectorID     ,"det/I");
  _rootEventTree->Branch("nhits"           ,&_rootNHits          ,"nhits/I");
  _rootEventTree->Branch("ngoodhits"       ,&_rootNGoodHits      ,"ngoodhits/I");
  
  // Get number of sensors
  int nSens = _detector.GetNSensors();
  
  std::string histoName;

  histoName = "hnhits";
  _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", nSens, 0, nSens);
  _histoMap[ histoName ]->SetXTitle("plane number"); 
  _histoMap[ histoName ]->SetYTitle("mean hits per event (before mask)");   
  
  histoName = "hnhits_mask";
  _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", nSens, 0, nSens);
  _histoMap[ histoName ]->SetXTitle("plane number"); 
  _histoMap[ histoName ]->SetYTitle("mean hits per event (after mask)");   
  
  histoName = "hnmasked";
  _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", nSens, 0, nSens);
  _histoMap[ histoName ]->SetXTitle("plane number"); 
  _histoMap[ histoName ]->SetYTitle("masked pixels");   
  
}

} // Namespace

