// HotStripKiller processor   
// 		
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "HotStripKiller.h"

// Include DEPFETTrackTools 
#include "DEPFET.h" 

// Include basic C
#include <iostream>
#include <iomanip>

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
HotStripKiller aHotStripKiller ;

//
// Constructor
//
HotStripKiller::HotStripKiller() : Processor("HotStripKiller")
{
   
   // Processor description
   _description = "HotStripKiller: Masking of hot strips";
   
   
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
                              "Maximum hit rate (hits/events) for normal strips",
                              _maxOccupancy, static_cast < float >(0.5));
   
   registerProcessorParameter ("EventsForMask",
                              "Events used to compute hot pixel mask",
                              _eventsForMask, static_cast < int >(20000));
   
   registerProcessorParameter ("ZSThresholdU",
                              "Zero suppression threshold for u strips [ADU]",
                              _offlineZSCutU, static_cast < float >(0));

   registerProcessorParameter ("ZSThresholdV",
                              "Zero suppression threshold for v strips [ADU]",
                              _offlineZSCutV, static_cast < float >(0));


}

//
// Method called at the beginning of data processing
//
void HotStripKiller::init() {
   
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
void HotStripKiller::processRunHeader(LCRunHeader * run)
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
void HotStripKiller::processEvent(LCEvent * evt)
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
     streamlog_out(MESSAGE4) << "Compute hot strip mask "
                             << (evt->getEventNumber())
                             << std::endl << std::endl;
      
     computeMask();
   }
   
}


//
// Method called after each event to check the data processed
//
void HotStripKiller::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void HotStripKiller::end()
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
    Det& Sensor = _detector.GetDet(ipl);      
      
    int nUCells = Sensor.GetNColumns(); 
    int nVCells = Sensor.GetNRows(); 
    int sensorID = Sensor.GetDAQID(); 
    
    TrackerRawDataImpl * statusMatrix   = new TrackerRawDataImpl;
    
    CellIDEncoder<TrackerRawDataImpl> idStatusEncoder(DEPFET::MATRIXDEFAULTENCODING, statusCollection);
  
    idStatusEncoder["sensorID"]   = sensorID;
    idStatusEncoder["uMax"]       = nUCells-1;
    idStatusEncoder["vMax"]       = nVCells-1;
    idStatusEncoder.setCellID(statusMatrix);

    // Concatenate u/v masks 
    ShortVec _status( nUCells + nVCells ,0);
    
    for (int cell = 0; cell < nUCells; cell++) {
      
      _status[cell] = _statusU[iDetector][cell];     
      if ( _status[cell] != 0) {
        streamlog_out(MESSAGE3) << "Masking uStrip " << cell << endl;
      }      
           
    }

    for (int cell = 0; cell < nVCells; cell++) {
      
      _status[cell+nUCells] = _statusU[iDetector][cell];    
      if ( _status[cell] != 0) {
        streamlog_out(MESSAGE3) << "Masking vStrip " << cell << endl;
      }      
           
    }
    
    statusMatrix->setADCValues( _status);
    statusCollection->push_back(statusMatrix);  
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
      
    int nUCells = Sensor.GetNColumns(); 
    int nVCells = Sensor.GetNRows(); 
    
    // Count masked channels 
    double nmaskedU = 0;  
    double nmaskedV = 0;  

    for (int cell = 0; cell < nUCells; cell++) {
      
      if ( _statusU[iDetector][cell] != 0 ) nmaskedU++;    
      
      histoName = "hmaskU_sensor"+to_string( ipl );
      _histoMap[histoName]->Fill(cell, _statusU[iDetector][cell]  ); 
      
      histoName = "hoccU_sensor"+to_string( ipl );
      _histoMap[histoName]->Fill(cell, _hitCounterU[iDetector][cell]/_nEvt  ); 
      
      histoName = "hoccU_histo_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill(_hitCounterU[iDetector][cell]/_nEvt); 
      
    }

    for (int cell = 0; cell < nVCells; cell++) {
      
      if ( _statusV[iDetector][cell] != 0 ) nmaskedV++;  
      
      histoName = "hmaskV_sensor"+to_string( ipl );
      _histoMap[histoName]->Fill(cell, _statusV[iDetector][cell]  ); 
      
      histoName = "hoccV_sensor"+to_string( ipl );
      _histoMap[histoName]->Fill(cell, _hitCounterV[iDetector][cell]/_nEvt  ); 
      
      histoName = "hoccU_histo_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill(_hitCounterV[iDetector][cell]/_nEvt); 
           
    }
    
    histoName = "hnhitsU_sensor"+to_string( ipl );
    double nhitsU = _histoMap[ histoName ]->GetMean();
      
    histoName = "hnhitsU_mask_sensor"+to_string( ipl );
    double ngoodhitsU = _histoMap[ histoName ]->GetMean();
    
    histoName = "hnhitsU";
    _histoMap[ histoName ]->SetBinContent(ipl+1, nhitsU);   
    
    histoName = "hnhitsU_mask";
    _histoMap[ histoName ]->SetBinContent(ipl+1, ngoodhitsU);
    
    histoName = "hnmaskedU";
    _histoMap[ histoName ]->SetBinContent(ipl+1, nmaskedU);

    histoName = "hnhitsV_sensor"+to_string( ipl );
    double nhitsV = _histoMap[ histoName ]->GetMean();
      
    histoName = "hnhitsV_mask_sensor"+to_string( ipl );
    double ngoodhitsV = _histoMap[ histoName ]->GetMean();
    
    histoName = "hnhitsV";
    _histoMap[ histoName ]->SetBinContent(ipl+1, nhitsV);  
    
    histoName = "hnhitsV_mask";
    _histoMap[ histoName ]->SetBinContent(ipl+1, ngoodhitsV);
    
    histoName = "hnmaskedV";
    _histoMap[ histoName ]->SetBinContent(ipl+1, nmaskedV);
    
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
void HotStripKiller::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "HotStripKiller Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   

}


 

//
// Method doing digit counting 
//
void HotStripKiller::accumulateHits(LCEvent * evt) {
  
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
      Det& Sensor = _detector.GetDet(ipl);
      
      int nUCells = Sensor.GetNColumns(); 
      int nVCells = Sensor.GetNRows();
          
      FloatVec sparseDigits = trackerData->getChargeValues();
      int nDigits = sparseDigits.size()/3; 
       
      streamlog_out(MESSAGE2) << "Module sensorID " << sensorID << " having digits " <<  nDigits << std::endl; 
      
      int nDigitsU = 0;  
      int nDigitsV = 0;  
      int nGoodDigitsU = 0;  
      int nGoodDigitsV = 0;            
      
      for ( int index=0; index<nDigits;  index++) { 
            
        int isU = static_cast<int> (sparseDigits[index * 3]);
        int cell = static_cast<int> (sparseDigits[index * 3 + 1]);
        float signal =  sparseDigits[index * 3 + 2]; 
        
        if (isU) {

          nDigitsU++;

          if ( signal < _offlineZSCutU ) {
            streamlog_out(MESSAGE2) << " Signal on uStrip below ZS cut. Skipping it." << std::endl; 
            continue;
          }
               
          if ( cell < 0 || cell >= nUCells  ) {
            streamlog_out(MESSAGE2) << "Digit on sensor " << sensorID 
                                    << " uCell: " << cell 
                                    << " in event " << evt->getEventNumber() 
                                    << " is out of range!! " << std::endl;    
            continue;  
          }
          
          // Count the digit
          _hitCounterU[ iDetector ][ cell ]++; 
        
          // Count number of good hits per frame  
          if (_statusU[iDetector][cell] == 0 ) ++nGoodDigitsU;
          
          // Print detailed strip summary, for testing/debugging only !!! 
          streamlog_out(MESSAGE1) << "Digit on sensor " << sensorID 
                                  << " uCell: " << cell 
                                  << ", charge: " << signal << std::endl;
        
        } else {

          nDigitsV++;

          if ( signal < _offlineZSCutV ) {
            streamlog_out(MESSAGE2) << " Signal on vStrip below ZS cut. Skipping it." << std::endl; 
            continue;
          }
               
          if ( cell < 0 || cell >= nVCells  ) {
            streamlog_out(MESSAGE2) << "Digit on sensor " << sensorID 
                                    << " vCell: " << cell 
                                    << " in event " << evt->getEventNumber() 
                                    << " is out of range!! " << std::endl;    
            continue;  
          }
          
          // Count the digit
          _hitCounterV[ iDetector ][ cell ]++; 
        
          // Count number of good hits per frame  
          if (_statusV[iDetector][cell] == 0 ) ++nGoodDigitsV;
          
          // Print detailed strip summary, for testing/debugging only !!! 
          streamlog_out(MESSAGE1) << "Digit on sensor " << sensorID 
                                  << " vCell: " << cell 
                                  << ", charge: " << signal << std::endl;
           
          
        }  // isU
         
      } // End digit loop
       
      std::string histoName; 
      
      histoName = "hnhitsU_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill(nDigitsU);
           
      if ( _nEvt > _eventsForMask ) {
        histoName = "hnhitsU_mask_sensor"+to_string( ipl );
        _histoMap[ histoName ]->Fill(nGoodDigitsU);
      }

      histoName = "hnhitsV_sensor"+to_string( ipl );
      _histoMap[ histoName ]->Fill(nDigitsV);
           
      if ( _nEvt > _eventsForMask ) {
        histoName = "hnhitsV_mask_sensor"+to_string( ipl );
        _histoMap[ histoName ]->Fill(nGoodDigitsV);
      }
            
    }  // End loop on detectors
    
    // increment the event number
    ++_nEvt;
    
  } catch (DataNotAvailableException& e) {
        streamlog_out ( WARNING2 ) << "No input collection " << _zeroSuppressedDataCollectionName << " is not available in the current event" << endl;
  }
  
}


void HotStripKiller::initializeAlgorithms(LCEvent * evt) {
  
  // Clear all internal vectors
  _statusU.clear();
  _hitCounterU.clear();
  _statusV.clear();
  _hitCounterV.clear();
   
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
      _planeNumbers.push_back(ipl); 

    
      Det& Sensor = _detector.GetDet(ipl);
             
      // Initialize all u strips as GOOD
      _statusU.push_back(ShortVec( Sensor.GetNColumns(), 0 ));
      
      // Initialize hit counter for all u strips to zero
      _hitCounterU.push_back(FloatVec( Sensor.GetNColumns(), 0.));

      // Initialize all v strips as GOOD
      _statusV.push_back(ShortVec( Sensor.GetNRows(), 0 ));
      
      // Initialize hit counter for all v strips to zero
      _hitCounterV.push_back(FloatVec( Sensor.GetNRows(), 0.));
          
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
      
      histoName = "hmaskU_sensor"+to_string( ipl );
      _histoMap[histoName] = new TH1D(histoName.c_str(), "" ,uBins, 0, uBins);
      _histoMap[histoName]->SetXTitle("uCell [cellID]"); 
      _histoMap[histoName]->SetYTitle("mask");      
      _histoMap[histoName]->SetStats( false );   
      
      histoName = "hmaskV_sensor"+to_string( ipl );
      _histoMap[histoName] = new TH1D(histoName.c_str(), "" ,vBins, 0, vBins);
      _histoMap[histoName]->SetXTitle("vCell [cellID]"); 
      _histoMap[histoName]->SetYTitle("mask");      
      _histoMap[histoName]->SetStats( false );   
      
      histoName = "hoccU_sensor"+to_string( ipl );
      _histoMap[histoName] = new TH1D(histoName.c_str(), "" ,uBins, 0, uBins);
      _histoMap[histoName]->SetXTitle("uCell [cellID]");  
      _histoMap[histoName]->SetYTitle("hit occupancy");       
      _histoMap[histoName]->SetStats( false );    

      histoName = "hoccV_sensor"+to_string( ipl );
      _histoMap[histoName] = new TH1D(histoName.c_str(), "" ,vBins, 0, vBins);
      _histoMap[histoName]->SetXTitle("vCell [cellID]");  
      _histoMap[histoName]->SetYTitle("hit occupancy");       
      _histoMap[histoName]->SetStats( false );   
      
      histoName = "hoccU_histo_sensor"+to_string( ipl );
      _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 10000000, 0, 0.01);
      _histoMap[ histoName ]->SetXTitle("hit occupancy"); 
      _histoMap[ histoName ]->SetYTitle("u strips");   

      histoName = "hoccV_histo_sensor"+to_string( ipl );
      _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 10000000, 0, 0.01);
      _histoMap[ histoName ]->SetXTitle("hit occupancy"); 
      _histoMap[ histoName ]->SetYTitle("v strips");   
    
      histoName = "hnhitsU_sensor"+to_string( ipl );
      _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, 0, 200);
      _histoMap[ histoName ]->SetXTitle("u hits per event"); 
      _histoMap[ histoName ]->SetYTitle("events");   

      histoName = "hnhitsV_sensor"+to_string( ipl );
      _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, 0, 200);
      _histoMap[ histoName ]->SetXTitle("v hits per event"); 
      _histoMap[ histoName ]->SetYTitle("events"); 
      
      histoName = "hnhitsU_mask_sensor"+to_string( ipl );
      _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, 0, 200);
      _histoMap[ histoName ]->SetXTitle("maked u hits per event"); 
      _histoMap[ histoName ]->SetYTitle("events");

      histoName = "hnhitsV_mask_sensor"+to_string( ipl );
      _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", 200, 0, 200);
      _histoMap[ histoName ]->SetXTitle("maked v hits per event"); 
      _histoMap[ histoName ]->SetYTitle("events");
            
    }
   
  } catch (DataNotAvailableException& e) {
        streamlog_out ( WARNING2 ) << "No input collection " << _zeroSuppressedDataCollectionName << " is not available in the current event" << endl;
  }
     
}

   



void HotStripKiller::computeMask() {
   
  for (  int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
              
    // Read geometry info for sensor
    int ipl = _planeNumbers[iDetector];       
    Det& Sensor = _detector.GetDet(ipl);      
     
    int nUCells = Sensor.GetNColumns(); 
    int nVCells = Sensor.GetNRows();
    int sensorID = Sensor.GetDAQID(); 
        
    // Loop over all vStrips  
    for (int cell = 0; cell < nVCells; cell++) {
          
      double occupancy =  _hitCounterV[ iDetector ][ cell ] / _nEvt;
      if ( occupancy  > _maxOccupancy ) {
           
        _statusV[iDetector][cell] = 1;
          
        streamlog_out(MESSAGE1) << "Mask vCell on sensor " << sensorID 
                                << " cellID: " << cell 
                                << " occupancy: " << occupancy  << endl;
                                  
      }      
           
    }
    
    // Loop over all uStrips  
    for (int cell = 0; cell < nUCells; cell++) {
          
      double occupancy =  _hitCounterU[ iDetector ][ cell ] / _nEvt;
      if ( occupancy  > _maxOccupancy ) {
           
        _statusU[iDetector][cell] = 1;
          
        streamlog_out(MESSAGE1) << "Mask uCell on sensor " << sensorID 
                                << " cellID: " << cell 
                                << " occupancy: " << occupancy  << endl;
                                  
      }      
           
    }
    
  } // iDetector loop   
  
}

void HotStripKiller::bookHistos()
{   
  
  _rootFile = new TFile(_rootFileName.c_str(),"recreate");
  _rootFile->cd("");
   
  // Get number of sensors
  int nSens = _detector.GetNSensors();
  
  std::string histoName;

  histoName = "hnhitsU";
  _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", nSens, 0, nSens);
  _histoMap[ histoName ]->SetXTitle("plane number"); 
  _histoMap[ histoName ]->SetYTitle("mean hits on u strips per event");   
  
  histoName = "hnhitsU_mask";
  _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", nSens, 0, nSens);
  _histoMap[ histoName ]->SetXTitle("plane number"); 
  _histoMap[ histoName ]->SetYTitle("mean hits on u strips per event (masked)");   
  
  histoName = "hnmaskedU";
  _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", nSens, 0, nSens);
  _histoMap[ histoName ]->SetXTitle("plane number"); 
  _histoMap[ histoName ]->SetYTitle("masked u strips");   

  histoName = "hnhitsV";
  _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", nSens, 0, nSens);
  _histoMap[ histoName ]->SetXTitle("plane number"); 
  _histoMap[ histoName ]->SetYTitle("mean hits on v strips per event");   
  
  histoName = "hnhitsV_mask";
  _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", nSens, 0, nSens);
  _histoMap[ histoName ]->SetXTitle("plane number"); 
  _histoMap[ histoName ]->SetYTitle("mean hits on v strips per event (masked)");   
  
  histoName = "hnmaskedV";
  _histoMap[ histoName ] = new TH1D(histoName.c_str(), "", nSens, 0, nSens);
  _histoMap[ histoName ]->SetXTitle("plane number"); 
  _histoMap[ histoName ]->SetYTitle("masked v strips");   
  
}

} // Namespace

