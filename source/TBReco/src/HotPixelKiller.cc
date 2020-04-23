// HotPixelKiller processor   
// 		
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "HotPixelKiller.h"

// Include TBTools 
#include "DEPFET.h" 
#include "TBDetector.h"

// Include basic C
#include <iostream>
#include <iomanip>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>

// Include ROOT classes
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>


// Used namespaces
using namespace std; 
using namespace lcio ;
using namespace marlin ;
using namespace std::string_literals;

namespace depfet {

  //
  // Instantiate this object
  //
  HotPixelKiller aHotPixelKiller ;
  
  //
  // Constructor
  //
  HotPixelKiller::HotPixelKiller() : Processor("HotPixelKiller"),_inputDecodeHelper("")
  {
      
    // Processor description
    _description = "HotPixelKiller: Masking of hot pixels";
    
    //   
    // First of all, we need to register the input/output collections
    
    registerInputCollection (LCIO::TRACKERDATA, "InputCollectionName",
                             "Name of zero suppressed data collection",
                             _zeroSuppressedDataCollectionName, string("zsdata"));
     
    registerProcessorParameter("NoiseDBFileName",
                               "This is the name of the noiseDB file",
                               _noiseDBFileName, string("NoiseDB.root"));
    
    registerProcessorParameter ("MaxOccupancy",
                                "Maximum hit rate (hits/events) for normal pixels",
                                _maxOccupancy, static_cast < float >(0.01));
    
    registerProcessorParameter ("MinOccupancy",
                                "Minimum hit rate (hits/events) for normal pixels. Use negative value to deactivate cut. Use 0 to mask pixel with zero hits in run data.",
                                _minOccupancy, static_cast < float >(-1));

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
               
    
    
    // Print set parameters
    printProcessorParams();
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
    
    
    try {
      
      // Open sparse pixel collection 
      LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection ( _zeroSuppressedDataCollectionName ));
          
      for ( size_t iDetector = 0; iDetector < collectionVec->size() ; iDetector++) {
           
        // Get the TrackerRawData object from the collection for this plane
        TrackerDataImpl *trackerData = dynamic_cast < TrackerDataImpl * >(collectionVec->getElementAt (iDetector)); 
        
        // Get DAQ ID of sensor 
        CellIDDecoder< TrackerDataImpl > DataDecoder( collectionVec,&_inputDecodeHelper);
        int sensorID =  DataDecoder( trackerData ) ["sensorID"s] ;
        
        // Read geometry info for sensor 
        int ipl = TBDetector::GetInstance().GetPlaneNumber(sensorID);      
        const Det& adet = TBDetector::Get(ipl);
   
        int minUCell = adet.GetMinUCell();
        int maxUCell = adet.GetMaxUCell();
        int minVCell = adet.GetMinVCell();
        int maxVCell = adet.GetMaxVCell();	
        int nUCells = maxUCell-minUCell+1;
        int nVCells = maxVCell-minVCell+1;
        
        // Register counter variables for new sensorID 
        if (_hitCounterMap.find(sensorID) == _hitCounterMap.end() ) {
          
          int nPixel = nUCells * nVCells ;
          
          // Initialize hit counter 
          _hitCounterMap[sensorID] = FloatVec( nPixel, 0.);    
        }
        
        FloatVec sparseDigits = trackerData->getChargeValues();
        int nDigits = sparseDigits.size()/3; 
            
        streamlog_out(MESSAGE2) << "Module sensorID " << sensorID << " having digits " <<  nDigits << std::endl; 
         
        for ( int index=0; index<nDigits;  index++) { 
            
          int iU = static_cast<int> (sparseDigits[index * 3]);
          int iV = static_cast<int> (sparseDigits[index * 3 + 1]);
          float chargeValue =  sparseDigits[index * 3 + 2]; 
           
          // Describe pixel by unique ID
          int uniqPixelID  = adet.encodePixelID(iV, iU);   
          
          // If digit signal below threshold, skip it
          if ( chargeValue < _offlineZSCut ) {
            streamlog_out(MESSAGE2) << "  Signal below ZS cut. Skipping it." << std::endl; 
            continue;
          }
          
          // Check if digit cellIDs are valid
          if ( iU < minUCell || iU > maxUCell || iV < minVCell || iV > maxVCell ) 
          {   
            streamlog_out(MESSAGE2) << "Digit on sensor " << sensorID 
                                      << "   iU:" << iU << ", iV:" << iV 
                                      << " in event " << evt->getEventNumber() 
                                      << " is out of range!! " << std::endl;    
            continue; 
          }
               
          // Check if digit is a duplicate   
          for ( int id=0; id<index;  id++) {
            
            int iU2 = static_cast<int> (sparseDigits[id * 3]);
            int iV2 = static_cast<int> (sparseDigits[id * 3 + 1]);
          
            if ( iU2 == iU && iV2==iV ) {
              streamlog_out(MESSAGE2) << "Digit on sensor " << sensorID 
                                      << "   iU:" << iU << ", iV:" << iV 
                                      << " in event " << evt->getEventNumber() 
                                      << " is duplicated!! " << std::endl; 
              continue; 
            }     
          }  
          
          // Estimate firing rate for pixel 
          _hitCounterMap[ sensorID ][ uniqPixelID ]++; 
          
        } // End digit loop
            
      }  // End loop on detectors
    
      // increment the event number
      ++_nEvt;
      
    } catch (DataNotAvailableException& e) {
      streamlog_out ( WARNING2 ) << "No input collection " << _zeroSuppressedDataCollectionName << " is not available in the current event" << endl;
    }
     
  }
  
  //
  // Method called after each event to check the data processed
  //
  void HotPixelKiller::check( LCEvent * )
  {
  }
  
  //
  // Method called after all data processing
  //
  void HotPixelKiller::end()
  {
    
   
    
    // CPU time end
    _timeCPU = clock()/1000 - _timeCPU;
    
    streamlog_out ( MESSAGE3 ) << "Writing the noise DB file" << endl;
    
    TFile * _rootFile = new TFile( _noiseDBFileName.c_str(),"recreate");
    _rootFile->cd("");

    std::map< std::string, TH2F *> _histoMap;
    std::map< std::string, TH2F *> _histoMapOcc;
    
    // Loop over all registered sensors 
    for(auto it = _hitCounterMap.begin(); it != _hitCounterMap.end(); it++) {
      auto sensorID = it->first;
      auto&& hitVec = it->second;
        
      // Read geometry info for sensor 
      int ipl = TBDetector::GetInstance().GetPlaneNumber(sensorID);      
      const Det& adet = TBDetector::Get(ipl);
      
      int minUCell = adet.GetMinUCell();
      int maxUCell = adet.GetMaxUCell();
      int minVCell = adet.GetMinVCell();
      int maxVCell = adet.GetMaxVCell();	
      int nUCells = maxUCell-minUCell+1;
      int nVCells = maxVCell-minVCell+1; 

      string histoName = "hDB_sensor"+to_string(sensorID) + "_mask";
      _histoMap[histoName] = new TH2F(histoName.c_str(), "" , nUCells, minUCell, maxUCell+1, nVCells, minVCell, maxVCell+1); // +1 because maxCell is the not included upper border of the last bin, but maxCell should have a bin
      _histoMap[histoName]->SetXTitle("uCell [cellID]"); 
      _histoMap[histoName]->SetYTitle("vCell [cellID]"); 
      _histoMap[histoName]->SetZTitle("mask");     
      _histoMap[histoName]->SetStats( false );

      string occhistoName = "hDB_sensor"+to_string(sensorID) + "_occupancy";
      _histoMapOcc[occhistoName] = new TH2F(occhistoName.c_str(), "" , nUCells, minUCell, maxUCell+1, nVCells, minVCell, maxVCell+1); // +1 because maxCell is the not included upper border of the last bin, but maxCell should have a bin
      _histoMapOcc[occhistoName]->SetXTitle("uCell [cellID]");
      _histoMapOcc[occhistoName]->SetYTitle("vCell [cellID]");
      _histoMapOcc[occhistoName]->SetZTitle("occupancy");
      _histoMapOcc[occhistoName]->SetStats( false );
      
      int nMasked = 0; 
       
      // Loop over all pixels / histogram bins 
      for (int iV = minVCell; iV < maxVCell; iV++) {
        for (int iU = minUCell; iU < maxUCell; iU++) {
            
          int uniqPixelID  = adet.encodePixelID(iV, iU);   
          
          // Mask pixel with very high hit frequency -> hot pixel killer 
          double occupancy =  hitVec[ uniqPixelID ] / _nEvt;
          _histoMapOcc[occhistoName]->SetBinContent(iU-minUCell+1,iV-minVCell+1, occupancy );

          if ( (occupancy  > _maxOccupancy) or (occupancy <= _minOccupancy ) ) {
             
            nMasked++;     
            streamlog_out(MESSAGE1) << "Mask pixel on sensorID " << sensorID 
                                    << "   iU: " << iU << ", iV: " << iV 
                                    << "   (" << occupancy <<  ")" << endl;
            
            _histoMap[histoName]->SetBinContent(iU-minUCell+1,iV-minVCell+1, 1 );

          } else {
            _histoMap[histoName]->SetBinContent(iU-minUCell+1,iV-minVCell+1, 0 );
          }      
        }
      } 
        
      streamlog_out(MESSAGE3) << "Number of masked pixels on sensorID " << sensorID 
                              << " is " << nMasked << endl;  
    }
    
    streamlog_out(MESSAGE3) << std::endl << " " << "Writing ROOT file ..." << std::endl;   
     
    _rootFile->Write(); 
    _rootFile->Close();
    delete _rootFile;    
    
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
  
 
 


} // Namespace

