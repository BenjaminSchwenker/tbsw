// DEPFET HotPixelKiller Tool   
// 		
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "DEPFETHotPixelKillerForHits.h"

// Include DEPFETTrackTools 
#include "DEPFET.h" 
#include "MatrixDecoder.h"
#include "HitFactory.h"
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
DEPFETHotPixelKillerForHits aDEPFETHotPixelKillerForHits ;

//
// Constructor
//
DEPFETHotPixelKillerForHits::DEPFETHotPixelKillerForHits() : Processor("DEPFETHotPixelKillerForHits")
{
   
   // Processor description
   _description = "DEPFETHotPixelKillerForHits: Masking of hot pixels";
   
    
   //   
   // First of all, we need to register the input/output collections
    
   vector< string > inputHitCollectionNameVecExample;
   inputHitCollectionNameVecExample.push_back( "hit" );
   
   registerInputCollections (LCIO::TRACKERHIT, "InputHitCollectionNameVec",
                            "Hit collection names",
                            _inputHitCollectionNameVec, inputHitCollectionNameVecExample );

   registerProcessorParameter("OutputRootFileName",
                              "This is the name of the output root file",
                              _rootFileName, string("NoiseDB.root"));

   registerProcessorParameter("NoiseDBFileName",
                              "This is the name of the noise (hotpixel) data base",
                              _noiseDBFileName, string("NoiseDB.slcio"));

   registerProcessorParameter ("MaxOccupancy",
                              "Maximum occupancy for good pixels",
                              _maxOccupancy, static_cast < float >(0.01));

   


}

//
// Method called at the beginning of data processing
//
void DEPFETHotPixelKillerForHits::init() {
   
   // Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
   _timeCPU = clock()/1000;
              
   // Read detector constants from gear file
   _detector.ReadGearConfiguration();    

   // Initialize pedestal, commom mode and quality algorithms
   initializeAlgorithms(); 
   
   // Print set parameters
   printProcessorParams();
   
   // ROOT Output 
   _rootFile = new TFile(_rootFileName.c_str(),"recreate");
   _rootFile->cd("");
     
   // Note: this tree contains snapshots for complete module, sampled every 5000 events
   _rootOccTree = new TTree("Occupancy","Occupancy info");
   _rootOccTree->Branch("det"             ,&_rootDetectorID     ,"det/I");
   _rootOccTree->Branch("ipl"             ,&_rootPlaneNumber    ,"ipl/I");
   _rootOccTree->Branch("px_x"            ,&_rootCol            ,"px_x/I");
   _rootOccTree->Branch("px_y"            ,&_rootRow            ,"px_y/I");
   _rootOccTree->Branch("status"          ,&_rootStatus         ,"status/I");
   _rootOccTree->Branch("hitFreq"         ,&_rootHitFrequency   ,"hitFreq/D");
   
   _rootEventTree = new TTree("Event","Event info");
   _rootEventTree->Branch("iEvt"            ,&_rootEventNumber    ,"iEvt/I");
   _rootEventTree->Branch("det"             ,&_rootDetectorID     ,"det/I");
   _rootEventTree->Branch("ipl"             ,&_rootPlaneNumber    ,"ipl/I");
   _rootEventTree->Branch("nhits"           ,&_rootNHits          ,"nhits/I");
   _rootEventTree->Branch("ngoodhits"       ,&_rootNGoodHits      ,"ngoodhits/I");
    
}

//
// Method called for each run
//
void DEPFETHotPixelKillerForHits::processRunHeader(LCRunHeader * run)
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
void DEPFETHotPixelKillerForHits::processEvent(LCEvent * evt)
{
   
   // Print event number
   if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                 << (evt->getEventNumber())
                                                                 << std::endl << std::endl;
   
   // More detailed event numbering for testing
   streamlog_out(MESSAGE2) << std::endl << "Starting with Event Number " << evt->getEventNumber()  << std::endl; 
    
   
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
void DEPFETHotPixelKillerForHits::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void DEPFETHotPixelKillerForHits::end()
{

  // CPU time end
  _timeCPU = clock()/1000 - _timeCPU;
  
  streamlog_out ( MESSAGE4 ) << "Writing the noise data base file" << endl;
  
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
  
  LCCollectionVec * statusCollection = new LCCollectionVec(LCIO::TRACKERRAWDATA);
  
  for ( int ipl = 0; ipl < _detector.GetNSensors(); ipl++) {
    
    // Read geometry info for sensor     
    Det& adet = _detector.GetDet(ipl);      
      
    int noOfXPixels = adet.GetNColumns(); 
    int noOfYPixels = adet.GetNRows(); 
    
    
    TrackerRawDataImpl * statusMatrix   = new TrackerRawDataImpl;
    
    CellIDEncoder<TrackerRawDataImpl> idStatusEncoder(DEPFET::MATRIXDEFAULTENCODING, statusCollection);
  
    idStatusEncoder["sensorID"]   = adet.GetDAQID();
    idStatusEncoder["uMax"]       = noOfXPixels-1;
    idStatusEncoder["vMax"]       = noOfYPixels-1;
    idStatusEncoder.setCellID(statusMatrix);

    statusMatrix->setADCValues(_status[ipl]);
    statusCollection->push_back(statusMatrix);
  }

  event->addCollection(statusCollection, "status");

  lcWriter->writeEvent(event);
  delete event;
  lcWriter->close();
   
  streamlog_out(MESSAGE3) << std::endl << " " << "Writing ROOT file ..." << std::endl;   
  
  // Compute hit pixel mask   
  computeMask(); 
     
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
void DEPFETHotPixelKillerForHits::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "DEPFETHotPixelKillerForHits Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   

}


 

//
// Method doing sparse pixel conversion 
//
void DEPFETHotPixelKillerForHits::readEventData(LCEvent * evt) {

  // Read input hit collections
  // ===========================
  
  // Hit factory sorts hits according to plane    
  HitFactory HitStore( _detector );   
    
  for ( size_t iCol = 0 ; iCol < _inputHitCollectionNameVec.size(); ++iCol ) {
     
    try {
       
       LCCollectionVec * hitcol = dynamic_cast < LCCollectionVec * > (evt->getCollection( _inputHitCollectionNameVec.at( iCol ) ) );
              
       for ( int ihit = 0 ; ihit < (int) hitcol->size() ; ++ihit ) {
             
         // Built a TBHit     
         TrackerHitImpl * lciohit = dynamic_cast< TrackerHitImpl* > ( hitcol->getElementAt( ihit ) );
         TBHit RecoHit ( lciohit ); 

         // Add hit to store
         HitStore.AddRecoHit(RecoHit);
                      
       } // End for hit loop      
          
    } catch (lcio::DataNotAvailableException& e) {
       streamlog_out ( MESSAGE2 ) << "Not able to get collection "
                                  << _inputHitCollectionNameVec.at( iCol )
                                  << "\nfrom event " << evt->getEventNumber()
                                  << " in run " << evt->getRunNumber()  << endl;
       
    }  
  }   
  
  //                        
  // Loop over detector planes 
     
  for(int ipl=0;ipl<_detector.GetNSensors();++ipl)  { 
    
    int nHits = HitStore.GetNHits(ipl); 
    
    Det& adet = _detector.GetDet(ipl);
      
    int noOfXPixels = adet.GetNColumns(); 
    int noOfYPixels = adet.GetNRows();
    MatrixDecoder matrixDecoder( noOfXPixels, noOfYPixels); 
            
    int nGoodPixels =  0; 

    // Loop over all hits on this detector 
    for (int ihit = 0; ihit < nHits; ++ihit ) {  

      streamlog_out(MESSAGE2) << "  hit: " <<  ihit << std::endl; 

      TBHit & anyhit = HitStore.GetRecoHitFromID(ihit, ipl);
           
      // Measured hit coordinates
      double um = anyhit.GetCoord()[0][0];
      double vm = anyhit.GetCoord()[1][0];   
      
      int xPixel = adet.GetColumnFromCoord( um, vm );  
      int yPixel = adet.GetRowFromCoord( um, vm );  
      int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);  

      streamlog_out(MESSAGE2) << " xpixel: " << xPixel  << " ypixel: " << yPixel  << std::endl; 

      if ( xPixel >= noOfXPixels  || xPixel < 0 || yPixel >= noOfYPixels  || yPixel < 0) {
        streamlog_out(MESSAGE2) << " skipping bad hit!!!!! ---------------- " << std::endl; 
        continue;
      }
             
      // Estimate firing rate for all pixels 
      _hitCounter[ ipl ][ iPixel ]++; 
        
      // Count number of good hits per frame
      if (_status[ipl][iPixel] == 0) ++nGoodPixels;
          
      
          
    } // End pixel loop
      
    // Time to review quality of event
    _rootEventNumber = evt->getEventNumber();
    _rootDetectorID = adet.GetDAQID();
    _rootPlaneNumber = ipl;
    _rootNHits = nHits; 
    _rootNGoodHits = nGoodPixels;
              
    _rootFile->cd("");
    _rootEventTree->Fill();
        
  }  // End loop on detectors
    
  // increment the event number
  ++_nEvt;
    
}


void DEPFETHotPixelKillerForHits::initializeAlgorithms() {
  
  // Clear all internal vectors
  _status.clear();
  _hitCounter.clear();
      
  for ( int ipl = 0 ; ipl < _detector.GetNSensors() ; ++ipl ) {
        
    Det& adet = _detector.GetDet(ipl);
       
    int nPixel = adet.GetNColumns() * adet.GetNRows() ;
       
    // Initialize all pixels as GOODPIXEL
    _status.push_back(ShortVec( nPixel, 0 ));
      
    // Initialize hit counter 
    _hitCounter.push_back(FloatVec( nPixel, 0.));
         
  } // End of detector loop
       
}

   
void DEPFETHotPixelKillerForHits::fillOccupancyTuple() {

  for ( int ipl = 0 ; ipl < _detector.GetNSensors() ; ++ipl ) {
        
    // Read geometry info for sensor
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
        _rootPlaneNumber = ipl ;  
        _rootCol = xPixel;                
        _rootRow = yPixel; 
        _rootStatus = _status[ipl][iPixel];     
        _rootHitFrequency = _hitCounter[ipl][iPixel]/_nEvt;  
        _rootFile->cd("");
        _rootOccTree->Fill();
    
      } // end loop on xPixel
    } // end loop on yPixel
  }  // end loop on detectors

}


void DEPFETHotPixelKillerForHits::computeMask() {

    
  for ( int ipl = 0 ; ipl < _detector.GetNSensors() ; ++ipl ) {
              
    // Read geometry info for sensor       
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
        double firingFreq =  _hitCounter[ ipl ][ iPixel ] / _nEvt;
        if ( firingFreq  > _maxOccupancy ) {
           
          
          _status[ipl][iPixel] = 1;
          streamlog_out( MESSAGE2 ) << "HotPixel masking of pixel number " << iPixel
                                      << " on detector " << sensorID
                                      << " (" << firingFreq <<  ")" << endl;  
        }      
           
      }
    } // 2x pixel loop      
  } // ipl   
  
}

} // Namespace

