// SimpleNoiseDBCreator - Marlin Processor   
//
// See SimpleNoiseDBCreator.h for full documentation of this processor. 
//		
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "SimpleNoiseDBCreator.h"

// Include DEPFETTrackTools header files
#include "DEPFET.h"

// Include basic C
#include <iostream>
#include <memory>
#include <iomanip>

// Include LCIO classes
#include <lcio.h>

// Used namespaces
using namespace std; 
using namespace lcio; 
using namespace marlin;

namespace depfet {


//
// Instantiate this object
//
SimpleNoiseDBCreator aSimpleNoiseDBCreator ;


//
// Constructor
//
SimpleNoiseDBCreator::SimpleNoiseDBCreator() : Processor("SimpleNoiseDBCreator"), 
                                                       _pedestalCollectionVec(NULL),
                                                       _noiseCollectionVec(NULL),
                                                       _statusCollectionVec(NULL) 
{

   // Processor description
   _description = "SimpleNoiseDBCreator: produces pedestal / noise / status with user provided values";
     
   registerOutputCollection (LCIO::TRACKERDATA, "PedestalCollectionName",
                            "Name of pedestal collection",
                            _pedestalCollectionName, string ("pedestal"));
   
   registerOutputCollection (LCIO::TRACKERDATA, "NoiseCollectionName",
                            "Name of noise collection",
                            _noiseCollectionName, string("noise"));
   
   registerOutputCollection (LCIO::TRACKERRAWDATA, "StatusCollectionName",
                            "Name of status collection",
                            _statusCollectionName, string("status"));
   
   
   //   
   // Define processor parameters 
   
   const size_t nDetectorExample = 6;   
   
   
   FloatVec initPedeExample(nDetectorExample, 0.);
   registerOptionalParameter("InitPedestalValue",
                            "The initial value of pedestal (one value for detector)",
                            _initPedestal, initPedeExample );
   
   FloatVec initNoiseExample(nDetectorExample, 1.);
   registerOptionalParameter("InitNoiseValue",
                            "The initial value of noise (one value for detector)",
                            _initNoise, initNoiseExample );
   
}

//
// Method called at the beginning of data processing
//
void SimpleNoiseDBCreator::init() {
  
  // Initialize variables
  _iRun = 0 ;
  _iEvt = 0 ;
  
  // Read detector constants from gear file
  _detector.ReadGearConfiguration();    
}

//
// Method called for each run
//
void SimpleNoiseDBCreator::processRunHeader(LCRunHeader * run)
{
     
  // Print run number
  streamlog_out(MESSAGE3) << "Processing run: "
                          << (run->getRunNumber())
                          << std::endl << std::endl;
  _iRun++ ;

 
 
}

//
// Method called for each event
//
void SimpleNoiseDBCreator::processEvent(LCEvent * evt)
{
   
  ++_iEvt;
   
  if (isFirstEvent()) {
    
    _pedestalCollectionVec = new LCCollectionVec(LCIO::TRACKERDATA);
    _noiseCollectionVec    = new LCCollectionVec(LCIO::TRACKERDATA);
    _statusCollectionVec   = new LCCollectionVec(LCIO::TRACKERRAWDATA);
    
    int nDetectors = _detector.GetNSensors();
    
    if ( (int)_initPedestal.size() != nDetectors ) {
      streamlog_out(ERROR) << "Wrong number of pedestal values, does not match number of sensors!" << endl;  
      exit(-1);  
    }

    if ( (int)_initNoise.size() != nDetectors ) {
      streamlog_out(ERROR) << "Wrong number of noise values, does not match number of sensors!" << endl; 
      exit(-1);    
    }
    
    for (int iDetector = 0; iDetector < nDetectors; iDetector++) {

      Det& adet = _detector.GetDet(iDetector);
      
      int nPixel = adet.GetNColumns() * adet.GetNRows() ;
       
      TrackerRawDataImpl * status   = new TrackerRawDataImpl;
      CellIDEncoder<TrackerRawDataImpl> statusEncoder(DEPFET::MATRIXDEFAULTENCODING, _statusCollectionVec);
      statusEncoder["sensorID"] = adet.GetDAQID();
      statusEncoder["xMin"]     = 0;
      statusEncoder["yMin"]     = 0;
      statusEncoder["xMax"]     = adet.GetNColumns()-1;
      statusEncoder["yMax"]     = adet.GetNRows()-1;
      statusEncoder.setCellID(status);
      ShortVec statusVec(nPixel, 0);
      status->setADCValues(statusVec);
      _statusCollectionVec->push_back(status);

      TrackerDataImpl * pedestal    = new TrackerDataImpl;
      CellIDEncoder<TrackerDataImpl>  pedestalEncoder(DEPFET::MATRIXDEFAULTENCODING, _pedestalCollectionVec);
      pedestalEncoder["sensorID"] = adet.GetDAQID();
      pedestalEncoder["xMin"]     = 0;
      pedestalEncoder["yMin"]     = 0;
      pedestalEncoder["xMax"]     = adet.GetNColumns()-1;
      pedestalEncoder["yMax"]     = adet.GetNRows()-1;
      pedestalEncoder.setCellID(pedestal);
      FloatVec pedestalVec(nPixel, _initPedestal[iDetector]);
      pedestal->setChargeValues(pedestalVec);
      _pedestalCollectionVec->push_back(pedestal);

      TrackerDataImpl * noise    = new TrackerDataImpl;
      CellIDEncoder<TrackerDataImpl>  noiseEncoder(DEPFET::MATRIXDEFAULTENCODING, _noiseCollectionVec);
      noiseEncoder["sensorID"] =  adet.GetDAQID();
      noiseEncoder["xMin"]     = 0;
      noiseEncoder["yMin"]     = 0;
      noiseEncoder["xMax"]     = adet.GetNColumns()-1;
      noiseEncoder["yMax"]     = adet.GetNRows()-1;
      noiseEncoder.setCellID(noise);
      FloatVec noiseVec(nPixel, _initNoise[iDetector]);
      noise->setChargeValues(noiseVec);
      _noiseCollectionVec->push_back(noise);
    }

    _isFirstEvent = false;
  }

  evt->addCollection(_pedestalCollectionVec, _pedestalCollectionName);
  evt->takeCollection(_pedestalCollectionName);
  evt->addCollection(_noiseCollectionVec, _noiseCollectionName);
  evt->takeCollection(_noiseCollectionName);
  evt->addCollection(_statusCollectionVec, _statusCollectionName);
  evt->takeCollection(_statusCollectionName); 
  
}


//
// Method called after each event to check the data processed
//
void SimpleNoiseDBCreator::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void SimpleNoiseDBCreator::end()
{
  
  delete _pedestalCollectionVec;
  delete _noiseCollectionVec;
  delete _statusCollectionVec;
    
  // Print message
  streamlog_out(MESSAGE2) << std::endl
                          << "Processor succesfully finished!"
                          << std::endl;

}


} // Namespace

