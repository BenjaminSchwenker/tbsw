// PixelChargeCalibrator implementation file
// 
// Author: Helge Christoph Beck, University of Göttingen 
// <mailto:helge-christoph.beck@phys.uni-goettingen.de>

#include "PixelChargeCalibrator.h"

// Include TBTools 
#include "DEPFET.h" 
#include "TBDetector.h"

// Include ROOT classes
#include <TFile.h>

#include <iomanip>
using namespace std::string_literals;
// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;

namespace depfet {

//
// Instantiate this object
//
PixelChargeCalibrator aPixelChargeCalibrator ;

//
// Constructor
//
PixelChargeCalibrator::PixelChargeCalibrator() : Processor("PixelChargeCalibrator"),_inputDecodeHelper(""),
    _calibOutputEncoderHelper(DEPFET::ZSDATADEFAULTENCODING) 
{

// Processor description
   _description = "PixelChargeCalibrator: Applying calibration function pixel specific to zs pixel data for one specific detector" ;

//   
// First of all, we need to register the input/output collections
   
   registerInputCollection (LCIO::TRACKERDATA, "SparseDataCollectionName",
                            "Name of input sparsified pixel data collection",
                            _sparseDataCollectionName, string("zsdata"));
    
   registerOutputCollection (LCIO::TRACKERDATA, "CalibratedCollectionName",
                            "Name of the output charge calibrated collection",
                            _calibratedCollectionName, string("zscalib"));
    
   registerProcessorParameter( "SparseZSCut","Threshold for zero suppression",
                               _sparseZSCut, static_cast<float > (0));

   registerProcessorParameter("NoiseDBFileName",
                               "This is the name of the ROOT file with the status mask (add .root)",
                               _noiseDBFileName, static_cast< string > ( "NoiseDB.root" ) ); 

   registerProcessorParameter("CalibrationDBFileName", 
		   	       "This is the name of the ROOT file with the calibration function and the parameters per pixel for the function (add .root)", 
			       _calibDBFileName, static_cast< string > ( "CalibDB.root" ) );
  
   registerProcessorParameter ("DUTPlane",
                              "Plane number of DUT along the beam line which gets calibrated",
                              _idut,  static_cast < int > (3));

   registerProcessorParameter("CalibFuncName",
		   	      "This is the name of the ROOT TF1 function used to calibrate the DUT charge response",
			      _calibFuncName, static_cast< string > ( "fCalibFunc" ) ); 

   registerProcessorParameter("CalibParaBaseName",
		   	      "This is the base name of the ROOT TH2F's that contain each the pixel by pixel values for one parameter of the calibration function, full name: base + '_' + parnumber",
			      _calibParaBaseName, static_cast< string > ( "para" ) ); 

}

//
// Method called at the beginning of data processing
//
void PixelChargeCalibrator::init() {
   
   // Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
      
   // Open clusterDB file 
   TFile * noiseDBFile = new TFile(_noiseDBFileName.c_str(), "READ");
    
   int sensorID = TBDetector::Get(_idut).GetSensorID();  
   string histoName = "hDB_sensor"+to_string(sensorID) + "_mask";
   if ( (TH2F *) noiseDBFile->Get(histoName.c_str()) != nullptr) {
     _DB_Map_Mask[sensorID] = (TH2F *) noiseDBFile->Get(histoName.c_str());  
     _DB_Map_Mask[sensorID]->SetDirectory(0);
   }
     
   // Close root  file
   noiseDBFile->Close();
   delete noiseDBFile;
   
   // Open calibrationDB file
   TFile * calibDBFile = new TFile(_calibDBFileName.c_str(), "READ");

   if ((TF1 *) calibDBFile->Get(_calibFuncName.c_str()) != nullptr){ 
     _calibFunc = (TF1 *) calibDBFile->Get(_calibFuncName.c_str()); 
     _nparFunc = _calibFunc->GetNpar();
   }

   for(int par=0;par<_nparFunc;par++){
     string histoNamePar = _calibParaBaseName + "_" + to_string(par);
     if ((TH2F *) calibDBFile->Get(histoNamePar.c_str()) != nullptr){
       _DB_Map_CalibPar[par] = (TH2F *) calibDBFile->Get(histoNamePar.c_str());
       _DB_Map_CalibPar[par]->SetDirectory(0);
     }
   }

   // Close calibrationDB file
   calibDBFile->Close();
   delete calibDBFile;
   
   // Print set parameters
   printProcessorParams();
   
   // CPU time start
   _timeCPU = clock()/1000;
   
}

//
// Method called for each run
//
void PixelChargeCalibrator::processRunHeader(LCRunHeader * run)
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
void PixelChargeCalibrator::processEvent(LCEvent * evt)
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
     
     // Output collection containing calibrated data // TRACKERPULSE?, TRACKERDATA from unpacker
     LCCollectionVec * calibCollection = new LCCollectionVec(LCIO::TRACKERDATA) ;
     
     // Calibrate hits pixel by pixel  
     calibrate( evt , calibCollection  ); 
     
     // Add clusterCollection to event
     evt->addCollection( calibCollection, _calibratedCollectionName );
             	    
   } catch(DataNotAvailableException &e){}  
   _nEvt ++ ;
}


//
// Method called after each event to check the data processed
//
void PixelChargeCalibrator::check( LCEvent * )
{
}

//
// Method called after all data processing
//
void PixelChargeCalibrator::end()
{
   
   // CPU time end
   _timeCPU = clock()/1000 - _timeCPU;
   
   // Print message
   streamlog_out(MESSAGE) << std::endl
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
void PixelChargeCalibrator::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "PixelChargeCalibrator Development Version, be carefull!!"
			    << std::endl
			    << "Calibration function name: " << _calibFuncName
			    << std::endl
			    << "with npars: " << _nparFunc
                            << " "
                            << std::endl  << std::endl;   


}

// Called by the processEvent()

void PixelChargeCalibrator::calibrate( LCEvent * evt , LCCollectionVec * calibCollection  ) 
{  
    
  // Open zero suppressed pixel data  
  LCCollectionVec * Pix_collection = dynamic_cast < LCCollectionVec * > (evt->getCollection(_sparseDataCollectionName)); 
  // Helper class for decoding pixel data 
  CellIDDecoder<TrackerDataImpl> PixelID( Pix_collection,&_inputDecodeHelper );
  
  CellIDEncoder<TrackerDataImpl> calibEncoder(DEPFET::ZSDATADEFAULTENCODING, calibCollection,&_calibOutputEncoderHelper ); // DEPFET param?

  // loop over pixel detectors to select the right one from collection
  for (unsigned int iDet = 0; iDet < Pix_collection->size(); iDet++) { 
     
    // Get zs pixels from next pixel detector   
    TrackerDataImpl * pixModule = dynamic_cast<TrackerDataImpl* > ( Pix_collection->getElementAt(iDet) );

    // Sensor ID for pixel detector
    int sensorID = PixelID( pixModule ) ["sensorID"s];
    
    // Read geometry info for sensor 
    int ipl = TBDetector::GetInstance().GetPlaneNumber(sensorID);
    if (ipl != _idut) // select the wanted DUT
      //delete pixModule; ?
      continue; 

    const Det& adet = TBDetector::Get(ipl);

    // Get min channel numbers
    int minUCell = adet.GetMinUCell();
    int minVCell = adet.GetMinVCell();    
    // Get max channel numbers 
    int maxUCell = adet.GetMaxUCell();   
    int maxVCell = adet.GetMaxVCell(); 
    
    // List of firing pixels. Each pixel has a iU, iV and charge 
    FloatVec pixVector = pixModule->getChargeValues();
    int npixels = pixVector.size()/3; 
    
    // Prepare a TrackerData to store the calibrated data
    TrackerDataImpl* calibData = new TrackerDataImpl ; 
     
    // Loop over zspixels and calibrate each
    for (int iPix = 0; iPix < npixels; iPix++) 
    {   
      
      int iU = static_cast<int> (pixVector[iPix * 3]);
      int iV = static_cast<int> (pixVector[iPix * 3 + 1]);
      float charge =  pixVector[iPix * 3 + 2];     
       
      // Try to get status code for pixel 
      float status = 0; 
      if ( _DB_Map_Mask.find(sensorID) != _DB_Map_Mask.end() ) {
        // Here we use the same numbering convention as in the HotPixelKiller processor to map iu, iv to a bin in the mask. 
        status = _DB_Map_Mask[sensorID]->GetBinContent(iU-minUCell+1, iV-minVCell+1); 
      }
      
      // Print detailed pixel summary, for testing/debugging only !!! 
      streamlog_out(MESSAGE1) << "Pixel Nr. " << iPix << " on sensor " << sensorID  
                              << std::endl;  
      streamlog_out(MESSAGE1) << "   iU:" << iU << ", iV:" << iV
                              << ", charge:" << charge
                              << std::endl;
      
      // If pixel signal below threshold, skip it in clusterization
      if ( charge < _sparseZSCut ) {
        streamlog_out(MESSAGE2) << "  Signal below ZS cut. Skipping it." << std::endl; 
        continue;
      }
       
      // If a pixel is out of range, skip it in clusterization
      if ( iU < minUCell || iU > maxUCell || iV < minVCell || iV > maxVCell ) {
        streamlog_out(MESSAGE2) << "  Invalid pixel address found. Skipping it." << std::endl; 
        continue;
      }
      
      // Skip masked pixels 
      if ( status  != 0 ) {
        streamlog_out(MESSAGE2) << "  Bad pixel found. Skipping it." << std::endl; 
        continue;
      }
      
      // Use calibration function and the parameters for the pixel 
      for (int par = 0; par < _nparFunc; par++){
        _calibFunc->SetParameter(par, _DB_Map_CalibPar[par]->GetBinContent(iU-minUCell+1, iV-minVCell+1));
      }
      float calibCharge = _calibFunc->Eval(charge);

      calibData->chargeValues().push_back( iU );
      calibData->chargeValues().push_back( iV );
      calibData->chargeValues().push_back( calibCharge );
      
      streamlog_out(MESSAGE2) << " On sensor " << sensorID << " having DAC unit charge " << charge << " calibrated to " << calibCharge << std::endl;
    }

    static auto idx_orig_sensorID=calibEncoder.index("sensorID"s); //find the address ONCE.
    static auto idx_orig_sparsePixelType=calibEncoder.index("sparsePixelType"s);

    calibEncoder[idx_orig_sensorID] = sensorID;
    calibEncoder[idx_orig_sparsePixelType] = static_cast<int> (kSimpleSparsePixel);
    calibEncoder.setCellID( calibData );
    calibCollection->push_back( calibData );
  
  } // Detector selection loop

}

} // Namespace

