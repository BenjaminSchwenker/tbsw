#include "DEPFETPedestalNoiseProcessor.h"

// Include DEPFETTrackTools 
#include "DEPFET.h" 
#include "Utilities.h"
#include "MatrixDecoder.h"

// Include basic C
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

// Include LCIO classes
#include <lcio.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCEventImpl.h>
#include <UTIL/LCTime.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerRawDataImpl.h>


// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;


namespace depfet {

//
// Instantiate this object
//
DEPFETPedestalNoiseProcessor aDEPFETPedestalNoiseProcessor ;

//
// Destructor - clean up root histos
// 
DEPFETPedestalNoiseProcessor::~DEPFETPedestalNoiseProcessor() 
{
}

//
// Constructor
//
DEPFETPedestalNoiseProcessor::DEPFETPedestalNoiseProcessor() : Processor("DEPFETPedestalNoiseProcessor")
{

// Processor description
   _description = "DEPFETPedestalNoiseProcessor: Preprocessing of DEPFET DCDB (or CURO) and TAKI raw data";


//   
// First of all, we need to register the input/output collections
   
   registerInputCollection (LCIO::TRACKERRAWDATA, "RawDataCollectionName",
                            "Name of input raw data collection",
                            _rawDataCollectionName, string("rawdata_dep"));
   
   registerOutputCollection (LCIO::TRACKERDATA, "DataCollectionName",
                            "Name of the calibrated data collection",
                            _calibratedDataCollectionName, string("data_dep"));

   registerOutputCollection (LCIO::TRACKERDATA, "ZSDataCollectionName",
                            "Name of the zero suppressed data collection",
                            _zsDataCollectionName, string("zsdata_dep"));
   
   registerProcessorParameter("ConfigBadPixels",
                              "This file contains a list of bad pixels",
                              _blockPixelFileName, string("cfg/h4.1.04_blockpixel.cfg"));
    
   registerProcessorParameter("OutputPedeDB",
                              "Name of output pedestal DB",
                              _outputNoiseDB, string("PedeDB-DEP-TB2011.slcio"));
   
   registerProcessorParameter ("PedestalAlgorithm", 
                               "Select algorithm for pedestal and common mode correction:\n"
                               "DHP: DHP data processing\n"
                               "TAKI: TAKI data processing",
                               _pedestalAlgorithm, string( "DHP" ) );
   
   registerProcessorParameter ("ApplyCommonMode",
                              "Use common mode correction in preprocessing",
                              _useCMC,  static_cast < bool > (true));
   
   registerProcessorParameter ("nPedestalEvents",
                               "Number of events between update of pedestal/noise data",
                               _nPedeEvents, static_cast < int >(200));
   
   registerProcessorParameter ("nStatusEvents",
                               "Number of events between update of status data",
                               _nStatusEvents, static_cast < int >(10000));
    
   registerProcessorParameter ("NFoldCM",
                               "N-Fold common mode for N-Fold matrix readout (2or4)",
                               _nFold, static_cast < int >(4));     
   
   registerProcessorParameter ("HitThreshold",
                               "Signal threshold for firing pixels in zero suppression",
                               _hitThresholdZS, static_cast < float >(3));
   
   registerProcessorParameter( "UseSNRCut","Iff true, cut on signal to noise ratio",
                               _useSNRCut, static_cast<bool > (false) );
   
   registerProcessorParameter ("MaxFiringFrequency",
                              "Maximum mean firing frequency for good pixels ",
                              _maxFiringFreq, static_cast < float >(0.5));
   
   registerProcessorParameter ("MaskUpperNoiseCut",
                              "Upper noise threshold for bad pixel identification",
                              _pixelMaskUpperNoise, static_cast < float >(0.8));
   
   registerProcessorParameter ("MaskLowerNoiseCut",
                              "Lower noise threshold for bad pixel identification",
                              _pixelMaskLowerNoise, static_cast < float >(0.2));
   
   registerProcessorParameter ("MaskUpperPedestal",
                              "Upper pedestal threshold for bad pixel identification",
                              _pixelMaskUpperPede, static_cast< float > ( 80 ) );
   
   registerProcessorParameter ("MaskLowerPedestal",
                              "Lower pedestal threshold for bad pixel identification",
                              _pixelMaskLowerPede, static_cast< float > ( -100 ) );
   
   registerProcessorParameter("RootFileName", "Output root file name",
                               _rootFileName, string("PedeDB-DEP-TB2011.root"));
   
   IntVec rowMonitorExample( 1 , 3);
   registerProcessorParameter( "RowMonitor", "Monitor pixels in selected rows",
                               _rowMonitor, rowMonitorExample );
   
   IntVec colMonitorExample( 1 , 15);
   registerProcessorParameter( "ColMonitor", "Monitor pixels in selected columns",
                               _colMonitor, colMonitorExample );
   
   
   // ROOT pointer stuff
   
   _rootFile = NULL;
   _rootPedeTree = NULL; 
   _rootPixelTree = NULL;
   _rootCommonModeTree = NULL;
   _rootEventTree = NULL;  
}

//
// Method called at the beginning of data processing
//
void DEPFETPedestalNoiseProcessor::init() {
   
   // Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
   _timeCPU = clock()/1000;
              
   // Print set parameters
   printProcessorParams();
   
   // No initialization of pedestal alogorithm 
   _isAlgoInit = false; 
   
   // Hit threshold in SNR for pre-processing
   _hitThresholdPre = 4.0; 
   
   // Sort row numbers to improve search speed
   sort (_rowMonitor.begin(), _rowMonitor.end());
      
   // Transform the algorithm string to small letters
   transform(_pedestalAlgorithm.begin(), _pedestalAlgorithm.end(), _pedestalAlgorithm.begin(), ::tolower);
     
   // ROOT Output 
   _rootFile = new TFile(_rootFileName.c_str(),"recreate");
    
   // Note: this tree contains snapshots for complete module
   _rootPedeTree = new TTree("Pedestal","Pedestal info");
   _rootPedeTree->Branch("iEvt"            ,&_rootEventNumber    ,"iEvt/I");
   _rootPedeTree->Branch("det"             ,&_rootDetectorID     ,"det/I");
   _rootPedeTree->Branch("cycle"           ,&_rootCycle          ,"cycle/I");
   _rootPedeTree->Branch("px_x"            ,&_rootCol            ,"px_x/I");
   _rootPedeTree->Branch("px_y"            ,&_rootRow            ,"px_y/I");
   _rootPedeTree->Branch("pedestal"        ,&_rootPedestal       ,"pedestal/D");
   _rootPedeTree->Branch("noise"           ,&_rootNoise          ,"noise/D");
   _rootPedeTree->Branch("status"          ,&_rootStatus         ,"status/I");
   _rootPedeTree->Branch("cmc"             ,&_rootCMC            ,"cmc/D");
   _rootPedeTree->Branch("adc"             ,&_rootADC            ,"adc/D");
   _rootPedeTree->Branch("data"            ,&_rootDATA           ,"data/D");
   _rootPedeTree->Branch("hitFreq"         ,&_rootHitFrequency   ,"hitFreq/D");
   
   // Note: this tree contains only monitored pixels; sampled every event
   _rootPixelTree = new TTree("Pixel","Pixel monitor info");
   _rootPixelTree->Branch("iEvt"            ,&_rootEventNumber    ,"iEvt/I");
   _rootPixelTree->Branch("det"             ,&_rootDetectorID     ,"det/I");
   _rootPixelTree->Branch("px_x"            ,&_rootCol            ,"px_x/I");
   _rootPixelTree->Branch("px_y"            ,&_rootRow            ,"px_y/I");
   _rootPixelTree->Branch("pedestal"        ,&_rootPedestal       ,"pedestal/D");
   _rootPixelTree->Branch("noise"           ,&_rootNoise          ,"noise/D");
   _rootPixelTree->Branch("status"          ,&_rootStatus         ,"status/I");
   _rootPixelTree->Branch("cmc"             ,&_rootCMC            ,"cmc/D");
   _rootPixelTree->Branch("adc"             ,&_rootADC            ,"adc/D");  
   _rootPixelTree->Branch("data"            ,&_rootDATA           ,"data/D");
   
   _rootCommonModeTree = new TTree("CommonMode","CommonMode info");
   _rootCommonModeTree->Branch("iEvt"            ,&_rootEventNumber    ,"iEvt/I");
   _rootCommonModeTree->Branch("det"             ,&_rootDetectorID     ,"det/I");
   _rootCommonModeTree->Branch("gate"            ,&_rootGate           ,"gate/I");
   _rootCommonModeTree->Branch("cmc"             ,&_rootCMC            ,"cmc/D");
   
   _rootEventTree = new TTree("Event","Event info");
   _rootEventTree->Branch("iEvt"            ,&_rootEventNumber    ,"iEvt/I");
   _rootEventTree->Branch("det"             ,&_rootDetectorID     ,"det/I");
   _rootEventTree->Branch("nhit"            ,&_rootNFiring        ,"nhit/I");
   
}

//
// Method called for each run
//
void DEPFETPedestalNoiseProcessor::processRunHeader(LCRunHeader * run)
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
void DEPFETPedestalNoiseProcessor::processEvent(LCEvent * evt)
{
   
   // Print event number
   if ( evt->getEventNumber()%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                 << (evt->getEventNumber())
                                                                 << std::endl << std::endl;
   
   if ( !_isAlgoInit ) {
        
      // Initialize pedestal, commom mode and quality algorithms
      _isAlgoInit = initializeAlgorithms(evt); 
       
      streamlog_out(MESSAGE3) << "Initialization successfull ..." << endl; 
      
      for (size_t iDet = 0; iDet < _sensorIDVec.size(); iDet++) 
      {  
        streamlog_out(MESSAGE3) << "SensorID " << _sensorIDVec[iDet]  << ": " << endl 
                                << "uMin: " << _minX[iDet] << " uMax: " <<  _maxX[iDet] << endl
                                << "vMin: " << _minY[iDet] << " vMax: " <<  _maxY[iDet] << endl << endl; 
      }
   }
   
   // Choose the pedestal/noise algorithm
   if ( _pedestalAlgorithm == "dhp" ) {
     streamlog_out ( MESSAGE1 ) << "Using DHP methods ..." << endl;
     processRawDataDHP(evt);  	
   } else if ( _pedestalAlgorithm == "taki"  ) {
     streamlog_out ( MESSAGE1 ) << "Using TAKI methods ..." << endl;
     processRawDataTAKI(evt);  
   } else if ( _pedestalAlgorithm == "dhp2" ) {
     streamlog_out ( MESSAGE1 ) << "Using DHP2 methods ..." << endl;
     processRawDataDHP2(evt);  
   } else {
     streamlog_out ( ERROR4 )  << "Choosen correction method not available. Please change selection."<< endl;
     exit(-1); 
   }
   
   // Store corrected frames in LCIO
   streamlog_out ( MESSAGE1 ) << "Calibrate event ..." << endl;
   calibrateEvent(evt);
   
   // Store 0-suppressed pixels in LCIO
   streamlog_out ( MESSAGE1 ) << "Sparsify event ..." << endl;
   sparsifyEvent(evt);
   
   // Fill pixel data monitor ntuple for DQM  
   streamlog_out ( MESSAGE1 ) << "Monitor pixel data ..." << endl;
   fillPixelTuple(evt);
   
   // Create map of bad pixels and monitor pedestals    
   if ( ( _nEvt>0 ) && ( _nEvt%_nStatusEvents == 0 ) ) {
     streamlog_out(MESSAGE3) << "Masking bad pixels ..." << endl;    
     maskBadPixel(evt); 
   }
   
   // Increment internal event number
   ++_nEvt;
   
}


//
// Method called after each event to check the data processed
//
void DEPFETPedestalNoiseProcessor::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void DEPFETPedestalNoiseProcessor::end()
{
  
  // CPU time end
  _timeCPU = clock()/1000 - _timeCPU;
  
  streamlog_out ( MESSAGE4 ) << "Writing the output condition file" << endl;
  
  try {
     
    // Create lcio output noiseDB 
    LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
    lcWriter->open( _outputNoiseDB, LCIO::WRITE_NEW );
    
    // Write an almost empty run header
    LCRunHeaderImpl * lcHeader  = new LCRunHeaderImpl();
    lcHeader->setRunNumber( 0 );
        
    lcWriter->writeRunHeader(lcHeader);
    delete lcHeader;
 
    LCEventImpl * event = new LCEventImpl;
    event->setRunNumber( _nRun );
    event->setEventNumber( 0 );
       
    LCTime * now = new LCTime;
    event->setTimeStamp( now->timeStamp() );
    delete now;
    
    LCCollectionVec * pedestalCollection = new LCCollectionVec(LCIO::TRACKERDATA);
    LCCollectionVec * noiseCollection    = new LCCollectionVec(LCIO::TRACKERDATA);
    LCCollectionVec * statusCollection   = new LCCollectionVec(LCIO::TRACKERRAWDATA);
    
    
    for ( size_t iDetector = 0; iDetector < _noOfDetector; iDetector++) {
      
      TrackerDataImpl    * pedestalMatrix = new TrackerDataImpl;
      TrackerDataImpl    * noiseMatrix    = new TrackerDataImpl;
      TrackerRawDataImpl * statusMatrix   = new TrackerRawDataImpl;
     
      
      CellIDEncoder<TrackerDataImpl>    idPedestalEncoder(DEPFET::MATRIXDEFAULTENCODING, pedestalCollection);
      CellIDEncoder<TrackerDataImpl>    idNoiseEncoder(DEPFET::MATRIXDEFAULTENCODING, noiseCollection);
      CellIDEncoder<TrackerRawDataImpl> idStatusEncoder(DEPFET::MATRIXDEFAULTENCODING, statusCollection);
     
      
      idPedestalEncoder["sensorID"] = _sensorIDVec.at( iDetector );
      idNoiseEncoder["sensorID"]    = _sensorIDVec.at( iDetector );
      idStatusEncoder["sensorID"]   = _sensorIDVec.at( iDetector );
      idPedestalEncoder["uMax"]     = _maxX[iDetector];
      idNoiseEncoder["uMax"]        = _maxX[iDetector];
      idStatusEncoder["uMax"]       = _maxX[iDetector];
      idPedestalEncoder["vMax"]     = _maxY[iDetector];
      idNoiseEncoder["vMax"]        = _maxY[iDetector];
      idStatusEncoder["vMax"]       = _maxY[iDetector];
 
      
      idPedestalEncoder.setCellID(pedestalMatrix);
      idNoiseEncoder.setCellID(noiseMatrix);
      idStatusEncoder.setCellID(statusMatrix);
      
      
      pedestalMatrix->setChargeValues(_pedestal[iDetector]);
      noiseMatrix->setChargeValues(_noise[iDetector]);
      statusMatrix->setADCValues(_status[iDetector]);
      
      pedestalCollection->push_back(pedestalMatrix);
      noiseCollection->push_back(noiseMatrix);
      statusCollection->push_back(statusMatrix);
      
       
    }   
    
    event->addCollection(pedestalCollection, "pedestal");
    event->addCollection(noiseCollection, "noise");
    event->addCollection(statusCollection, "status");
    
    lcWriter->writeEvent(event);
    delete event;
    lcWriter->close();
    
  }
  catch ( IOException& e ) 
  {
    streamlog_out ( ERROR4 ) << e.what() << std::endl
                             << "Problem creating noiseDB. Sorry for quitting. " << std::endl;
    exit(-1);
  }
  
  // Create pedestal/noise/status maps in root 
  
  for ( size_t iDetector = 0; iDetector < _noOfDetector ; iDetector++) {
            
    int noOfXPixels = _maxX[iDetector] + 1; 
    int noOfYPixels = _maxY[iDetector] + 1;
    int sensorID =  _sensorIDVec[iDetector];   
     
    // Write status/noise/pedestal maps
    string histoName, histoTitle; 
    
    histoName = "noiseMap_ModID" + to_string(sensorID) ; 
    histoTitle = "noiseMap Mod" + to_string(sensorID);
    TH2D noiseMap(histoName.c_str(),histoTitle.c_str(),noOfXPixels,-0.5,noOfXPixels-0.5,noOfYPixels,-0.5,noOfYPixels-0.5);
    noiseMap.SetXTitle("X Axis"); 
    noiseMap.SetYTitle("Y Axis");
      
    histoName = "pedestalMap_ModID" + to_string(sensorID) ; 
    histoTitle = "pedestalMap Mod" + to_string(sensorID);
    TH2D  pedestalMap(histoName.c_str(),histoTitle.c_str(),noOfXPixels,-0.5,noOfXPixels-0.5,noOfYPixels,-0.5,noOfYPixels-0.5);
    pedestalMap.SetXTitle("X Axis"); 
    pedestalMap.SetYTitle("Y Axis");
         
    histoName = "statusMap_ModID" + to_string(sensorID ) ; 
    histoTitle = "statusMap Mod" + to_string(sensorID);
    TH2D  statusMap(histoName.c_str(),histoTitle.c_str(),noOfXPixels,-0.5,noOfXPixels-0.5,noOfYPixels,-0.5,noOfYPixels-0.5);
    statusMap.SetXTitle("X Axis"); 
    statusMap.SetYTitle("Y Axis");
    
    
      
    MatrixDecoder matrixDecoder(noOfXPixels, noOfYPixels); 
      
    for (int yPixel = 0; yPixel < noOfYPixels; yPixel++) {
      for (int xPixel = 0; xPixel < noOfXPixels; xPixel++) {
   
        int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);
        noiseMap.Fill(xPixel,yPixel, _noise[iDetector][iPixel] );
        pedestalMap.Fill(xPixel,yPixel, _pedestal[iDetector][iPixel]  );
        statusMap.Fill(xPixel,yPixel, _status[iDetector][iPixel]  );     
         
      } 
    }
        
    _rootFile->cd(""); 
    noiseMap.Write(); 
    statusMap.Write();
    pedestalMap.Write(); 
     
  } 
  
  _rootPedeTree->Write(); 
  _rootPixelTree->Write();
  _rootCommonModeTree->Write();
  _rootEventTree->Write();   
  
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
void DEPFETPedestalNoiseProcessor::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "DEPFETPedestalNoiseProcessor Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   

}


void DEPFETPedestalNoiseProcessor::finalizePedestalsDHP( ) {
  
  // A pede loop on events is over. So we need to move temporary pedestal
  // noise estimates to buffers for next pede loop. 
  for ( size_t iDetector = 0; iDetector < _noOfDetector; iDetector++) {
    for ( size_t iPixel = 0; iPixel <  _noise[iDetector].size(); iPixel++ ) { 
      
      // Calculate online noise estimate 
      _tmpNoise[iDetector][iPixel] = sqrt(_ttmpNoise[iDetector][iPixel]/(_tmpEntries[iDetector][iPixel]-1));
      
      // Calculate online pedestal estimate
      _tmpPede[iDetector][iPixel]  = _ttmpPede[iDetector][iPixel];
      
      // Reset variables for next loop
      _tmpEntries[iDetector][iPixel] = 0;       
      _ttmpNoise[iDetector][iPixel] = 0; 
      _ttmpPede[iDetector][iPixel] = 0; 
      
    }
  }
  
  // The calculation of pedestal/noise is finished at odd pede loops.
  // Upload pedestal/noise at odd cycles  
  if ( _iPedeLoop%2 == 1 ) { 
     streamlog_out(MESSAGE2) << "Upload new calibration constants ..." << endl;
    _pedestal = _tmpPede;
    _noise    = _tmpNoise; 
  }
  
  // Check if data processing can be started. This needs stable pedestal 
  // noise estimates.   
  if ( _iPedeLoop == 3 ) {
    streamlog_out(MESSAGE3) << "Start calibration of rawdata ..." << endl;
    _isProcessData = true;
  }
  
  // Increment the loop counter
  streamlog_out(MESSAGE2) << "Pedestal loop: " << _iPedeLoop << endl;
  ++_iPedeLoop;
    
}


//
// Method doing all raw data processing in DHP style
//
void DEPFETPedestalNoiseProcessor::processRawDataDHP(LCEvent * evt) {
  
  if ( ( _nEvt > 0 ) && (_nEvt % _nPedeEvents == 0) ) finalizePedestalsDHP( );
  
  // If true, us all pixels in common 
  // mode correction 
  bool no_cmmask = true;
  
  try {
      
    LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionName ));
    
    for ( size_t iDetector = 0; iDetector < collectionVec->size() ; iDetector++) {
        
      // Get the TrackerRawData object from the collection for this plane
      TrackerRawDataImpl *trackerRawData = dynamic_cast < TrackerRawDataImpl * >(collectionVec->getElementAt (iDetector));
      ShortVec adcValues = trackerRawData->getADCValues ();
        
      // Get number of pixels of current module 
      CellIDDecoder< TrackerRawDataImpl > rawDataDecoder( collectionVec );
      MatrixDecoder matrixDecoder( rawDataDecoder, trackerRawData);  
      int noOfXPixels = rawDataDecoder( trackerRawData ) ["uMax"] + 1; 
      int noOfYPixels = rawDataDecoder( trackerRawData ) ["vMax"] + 1;
      int sensorID = rawDataDecoder( trackerRawData ) ["sensorID"]; 
      
      // Veto for data frame that should not be corrected and saved to lcio 
      _isFrameValid[iDetector] = true;
      
      // No calibrated data in warm up     
      if  ( _isWarmUp ) _isFrameValid[iDetector] = false;   
      
      // Loop over all DEPFET gates - pixels read out at the same time 
      // Note: Do all preprocessing of signals gate wise 
      for (int iGate = 0; iGate < noOfYPixels/_nFold ; iGate ++) {
              
        // Calibrate DEPEFT pixel signals 
        if ( _isProcessData ) {
        
          // Common mode estimate for data processing         
          float CommonMode = 0.;
          float pixelSum = 0; 
          int nGoodPixels = 0; 
          
          // 1st loop over all drains: 1st pass common mode correction
          for (int iDrain = 0; iDrain < _nFold*noOfXPixels; iDrain++) {
            
            // Decode pixel row/col
            int xPixel = iDrain/_nFold;
            int yPixel = ( _nFold*iGate ) + ( iDrain%_nFold );
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);
            
            // Subtract pedestal from rawdata     
            float I = adcValues[iPixel] - _pedestal[iDetector][iPixel];    
            
            // Only use good pixel in common mode 
            bool isGood = ( _status[iDetector][iPixel] == 0 );       

            // Only use good pixels for common mode 
            if ( no_cmmask ||  isGood ) { 
              nGoodPixels++; 
              pixelSum += I;
            }
                         
          } 
          
          //std::cout << "nGoodPixels (1loop): " << nGoodPixels << " frame " << iDetector << std::endl;  
          
          // Common mode estimate - status filtered 
          CommonMode = pixelSum / nGoodPixels;
          pixelSum = 0; 
          nGoodPixels = 0;  
            
          if (!_useCMC)  CommonMode = 0; 
               
	  // 2nd loop over all drains: 2nd pass common mode correction  
          for (int iDrain = 0; iDrain < _nFold*noOfXPixels; iDrain++) {
                     
            int xPixel = iDrain/_nFold;
            int yPixel = ( _nFold*iGate ) + ( iDrain%_nFold );      
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);  
            
            // Subtract pedestal and common mode from rawdata  
            float I = adcValues[iPixel] - _pedestal[iDetector][iPixel] -  CommonMode;     
            
            // Firing pixel threshold   
            float threshold = _hitThresholdZS; 
            if ( _useSNRCut ) threshold *= _noise[iDetector][iPixel]; 
            
            // Not use firing pixels in common mode
            bool isHit = std::abs( I ) > threshold;
                        
            // Not use bad pixels in common mode 
            bool isGood = ( _status[iDetector][iPixel] == 0 );              
            
            // Only use good pixels for common mode 
            if ( no_cmmask || (isGood && !isHit)  ) { 
              nGoodPixels++; 
              pixelSum += adcValues[iPixel] - _pedestal[iDetector][iPixel];
            }                 
          } 
          
          //std::cout << "nGoodPixels (2loop): " << nGoodPixels << " frame " << iDetector << std::endl;  
           
          // Common mode estimate - status and hit filtered  
          CommonMode = pixelSum / nGoodPixels ; 
          pixelSum = 0; 
          nGoodPixels = 0;  
           
          if (!_useCMC)  CommonMode = 0; 
           
          // Final loop over all drains
          for (int iDrain = 0; iDrain < _nFold*noOfXPixels; iDrain++) {
          
            int xPixel = iDrain/_nFold;
            int yPixel = ( _nFold*iGate ) + ( iDrain%_nFold );      
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);  
                      
            _commonMode[iDetector][iPixel] = CommonMode;
            _correctedData[iDetector][iPixel] = adcValues[iPixel] - _pedestal[iDetector][iPixel] - CommonMode; 
            
            
            
            // Firing pixel threshold   
            float threshold = _hitThresholdZS; 
            if ( _useSNRCut ) threshold *= _noise[iDetector][iPixel]; 
                
            if ( _correctedData[iDetector][iPixel] > threshold ) {  
              _hitCounter[iDetector][iPixel]++;  
            }
                 
          }
          
          // Time to review the quality of common mode 
          _rootEventNumber = evt->getEventNumber();
          _rootDetectorID = sensorID;
          _rootGate = iGate;
          _rootCMC = CommonMode;
             
          _rootCommonModeTree->Fill(); 
              
        }
        
        // Calculate calibration constants
        
        if ( _iPedeLoop%2 == 0 ) { 
          
          // Initially, we do not know pixel noise or pedestal value. This means, 
          // we have no common mode correction or hit filtering. 
          // It means that noise will be too large (~50%), since no common mode 
          // correction is applied. Pedestals are hit biased.        
          
          for (int iDrain = 0; iDrain < _nFold*noOfXPixels; iDrain++) {
            
            // Decode pixel row/col
            int xPixel = iDrain/_nFold;
            int yPixel = ( _nFold*iGate ) + ( iDrain%_nFold );
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);  
             
            // Online calculation of pedestal and nosie as mean and standard deviation 
            // of raw adc codes.
            
            // For reference, see http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance         
            _tmpEntries[iDetector][iPixel]++;
            double delta = adcValues[iPixel] - _ttmpPede[iDetector][iPixel];
            _ttmpPede[iDetector][iPixel] += delta/_tmpEntries[iDetector][iPixel];
            _ttmpNoise[iDetector][iPixel] += delta*(adcValues[iPixel] - _ttmpPede[iDetector][iPixel]);
               
          } 
        } else { 
          
          // Now, we have intermediate(!) estimates for pixel pedestal and noise values from prev. loop.
          // This time we may use hit filtering and common mode correction :)   
            
          float CommonMode = 0;
          float pixelSum = 0; 
          int nGoodPixels = 0;   
          
          // 1st loop over all drains: 1st pass common mode correction
          for (int iDrain = 0; iDrain < _nFold*noOfXPixels; iDrain++) {
          
            // Decode pixel row/col
            int xPixel = iDrain/_nFold;
            int yPixel = ( _nFold*iGate ) + ( iDrain%_nFold );
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);  
                 
            // Only use good pixel in common mode 
            bool isGood = ( _status[iDetector][iPixel] == 0 );      
            
           
            if ( no_cmmask || isGood ) { 
              nGoodPixels++; 
              pixelSum += adcValues[iPixel] - _tmpPede[iDetector][iPixel];
            }                 
          } 
          
          // Common mode estimate - status filtered 
          CommonMode = pixelSum / nGoodPixels;
          pixelSum = 0; 
          nGoodPixels = 0;  
          
          if (!_useCMC)  CommonMode = 0; 
           
          // 2nd loop over all drains: 2nd pass common mode correction  
          for (int iDrain = 0; iDrain < _nFold*noOfXPixels; iDrain++) {
            
            int xPixel = iDrain/_nFold;
            int yPixel = ( _nFold*iGate ) + ( iDrain%_nFold );        
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);  
                
            // Filter firing pixels 
            float I = adcValues[iPixel] - _tmpPede[iDetector][iPixel] -  CommonMode;
              
            float threshold = _hitThresholdZS; 
            if ( _useSNRCut ) threshold *= _tmpNoise[iDetector][iPixel]; 
                 
            // Not use firing pixels for common mode
            bool isHit = std::abs( I ) >  threshold;
            
            // Not use bad pixels for common mode
            bool isGood = ( _status[iDetector][iPixel] == 0 );              
             
            if ( no_cmmask || (isGood && !isHit)   ) { 
              nGoodPixels++; 
              pixelSum += adcValues[iPixel] - _tmpPede[iDetector][iPixel];
            }               
          } 
          
          // Common mode estimate - status and hit filtered  
          CommonMode = pixelSum / nGoodPixels ;
          pixelSum = 0; 
          nGoodPixels = 0;  
          
          if (!_useCMC)  CommonMode = 0; 
           
          // Final loop over all drains
          for (int iDrain = 0; iDrain < _nFold*noOfXPixels; iDrain++) {
            
            int xPixel = iDrain/_nFold;
            int yPixel = ( _nFold*iGate ) + ( iDrain%_nFold );        
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);  
            
            // Filter firing pixels 
            float I = adcValues[iPixel] - _tmpPede[iDetector][iPixel] - CommonMode;  
          
            float threshold = _hitThresholdZS; 
            if ( _useSNRCut ) threshold *= _tmpNoise[iDetector][iPixel]; 
                 
            // Not use firing pixels for common mode
            bool isHit = std::abs( I ) >  threshold; 
             
            // Only use good pixels for common mode 
            if ( !isHit ) { 
              
              // Use measurment for pedestal calculations      
              _tmpEntries[iDetector][iPixel]++;
               
              // Pedestal is simple average (but using a hit filter)  
              double delta = adcValues[iPixel] - _ttmpPede[iDetector][iPixel];
              _ttmpPede[iDetector][iPixel] += delta/_tmpEntries[iDetector][iPixel];
              
              // Noise is rms of baseline fluctuations (after common mode correction)
              _ttmpNoise[iDetector][iPixel] += (_tmpEntries[iDetector][iPixel] - 1)*pow( I, 2)
                                              / _tmpEntries[iDetector][iPixel];
              
            } 
          }
        }
             
      } // End gate loop  
             
    } // End loop on detectors
    
  } catch (DataNotAvailableException& e) {
    streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionName << " is not available in the current event" << endl;
  }
  
}

void DEPFETPedestalNoiseProcessor::finalizePedestalsDHP2( ) {
   
  // A pede loop on events is over. So we need to move temporary pedestal
  // noise estimates to buffers for next pede loop. 
  for ( size_t iDetector = 0; iDetector < _noOfDetector; iDetector++) {
    for ( size_t iPixel = 0; iPixel <  _noise[iDetector].size(); iPixel++ ) { 
      
      // Calculate online noise estimate 
      _tmpNoise[iDetector][iPixel] = sqrt(_ttmpNoise[iDetector][iPixel]/(_tmpEntries[iDetector][iPixel]-1));
      
      // Calculate online pedestal estimate
      _tmpPede[iDetector][iPixel]  = _ttmpPede[iDetector][iPixel];
      
      // Reset variables for next loop
      _tmpEntries[iDetector][iPixel] = 0;       
      _ttmpNoise[iDetector][iPixel] = 0; 
      _ttmpPede[iDetector][iPixel] = 0; 
      
    }
  }
  
  // The calculation of pedestal/noise is finished at odd pede loops.
  // Upload pedestal/noise at odd cycles  
  if ( _iPedeLoop%2 == 1 ) { 
     streamlog_out(MESSAGE2) << "Upload new calibration constants ..." << endl;
    _pedestal = _tmpPede;
    _noise    = _tmpNoise; 
  }
  
  // Check if data processing can be started. This needs stable pedestal 
  // noise estimates.   
  if ( _iPedeLoop == 3 ) {
    streamlog_out(MESSAGE3) << "Start calibration of rawdata ..." << endl;
    _isProcessData = true;
  }
  
  // Increment the loop counter
  ++_iPedeLoop;
    
}

//
// Method doing all raw data processing in DHP2 style
//
void DEPFETPedestalNoiseProcessor::processRawDataDHP2(LCEvent * evt) {
  
  if ( ( _nEvt > 0 ) && (_nEvt % _nPedeEvents == 0) ) finalizePedestalsDHP2( ); 
  
  try {
      
    LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionName ));
    
    for ( size_t iDetector = 0; iDetector < collectionVec->size() ; iDetector++) {
        
      // Get the TrackerRawData object from the collection for this plane
      TrackerRawDataImpl *trackerRawData = dynamic_cast < TrackerRawDataImpl * >(collectionVec->getElementAt (iDetector));
      ShortVec adcValues = trackerRawData->getADCValues ();
        
      // Get number of pixels of current module 
      CellIDDecoder< TrackerRawDataImpl > rawDataDecoder( collectionVec );
      MatrixDecoder matrixDecoder( rawDataDecoder, trackerRawData);  
      int noOfXPixels = rawDataDecoder( trackerRawData ) ["uMax"] + 1; 
      int noOfYPixels = rawDataDecoder( trackerRawData ) ["vMax"] + 1;
      int sensorID = rawDataDecoder( trackerRawData ) ["sensorID"]; 
      
      // Veto for data frame that should not be corrected and saved to lcio 
      _isFrameValid[iDetector] = true;
      
      // No calibrated data in warm up     
      if  ( _isWarmUp ) _isFrameValid[iDetector] = false;   
      
      // Loop over all DEPFET gates - pixels read out at the same time 
      // Note: Do all preprocessing of signals gate wise 
      for (int iGate = 0; iGate < noOfYPixels/_nFold ; iGate ++) {
          
        bool GateOK = true;  
        
        // Calibrate DEPEFT pixel signals 
        if ( _isProcessData ) {
          
          // 
          // Common mode loop
              
          float CommonMode = 0.;
          float pixelSum = 0; 
          int nGoodPixels = 0; 
          
          for (int iDrain = 0; iDrain < _nFold*noOfXPixels; iDrain++) {
              
            // Decode pixel row/col
            int xPixel = iDrain/_nFold;
            int yPixel = ( _nFold*iGate ) + ( iDrain%_nFold );
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);
            
            // Mean based common mode
            pixelSum += adcValues[iPixel];    
            ++nGoodPixels;                
          } 
           
          CommonMode = pixelSum / nGoodPixels;
          pixelSum = 0; 
          nGoodPixels = 0;
           
          // 
          // Veto loop 
           
          GateOK = true; 
          int HitPerGate = 0;
          
          for (int iDrain = 0; iDrain < _nFold*noOfXPixels; iDrain++) {
              
            // Decode pixel row/col
            int xPixel = iDrain/_nFold;
            int yPixel = ( _nFold*iGate ) + ( iDrain%_nFold );
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);
            
            // Count number of hits (outliers) per gate 
            double cdata = std::abs( adcValues[iPixel] - _pedestal[iDetector][iPixel] - CommonMode ); 
            if ( cdata > _hitThresholdPre * _noise[iDetector][iPixel] )  {
              HitPerGate++; 
            }        
          }
          
          if ( HitPerGate > 10) {
            streamlog_out(MESSAGE1) << "Veto in event " << evt->getEventNumber() 
                                    << " in gate " << iGate  
                                    << " in frame " << iDetector
                                    << std::endl;
            GateOK = false;   
          }      
               
          // 
          // Final calibration loop  
	     
          for (int iDrain = 0; iDrain < _nFold*noOfXPixels; iDrain++) {
                     
            int xPixel = iDrain/_nFold;
            int yPixel = ( _nFold*iGate ) + ( iDrain%_nFold );      
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);  
                
            _commonMode[iDetector][iPixel] = CommonMode;
            
            if ( GateOK ) { 
              _correctedData[iDetector][iPixel] = adcValues[iPixel] - _pedestal[iDetector][iPixel] - CommonMode;  
              if ( _correctedData[iDetector][iPixel] > _hitThresholdPre * _noise[iDetector][iPixel] )  {
                _hitCounter[iDetector][iPixel]++;   
              }
            } else { 
              _correctedData[iDetector][iPixel] = -9999; // neg. value not in ADC range!!
            }
          } 
          
          // Time to review the quality of common mode 
          _rootEventNumber = evt->getEventNumber();
          _rootDetectorID = sensorID;
          _rootGate = iGate;
          _rootCMC = CommonMode;
          _rootCommonModeTree->Fill(); 
              
        }
        
        // Calculate pedestal + noise constants
        if ( GateOK ) {   // Do not use bad gates here!! 
        if ( _iPedeLoop%2 == 0 ) { 
          
          // Initially, we do not know pixel noise or pedestal value. We have to 
          // to find rough estimates for the pedestals.
          
          float CommonMode = 0;
          float pixelSum = 0; 
          int nGoodPixels = 0;   
                
          for (int iDrain = 0; iDrain < _nFold*noOfXPixels; iDrain++) {
            
            // Decode pixel row/col
            int xPixel = iDrain/_nFold;
            int yPixel = ( _nFold*iGate ) + ( iDrain%_nFold );
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);  
            
            // Mean based common mode
            pixelSum += adcValues[iPixel];    
            ++nGoodPixels; 
          } 
          
          // Common mode estimate 
          CommonMode = pixelSum / nGoodPixels;
          pixelSum = 0; 
          nGoodPixels = 0;
          
          for (int iDrain = 0; iDrain < _nFold*noOfXPixels; iDrain++) {     
              
            // Decode pixel row/col
            int xPixel = iDrain/_nFold;
            int yPixel = ( _nFold*iGate ) + ( iDrain%_nFold );
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);  
            
            // Subtract common mode offset: rest is pedestal + noise + signal  
            float I  = adcValues[iPixel] - CommonMode;
               
            // Online calculation of pedestal and nosie as mean and standard deviation 
            // of common mode corrected signals. 
             
            // For reference, see http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance         
            _tmpEntries[iDetector][iPixel]++;
            double delta = I - _ttmpPede[iDetector][iPixel];
            _ttmpPede[iDetector][iPixel] += delta/_tmpEntries[iDetector][iPixel];
            _ttmpNoise[iDetector][iPixel] += delta*(I - _ttmpPede[iDetector][iPixel]);   
          } 
          
           
        } else { 
          
          // Now, we have intermediate(!) estimates for pixel pedestal and noise values from prev. loop.
          // This time we may use hit filtering :)   
            
          float CommonMode = 0;
          float pixelSum = 0; 
          int nGoodPixels = 0;   
          
          for (int iDrain = 0; iDrain < _nFold*noOfXPixels; iDrain++) {
          
            // Decode pixel row/col
            int xPixel = iDrain/_nFold;
            int yPixel = ( _nFold*iGate ) + ( iDrain%_nFold );
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);  
                 
            // Mean based common mode
            pixelSum += adcValues[iPixel];    
            ++nGoodPixels;     
          } 
          
          // Common mode estimate - status filtered 
          CommonMode = pixelSum / nGoodPixels;
          pixelSum = 0; 
          nGoodPixels = 0;  
           
          for (int iDrain = 0; iDrain < _nFold*noOfXPixels; iDrain++) {
            
            int xPixel = iDrain/_nFold;
            int yPixel = ( _nFold*iGate ) + ( iDrain%_nFold );        
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);  
                    
            // Rough signal estimate using intermediate pedestal 
            float Isig = adcValues[iPixel] - _tmpPede[iDetector][iPixel] - CommonMode;
            
            // Hit filter for final noise + pedestal estimation  
            bool isHit = std::abs( Isig ) > _hitThresholdPre *_tmpNoise[iDetector][iPixel];           
             
            //if ( !isHit ) {   
            if ( true ) {   
               
              // Subtract common mode offset: rest is pedestal + noise + signal  
              float I  = adcValues[iPixel] - CommonMode;
              
              // Online calculation of pedestal and nosie as mean and standard deviation 
              // of common mode corrected signals. 
               
              // For reference, see http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance        
              _tmpEntries[iDetector][iPixel]++;  
              double delta = I - _ttmpPede[iDetector][iPixel];
              _ttmpPede[iDetector][iPixel] += delta/_tmpEntries[iDetector][iPixel];
              _ttmpNoise[iDetector][iPixel] += delta*(I - _ttmpPede[iDetector][iPixel]);   
                             
              /* This calculated noise as rms around old pedestal 
              double Ip = adcValues[iPixel] - _tmpPede[iDetector][iPixel] - CommonMode; 
              _ttmpNoise[iDetector][iPixel] += (_tmpEntries[iDetector][iPixel] - 1)*pow( Ip, 2)
                                              / _tmpEntries[iDetector][iPixel];
              */
            }  
                
          } 
            
        }
        } // GateOK     
      } // End gate loop  
             
    } // End loop on detectors
    
  } catch (DataNotAvailableException& e) {
    streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionName << " is not available in the current event" << endl;
  }
  
}


void DEPFETPedestalNoiseProcessor::finalizePedestalsTAKI( ) {
  
  // A pede loop on events is over. So we need to move temporary pedestal
  // noise estimates to buffers for next pede loop. 
  
  for ( size_t iDetector = 0; iDetector < _noOfDetector; iDetector++) {
    for ( size_t iPixel = 0; iPixel <  _noise[iDetector].size(); iPixel++ ) { 
      
      // Calculate online noise estimate 
      _tmpNoise[iDetector][iPixel] = sqrt( _ttmpNoise[iDetector][iPixel]/(_tmpEntries[iDetector][iPixel]-1) ) ;
      
      // Calculate online pedestal estimate
      _tmpPede[iDetector][iPixel]  = _ttmpPede[iDetector][iPixel];
      
      // Reset variables for next loop
      _tmpEntries[iDetector][iPixel] = 0;       
      _ttmpNoise[iDetector][iPixel] = 0; 
      _ttmpPede[iDetector][iPixel] = 0; 
      
    }
  }
  
  // The calculation of pedestal/noise is finished.
  // Upload pedestal/noise data   
  streamlog_out(MESSAGE2) << "Upload new calibration constants ..." << endl;
  _pedestal = _tmpPede;
  _noise    = _tmpNoise; 
  
  // Check if data processing can be started. This needs stable pedestal 
  // noise estimates.   
  if ( _iPedeLoop == 1 ) {
    streamlog_out(MESSAGE3) << "Start calibration of rawdata ..." << endl;
    _isProcessData = true;

    streamlog_out(MESSAGE3) << "Taki end of warm up period ..." << endl;
    _isWarmUp = false; 
  }
  
  // Increment the loop counter
  streamlog_out(MESSAGE2) << "Pedestal loop: " << _iPedeLoop << endl;
  ++_iPedeLoop;
    
}

//
// Method doing all raw data processing for TAKI data
//
void DEPFETPedestalNoiseProcessor::processRawDataTAKI(LCEvent * evt) {
  
  if ( ( _nEvt > 0 ) && (_nEvt % _nPedeEvents == 0) ) finalizePedestalsTAKI( ); 
  
  try {
    
    LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionName ));
        
    for ( size_t iDetector = 0; iDetector < collectionVec->size() ; iDetector++) {
         
      // Get the TrackerRawData object from the collection for this plane
      TrackerRawDataImpl *trackerRawData = dynamic_cast < TrackerRawDataImpl * >(collectionVec->getElementAt (iDetector));
      ShortVec adcValues = trackerRawData->getADCValues ();
        
      // Get number of pixels of current module 
      CellIDDecoder< TrackerRawDataImpl > rawDataDecoder( collectionVec );
      MatrixDecoder matrixDecoder( rawDataDecoder, trackerRawData);  
      int noOfXPixels = rawDataDecoder( trackerRawData ) ["uMax"] + 1; 
      int noOfYPixels = rawDataDecoder( trackerRawData ) ["vMax"] + 1;
      int sensorID = rawDataDecoder( trackerRawData ) ["sensorID"]; 
      
      // Veto for data frame that should not be corrected and saved to lcio 
      _isFrameValid[iDetector] = true;
      
      // No calibrated data in warm up     
      if  ( _isWarmUp ) _isFrameValid[iDetector] = false;       
      
      // Loop over rows of sensor matrix  
      for (int iRow = 0; iRow < noOfYPixels ; iRow ++) {
        
        // Calibrate TAKI pixel signals 
        if ( _isProcessData ) {
           
          // Common mode corrections   
          float CommonMode = 0.;
          float pixelSum = 0.;
          int nGoodPixels = 0;
          FloatVec correctedRow;
          
          for (int iCol = 0; iCol < noOfXPixels; iCol++) {
                   
            int iPixel = matrixDecoder.getIndexFromXY(iCol, iRow); 
            
            // Subtract pedestal from rawdata     
            float I = adcValues[iPixel] - _pedestal[iDetector][iPixel];  
            
            // Only use good pixel in common mode  
            if ( _status[iDetector][iPixel] == 0 ) 
            {   
              correctedRow.push_back( I );
              pixelSum += I;
              ++nGoodPixels;
            }      
          }
                      
          // So let's calculate the common mode 
          if ( nGoodPixels != 0 ) {
	    sort (correctedRow.begin(), correctedRow.end());
	    if (nGoodPixels%2 == 0) CommonMode= (correctedRow[(nGoodPixels/2)-1] + correctedRow[(nGoodPixels/2)])/2;
	    else CommonMode=correctedRow[((nGoodPixels-1)/2)];
          }
          
          if (!_useCMC)  CommonMode = 0;      
          
          // 2nd loop, subract common mode offset
          for (int iCol = 0; iCol < noOfXPixels; iCol++) {
                  
            int iPixel = matrixDecoder.getIndexFromXY(iCol, iRow);
            _commonMode[iDetector][iPixel] = CommonMode;  
            _correctedData[iDetector][iPixel] = adcValues[iPixel] - _pedestal[iDetector][iPixel] - CommonMode;
              
            if ( _correctedData[iDetector][iPixel] > _hitThresholdPre * _noise[iDetector][iPixel] )  {
              _hitCounter[iDetector][iPixel]++;  
            }  
          }
        }
         
        // Calculate calibration constants
        
        // Common mode variables  
        float CommonMode = 0;    
        float pixelSum = 0.;
        int nGoodPixels = 0;
        FloatVec correctedRow;
        
        for (int iCol = 0; iCol < noOfXPixels; iCol++) {
                   
          int iPixel = matrixDecoder.getIndexFromXY(iCol, iRow);       
          float I = adcValues[iPixel] - _pedestal[iDetector][iPixel];  
               
          // Only use good pixel in common mode  
          if ( _status[iDetector][iPixel] == 0 ) 
          {   
            correctedRow.push_back( I );
            pixelSum += I;
            ++nGoodPixels;
          }
                
        }         
        
        // So let's calculate the common mode 
        if ( nGoodPixels != 0 ) {
	  sort (correctedRow.begin(), correctedRow.end());
	  if (nGoodPixels%2 == 0) CommonMode= (correctedRow[(nGoodPixels/2)-1] + correctedRow[(nGoodPixels/2)])/2;
	  else CommonMode=correctedRow[((nGoodPixels-1)/2)];
        }
          
        if (!_useCMC)  CommonMode = 0;     
        
        // 2nd loop, subract common mode offset
        for (int iCol = 0; iCol < noOfXPixels; iCol++) {
                  
          int iPixel = matrixDecoder.getIndexFromXY(iCol, iRow);  
          float I = adcValues[iPixel] - _pedestal[iDetector][iPixel] - CommonMode;  
          
          // Check for firing pixel iff noise is available 
          bool isHit = false;
          if ( _iPedeLoop > 1 ) {
            isHit = std::abs( I ) > _hitThresholdPre * _noise[iDetector][iPixel];
          } 
          
          if ( !isHit ) {
            
            // Use this measurment for pedestal calculations    
            _tmpEntries[iDetector][iPixel]++;   
            
            // Pedestal is simple average of raw data 
            double delta = adcValues[iPixel] - _ttmpPede[iDetector][iPixel];
            _ttmpPede[iDetector][iPixel] += delta/_tmpEntries[iDetector][iPixel];            
            
            // Noise is rms of baseline fluctuations (after common mode correction)
            _ttmpNoise[iDetector][iPixel] += (_tmpEntries[iDetector][iPixel] - 1)*pow( I, 2)
                                              / _tmpEntries[iDetector][iPixel];
              
          }     
        } 
        
        // Time to review the quality of common mode 
        _rootEventNumber = evt->getEventNumber();
        _rootDetectorID = sensorID;
        _rootGate = iRow;
        _rootCMC = CommonMode; 
        
        _rootCommonModeTree->Fill();  

      } // End row loop 
           
    } // End loop on detectors

  } catch (DataNotAvailableException& e) {
        streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionName << " is not available in the current event" << endl;
  }
  
}

bool DEPFETPedestalNoiseProcessor::initializeAlgorithms(LCEvent * evt) {
  
  // Return value
  bool success = true;  
   
  // Calibration constants not available, suspend data processing 
  _isProcessData = false; 
  
  // Do not save data in warm up period, wait until process is stable 
  _isWarmUp = true; 
  
  _iPedeLoop = 0; 
  _iStatusLoop = 0; 
  
  try {
      
    LCCollectionVec * collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection ( _rawDataCollectionName ));
    _noOfDetector = collectionVec->size();
      
    for ( int iDetector = 0 ; iDetector < (int) collectionVec->size() ; ++iDetector ) {
      
      // Get the TrackerRawData object from the collection for this detector
      TrackerRawData *trackerRawData = dynamic_cast < TrackerRawData * >(collectionVec->getElementAt (iDetector));
      CellIDDecoder< TrackerRawData > rawDataDecoder( collectionVec );
      
      // Read geometry info for sensor from cellid
      _minX.push_back( 0 ) ;
      _maxX.push_back( rawDataDecoder( trackerRawData ) ["uMax"] ) ;
      _minY.push_back( 0 ) ;
      _maxY.push_back( rawDataDecoder( trackerRawData ) ["vMax"] ) ;
      _sensorIDVec.push_back( rawDataDecoder( trackerRawData ) ["sensorID"] );
      MatrixDecoder matrixDecoder( rawDataDecoder, trackerRawData);  
       
      int npixel = (int) trackerRawData->getADCValues().size();
      int noOfXPixels = rawDataDecoder( trackerRawData ) ["uMax"] + 1; 
      int noOfYPixels = rawDataDecoder( trackerRawData ) ["vMax"] + 1;
      
      if ( npixel !=  noOfXPixels* noOfYPixels) {
        cout << "Inconsistent data stream. Quit!" << endl; 
        exit(-1); 
      }   
      
      // Initialize _pedestals with zeros 
      _pedestal.push_back(FloatVec(npixel, 0.));
         
      // Initialize _noise with zeros
      _noise.push_back(FloatVec(npixel, 0.));
      
      // Initialize all pixels as GOODPIXEL
      _status.push_back(ShortVec(npixel, 0 ));
      
       
      // Open block pixel config file
      TEnv blockEnv(_blockPixelFileName.c_str());
            
      // Bad column masking 
      for (int xPixel=0; xPixel < noOfXPixels; xPixel++){         
        if (blockEnv.GetValue(Form("block_all_Column%d", xPixel ), false)){ 
          for (int yPixel = 0; yPixel < noOfYPixels; yPixel++) {
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);
            _status[iDetector][iPixel] = 1;          
          }
        }
      }  
      
      // Bad row masking 
      for (int yPixel=0; yPixel < noOfYPixels; yPixel++){         
        if (blockEnv.GetValue(Form("block_all_Row%d", yPixel ), false)){ 
          for (int xPixel = 0; xPixel < noOfXPixels; xPixel++) {
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);
            _status[iDetector][iPixel] = 1;          
          }
        }
      }  
      
      // Bad pixel masking 
      for (int xPixel=0; xPixel < noOfXPixels; xPixel++){
        for (int yPixel=0; yPixel < noOfYPixels; yPixel++){
          if (blockEnv.GetValue(Form("block_Column%d_Row%d", xPixel, yPixel ), false)){          
            int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);
            _status[iDetector][iPixel] = 1;      
          }
        }
      }
        
      _isFrameValid.push_back(true); 
      _commonMode.push_back(FloatVec(npixel, 0.));
      _correctedData.push_back(FloatVec(npixel, 0.));
      
      // Initialize temporary buffers for pedstal + noise 
      _tmpNoise.push_back(FloatVec(npixel, 0.));
      _ttmpNoise.push_back(FloatVec(npixel, 0.));
      _tmpPede.push_back(FloatVec(npixel, 0.));
      _ttmpPede.push_back(FloatVec(npixel, 0.));
      _tmpEntries.push_back(IntVec(npixel, 0));
      
      // Count hits to estimate firing frequency 
      _hitCounter.push_back(FloatVec(npixel, 0.));
      
  
    } // End of detector loop
    
  } catch (DataNotAvailableException& e) { success = false;  }
  
  return success;
     
}

   
void DEPFETPedestalNoiseProcessor::fillPedeTuple(LCEvent * evt) {

  try {
    
    LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionName ));
    for ( size_t iDetector = 0; iDetector < collectionVec->size() ; iDetector++) {
        
      // Get the TrackerRawData object from the collection for this plane
      TrackerRawDataImpl *trackerRawData = dynamic_cast < TrackerRawDataImpl * >(collectionVec->getElementAt (iDetector));
      ShortVec adcValues = trackerRawData->getADCValues ();
        
      // Get number of pixels of current module 
      CellIDDecoder< TrackerRawDataImpl > rawDataDecoder( collectionVec );
      MatrixDecoder matrixDecoder( rawDataDecoder, trackerRawData);  
      int noOfXPixels = rawDataDecoder( trackerRawData ) ["uMax"] + 1; 
      int noOfYPixels = rawDataDecoder( trackerRawData ) ["vMax"] + 1;
      int sensorID = rawDataDecoder( trackerRawData ) ["sensorID"]; 
            
      // Start looping on all pixels
      for (int yPixel = 0; yPixel < noOfYPixels; yPixel++) {
        for (int xPixel = 0; xPixel < noOfXPixels; xPixel++) {
          
          int iPixel = matrixDecoder.getIndexFromXY(xPixel, yPixel);
        
          // ROOT Output
          _rootEventNumber = evt->getEventNumber(); 
          _rootDetectorID = sensorID ;  
          _rootCol = xPixel;                
          _rootRow = yPixel; 
   
          _rootPedestal = _pedestal[iDetector][iPixel] ; 
          _rootNoise = _noise[iDetector][iPixel] ;    
          _rootStatus = _status[iDetector][iPixel];     
          _rootCMC = _commonMode[iDetector][matrixDecoder.getIndexFromXY(xPixel, yPixel)];
          _rootADC = adcValues[ matrixDecoder.getIndexFromXY(xPixel, yPixel) ]; 
          _rootDATA = _correctedData[iDetector][matrixDecoder.getIndexFromXY(xPixel, yPixel)];
          _rootHitFrequency = _hitCounter[iDetector][matrixDecoder.getIndexFromXY(xPixel, yPixel)]/_nStatusEvents;  
          _rootCycle = _iStatusLoop;
          
          _rootPedeTree->Fill();
          
        } // end loop on xPixel
      } // end loop on yPixel
      
    }  // end loop on detectors

  } catch (DataNotAvailableException& e) {
      streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionName << " is not available in the current event" << endl;
  }
    

}


void DEPFETPedestalNoiseProcessor::fillPixelTuple(LCEvent * evt) {

  try {
  
    LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionName ));
    for ( size_t iDetector = 0; iDetector < collectionVec->size() ; iDetector++) {
        
      // Get the TrackerRawData object from the collection for this plane
      TrackerRawDataImpl *trackerRawData = dynamic_cast < TrackerRawDataImpl * >(collectionVec->getElementAt (iDetector));
      ShortVec adcValues = trackerRawData->getADCValues ();
        
      // Get number of pixels of current module 
      CellIDDecoder< TrackerRawDataImpl > rawDataDecoder( collectionVec );
      int noOfXPixels = rawDataDecoder( trackerRawData ) ["uMax"] + 1; 
      int noOfYPixels = rawDataDecoder( trackerRawData ) ["vMax"] + 1;
      int sensorID = rawDataDecoder( trackerRawData ) ["sensorID"]; 
        
      // Use standard matrix encoding 
      MatrixDecoder matrixDecoder( rawDataDecoder, trackerRawData); 
      
      // Loop over monitored rows
      for (unsigned int iRow = 0; iRow < _rowMonitor.size() ; iRow++) {
        
        int yPixel = _rowMonitor[iRow]; 
        for (int xPixel = 0; xPixel < noOfXPixels; xPixel++) {
          
          _rootEventNumber = evt->getEventNumber(); 
          _rootDetectorID = sensorID ;    
          _rootCol = xPixel;                
          _rootRow = yPixel; 

          _rootPedestal = _pedestal[iDetector][matrixDecoder.getIndexFromXY(xPixel, yPixel)] ; 
          _rootNoise = _noise[iDetector][matrixDecoder.getIndexFromXY(xPixel, yPixel)] ; 
          _rootStatus = _status[iDetector][matrixDecoder.getIndexFromXY(xPixel, yPixel)] ; 
           
          _rootCMC = _commonMode[iDetector][matrixDecoder.getIndexFromXY(xPixel, yPixel)];
          _rootADC = adcValues[ matrixDecoder.getIndexFromXY(xPixel, yPixel) ];
          _rootDATA = _correctedData[iDetector][matrixDecoder.getIndexFromXY(xPixel, yPixel)];           
         
          _rootPixelTree->Fill();

        }
      }
      
      // loop over monitored columns
      for (unsigned int iCol = 0; iCol < _colMonitor.size() ; iCol++) {
        
        int xPixel = _colMonitor[iCol]; 
        for (int yPixel = 0; yPixel < noOfYPixels; yPixel++) {
          
          // check if pixel in row monitor
          if ( !binary_search (_rowMonitor.begin(), _rowMonitor.end(), yPixel) ) {
            _rootEventNumber = evt->getEventNumber(); 
            _rootDetectorID = sensorID ;    
            _rootCol = xPixel;                
            _rootRow = yPixel; 
              
            _rootPedestal = _pedestal[iDetector][matrixDecoder.getIndexFromXY(xPixel, yPixel)] ; 
            _rootNoise = _noise[iDetector][matrixDecoder.getIndexFromXY(xPixel, yPixel)] ; 
            _rootStatus = _status[iDetector][matrixDecoder.getIndexFromXY(xPixel, yPixel)] ;  
             
            _rootCMC = _commonMode[iDetector][matrixDecoder.getIndexFromXY(xPixel, yPixel)];
            _rootADC = adcValues[ matrixDecoder.getIndexFromXY(xPixel, yPixel) ];
            _rootPedestal = _pedestal[iDetector][matrixDecoder.getIndexFromXY(xPixel, yPixel)] ;   
                       
            _rootPixelTree->Fill();
          }  
        }
      }

    }  // End loop on detectors

  } catch (DataNotAvailableException& e) {
      streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionName << " is not available in the current event" << endl;
  }
    

}

void DEPFETPedestalNoiseProcessor::maskBadPixel(LCEvent * evt) {

  try {
    
    LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionName ));
    for ( unsigned int iDetector = 0; iDetector < collectionVec->size(); iDetector++) {
       
      // Get the TrackerRawData object from the collection for this plane
      TrackerRawDataImpl *trackerRawData = dynamic_cast < TrackerRawDataImpl * >(collectionVec->getElementAt (iDetector));
        
      // Get number of pixels of current module 
      CellIDDecoder< TrackerRawDataImpl > rawDataDecoder( collectionVec );
      
      
      // Use standard matrix encoding 
      MatrixDecoder matrixDecoder( rawDataDecoder, trackerRawData);  
      
      // Count masked pixels 
      int nMasked = 0; 

      // Loop over all pixels
      for ( int iPixel = 0; iPixel < (int) _status[iDetector].size(); iPixel++ ) { 
           
          // Mask pixel with very low noise  
          if (  _noise[iDetector][iPixel] < _pixelMaskLowerNoise ) {
          
            
            _status[iDetector][iPixel] = 1;
            nMasked++; 
            streamlog_out ( MESSAGE3 ) <<  "DeadPixel masking of pixel number " << iPixel
                                       << " on detector " << _sensorIDVec.at( iDetector )
                                       << " (" << _noise[iDetector][iPixel] << ")" << endl;
          }
           
          // Mask pixel with very high noise  
          if ( _noise[iDetector][iPixel] > _pixelMaskUpperNoise ) {
           
            
            _status[iDetector][iPixel] = 1;
            nMasked++;
            streamlog_out ( MESSAGE3 ) <<  "NoisyPixel masking of pixel number " << iPixel
                                       << " on detector " << _sensorIDVec.at( iDetector )
                                       << " (" << _noise[iDetector][iPixel] << ")" << endl;
          }
                     
          // Mask pixel with too high or low pedestals 
          if ( ( _pedestal[iDetector][iPixel] > _pixelMaskUpperPede) ||
               ( _pedestal[iDetector][iPixel] < _pixelMaskLowerPede ) ) {
           
            
            _status[iDetector][iPixel] = 1;
            nMasked++;
            streamlog_out( MESSAGE3 ) << "BrigthPixel masking of pixel number " << iPixel
                                      << " on detector " <<  _sensorIDVec.at( iDetector )
                                      << " (" << _pedestal[iDetector][iPixel] <<  ")" << endl;
          }
           
          // Mask pixel with very high hit frequency -> hot pixel killer 
          double firingFreq =  _hitCounter[ iDetector ][ iPixel ] / _nStatusEvents;
          if ( firingFreq  > _maxFiringFreq ) {
           
            
            _status[iDetector][iPixel] = 1;
            nMasked++;
            streamlog_out( MESSAGE3 ) << "HotPixel masking of pixel number " << iPixel
                                      << " on detector " << _sensorIDVec.at( iDetector ) 
                                      << " (" << firingFreq <<  ")" << endl;  
          }
                     
      } // Pixel loop
      
      streamlog_out( MESSAGE3 ) << endl
                                << "Total of masked pixels is " << nMasked
                                << " on detector " << _sensorIDVec.at( iDetector ) 
                                << endl << endl;      
      
    } // Detector loop
    
    // Sample pedestals for monitoring 
    if ( _iStatusLoop < 6 ) { 
      streamlog_out(MESSAGE2) << "Monitor pedestal data ..." << endl;
      fillPedeTuple(evt);
    }
    
    // Reset counters for next loop
    _hitCounter.clear();
    for ( size_t iDetector = 0; iDetector < _noOfDetector; iDetector++) {
      _hitCounter.push_back(FloatVec(_noise[iDetector].size(), 0.));
    }     
    
    // Now, full set of calibration constants is available :) 
    _isWarmUp = false; 
        
  } catch (DataNotAvailableException& e) {
      streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionName << " is not available in the current event" << endl;
  }  
  
  // Increment the loop counter
  streamlog_out(MESSAGE2) << "Status loop: " << _iStatusLoop << endl;
  ++_iStatusLoop; 
  
}


void DEPFETPedestalNoiseProcessor::sparsifyEvent(LCEvent * evt) {
   
  try {
    
    LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionName ));
    
    // Initialize sparse data collection
    LCCollectionVec * sparseDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);
    CellIDEncoder<TrackerDataImpl> sparseDataEncoder( "sensorID:6,sparsePixelType:5" , sparseDataCollection );  
    
    for ( size_t iDetector = 0; iDetector < collectionVec->size() ; iDetector++) {
         
      // Get the TrackerRawData object 
      TrackerRawDataImpl * matrix = dynamic_cast < TrackerRawDataImpl * >(collectionVec->getElementAt (iDetector));
      CellIDDecoder< TrackerRawDataImpl > rawDataDecoder( collectionVec );
      
      MatrixDecoder matrixDecoder( rawDataDecoder, matrix);  
      int sensorID = static_cast<int > (rawDataDecoder(matrix)["sensorID"]); 
       
      if ( _sensorIDVec.at( iDetector ) !=  sensorID   ) {
        streamlog_out ( ERROR4 ) << "Inconsistent Sensor ID's" << endl;
        exit(-1); 
      } 
      
      // Prepare a TrackerData to store zs pixels
      TrackerDataImpl* zspixels = new TrackerDataImpl;
      
      // Set description for zspixels 
      sparseDataEncoder["sensorID"] = sensorID;
      sparseDataEncoder["sparsePixelType"] = 0;
      sparseDataEncoder.setCellID( zspixels );
      
      // Cound firing pixels 
      int nFiring = 0; 

      // Check if frame valid
      if ( _isFrameValid[iDetector] ) {
         
        for ( int iPixel = 0; iPixel < (int) _correctedData[iDetector].size(); iPixel++ ) {
            
          if (  _status[iDetector][iPixel]  == 0 ) {
            
            float data  = _correctedData[iDetector][iPixel];
                 
            // Use a global ADU threshold 
            float threshold = _hitThresholdZS; 
            
            // Use a global SNR threshold
            if ( _useSNRCut ) threshold *= _noise[iDetector][iPixel]; 
             
            if ( data >= threshold  ) {
              
              int xPixel,yPixel; 
              matrixDecoder.getXYFromIndex(iPixel,  xPixel, yPixel);
              
              // Store pixel data int EUTelescope format 
              zspixels->chargeValues().push_back( xPixel );
              zspixels->chargeValues().push_back( yPixel );
              zspixels->chargeValues().push_back( data );    
               
              nFiring++;
               
              // Print detailed pixel summary 
              streamlog_out(MESSAGE1) << "Found pixel Nr. " << iPixel << " on sensor " << sensorID << endl
                                      << "   x:" << xPixel << ", y:" << yPixel << ", charge:" << data
                                      << ", quality:" <<   _status[iDetector][iPixel] << endl;
               
              
               
            }  // Endif threshold 
          }  // Endif status
        } // Endfor ipixel
      } // Endif valid frame
      
      // Add zs pixels to collection 
      sparseDataCollection->push_back( zspixels );
      
      // ROOT OUTPUT
      _rootEventNumber = evt->getEventNumber();
      _rootDetectorID = sensorID;
      _rootNFiring = nFiring;               

      _rootEventTree->Fill();
       
    }  // Endfor detectors
    
    evt->addCollection(sparseDataCollection, _zsDataCollectionName );
    
  } catch (DataNotAvailableException& e) {
      streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionName << " is not available in the current event" << endl;
  }

}

void DEPFETPedestalNoiseProcessor::calibrateEvent(LCEvent * evt) {
  
  try {
    
    LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionName ));
    LCCollectionVec * correctedDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);
    
    for ( size_t iDetector = 0; iDetector < collectionVec->size() ; iDetector++) {
         
      // Get the TrackerRawData object from the collection for this plane
      TrackerRawDataImpl *trackerRawData = dynamic_cast < TrackerRawDataImpl * >(collectionVec->getElementAt (iDetector));
      CellIDDecoder< TrackerRawDataImpl > rawDataDecoder( collectionVec );
      MatrixDecoder matrixDecoder( rawDataDecoder, trackerRawData);  
      int sensorID = static_cast<int > (rawDataDecoder(trackerRawData)["sensorID"]); 
      
      TrackerDataImpl  * corrected = new TrackerDataImpl;
      CellIDEncoder<TrackerDataImpl> idDataEncoder(DEPFET::MATRIXDEFAULTENCODING, correctedDataCollection);
      idDataEncoder["sensorID"] =  sensorID;
      idDataEncoder["uMax"]     = static_cast<int > (rawDataDecoder(trackerRawData)["uMax"]);
      idDataEncoder["vMax"]     = static_cast<int > (rawDataDecoder(trackerRawData)["vMax"]);
      idDataEncoder.setCellID(corrected);
          
      if ( _sensorIDVec.at( iDetector ) !=  sensorID   ) {
        streamlog_out ( ERROR4 ) << "Inconsistent Sensor ID's" << endl;
        exit(-1); 
      } 
      
      if ( _isFrameValid[iDetector] ) corrected->setChargeValues(_correctedData[iDetector]);
      else corrected->setChargeValues( FloatVec(_correctedData[iDetector].size(), 0) );
      
      correctedDataCollection->push_back(corrected);
      
    }  // End loop on detectors
    
    evt->addCollection(correctedDataCollection, _calibratedDataCollectionName);
    
  } catch (DataNotAvailableException& e) {
      streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionName << " is not available in the current event" << endl;
  }
  
}

} // Namespace

