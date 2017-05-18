// ///////////////////////////////////////////////////////////////////////////////////////     //
//                                                                                             //
//    DEPFETPedestalNoiseProcessor - Marlin Processor - Preprocessing for DEPFET-DCDB raw data //
// ///////////////////////////////////////////////////////////////////////////////////////     //


#ifndef DEPFETPedestalNoiseProcessor_H
#define DEPFETPedestalNoiseProcessor_H 1

// Include LCIO classes
#include <EVENT/LCParameters.h>
#include <IMPL/LCCollectionVec.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include <marlin/Exceptions.h>

// Include ROOT classes
#include <TEnv.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>

// Include basic C
#include <string>
#include <map>
#include <vector>



namespace depfet {

  //! Pedestal and noise processor for DEPFET pixel sensors
  /*! This processor performs the pre-processing of DEPFET-DCDB or S3B
   *  or TAKI raw data. In particular, the DEPFETPedestalNoiseProcessor
   *  serves as a software simulation of the Data Handling Hybrid (DHP) 
   *  planned for the DEPFET based Belle II pixel detector. The pre-
   *  processing consists of pedestal subtraction, common mode correction
   *  and zero suppression of full analog raw data from ADC chips. 
   *  The required calibration constants (pedestals, noise, status)
   *  are directly calculated from the incoming raw data stream.  
   *  
   *  The input data is organized in a collection of type TrackerRawData 
   *  named 'rawdata'. This collection has as many elements as DEPFET 
   *  modules in the experimental setup. Each element stores the full frame 
   *  'raw' ADC codes from the DCDB of all DEPFET pixels on the module.
   *  
   *  The output of the processor is a list of firing pixels for each 
   *  sensor module organized in a collection of type TrackerData. 
   *  Optionally, the processor can provide full frames of corrected 
   *  pixel signals (after pedestal and common mode corrections). The 
   *  processor produces a seperate lcio file with status, pedestal and 
   *  noise maps for each module.  
   *  
   *  The user can select multiple preprocessing methods. The most important
   *  method is labeled 'DHP'. Here a brief description: The DHP performs a 
   *  'static' pedestal subtraction for each pixel. Pedestal and noise values 
   *  are re-calculated every 'nPedeEvents' events from raw data. Afterwards, 
   *  the DHP performs a 2-pass common mode correction. The 1st pass common mode 
   *  is a mean from all pixels in a DEPFET gate. The second pass blocks firing 
   *  pixels. Finally, a hit threshold (~4xnoise) is applied to the corrected 
   *  signal to find firing pixels. Bad pixels are suppressed. 
   *  
   *  The status matrix is calculated after 'nStatusEvent' events. The criteria 
   *  to identify bad pixels are: 
   *  
   *  A) Noise: Noise value (rms around baseline) is too high or too low compared
   *  		with module wise noise mean. 		
   *  
   *  B) Firing Frequency: The pixels mean firing rate is higher than expected from 
   *  		pixel noise (and tracks).  		        		
   *  
   *  C) Pedestal: The pixel pedestal is too low (dead transistor) or too high (dynamic
   *  		range of ADC) for proper working.
   *  
   *  In addition, the user can supply a list of bad pixels in a seperate 
   *  config file.
   *  
   *  A seperate root file is created to monitor the quality of preprocessing and 
   *  stability of the system over time. 
   *  
   *  Author: Benjamin Schwenker, Göttingen University 
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de> 
   *  
   */


class DEPFETPedestalNoiseProcessor : public marlin::Processor {
   
 public:
   
//!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new DEPFETPedestalNoiseProcessor ; }
   
//!Constructor - set processor description and register processor parameters
   DEPFETPedestalNoiseProcessor();
   
//!Destructor - clean up root histos
   ~DEPFETPedestalNoiseProcessor();
   
//!Method called at the beginning of data processing - used for initialization
   virtual void init();
   
//!Method called for each run - used for run header processing
   virtual void processRunHeader(LCRunHeader * run);
   
//!Method called for each event - used for event data processing
   virtual void processEvent(LCEvent * evt);
   
//!Method called after each event - used for data checking
   virtual void check(LCEvent * evt);
   
//!Method called after all data processing
   virtual void end();
   
protected:
   
//! Process raw data from DCDB
/*! This method is called in each event to correct pixel raw data 
 *  coming from DCDB (or CURO), i.e. apply pedestal and common mode corrections
 *  in DHP style.  
 */
   void processRawDataDHP(LCEvent * evt);
   
//! Finalize pedestal/noise DHP 
/*! This method is called at the end of every pedestal/noise loop for DHP. It  
 *  provides new pixel pedestal and noise data. 
 */
   void finalizePedestalsDHP();

//! Process raw data from DCDB
/*! This is a variation of DHP preprcessing. It performs the 
 *  common mode correction first. 
 */
   void processRawDataDHP2(LCEvent * evt);
   
//! Finalize pedestal/noise DHP 
/*! This is a variation of DHP preprcessing. It performs the 
 *  common mode correction first. 
 */
   void finalizePedestalsDHP2();
   
//! Process raw data from TAKI 
/*! This method is called in each event to correct pixel raw data 
 *  coming from TAKI, i.e. apply static pedestal and common mode 
 *  corrections in standard test beam style.  
 */
   void processRawDataTAKI(LCEvent * evt);

//! Finalize pedestal/noise TAKI
/*! This method is called at the end of every pedestal/noise loop for TAKI. It  
 *  provides new pixel pedestal and noise data. 
 */
   void finalizePedestalsTAKI();
   
//! Calibrate event 
/*! Add calibrated data collection to current LCIO::Event 
 */
   void calibrateEvent(LCEvent * evt);
   
//! Sparsify event 
/*! Add 0-suppressed pixel signals to current LCIO::Event 
 */
   void sparsifyEvent(LCEvent * evt);
   
//! Mask Pixel Quality
/*! This method is called at the end of every status cycle and is used to
 *  update pixel status. 
 */
   void maskBadPixel(LCEvent * evt);
   
//! This method is used to fill pedestal ntuple   
   void fillPedeTuple(LCEvent * evt);
   
//! This method is used to fill pixel ntuple   
   void fillPixelTuple(LCEvent * evt); 
   
//! This method is used to initialize algorithms  
   bool initializeAlgorithms(LCEvent * evt); 
   
//! Method printing processor parameters
   void printProcessorParams() const;
   
// Processor Parameters 
   
//! Input rawdata collection name
   std::string _rawDataCollectionName; 
   
//! Output calibrated data collection name
   std::string _calibratedDataCollectionName;  
   
//! Output calibrated data collection name
   std::string _zsDataCollectionName; 
   
//! Name of block pixel config file
   std::string _blockPixelFileName;
   
//! Name of output noiseDB
   std::string _outputNoiseDB; 
   
//! Algorithm for pedestal and common mode correction
   std::string _pedestalAlgorithm;

//! Apply common mode correction  
   bool _useCMC;     

//! Update pedestal/noise map after n events 
   int _nPedeEvents; 
   
//! Update status map after n events 
   int _nStatusEvents; 
   
//! Hit threshold for zero suppression  
   float _hitThresholdZS;

//! Hit threshold for preprocessing   
   float _hitThresholdPre;
   
//! Use SNR cut in zero suppression
   bool _useSNRCut;
   
// 
// BAD PIXEL MASKING 
//    

//! Upper noise cut in absolute ADC values
   float _pixelMaskUpperNoise;
   
//! Lower noise cut in absolute ADC values
   float _pixelMaskLowerNoise;
   
//! Upper pedestal cut in absolute ADC values
   float _pixelMaskUpperPede;
   
//! Lower pedestal cut in absolute ADC value
   float _pixelMaskLowerPede;
    
//! Maximum firing frequency
/*! This is the maximum allowed frequency (in %) at which a good 
 *  pixel can fire
 */
   float _maxFiringFreq;
   
//! N-fold readout means Nfold common mode correction 
   int _nFold;
   
//! ROOT output file name  
   std::string _rootFileName;  
   
//! Monitoring of raw ADC signals from pixels in selected rows
   IntVec _rowMonitor;        
   
//! Monitoring of raw ADC signals from pixels in selected cols
   IntVec _colMonitor;          
   
 private:
   
   //! ROOT_Output
   TFile * _rootFile;
   TTree * _rootPedeTree; 
   TTree * _rootPixelTree;
   TTree * _rootCommonModeTree;
   TTree * _rootEventTree;  
   
   int _rootEventNumber; 
   int _rootDetectorID; 
   int _rootGate; 
   int _rootCol;                
   int _rootRow;  
   double _rootPedestal;
   double _rootNoise;    
   int _rootStatus;     
   double _rootCMC;
   double _rootADC;
   double _rootDATA;
   int _rootCycle; 
   double _rootHitFrequency;
   int _rootNFiring;  
   
   // Initialization successfull? 
   bool _isAlgoInit; 
   
   // Calibration constants ready? 
   bool _isProcessData; 
   
   // Do not save data in warm up
   bool _isWarmUp; 
      
   // Number of DEPFET modules to be processed 
   size_t _noOfDetector; 
   
   //! The sensorID vector
   std::vector< int > _sensorIDVec;
   
   //! First pixel along X
   IntVec _minX; 
   
   //! Last pixel along X
   IntVec _maxX; 
   
   //! First pixel along Y
   IntVec _minY; 
   
   //! Last pixel along Y
   IntVec _maxY;  
   
   // Final pedestal, noise and status for data processing 
   std::vector < FloatVec > _pedestal;
   std::vector < FloatVec > _noise;
   std::vector < ShortVec > _status;
      
   // Store corrected data to lcio file 
   std::vector < bool > _isFrameValid;  
   
   // Current common mode offsets   
   std::vector < FloatVec > _commonMode;
  
  // Pedestal and commom mode corrected signals 
  std::vector < FloatVec > _correctedData;
  
  // These variables are used to calculate the 
  // calibration constants 
  
  // Count pixel hits (-> mean firing frequency)
  std::vector < FloatVec > _hitCounter;
   
  // Online mean & variance calculations 
  std::vector < FloatVec > _tmpNoise;
  std::vector < FloatVec > _ttmpNoise;
  std::vector < FloatVec > _tmpPede;
  std::vector < FloatVec > _ttmpPede;
  std::vector < IntVec > _tmpEntries;
   
  int _iPedeLoop; 
  int _iStatusLoop;
  
  double _timeCPU; //!< CPU time
  int    _nRun ;   //!< Run number
  int    _nEvt ;   //!< Event number
       
}; // Class

} // Namespace

#endif 

