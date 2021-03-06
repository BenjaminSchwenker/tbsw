/////////////////////////////////////////////////////////  //
//                                                         //
//    PixelChargeCalibrator - Marlin Processor             //
/////////////////////////////////////////////////////////  //

#ifndef PixelChargeCalibrator_H
#define PixelChargeCalibrator_H 1



// Include LCIO classes
#include <lcio.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include basic C
#include <string>


// Include ROOT classes
#include <TH2F.h>
#include <TF1.h>

namespace depfet {

  /*
   * PixelChargeCalibrator:
   * Applies pixel specific charge calibration function to 
   * sparsified data for all sensors in collection. 
   *
   * Optional: NoiseDB file generated by i.e. HotPixelKiller 
   * because calibration of noise is not meaningfull. 
   *
   * Needs a root file containing a TF1 function per sensor 
   * and parameters for the funtions stored in TH2F histograms 
   * binned in pixel IDs (i.e. column, row).
   *   Base function name is a processor parameter, full name:
   *     BaseFuncName + _ + sensorID
   *  
   *   The name of the TH2F parameter histograms is expected to be: 
   *     CalibParaBaseName + _ + sensorID + _ + X where X is the number of the parameter. 
   * 
   * Plusses and spaces only for simpler display, not in name! 
   *
   * A threshold can be applied to the uncalibrated charge. 
   * Usefull if threshold in detector specific units is easier then calibrated values. 
   * 
   * ! source/tbsw/calibFuncExample.py
   * !   Script generates a root file with the expected structure 
   * !   and names of the function and parameters in the gainCalibFile. 
   *
   *
   * Author: Helge Christoph Beck
   * Email: helge-christoph.beck@phys.uni-goettingen.de
   *
   */
  
  class PixelChargeCalibrator : public marlin::Processor {
   
   public:
    
    //!Method that returns a new instance of this processor
    virtual Processor*  newProcessor() { return new PixelChargeCalibrator ; }

    //!Constructor - set processor description and register processor parameters
    PixelChargeCalibrator();
    
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
      
    // Called by the processEvent() once for the calibration
    // of sparsified pixels. 
    void calibrate( LCEvent * evt , LCCollectionVec * calibCollection );  
               
    //!Method printing processor parameters
    void printProcessorParams() const;
      
    // PROCESSOR PARAMETERS
    
    //! Input sparsified data collection name
    std::string _sparseDataCollectionName;
    
    //! Output hit calibrated collection name
    std::string _calibratedCollectionName;
    
    //! Name of clusterDB file 
    std::string  _noiseDBFileName;

    //! Name of calibrationDB file
    std::string _gainCalibDBFileName;

    //! Name of calibration function
    std::string _baseCalibFuncName;

    //! Base name of the parameters of the calibration function
    std::string _calibParaBaseName;

    //! Minimimum signal for zero suppression, applied before calibration
    float _sparseZSCut;
       	       
  private: 
   
   //! internally used as storage for input decoding
   UTIL::BitField64 _inputDecodeHelper;
   CellIDEncodeConstructHelper _calibOutputEncoderHelper;
   // Pixel mask to filter brocken (hot) channels  
   // Key is sensorID 
   std::map< int, TH2F * > _DB_Map_Mask;

   // Calibration functions
   // Key is sensorID
   std::map< int, TF1 * > _DB_Map_CalibFunc;

   // Number of parameters for calibration functions
   // Key is sensorID
   std::map< int, int > _DB_Map_NparFunc;

   // Parameter for calibration function
   // First key is sensorID, Second key is parameter number
   std::map< int, std::map< int, TH2F * > > _DB_Map_CalibPar;

   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
   
}; // Class

} // Namespace

#endif 

