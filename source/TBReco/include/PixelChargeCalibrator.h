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
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerPulseImpl.h>

// Include Marlin classes
#include <marlin/Global.h>
#include <marlin/Processor.h>
#include "marlin/Exceptions.h"

// Include basic C
#include <vector>
#include <string>


// Include ROOT classes
#include <TH2F.h>
#include <TF1.h>

namespace depfet {

  //typedef std::vector<FloatVec> Pix_GroupVector; not needed?
  
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
    
    //! Output hit cluster collection name
    std::string _calibratedCollectionName;
    
    //! Output original pixel data collection name // needed?
    //std::string _dummyCollectionName;
    
    //! Name of clusterDB file 
    std::string  _noiseDBFileName;

    //! Name of calibrationDB file
    std::string _calibDBFileName;

    //! Name of calibration function
    std::string _calibFuncName;

    //! Base name of the parameters of the calibration function
    std::string _calibParaBaseName;

    // Are cuts needed, clusterizer needs appropriate calibrated cut!
    //! Minimimum signal for zero suppression
    float _sparseZSCut;
    
    //! Minimimum signal for seed pixel in clusters
    //float _sparseSeedCut;
    
    //! Minimum signal for total cluster charge
    //float _sparseClusterCut;

    //! Plane number of detector to be calibrated
    int _idut;
       	       
  private: 
   
   //! internally used as storage for input decoding
   UTIL::BitField64 _inputDecodeHelper;
   //CellIDEncodeConstructHelper _orginalOutputEncoderHelper;
   CellIDEncodeConstructHelper _calibOutputEncoderHelper;
   // Pixel mask to filter brocken (hot) channels  
   // Key is sensorID 
   std::map< int, TH2F *> _DB_Map_Mask;

   // Calibration function
   TF1 * _calibFunc;

   // Number of parameters for calibration function
   int _nparFunc;

   // Parameter for calibration function
   // Key is parameter number
   std::map< int, TH2F *> _DB_MAP_CalibPar;

   double _timeCPU; //!< CPU time
   int    _nRun ;   //!< Run number
   int    _nEvt ;   //!< Event number
   
}; // Class

} // Namespace

#endif 

