// /////////////////////////////////////////////////////////////////////////////////////// //
//                                                                                         //
//    SimpleNoiseDBCreator - Marlin Processor                                           //
// /////////////////////////////////////////////////////////////////////////////////////// //


#ifndef SimpleNoiseDBCreator_H
#define SimpleNoiseDBCreator_H 1

#include "TBDetector.h"

// Include LCIO classes
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <UTIL/CellIDEncoder.h>

// Include Marlin classes
#include <marlin/Processor.h>
#include <marlin/Exceptions.h>

namespace depfet {

//! SimpleNoiseDBCreator Processor
  /*! Create a simple NoiseDB from steering parameters
   *  
   *  This processor is intended to create prepare status 
   *  pedestal and noise collections for simulation data. 
   * 
   *  For real data, prefer data driven processors like 
   *  the HotPixelKiller. 
   *
   *  Author: B.Schwenker, Universität Göttingen
   *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
   */

class SimpleNoiseDBCreator : public marlin::Processor {

 public:

   //!Method that returns a new instance of this processor
   virtual Processor*  newProcessor() { return new SimpleNoiseDBCreator ; }

   //!Constructor - set processor description and register processor parameters
   SimpleNoiseDBCreator();
   
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
   
   //! Pedestal collection name
   /*! Input pedestal collection name. Default value is pedestal.
    */
   std::string _pedestalCollectionName;
   
   //! Noise collection name
   /*! Input noise collection name. Default value is noise.
    */
   std::string _noiseCollectionName;
   
   //! Status collection name
   /*! Input status collection name. Default value is status.
    */
   std::string _statusCollectionName;
   
   //! Initial pedestal value
   /*! This vector of floats is used to store the initial values of
    *  pedestal. The size of this vector should match the number of
    *  detectors in the telescope setup. This vector is provided by
    *  the user from steering file.
    */
   FloatVec _initPedestal;
   
   //! Initial noise value
   /*! This vector of floats is used to store the initial values of
    *  noise. The size of this vector should match the number of
    *  detectors in the telescope setup. This vector is provided by
    *  the user from steering file.
    */
   FloatVec _initNoise;
   
       
 private:
    
   //! Run number
   int _iRun;
     
   //! Event number
   int _iEvt;
   
   //! Pedestal collection
   IMPL::LCCollectionVec * _pedestalCollectionVec;
   
   //! Noise collection
   IMPL::LCCollectionVec * _noiseCollectionVec;
   
   //! Status collection
   IMPL::LCCollectionVec * _statusCollectionVec;
   
   // Handle to detector data sheets 
   TBDetector _detector;       
     
}; // Class

} // Namespace

#endif 



